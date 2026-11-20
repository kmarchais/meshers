# Copyright (c) 2026 Kevin Marchais

"""Scalar expression tracing and first derivatives for the embedded native compiler."""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING, Literal, Never, TypedDict

import numpy as np
from numpy.lib.mixins import NDArrayOperatorsMixin

from . import _meshers
from ._meshers import compile_graph
from ._options import check_keywords

if TYPE_CHECKING:
    from collections.abc import Callable

_Node = tuple[str, tuple[int, ...], float]
_MIN_CONSTRAINTS = 2
_MAX_NODES = 20_000
_AXES = "xyz"
_INTERP_NODES_PER_KNOT = 8
_TERNARY_ARITY = 3


_UFUNCS = {
    np.arctan2: "atan2",
    np.arccos: "acos",
    np.arcsin: "asin",
    np.arctan: "atan",
    np.sin: "sin",
    np.cos: "cos",
    np.tan: "tan",
    np.exp: "exp",
    np.log: "ln",
    np.sqrt: "sqrt",
    np.tanh: "tanh",
    np.sinh: "sinh",
    np.cosh: "cosh",
    np.absolute: "abs",
    np.negative: "neg",
    np.add: "add",
    np.subtract: "sub",
    np.multiply: "mul",
    np.divide: "div",
    np.power: "pow",
    np.minimum: "minimum",
    np.maximum: "maximum",
    np.less: "lt",
    np.less_equal: "le",
    np.greater: "gt",
    np.greater_equal: "ge",
    np.equal: "eq",
    np.not_equal: "ne",
}


class LegacyCompileOption(TypedDict, total=False):
    """Original keyword accepted by released versions of the mesh generators."""

    compile: bool


def resolve_compile_option(
    legacy: LegacyCompileOption,
    *,
    setting: bool | None,
    stacklevel: int = 3,
    schema: type = LegacyCompileOption,
) -> bool:
    """Resolve the automatic-compilation setting for a mesh generator.

    Returns:
        Whether to trace callable fields.

    Raises:
        TypeError: An unknown option or both spellings were supplied.
    """
    check_keywords(legacy, schema)
    if "compile" not in legacy:
        return True if setting is None else setting
    if setting is not None:
        msg = "Specify only one of compile_expression and compile"
        raise TypeError(msg)
    warnings.warn(
        "compile is deprecated; use compile_expression",
        DeprecationWarning,
        stacklevel=stacklevel,
    )
    return legacy["compile"]


class TraceError(TypeError):
    """A callable uses operations outside the supported expression subset."""


class _Expression(NDArrayOperatorsMixin):
    __slots__ = ("graph", "index")
    __array_priority__ = 1000
    __hash__ = None

    def __init__(self, graph: _Graph, index: int) -> None:
        self.graph, self.index = graph, index

    def _binary(self, op: str, other: object) -> _Expression:
        return self.graph.binary(op, self, self.graph.convert(other))

    def __add__(self, other: object) -> _Expression:
        return self._binary("add", other)

    def __radd__(self, other: object) -> _Expression:
        return self.graph.convert(other)._binary("add", self)

    def __sub__(self, other: object) -> _Expression:
        return self._binary("sub", other)

    def __rsub__(self, other: object) -> _Expression:
        return self.graph.convert(other)._binary("sub", self)

    def __mul__(self, other: object) -> _Expression:
        return self._binary("mul", other)

    def __rmul__(self, other: object) -> _Expression:
        return self.graph.convert(other)._binary("mul", self)

    def __truediv__(self, other: object) -> _Expression:
        return self._binary("div", other)

    def __pow__(self, other: object) -> _Expression:
        return self._binary("pow", other)

    def __pos__(self) -> _Expression:
        return self

    def __bool__(self) -> Never:
        msg = "data-dependent Python branches cannot be traced; use the NumPy callback path"
        raise TraceError(msg)

    # Traced values are immutable. Augmented assignments rebind the caller's variable.
    __iadd__ = __add__
    __isub__ = __sub__
    __imul__ = __mul__
    __itruediv__ = __truediv__
    __ipow__ = __pow__

    def __float__(self) -> Never:
        msg = "scalar conversion cannot be traced"
        raise TraceError(msg)

    def __int__(self) -> Never:
        msg = "integer conversion cannot be traced"
        raise TraceError(msg)

    def __len__(self) -> Never:
        msg = "array lengths require the NumPy callback path"
        raise TraceError(msg)

    def __getitem__(self, key: object) -> Never:
        msg = "array indexing requires the NumPy callback path"
        raise TraceError(msg)

    def __getattr__(self, name: str) -> Never:
        if name.startswith("__"):
            raise AttributeError(name)
        msg = f"array attribute {name} requires the NumPy callback path"
        raise TraceError(msg)

    def __array_function__(
        self,
        func: Callable,
        types: tuple[type, ...],
        args: tuple,
        kwargs: dict,
    ) -> _Expression:
        if func is np.interp:
            return _Interpolation(*args, **kwargs).trace(self.graph)
        if func is np.where and len(args) == _TERNARY_ARITY and not kwargs:
            condition, yes, no = args
            return self.graph.select(self.graph.convert(condition), yes, no)
        if func is np.clip and len(args) == _TERNARY_ARITY and not kwargs:
            value, lower, upper = args
            result = self.graph.convert(value)
            if lower is not None:
                result = self.graph.binary("maximum", result, self.graph.convert(lower))
            if upper is not None:
                result = self.graph.binary("minimum", result, self.graph.convert(upper))
            return result
        name = getattr(func, "__name__", type(func).__name__)
        msg = f"unsupported traced NumPy function: {name}"
        raise TraceError(msg)

    def __array__(self, *args: object, **kwargs: object) -> Never:
        msg = "array construction is unsupported by the Rust scalar tracer"
        raise TraceError(msg)

    def __array_ufunc__(
        self,
        ufunc: np.ufunc,
        method: str,
        *args: object,
        **kwargs: object,
    ) -> _Expression:
        if method != "__call__" or kwargs or ufunc not in _UFUNCS:
            msg = f"unsupported traced operation: {ufunc.__name__}.{method}"
            raise TraceError(msg)
        op = _UFUNCS[ufunc]
        first = self.graph.convert(args[0])
        if ufunc.nin == 1:
            return self.graph.unary(op, first)
        return self.graph.binary(op, first, self.graph.convert(args[1]))


@dataclass(frozen=True, slots=True)
class _Interpolation:
    """Arguments of a NumPy interpolation call, in NumPy's positional order."""

    x: object
    xp: object
    fp: object
    left: float | None = None
    right: float | None = None
    period: float | None = None

    def trace(self, graph: _Graph) -> _Expression:
        # Match NumPy's signature; only finite, static, increasing real tables trace.
        if self.period is not None:
            msg = "np.interp with period requires the NumPy callback path"
            raise TraceError(msg)
        try:
            knots, values = np.asarray(self.xp), np.asarray(self.fp)
        except (TypeError, ValueError) as error:
            msg = "np.interp requires constant real tables"
            raise TraceError(msg) from error
        if (
            any(table.ndim != 1 or table.dtype.kind not in "biuf" for table in (knots, values))
            or not len(knots)
            or len(knots) != len(values)
        ):
            msg = "np.interp requires equal nonempty one-dimensional real tables"
            raise TraceError(msg)
        knots, values = knots.astype(np.float64), values.astype(np.float64)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            spans = np.diff(knots)
            slopes = np.diff(values) / spans
        if (
            not np.isfinite(knots).all()
            or not np.isfinite(values).all()
            or not (spans > 0).all()
            or not np.isfinite(slopes).all()
            or any(
                v is not None
                and (
                    not isinstance(v, (int, float, np.integer, np.floating))
                    or not math.isfinite(float(v))
                )
                for v in (self.left, self.right)
            )
        ):
            msg = "np.interp requires finite constants and strictly increasing knots"
            raise TraceError(msg)
        if len(graph.nodes) + _INTERP_NODES_PER_KNOT * len(knots) > _MAX_NODES:
            msg = "np.interp table exceeds the native expression graph limit"
            raise TraceError(msg)
        coordinate = graph.convert(self.x)

        def branch(lo: int, hi: int) -> float | _Expression:
            if lo == hi:
                return float(values[lo])
            if hi - lo == 1:
                return float(values[lo]) + float(slopes[lo]) * (coordinate - float(knots[lo]))
            middle = (lo + hi) // 2
            return graph.select(
                coordinate < float(knots[middle]),
                branch(lo, middle),
                branch(middle, hi),
            )

        # ponytail: eager selection costs O(table size); add a lookup kernel only
        # if profiles justify it. Existing native operations preserve the table.
        value = branch(0, len(knots) - 1)
        first, last = float(knots[0]), float(knots[-1])
        # At a knot NumPy returns the tabulated value, including when a fill is supplied.
        lower = float(values[0])
        upper = float(values[-1])
        return graph.select(
            coordinate <= first,
            lower
            if self.left is None
            else graph.select(coordinate < first, float(self.left), lower),
            graph.select(
                coordinate >= last,
                upper
                if self.right is None
                else graph.select(coordinate > last, float(self.right), upper),
                value,
            ),
        )


class _Graph:
    def __init__(self) -> None:
        self.nodes: list[_Node] = []
        self._indices: dict[tuple[str, tuple[int, ...]], int] = {}

    def _record(self, op: str, operands: tuple[int, ...], value: float = 0.0) -> _Expression:
        # Hex keys keep +0.0 and -0.0 distinct in the interning table.
        key = (float(value).hex() if op == "constant" else op, operands)
        index = self._indices.get(key)
        if index is None:
            index = len(self.nodes)
            self._indices[key] = index
            self.nodes.append((op, operands, value))
        return _Expression(self, index)

    def input(self, axis: int) -> _Expression:
        return self._record("input", (axis,))

    def binary(self, op: str, left: _Expression, right: _Expression) -> _Expression:
        return self._record(op, (left.index, right.index))

    def unary(self, op: str, arg: _Expression) -> _Expression:
        return self._record(op, (arg.index,))

    def select(self, condition: _Expression, yes: object, no: object) -> _Expression:
        return self._record(
            "select",
            (condition.index, self.convert(yes).index, self.convert(no).index),
        )

    def convert(self, value: object) -> _Expression:
        if isinstance(value, _Expression):
            if value.graph is not self:
                msg = "cannot mix expression traces"
                raise TraceError(msg)
            return value
        if not isinstance(value, (int, float, np.integer, np.floating)):
            msg = "Rust fields must return a scalar expression"
            raise TraceError(msg)
        scalar = float(value)
        if not math.isfinite(scalar):
            msg = "nonfinite constants are unsupported"
            raise ValueError(msg)
        return self._record("constant", (), scalar)

    def derivatives(self, output: _Expression) -> list[_Expression]:
        derivative = _Derivatives(self)
        return [derivative.at(output.index, axis) for axis in range(len(_AXES))]


class _Derivatives:
    def __init__(self, graph: _Graph) -> None:
        self.graph = graph
        self.zero = graph.convert(0.0)
        self.one = graph.convert(1.0)
        self.cache: dict[tuple[int, int], _Expression] = {}

    def is_constant(self, expression: _Expression, literal: float) -> bool:
        # A small nonzero coefficient still contributes to the derivative.
        op, _, value = self.graph.nodes[expression.index]
        return op == "constant" and value == literal

    def mul(self, a: _Expression, b: _Expression) -> _Expression:
        if self.is_constant(a, 0.0) or self.is_constant(b, 0.0):
            return self.zero
        if self.is_constant(a, 1.0):
            return b
        if self.is_constant(b, 1.0):
            return a
        return a * b

    def add(self, a: _Expression, b: _Expression) -> _Expression:
        if self.is_constant(a, 0.0):
            return b
        if self.is_constant(b, 0.0):
            return a
        return a + b

    def at(self, index: int, axis: int) -> _Expression:
        key = (index, axis)
        if key not in self.cache:
            self.cache[key] = self._differentiate(index, axis)
        return self.cache[key]

    def _differentiate(self, index: int, axis: int) -> _Expression:
        op, ids, _ = self.graph.nodes[index]
        if op == "input":
            return self.one if ids[0] == axis else self.zero
        if op in {"constant", "lt", "le", "gt", "ge", "eq", "ne"}:
            return self.zero
        a = _Expression(self.graph, ids[0])
        if op == "select":
            return self.graph.select(a, self.at(ids[1], axis), self.at(ids[2], axis))
        da = self.at(ids[0], axis)
        if len(ids) == 1:
            return self._unary(op, a, da)
        b = _Expression(self.graph, ids[1])
        return self._binary(op, a, b, da, self.at(ids[1], axis))

    def _binary(
        self,
        op: str,
        a: _Expression,
        b: _Expression,
        da: _Expression,
        db: _Expression,
    ) -> _Expression:
        if op == "atan2":
            denominator = a * a + b * b
            return self.graph.select(
                denominator > self.zero,
                (self.mul(b, da) - self.mul(a, db)) / denominator,
                self.zero,
            )
        if op in {"minimum", "maximum"}:
            half = self.graph.convert(0.5)
            tied = self.add(self.mul(half, da), self.mul(half, db))
            return self.graph.select(
                a == b,
                tied,
                self.graph.select(a <= b if op == "minimum" else a >= b, da, db),
            )
        return self._arithmetic(op, a, b, da, db)

    def _arithmetic(
        self,
        op: str,
        a: _Expression,
        b: _Expression,
        da: _Expression,
        db: _Expression,
    ) -> _Expression:
        if op == "add":
            return self.add(da, db)
        if op == "sub":
            return da if self.is_constant(db, 0.0) else self.add(da, -db)
        if op == "mul":
            return self.add(self.mul(da, b), self.mul(a, db))
        if op == "pow":
            return self._power(a, b, da, db)
        if op != "div":
            msg = f"no derivative rule for {op}"
            raise TraceError(msg)
        if self.is_constant(db, 0.0):
            return self.zero if self.is_constant(da, 0.0) else da / b
        return (self.mul(da, b) - self.mul(a, db)) / (b * b)

    def _power(
        self,
        a: _Expression,
        b: _Expression,
        da: _Expression,
        db: _Expression,
    ) -> _Expression:
        if self.is_constant(b, 0.0):
            return self.zero
        if self.is_constant(b, 1.0):
            return da
        if self.is_constant(db, 0.0):
            return self.mul(self.mul(b, a ** (b - self.one)), da)
        return self.mul(
            a**b,
            self.add(self.mul(db, self.graph.unary("ln", a)), self.mul(b, da / a)),
        )

    def _unary(self, op: str, a: _Expression, da: _Expression) -> _Expression:
        if self.is_constant(da, 0.0):
            return self.zero
        if op == "neg":
            return -da
        if op == "abs":
            return self.graph.select(
                a > self.zero,
                da,
                self.graph.select(a < self.zero, -da, self.zero),
            )
        if op == "sqrt":
            # Keep the zero subgradient at a norm's origin and the singular slope of sqrt(x).
            quotient = da / (self.graph.convert(2.0) * self.graph.unary("sqrt", a))
            return self.graph.select(
                a == self.zero,
                self.graph.select(da == self.zero, self.zero, quotient),
                quotient,
            )
        slope = _UNARY_SLOPES.get(op)
        if slope is None:
            msg = f"no derivative rule for {op}"
            raise TraceError(msg)
        return self.mul(slope(self.graph, self.one, a), da)


_UNARY_SLOPES: dict[str, Callable[[_Graph, _Expression, _Expression], _Expression]] = {
    "asin": lambda graph, one, a: one / graph.unary("sqrt", one - a * a),
    "acos": lambda graph, one, a: -one / graph.unary("sqrt", one - a * a),
    "atan": lambda _graph, one, a: one / (one + a * a),
    "sin": lambda graph, _one, a: graph.unary("cos", a),
    "cos": lambda graph, _one, a: -graph.unary("sin", a),
    "tan": lambda graph, one, a: one / (graph.unary("cos", a) ** graph.convert(2.0)),
    "exp": lambda graph, _one, a: graph.unary("exp", a),
    "ln": lambda _graph, one, a: one / a,
    "sinh": lambda graph, _one, a: graph.unary("cosh", a),
    "cosh": lambda graph, _one, a: graph.unary("sinh", a),
    "tanh": lambda graph, one, a: one - graph.unary("tanh", a) ** graph.convert(2.0),
}


@lru_cache(maxsize=32)
def _compile_cached(
    nodes: tuple[_Node, ...],
    outputs: tuple[int, ...],
    constraints: tuple[tuple[int, ...], ...] = (),
    *,
    retain_sdf: bool = False,
) -> _meshers.CompiledField:
    if constraints:
        return compile_graph(nodes, outputs, constraints)
    if retain_sdf:
        return compile_graph(nodes, outputs, retain_sdf=True)
    return compile_graph(nodes, outputs)


def _compile_with_intersection(
    graph: _Graph,
    value: _Expression,
    outputs: tuple[int, ...],
) -> _meshers.CompiledField:
    original = tuple(graph.nodes)
    if not getattr(_meshers, "RETAINED_INTERSECTIONS", False):
        return _compile_cached(original, outputs)
    if graph.nodes[value.index][0] != "maximum":
        return _compile_cached(
            original,
            outputs,
            retain_sdf=getattr(_meshers, "SDF_INTERSECTIONS", False),
        )
    roots, pending, seen = [], [value.index], set()
    while pending:
        index = pending.pop()
        if index in seen:
            continue
        seen.add(index)
        op, args, _ = graph.nodes[index]
        if op == "maximum":
            pending.extend(reversed(args))
        else:
            roots.append(index)
    if not _MIN_CONSTRAINTS <= len(roots) <= _meshers.MAX_CONSTRAINTS or any(
        graph.nodes[index][0] == "constant" for index in roots
    ):
        return _compile_cached(original, outputs)
    constraints = tuple(
        (index, *(e.index for e in graph.derivatives(_Expression(graph, index)))) for index in roots
    )
    if len(graph.nodes) > _MAX_NODES:
        return _compile_cached(original, outputs)
    return _compile_cached(tuple(graph.nodes), outputs, constraints)


def compile_field(
    function: Callable[..., object],
    *,
    gradient: Callable[..., object] | Literal[False] | None = None,
    backend: str = "rust",
) -> _meshers.CompiledField:
    """Prepare an immutable native evaluator. Usually called automatically.

    Args:
        function: Pointwise NumPy-style function of x, y, z.
        gradient: None differentiates the expression graph automatically. A callable
            returns three scalar expressions; False retains finite differences.
        backend: ``rust`` or its alias ``jit``. Both use the embedded compiler.

    Returns:
        CompiledField: A cached float64 evaluator owning its executable memory.
            Captured constants are frozen; preparing again observes changed parameters.

    Raises:
        TraceError: The callable uses unsupported operations or exceeds the graph limit.
        ValueError: The backend name or a constant is invalid.

    Note:
        No Rust installation, Numba, SymPy or subprocess is used at runtime.
        Gradients at abs(0) are zero. Exact min/max ties average the operand
        derivatives. This is a subgradient choice at a true crease, not a
        branch constraint. Eligible meshing modes can additionally retain
        outer-maximum branches. Adaptive Tet4 meshing can also retain the
        zero-sublevel intersection of a standard clamped-norm signed-distance
        graph with unit-distance affine/radial branches, after numerical checks.
        Supply a gradient callable for another policy.
    """
    if backend not in {"rust", "jit"}:
        msg = "only the embedded Rust compiler is supported"
        raise ValueError(msg)
    graph = _Graph()
    xyz = [graph.input(axis) for axis in range(len(_AXES))]
    value = graph.convert(function(*xyz))
    outputs = [value]
    if gradient is None:
        outputs += graph.derivatives(value)
    elif gradient is not False:
        g = gradient(*xyz)
        if not isinstance(g, (tuple, list)) or len(g) != len(_AXES):
            msg = "gradient must return three scalar expressions"
            raise TraceError(msg)
        outputs.extend(graph.convert(v) for v in g)
    if len(graph.nodes) > _MAX_NODES:
        msg = "expression is too large to compile"
        raise TraceError(msg)
    indices = tuple(e.index for e in outputs)
    return (
        _compile_with_intersection(graph, value, indices)
        if gradient is None
        else _compile_cached(tuple(graph.nodes), indices)
    )
