"""Scalar expression tracing and first derivatives for the embedded native compiler."""

import math
from collections.abc import Callable
from functools import lru_cache
from typing import Literal

import numpy as np

from . import _meshers


class TraceError(TypeError):
    """A callable uses operations outside the supported expression subset."""


class _Expression:
    __array_priority__ = 1000

    def __init__(self, graph, index):
        self.graph, self.index = graph, index

    def _binary(self, op, other):
        return self.graph.node(op, self, self.graph.convert(other))

    def __add__(self, o):
        return self._binary("add", o)

    def __radd__(self, o):
        return self.graph.convert(o)._binary("add", self)

    def __sub__(self, o):
        return self._binary("sub", o)

    def __rsub__(self, o):
        return self.graph.convert(o)._binary("sub", self)

    def __mul__(self, o):
        return self._binary("mul", o)

    def __rmul__(self, o):
        return self.graph.convert(o)._binary("mul", self)

    def __truediv__(self, o):
        return self._binary("div", o)

    def __rtruediv__(self, o):
        return self.graph.convert(o)._binary("div", self)

    def __pow__(self, o):
        return self._binary("pow", o)

    def __rpow__(self, o):
        return self.graph.convert(o)._binary("pow", self)

    def __neg__(self):
        return self.graph.node("neg", self)

    def __pos__(self):
        return self

    def __abs__(self):
        return self.graph.node("abs", self)

    def __bool__(self):
        raise TraceError(
            "data-dependent Python branches cannot be traced; use the NumPy callback path"
        )

    def __lt__(self, o):
        return self._binary("lt", o)

    def __le__(self, o):
        return self._binary("le", o)

    def __gt__(self, o):
        return self._binary("gt", o)

    def __ge__(self, o):
        return self._binary("ge", o)

    def __eq__(self, o):
        return self._binary("eq", o)

    def __ne__(self, o):
        return self._binary("ne", o)

    def __float__(self):
        raise TraceError("scalar conversion cannot be traced")

    def __int__(self):
        raise TraceError("integer conversion cannot be traced")

    def __len__(self):
        raise TraceError("array lengths require the NumPy callback path")

    def __getitem__(self, key):
        raise TraceError("array indexing requires the NumPy callback path")

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        raise TraceError(f"array attribute {name} requires the NumPy callback path")

    def __array_function__(self, func, types, args, kwargs):
        if func is np.where and len(args) == 3 and not kwargs:
            return self.graph.node("select", *(self.graph.convert(a) for a in args))
        if func is np.clip and len(args) == 3 and not kwargs:
            value, lower, upper = args
            result = self.graph.convert(value)
            if lower is not None:
                result = self.graph.node("maximum", result, self.graph.convert(lower))
            if upper is not None:
                result = self.graph.node("minimum", result, self.graph.convert(upper))
            return result
        raise TraceError(f"unsupported traced NumPy function: {func.__name__}")

    def __array__(self, *args, **kwargs):
        raise TraceError("array construction is unsupported by the Rust scalar tracer")

    def __array_ufunc__(self, ufunc, method, *args, **kwargs):
        ops = {
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
        if method != "__call__" or kwargs or ufunc not in ops:
            raise TraceError(f"unsupported traced operation: {ufunc.__name__}.{method}")
        return self.graph.node(ops[ufunc], *(self.graph.convert(a) for a in args))


class _Graph:
    def __init__(self):
        self.nodes = []
        self.intern = {}

    def node(self, op, *args):
        key = (op, *(a.index for a in args))
        return self._intern(key)

    def _intern(self, key):
        if key not in self.intern:
            self.intern[key] = len(self.nodes)
            self.nodes.append(key)
        return _Expression(self, self.intern[key])

    def convert(self, value):
        if isinstance(value, _Expression):
            if value.graph is not self:
                raise TraceError("cannot mix expression traces")
            return value
        if not isinstance(value, (int, float, np.integer, np.floating)):
            raise TraceError("Rust fields must return a scalar expression")
        v = float(value)
        if not math.isfinite(v):
            raise ValueError("nonfinite constants are unsupported")
        # Hex identity preserves signed zero and avoids float-key equality surprises.
        return self._intern(("constant", v.hex()))

    def derivatives(self, output):
        zero = self.convert(0.0)
        one = self.convert(1.0)

        def constant(e):
            n = self.nodes[e.index]
            return float.fromhex(n[1]) if n[0] == "constant" else None

        def mul(a, b):
            if constant(a) == 0.0 or constant(b) == 0.0:
                return zero
            if constant(a) == 1.0:
                return b
            if constant(b) == 1.0:
                return a
            return a * b

        def add(a, b):
            if constant(a) == 0.0:
                return b
            if constant(b) == 0.0:
                return a
            return a + b

        def choose(c, a, b):
            return self.node("select", c, a, b)

        @lru_cache(None)
        def d(index, axis):
            op, *ids = self.nodes[index]
            if op == "input":
                return one if ids[0] == axis else zero
            if op == "constant" or op in {"lt", "le", "gt", "ge", "eq", "ne"}:
                return zero
            args = [_Expression(self, i) for i in ids]
            if op == "select":
                return choose(args[0], d(ids[1], axis), d(ids[2], axis))
            a = args[0]
            da = d(ids[0], axis)
            if op == "atan2":
                b = args[1]
                db = d(ids[1], axis)
                denominator = a * a + b * b
                return choose(denominator > zero, (mul(b, da) - mul(a, db)) / denominator, zero)
            if op in {"minimum", "maximum"}:
                b = args[1]
                db = d(ids[1], axis)
                return choose(a <= b if op == "minimum" else a >= b, da, db)
            if op in {"add", "sub", "mul", "div", "pow"}:
                b = args[1]
                db = d(ids[1], axis)
                if op == "add":
                    return add(da, db)
                if op == "sub":
                    return add(da, -db) if constant(db) != 0.0 else da
                if op == "mul":
                    return add(mul(da, b), mul(a, db))
                if op == "div":
                    if constant(da) == 0.0 and constant(db) == 0.0:
                        return zero
                    if constant(db) == 0.0:
                        return da / b
                    return (mul(da, b) - mul(a, db)) / (b * b)
                exponent = constant(b)
                if exponent == 0.0:
                    return zero
                if exponent == 1.0:
                    return da
                if constant(db) == 0.0:
                    return mul(mul(b, a ** (b - one)), da)
                # Real-valued variable powers require a positive base.
                return mul(a**b, add(mul(db, self.node("ln", a)), mul(b, da / a)))
            if constant(da) == 0.0:
                return zero
            if op == "neg":
                return -da
            if op == "abs":
                return choose(a > zero, da, choose(a < zero, -da, zero))
            if op == "sqrt":
                # A norm at its origin has no unique gradient. Choose its zero
                # subgradient when both the radicand and its derivative vanish.
                # Preserve infinities for sqrt(x) at zero and invalid domains.
                quotient = da / (self.convert(2.0) * self.node("sqrt", a))
                return choose(a == zero, choose(da == zero, zero, quotient), quotient)
            local = {
                "asin": lambda: one / self.node("sqrt", one - a * a),
                "acos": lambda: -one / self.node("sqrt", one - a * a),
                "atan": lambda: one / (one + a * a),
                "sin": lambda: self.node("cos", a),
                "cos": lambda: -self.node("sin", a),
                "tan": lambda: one / (self.node("cos", a) ** self.convert(2.0)),
                "exp": lambda: self.node("exp", a),
                "ln": lambda: one / a,
                "sqrt": lambda: one / (self.convert(2.0) * self.node("sqrt", a)),
                "sinh": lambda: self.node("cosh", a),
                "cosh": lambda: self.node("sinh", a),
                "tanh": lambda: one - self.node("tanh", a) ** self.convert(2.0),
            }
            if op not in local:
                raise TraceError(f"no derivative rule for {op}")
            return mul(local[op](), da)

        return [d(output.index, axis) for axis in ("x", "y", "z")]

    def encoded(self):
        encoded = []
        for op, *args in self.nodes:
            if op == "input":
                encoded.append((op, ("xyz".index(args[0]),), 0.0))
            elif op == "constant":
                encoded.append((op, (), float.fromhex(args[0])))
            else:
                encoded.append((op, tuple(args), 0.0))
        return tuple(encoded)


@lru_cache(maxsize=32)
def _compile_cached(nodes, outputs):
    return _meshers._compile_expression(nodes, outputs)


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
        Gradients at abs(0) are zero; min/max ties use the first branch.
    """
    if backend not in {"rust", "jit"}:
        raise ValueError("only the embedded Rust compiler is supported")
    graph = _Graph()
    xyz = [graph._intern(("input", a)) for a in ("x", "y", "z")]
    value = graph.convert(function(*xyz))
    outputs = [value]
    if gradient is None:
        outputs += graph.derivatives(value)
    elif gradient is not False:
        g = gradient(*xyz)
        if not isinstance(g, (tuple, list)) or len(g) != 3:
            raise TraceError("gradient must return three scalar expressions")
        outputs.extend(graph.convert(v) for v in g)
    if len(graph.nodes) > 20000:
        raise TraceError("expression is too large to compile")
    return _compile_cached(graph.encoded(), tuple(e.index for e in outputs))
