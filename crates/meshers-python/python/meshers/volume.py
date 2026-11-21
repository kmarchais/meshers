# Copyright (c) 2026 Kevin Marchais

"""Generate tetrahedral meshes from scalar fields."""

import time
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import TypedDict, Unpack

import numpy as np
import numpy.typing as npt

from . import _meshers
from ._meshers import CompiledField, MeshingError
from ._options import MeshingOptions
from .expression import LegacyCompileOption, TraceError, compile_field, resolve_compile_option
from .intersection import (
    MAX_INTERSECTION_PASSES,
    _affine_planes,
    _check_periodic_axes,
    _periodic_affine,
    _retained_constraints,
    generate_intersection,
)
from .mesh import Mesh

_QUADRATIC_ORDER = 2

Field = Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
    npt.NDArray[np.float64],
]


def _check_components(mesh: Mesh, expected_components: int | None) -> Mesh:
    if expected_components is not None:
        count = _meshers.face_components(mesh.tetrahedra)
        if count != expected_components:
            msg = f"Expected {expected_components} face-connected components, found {count}"
            raise MeshingError(msg)
        mesh.diagnostics["face_connected_components"] = count
    return mesh


def _prepare_volume_field(
    field: str | Field | CompiledField,
    gradient: Field | None,
    *,
    compile_expression: bool,
) -> tuple[str | Field | CompiledField, Field | None, str]:
    evaluator = (
        "native"
        if isinstance(field, str)
        else "compiled"
        if isinstance(field, CompiledField)
        else "python"
    )
    if callable(field) and not isinstance(field, CompiledField) and compile_expression:
        try:
            field = compile_field(field, gradient=gradient)
            gradient = None
            evaluator = "compiled"
        except TraceError as error:
            warnings.warn(f"Using NumPy callbacks: {error}", RuntimeWarning, stacklevel=3)
    return field, gradient, evaluator


class _NativeConstraints(TypedDict, total=False):
    halfspaces: Sequence[Sequence[float]]
    regions: Sequence[int]
    constraints: list[Callable | CompiledField | str]


def _native_constraints(
    field: str | Field | CompiledField,
    planes: dict[str, Callable | CompiledField | str] | None,
    *,
    adapt: bool,
    scalar_mode: bool,
) -> _NativeConstraints:
    if not isinstance(field, CompiledField):
        return {}
    if not scalar_mode and getattr(_meshers, "POLYHEDRAL_REGIONS", False):
        regions = field.polyhedral_regions
        if regions is not None:
            return {"halfspaces": regions[0], "regions": regions[1]}
    if adapt and planes is not None:
        halfspaces = field.halfspaces
        if halfspaces is None:
            return {"constraints": list(planes.values())}
        if getattr(_meshers, "ADAPTIVE_POLYHEDRAL", False):
            return {"halfspaces": halfspaces}
    return {}


class _VolumeOptions(MeshingOptions, LegacyCompileOption, total=False):
    """Tetrahedral meshing, adaptation and admission keywords."""

    gradient: Field | None
    periodic: Sequence[bool]
    band: tuple[float, float] | None
    minimum_quality: float
    max_tetrahedra: int
    optimize_passes: int
    interior_refine_nodes: int
    snap: float
    threads: int | None
    adapt: bool
    layers: float
    adapt_iterations: int
    quality_guard: float
    adapt_tolerance: float | None
    order: int
    size_field: Field | None
    curved_fit: bool
    expected_components: int | None


@dataclass(frozen=True, slots=True)
class _VolumeRoute:
    periodic: Sequence[bool]
    threads: int
    passes: int
    budget: int
    interior_nodes: int
    order: int
    adapt: bool

    @classmethod
    def from_options(cls, options: _VolumeOptions) -> "_VolumeRoute":
        periodic = options.get("periodic", (False, False, False))
        _check_periodic_axes(periodic)
        threads = options.get("threads", 1)
        return cls(
            periodic=periodic,
            threads=_meshers.available_threads() if threads is None else threads,
            passes=options.get("optimize_passes", 4),
            budget=options.get("max_tetrahedra", 2_000_000),
            interior_nodes=options.get("interior_refine_nodes", 0),
            order=options.get("order", 1),
            adapt=options.get("adapt", False),
        )

    def constraints(
        self,
        field: str | Field | CompiledField,
        bounds: npt.NDArray[np.float64],
        options: _VolumeOptions,
    ) -> tuple[dict[str, Callable | CompiledField | str] | None, bool]:
        scalar_mode = (
            self.order != 1,
            options.get("band") is not None,
            any(self.periodic)
            and not _periodic_affine(field, self.periodic, bounds if self.adapt else None),
            self.threads != 1,
            self.interior_nodes != 0,
            options.get("gradient") is not None,
            # Nonadaptive SDF retention would bypass the original parent's range checks.
            not self.adapt and isinstance(field, CompiledField) and field.sdf_intersection,
            not 0 <= self.passes <= MAX_INTERSECTION_PASSES,
            self.budget <= 0,
            not all(
                isinstance(value, (int, np.integer))
                for value in (self.threads, self.interior_nodes, self.passes, self.budget)
            ),
        )
        planes = _affine_planes(field, scalar_mode)
        scalar = any(scalar_mode)
        if (
            planes is None
            and not scalar
            and (not self.adapt or getattr(_meshers, "ADAPTIVE_RETAINED_INTERSECTIONS", False))
        ):
            planes = _retained_constraints(field, bounds)
        return planes, scalar


def _adaptation(
    options: _VolumeOptions,
    bounds: npt.NDArray[np.float64],
    tolerance: float,
) -> tuple[float, float, float, float, int, float, float] | None:
    if not options.get("adapt"):
        return None
    minimum, maximum = options.get("min_length"), options.get("max_length")
    target = options.get("adapt_tolerance")
    return (
        tolerance if minimum is None else float(minimum),
        float(np.min(bounds[1::2] - bounds[::2])) / 4 if maximum is None else float(maximum),
        float(options.get("gradation", 0.3)),
        float(options.get("layers", 0.0)),
        int(options.get("adapt_iterations", 4)),
        float(options.get("quality_guard", 0.3)),
        tolerance if target is None else float(target),
    )


def generate(
    field: str | Field | CompiledField,
    *,
    bounds: Sequence[float] = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
    cells: int | Sequence[int] = 24,
    geometry_tolerance: float = 0.01,
    compile_expression: bool | None = None,
    **options: Unpack[_VolumeOptions],
) -> Mesh:
    """Generate a tetrahedral mesh directly from a bounded implicit field.

    Geometry and topology checks raise `MeshingError` on failure. Cancellation raises
    `CancelledError`.

    Args:
        field: Pointwise function of physical x, y, z coordinates, or a built-in
            name: ``gyroid``, ``schwarz_p``, ``schwarz_d``. Built-ins have period 1.
        bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)`` in physical units.
        cells: Background lattice intervals, 4..128 per axis. An int applies to all axes.
        geometry_tolerance: Maximum sampled boundary-distance estimate, in physical units.
        compile_expression: Compile supported expressions and their first derivatives
            automatically. False selects NumPy callbacks. Unsupported operations
            emit a warning and use callbacks, which must return float64 (N,) arrays.
            Eligible Tet4 convex polyhedra retain exact half-space constraints;
            requested periodic axes must leave every plane invariant. Diagnostics
            report ``polyhedral_constraints``. Bounded min/max combinations of
            convex affine fields retain nonconvex regions on the nonperiodic path.
            Eligible nonperiodic outer-max intersections retain analytic branches;
            diagnostics report ``retained_constraints``.
        **options: Additional meshing keywords listed below.

    Other Parameters:
        gradient (Field | None): Default ``None``. Optional gradient of the complete field,
            including grading. Traced gradients return three expressions; NumPy callbacks return
            (N,3).
        periodic (Sequence[bool]): Default ``(False, False, False)``. Explicit matching constraints
            for x, y, z. The field must match across requested faces; periodicity is not inferred.
        band (tuple[float, float] | None): Default ``None``. Keep ``lower < field < upper``.
            Without a band, keep the negative side.
        minimum_quality (float): Default ``0.0``. Required minimum MMG tetrahedron quality, 0..1.
            Zero disables the requirement. Adaptive splitting and wall projection honor a positive
            requested minimum. Nonempty meshes below it are rejected, and requested geometry limits
            are not relaxed. Empty meshes remain valid.
        max_tetrahedra (int): Default ``2000000``. Budget for generated tetrahedra. Adaptation
            raises MeshingError when a local topology operation exceeds this budget.
        optimize_passes (int): Default ``4``. Number of quality-optimization passes, 0..100.
        interior_refine_nodes (int): Default ``0``. Maximum new nodes in one length-ranked interior
            bisection pass after scalar surface generation or adaptation. Zero disables it. For
            order=2, this counts new corner nodes; quadratic edge nodes are additional. Boundary
            vertices and triangles stay fixed. Disjoint cavities, the tetrahedron budget and a
            quality floor of max(minimum_quality, 0.1) may limit the number inserted.
        snap (float): Default ``0.2``. Maximum lattice snap as a fraction of minimum spacing,
            0..0.2.
        threads (int | None): Native worker limit, default one. None uses available CPU
            parallelism, capped at 256. One uses one worker; zero retains legacy serial update
            ordering.
        batch_size (int): Default ``4096``. Maximum NumPy callback batch size.
        cancel (CancellationToken | None): Default ``None``. Optional cooperative cancellation
            token.
        adapt (bool): Default ``False``. Experimental. Replace fixed surface refinement by
            size-driven adaptation: boundary edges follow the curvature needed for
            ``geometry_tolerance``, the interior grades towards ``max_length``. ``cells`` then only
            sets the starting lattice, which must resolve topology. Compiled convex polyhedra
            preserve their exact planes and retain existing crease vertices. Face/interior
            coarsening, splitting, spatial sizing and constrained quality improvement remain
            active.
        min_length (float | None): Smallest adapted edge length; default ``geometry_tolerance``.
        max_length (float | None): Largest adapted edge length; default a quarter of the smallest
            box side.
        gradation (float): Default ``0.3``. Maximum length growth per unit distance (MMG ``hgrad -
            1``).
        layers (float): Default ``0.0``. For bands, minimum element count across the local sheet
            thickness.
        adapt_iterations (int): Default ``4``. Split/collapse/swap/relocate rounds.
        quality_guard (float): Default ``0.3``. Operations may not create MMG quality below this
            value unless the affected elements were already worse.
        adapt_tolerance (float | None): Distance the adapter aims for; default
            ``geometry_tolerance``, which remains the acceptance check.
        order (int): Default ``1``. 1 for linear tetrahedra, 2 for 10-node tetrahedra
            (VTK/Gmsh/fedoo node order) with wall mid-edge nodes projected onto the field where
            sampled validity permits. Invalid corrections are reduced or straightened;
            ``diagnostics["quadratic"]`` records both counts. Quadratic meshes report
            ``tetrahedra`` of shape (M, 10) and ``surface`` of shape (K, 6); ``periodic_pairs``
            then include cap midpoint nodes and preserve complete constraint classes for existing
            periodic node images.
        size_field (Field | None): Default ``None``. Optional vectorized callback returning
            positive physical target edge lengths, shape (N,), for ``adapt=True``. Targets combine
            with curvature sizing by taking the smaller length, then applying ``min_length``,
            ``max_length`` and gradation. Opposite periodic traces must match. This callback is
            evaluated in batches without tracing.
        curved_fit (bool): Default ``False``. Experimental bounded fitting of the final quadratic
            wall faces. Requires order=2. Corners and connectivity stay fixed; wall edge nodes may
            move off the exact implicit surface to improve the complete patch. Cap coordinates and
            periodic aliases stay coupled. Sampled validity guards each accepted move. This changes
            the discrete geometry and FE space; it does not transfer an old solution or certify
            convergence.
        expected_components (int | None): Default ``None``. Require this many face-connected solid
            components in the cut box, before pairing periodic faces. Zero requires sampled-empty
            output.

    Returns:
        Mesh: Owned arrays and diagnostics. Empty geometries return empty arrays.
            Timing excludes field preparation, which is reported separately.

    Raises:
        ValueError: Options or callback outputs are invalid.

    Note:
        Fields must be deterministic and independent of batch contents. Tracing
        executes a function once per preparation, so avoid side effects. Finite
        differences need the field slightly beyond the box. Sampling cannot
        certify arbitrarily small features. Native threads release the Python GIL;
        Python fallback callbacks execute on the calling thread.
    """
    compile_expression = resolve_compile_option(
        options,
        schema=_VolumeOptions,
        setting=compile_expression,
    )
    expected_components = options.get("expected_components")
    if expected_components is not None and (
        isinstance(expected_components, (bool, np.bool_))
        or not isinstance(expected_components, (int, np.integer))
        or expected_components < 0
    ):
        msg = "expected_components must be a nonnegative integer"
        raise ValueError(msg)

    b = np.asarray(bounds, dtype=np.float64)
    if options.get("curved_fit") and options.get("order", 1) != _QUADRATIC_ORDER:
        msg = "curved_fit requires order 2"
        raise ValueError(msg)
    if options.get("size_field") is not None and not options.get("adapt"):
        msg = "size_field requires adapt=True"
        raise ValueError(msg)
    if b.shape != (6,):
        msg = "bounds must contain xmin,xmax,ymin,ymax,zmin,zmax"
        raise ValueError(msg)
    if isinstance(cells, (int, np.integer)):
        cells = (int(cells),) * 3
    route = _VolumeRoute.from_options(options)
    preparation_start = time.perf_counter()
    field, gradient, evaluator = _prepare_volume_field(
        field,
        options.get("gradient"),
        compile_expression=compile_expression,
    )
    planes, scalar_mode = route.constraints(field, b, options)
    preparation_seconds = time.perf_counter() - preparation_start
    if planes is not None and not route.adapt:
        mesh = generate_intersection(
            planes,
            periodic=route.periodic,
            bounds=b.tolist(),
            cells=cells,
            geometry_tolerance=geometry_tolerance,
            minimum_quality=options.get("minimum_quality", 0.0),
            max_tetrahedra=int(route.budget),
            optimize_passes=int(route.passes),
            snap=options.get("snap", 0.2),
            batch_size=options.get("batch_size", 4096),
            cancel=options.get("cancel"),
        )
        mesh.diagnostics.update(
            threads=1,
            evaluator="compiled",
            inserted_interior_nodes=0,
            preparation_seconds=preparation_seconds
            + mesh.diagnostics.get("preparation_seconds", 0.0),
        )
        mesh.diagnostics[
            "polyhedral_constraints"
            if isinstance(field, CompiledField) and field.halfspaces is not None
            else "retained_constraints"
        ] = len(planes)
        return _check_components(mesh, expected_components)
    result = _meshers.generate(
        field,
        gradient,
        [b[::2].tolist(), b[1::2].tolist()],
        cells,
        tuple(bool(v) for v in route.periodic),
        options.get("band"),
        geometry_tolerance,
        options.get("minimum_quality", 0.0),
        route.budget,
        route.passes,
        route.interior_nodes,
        options.get("snap", 0.2),
        route.threads,
        options.get("batch_size", 4096),
        options.get("cancel"),
        _adaptation(options, b, geometry_tolerance),
        int(route.order),
        options.get("size_field"),
        curved_fit=options.get("curved_fit", False),
        **_native_constraints(field, planes, adapt=route.adapt, scalar_mode=scalar_mode),
    )
    result["diagnostics"]["threads"] = route.threads
    result["diagnostics"]["evaluator"] = evaluator
    result["diagnostics"]["preparation_seconds"] = preparation_seconds
    _warn_defeaturing(result["diagnostics"], options.get("minimum_quality", 0.0))
    if result["diagnostics"].get("adaptation_abandoned"):
        warnings.warn(
            "The final checks refused the adapted mesh; returned the uniform lattice mesh "
            "of the same field instead. A finer starting lattice may let adaptation pass.",
            UserWarning,
            stacklevel=2,
        )
    return _check_components(_volume_mesh(result), expected_components)


def _warn_defeaturing(diagnostics: dict, minimum_quality: float) -> None:
    """Warn when sliver recovery moved wall vertices off the implicit surface."""
    moved = diagnostics.get("defeaturing_moves", 0)
    if moved:
        warnings.warn(
            f"Moved {moved} surface vertices off the implicit surface to reach "
            f"minimum_quality={minimum_quality}; the geometry changed there, within "
            "geometry_tolerance. A finer lattice may avoid it.",
            UserWarning,
            stacklevel=3,
        )


def _volume_mesh(result: "_meshers._GeneratedVolume") -> Mesh:
    """Return a mesh with the connectivity and periodic pairs of the requested order."""
    quadratic = result.get("quadratic")
    geometry = quadratic if quadratic is not None else result
    diagnostics = result["diagnostics"]
    if quadratic is not None:
        diagnostics["quadratic"] = {
            key: value
            for key, value in quadratic.items()
            if key not in {"points", "tetrahedra", "surface", "periodic_pairs"}
        }
    return Mesh(
        points=geometry["points"],
        tetrahedra=geometry["tetrahedra"],
        surface=geometry["surface"],
        periodic_pairs=tuple(geometry["periodic_pairs"]),
        boundary_tags=result["boundary_tags"],
        diagnostics=diagnostics,
    )
