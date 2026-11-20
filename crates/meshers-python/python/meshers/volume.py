# Copyright (c) 2026 Kevin Marchais

"""Generate tetrahedral meshes from scalar fields."""

import time
import warnings
from collections.abc import Callable, Sequence

import numpy as np
import numpy.typing as npt

from . import _meshers
from ._meshers import CancellationToken, CompiledField, MeshingError
from .expression import TraceError, compile_field
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


# Existing keyword arguments are part of the published Python API.
def generate(  # ruff: ignore[too-many-arguments]
    field: str | Field | CompiledField,
    *,
    gradient: Field | None = None,
    bounds: Sequence[float] = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
    cells: int | Sequence[int] = 24,
    periodic: Sequence[bool] = (False, False, False),
    band: tuple[float, float] | None = None,
    geometry_tolerance: float = 0.01,
    minimum_quality: float = 0.0,
    max_tetrahedra: int = 2_000_000,
    optimize_passes: int = 4,
    interior_refine_nodes: int = 0,
    snap: float = 0.2,
    threads: int | None = 1,
    batch_size: int = 4096,
    cancel: CancellationToken | None = None,
    compile: bool = True,  # ruff: ignore[builtin-argument-shadowing] - public keyword.
    adapt: bool = False,
    min_length: float | None = None,
    max_length: float | None = None,
    gradation: float = 0.3,
    layers: float = 0.0,
    adapt_iterations: int = 4,
    quality_guard: float = 0.3,
    adapt_tolerance: float | None = None,
    order: int = 1,
    size_field: Field | None = None,
    curved_fit: bool = False,
    expected_components: int | None = None,
) -> Mesh:
    """Generate a tetrahedral mesh directly from a bounded implicit field.

    Args:
        expected_components: Require this many face-connected solid components in the cut box,
            before pairing periodic faces. Zero requires sampled-empty output.
        field: Pointwise function of physical x, y, z coordinates, or a built-in
            name: ``gyroid``, ``schwarz_p``, ``schwarz_d``. Built-ins have period 1.
        gradient: Optional gradient of the complete field, including grading.
            Traced gradients return three expressions; NumPy callbacks return (N,3).
        bounds: ``(xmin, xmax, ymin, ymax, zmin, zmax)`` in physical units.
        cells: Background lattice intervals, 4..128 per axis. An int applies to all axes.
        periodic: Explicit matching constraints for x, y, z. The field must match
            across requested faces; periodicity is not inferred.
        band: Keep ``lower < field < upper``. Without a band, keep the negative side.
        geometry_tolerance: Maximum sampled boundary-distance estimate, in physical units.
        minimum_quality: Required minimum MMG tetrahedron quality, 0..1. Zero disables
            the requirement. Adaptive splitting and wall projection honor a positive
            requested minimum. Nonempty meshes below it are rejected, and requested
            geometry limits are not relaxed. Empty meshes remain valid.
        max_tetrahedra: Budget for generated tetrahedra. Adaptation raises
            MeshingError when a local topology operation exceeds this budget.
        optimize_passes: Number of quality-optimization passes, 0..100.
        interior_refine_nodes: Maximum new nodes in one length-ranked interior
            bisection pass after scalar surface generation or adaptation. Zero
            disables it. For order=2, this counts new corner nodes; quadratic
            edge nodes are additional. Boundary vertices and triangles stay
            fixed. Disjoint cavities, the tetrahedron budget and a quality floor
            of max(minimum_quality, 0.1) may limit the number inserted.
        snap: Maximum lattice snap as a fraction of minimum spacing, 0..0.2.
        threads: Native worker limit, default one. None uses available CPU parallelism, capped
            at 256. One uses one worker; zero retains legacy serial update ordering.
        batch_size: Maximum NumPy callback batch size.
        cancel: Optional cooperative cancellation token.
        compile: Compile supported expressions and their first derivatives
            automatically. False selects NumPy callbacks. Unsupported operations
            emit a warning and use callbacks, which must return float64 (N,) arrays.
            Eligible Tet4 convex polyhedra retain exact half-space constraints;
            requested periodic axes must leave every plane invariant. Diagnostics
            report ``polyhedral_constraints``. Bounded min/max combinations of
            convex affine fields retain nonconvex regions on the nonperiodic path.
            Eligible nonperiodic outer-max intersections retain analytic branches;
            diagnostics report ``retained_constraints``.
        adapt: Experimental. Replace fixed surface refinement by size-driven
            adaptation: boundary edges follow the curvature needed for
            ``geometry_tolerance``, the interior grades towards ``max_length``.
            ``cells`` then only sets the starting lattice, which must resolve topology.
            Compiled convex polyhedra preserve their exact planes and retain
            existing crease vertices. Face/interior coarsening, splitting,
            spatial sizing and constrained quality improvement remain active.
        min_length: Smallest adapted edge length; default ``geometry_tolerance``.
        max_length: Largest adapted edge length; default a quarter of the smallest box side.
        gradation: Maximum length growth per unit distance (MMG ``hgrad - 1``).
        layers: For bands, minimum element count across the local sheet thickness.
        adapt_iterations: Split/collapse/swap/relocate rounds.
        quality_guard: Operations may not create MMG quality below this value
            unless the affected elements were already worse.
        adapt_tolerance: Distance the adapter aims for; default ``geometry_tolerance``,
            which remains the acceptance check.
        order: 1 for linear tetrahedra, 2 for 10-node tetrahedra (VTK/Gmsh/fedoo
            node order) with wall mid-edge nodes projected onto the field where
            sampled validity permits. Invalid corrections are reduced or straightened;
            ``diagnostics["quadratic"]`` records both counts. Quadratic meshes
            report ``tetrahedra`` of shape (M, 10) and ``surface`` of shape (K, 6);
            ``periodic_pairs`` then include cap midpoint nodes and preserve
            complete constraint classes for existing periodic node images.
        size_field: Optional vectorized callback returning positive physical target
            edge lengths, shape (N,), for ``adapt=True``. Targets combine with
            curvature sizing by taking the smaller length, then applying
            ``min_length``, ``max_length`` and gradation. Opposite periodic traces
            must match. This callback is evaluated in batches without tracing.
        curved_fit: Experimental bounded fitting of the final quadratic wall faces.
            Requires order=2. Corners and connectivity stay fixed; wall edge nodes
            may move off the exact implicit surface to improve the complete patch.
            Cap coordinates and periodic aliases stay coupled. Sampled validity
            guards each accepted move. This changes the discrete geometry and FE
            space; it does not transfer an old solution or certify convergence.

    Returns:
        Mesh: Owned arrays and diagnostics. Empty geometries return empty arrays.
            Timing excludes field preparation, which is reported separately.

    Raises:
        ValueError: Options or callback outputs are invalid.
        MeshingError: Meshing fails its geometry or topology checks.
        CancelledError: The cancellation token is triggered.

    Note:
        Fields must be deterministic and independent of batch contents. Tracing
        executes a function once per preparation, so avoid side effects. Finite
        differences need the field slightly beyond the box. Sampling cannot
        certify arbitrarily small features. Native threads release the Python GIL;
        Python fallback callbacks execute on the calling thread.
    """  # ruff: ignore[docstring-extraneous-exception] - includes native errors.
    if expected_components is not None and (
        isinstance(expected_components, (bool, np.bool_))
        or not isinstance(expected_components, (int, np.integer))
        or expected_components < 0
    ):
        msg = "expected_components must be a nonnegative integer"
        raise ValueError(msg)

    b = np.asarray(bounds, dtype=np.float64)
    if curved_fit and order != _QUADRATIC_ORDER:
        msg = "curved_fit requires order 2"
        raise ValueError(msg)
    if size_field is not None and not adapt:
        msg = "size_field requires adapt=True"
        raise ValueError(msg)
    if b.shape != (6,):
        msg = "bounds must contain xmin,xmax,ymin,ymax,zmin,zmax"
        raise ValueError(msg)
    if isinstance(cells, (int, np.integer)):
        cells = (int(cells),) * 3
    _check_periodic_axes(periodic)
    if threads is None:
        threads = _meshers.available_threads()
    preparation_start = time.perf_counter()
    gradient_override = gradient is not None
    field, gradient, evaluator = _prepare_volume_field(field, gradient, compile_expression=compile)
    scalar_mode = (
        order != 1,
        band is not None,
        any(periodic) and not _periodic_affine(field, periodic, b if adapt else None),
        threads != 1,
        interior_refine_nodes != 0,
        gradient_override,
        # Nonadaptive SDF retention would bypass the original parent's range checks.
        not adapt and getattr(field, "_sdf_intersection", False),
        not 0 <= optimize_passes <= MAX_INTERSECTION_PASSES,
        max_tetrahedra <= 0,
        not all(
            isinstance(value, (int, np.integer))
            for value in (threads, interior_refine_nodes, optimize_passes, max_tetrahedra)
        ),
    )
    planes = _affine_planes(field, scalar_mode)
    if (
        planes is None
        and (not adapt or getattr(_meshers, "ADAPTIVE_RETAINED_INTERSECTIONS", False))
        and not any(scalar_mode)
    ):
        planes = _retained_constraints(field, b)
    preparation_seconds = time.perf_counter() - preparation_start
    if planes is not None and not adapt:
        mesh = generate_intersection(
            planes,
            periodic=periodic,
            bounds=b.tolist(),
            cells=cells,
            geometry_tolerance=geometry_tolerance,
            minimum_quality=minimum_quality,
            max_tetrahedra=int(max_tetrahedra),
            optimize_passes=int(optimize_passes),
            snap=snap,
            batch_size=batch_size,
            cancel=cancel,
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
            if getattr(field, "_halfspaces", None) is not None
            else "retained_constraints"
        ] = len(planes)
        return _check_components(mesh, expected_components)
    feature_options = (
        {"halfspaces": getattr(field, "_halfspaces", None)}
        if planes is not None and adapt and getattr(_meshers, "ADAPTIVE_POLYHEDRAL", False)
        else {}
    )
    feature_options.update(
        {"constraints": list(planes.values())}
        if planes is not None and adapt and getattr(field, "_halfspaces", None) is None
        else {},
    )
    planes = (
        getattr(field, "_polyhedral_regions", None)
        if isinstance(field, CompiledField)
        and not any(scalar_mode)
        and getattr(_meshers, "POLYHEDRAL_REGIONS", False)
        else None
    )
    if planes is not None:
        feature_options = {"halfspaces": planes[0], "regions": planes[1]}
    result = _meshers.generate(
        field,
        gradient,
        [b[::2].tolist(), b[1::2].tolist()],
        cells,
        tuple(bool(v) for v in periodic),
        band,
        geometry_tolerance,
        minimum_quality,
        max_tetrahedra,
        optimize_passes,
        interior_refine_nodes,
        snap,
        threads,
        batch_size,
        cancel,
        (
            geometry_tolerance if min_length is None else float(min_length),
            float(np.min(b[1::2] - b[::2])) / 4 if max_length is None else float(max_length),
            float(gradation),
            float(layers),
            int(adapt_iterations),
            float(quality_guard),
            geometry_tolerance if adapt_tolerance is None else float(adapt_tolerance),
        )
        if adapt
        else None,
        int(order),
        size_field,
        curved_fit=curved_fit,
        **feature_options,
    )
    result["diagnostics"]["threads"] = threads
    result["diagnostics"]["evaluator"] = evaluator
    result["diagnostics"]["preparation_seconds"] = preparation_seconds
    return _check_components(_volume_mesh(result), expected_components)


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
