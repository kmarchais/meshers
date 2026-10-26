"""Direct meshing of bounded implicit solids. CPU float64; Python 3.10+."""

import json
import time
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from . import _meshers
from ._meshers import CancellationToken, CancelledError, CompiledField, MeshingError
from .expression import TraceError, compile_field
from .intersection import (
    _affine_planes,
    _periodic_affine,
    _retained_constraints,
    generate_intersection,
)

if TYPE_CHECKING:
    import fedoo

Field = Callable[
    [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
    npt.NDArray[np.float64],
]


@dataclass(frozen=True)
class Mesh:
    """Owned arrays and diagnostics for a tetrahedral mesh.

    Attributes:
        points: Float64 coordinates, shape ``(N, 3)``.
        tetrahedra: Int64 connectivity, shape ``(M, 4)`` or ``(M, 10)``.
        surface: Boundary triangles, shape ``(K, 3)`` or ``(K, 6)``.
        boundary_tags: Zero for implicit faces, 1 and 2 for x-/x+, then y and z.
        periodic_pairs: Three arrays of paired point IDs, one for each box axis.
            For scalar fields, every material-cap node is matched. Existing
            edge/corner images form complete constraint classes. A wall-only
            contact can lack an opposite image when that image has no
            tetrahedral material support.
        diagnostics: Quality, timing, sampled error, evaluator and worker count.
        parameters: Parameter coordinates for mapped intersections, otherwise None.
        constraint_masks: Bit j marks membership in constraint j, otherwise None.
        constraint_names: Names corresponding to the constraint-mask bits.
        periodic_transforms: Three low-to-high rigid 4x4 transforms for intersections.
            For mapped meshes, use these to rotate vector degrees of freedom; the
            exported PeriodicShift alone only describes point displacement.

    Arrays remain mutable even though the container is frozen. Editing coordinates
    or connectivity can invalidate quality values and periodic metadata.
    """

    points: npt.NDArray[np.float64]
    tetrahedra: npt.NDArray[np.int64]
    surface: npt.NDArray[np.int64]
    boundary_tags: npt.NDArray[np.uint8]
    periodic_pairs: tuple[npt.NDArray[np.int64], ...]
    diagnostics: dict
    parameters: npt.NDArray[np.float64] | None = None
    constraint_masks: npt.NDArray[np.uint8] | npt.NDArray[np.uint32] | None = None
    constraint_names: tuple[str, ...] = ()
    periodic_transforms: npt.NDArray[np.float64] | None = None

    def to_fedoo(self) -> "fedoo.Mesh":
        """Copy linear or curved quadratic tetrahedra into a fedoo mesh.

        Requires the optional fedoo package. Coordinates and the VTK/fedoo node
        order are preserved. Boundary tags and periodic constraints are not
        transferred; configure constraints on the returned mesh before solving.
        Conversion does not check element validity or periodic conformity.

        Returns:
            A fedoo ``tet4`` or ``tet10`` mesh with independent arrays.

        Raises:
            ValueError: Connectivity does not have four or ten nodes per cell.
        """
        import fedoo

        if self.tetrahedra.ndim != 2 or self.tetrahedra.shape[1] not in {4, 10}:
            raise ValueError("fedoo conversion requires Tet4 or Tet10 connectivity")
        element_type = "tet4" if self.tetrahedra.shape[1] == 4 else "tet10"
        return fedoo.Mesh(self.points.copy(), self.tetrahedra.copy(), element_type)

    def surface_for(self, name: str) -> npt.NDArray[np.int64]:
        """Return boundary triangles assigned to a named intersection constraint.

        Returns:
            Triangle connectivity referencing ``points``.

        Raises:
            ValueError: The mesh has no constraint with this name.
        """
        if name not in self.constraint_names or self.constraint_masks is None:
            raise ValueError(f"Unknown constraint: {name}")
        bit = 1 << self.constraint_names.index(name)
        masks = self.constraint_masks[self.surface]
        labels = masks[:, 0] & masks[:, 1] & masks[:, 2]
        return self.surface[(labels & bit) != 0]

    @property
    def feature_edges(self) -> npt.NDArray[np.int64]:
        """Boundary edges shared by at least two named constraints."""
        if self.constraint_masks is None:
            return np.empty((0, 2), dtype=np.int64)
        edges = np.unique(
            np.sort(
                np.concatenate([
                    self.surface[:, [0, 1]],
                    self.surface[:, [1, 2]],
                    self.surface[:, [2, 0]],
                ]),
                axis=1,
            ),
            axis=0,
        )
        common = self.constraint_masks[edges[:, 0]] & self.constraint_masks[edges[:, 1]]
        return edges[(common & (common - np.uint8(1))) != 0]

    def write_vtkhdf(self, path: str | Path, *, surface: bool = False) -> None:  # ruff: ignore[complex-structure] - optional mesh metadata share one file lifecycle.
        """Write VTKHDF 2.0, with physical periodic shifts. Requires meshers[io]."""
        import h5py

        cells = self.surface if surface else self.tetrahedra
        points = self.points
        parent = np.arange(len(points), dtype=np.int64)

        def root(i: int) -> int:
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for pairs in self.periodic_pairs:
            for a, b in pairs:
                ra, rb = root(a), root(b)
                parent[max(ra, rb)] = min(ra, rb)
        masters = np.array([root(i) for i in range(len(points))], dtype=np.int64)
        # Retain full points for both outputs so master IDs retain the same meaning.
        with h5py.File(path, "w") as file:
            g = file.create_group("VTKHDF")
            g.attrs["Version"] = np.array([2, 0], dtype=np.int64)
            g.attrs["Type"] = np.bytes_("UnstructuredGrid")
            g["NumberOfPoints"] = np.array([len(points)], dtype=np.int64)
            g["NumberOfCells"] = np.array([len(cells)], dtype=np.int64)
            g["NumberOfConnectivityIds"] = np.array([cells.size], dtype=np.int64)
            g["Points"] = points
            g["Connectivity"] = cells.ravel()
            g["Offsets"] = np.arange(len(cells) + 1, dtype=np.int64) * cells.shape[1]
            # Linear or quadratic triangle/tetrahedron; corners come first in both.
            vtk_type = {3: 5, 4: 10, 6: 22, 10: 24}[cells.shape[1]]
            g["Types"] = np.full(len(cells), vtk_type, dtype=np.uint8)
            data = g.create_group("CellData")
            if surface:
                data["BoundaryTag"] = self.boundary_tags
            else:
                p = points[cells]
                determinant = np.einsum(
                    "ij,ij->i",
                    p[:, 1] - p[:, 0],
                    np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]),
                )
                edge_sum = sum(
                    np.sum((p[:, i] - p[:, j]) ** 2, axis=1)
                    for i in range(4)
                    for j in range(i + 1, 4)
                )
                data["Volume"] = determinant / 6
                data["MMGQuality"] = np.sqrt(np.maximum(432 * determinant**2 / edge_sum**3, 0))
            data = g.create_group("PointData")
            data["PeriodicMasterId"] = masters
            data["PeriodicShift"] = points - points[masters]
            if self.parameters is not None:
                data["ParameterCoordinates"] = self.parameters
            if self.constraint_masks is not None:
                data["SurfaceConstraints"] = self.constraint_masks
            metadata = g.create_group("FieldData")
            if self.periodic_transforms is not None:
                metadata["PeriodicTransforms"] = self.periodic_transforms.reshape(3, 16)
            for axis, pairs in enumerate(self.periodic_pairs):
                metadata[f"PeriodicPairs{axis}"] = pairs
            if self.constraint_names:
                metadata["ConstraintNamesUTF8"] = np.frombuffer(
                    json.dumps(self.constraint_names, ensure_ascii=False).encode("utf-8"),
                    dtype=np.uint8,
                )
                metadata["FeatureEdges"] = self.feature_edges


# Keep public option validation and pipeline dispatch together.
def generate(  # ruff: ignore[complex-structure]
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
    compile: bool = True,
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
) -> Mesh:
    """Generate a tetrahedral mesh directly from a bounded implicit field.

    Args:
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
    """
    b = np.asarray(bounds, dtype=np.float64)
    if size_field is not None and not adapt:
        raise ValueError("size_field requires adapt=True")
    if b.shape != (6,):
        raise ValueError("bounds must contain xmin,xmax,ymin,ymax,zmin,zmax")
    if isinstance(cells, (int, np.integer)):
        cells = (int(cells),) * 3
    if len(periodic) != 3 or any(not isinstance(v, (bool, np.bool_)) for v in periodic):
        raise ValueError("periodic must contain three booleans")
    if threads is None:
        threads = _meshers.available_threads()
    preparation_start = time.perf_counter()
    gradient_override = gradient is not None
    evaluator = (
        "native"
        if isinstance(field, str)
        else "compiled"
        if isinstance(field, CompiledField)
        else "python"
    )
    if callable(field) and not isinstance(field, CompiledField) and compile:
        try:
            field = compile_field(field, gradient=gradient)
            gradient = None
            evaluator = "compiled"
        except TraceError as error:
            warnings.warn(f"Using NumPy callbacks: {error}", RuntimeWarning, stacklevel=2)
    scalar_mode = (
        order != 1,
        band is not None,
        any(periodic) and not _periodic_affine(field, periodic, b if adapt else None),
        threads != 1,
        interior_refine_nodes != 0,
        gradient_override,
        not 0 <= optimize_passes <= 20,
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
        planes = _retained_constraints(field)
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
        return mesh
    feature_options = (
        {"halfspaces": getattr(field, "_halfspaces", None)}
        if planes is not None and adapt and getattr(_meshers, "ADAPTIVE_POLYHEDRAL", False)
        else {}
    )
    feature_options.update(
        {"constraints": list(planes.values())}
        if planes is not None and adapt and getattr(field, "_halfspaces", None) is None
        else {}
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
        **feature_options,
    )
    quadratic = result.pop("quadratic", None)
    if quadratic is not None:
        for name in ("points", "tetrahedra", "surface", "periodic_pairs"):
            result[name] = quadratic.pop(name)
        result["diagnostics"]["quadratic"] = quadratic
    result["diagnostics"]["threads"] = threads
    result["diagnostics"]["evaluator"] = evaluator
    result["diagnostics"]["preparation_seconds"] = preparation_seconds
    result["periodic_pairs"] = tuple(result["periodic_pairs"])
    return Mesh(**result)


__all__ = [
    "CancellationToken",
    "CancelledError",
    "CompiledField",
    "Field",
    "Mesh",
    "MeshingError",
    "compile_field",
    "generate",
    "generate_intersection",
]

if hasattr(_meshers, "generate_surface"):
    from .surface import SurfaceMesh, generate_surface

    __all__ += ["SurfaceMesh", "generate_surface"]
