# Copyright (c) 2026 Kevin Marchais

"""Named implicit intersections with retained curves and mapped periodic ends."""

from __future__ import annotations

import time
import warnings
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import TypedDict, cast

import numpy as np
import numpy.typing as npt

from . import _meshers
from .expression import TraceError, compile_field
from .mesh import Mesh

CLOSED_MAPPED_SEAMS = True
ADAPTIVE_MAPPED_INTERSECTIONS = True
SPATIAL_DIMENSIONS = 3
MAX_INTERSECTION_PASSES = 20
_MIN_CONSTRAINTS = 2
_MIN_GRID_CELLS = 4
_MAX_GRID_CELLS = 128
_MAX_SNAP = 0.2
_TET4_NODES = 4
_LOW_QUALITY = 0.1
_MAX_SEED_ATTEMPTS = 3

Scalar = Callable | _meshers.CompiledField | str
Points = npt.NDArray[np.float64]
Indices = npt.NDArray[np.int64]


def _periodic_affine(
    field: Scalar,
    periodic: Sequence[bool],
    bounds: np.ndarray | None = None,
) -> bool:
    if not isinstance(field, _meshers.CompiledField) or not getattr(
        _meshers,
        "PERIODIC_POLYHEDRAL",
        False,
    ):
        return False
    planes = getattr(field, "_halfspaces", None)
    certificate = getattr(_meshers, "_periodic_affine_traces", None)
    groups = getattr(field, "_polyhedral_regions", None)
    if (
        planes is None
        and groups is not None
        and bounds is not None
        and certificate is not None
        and getattr(_meshers, "PERIODIC_AFFINE_REGIONS", False)
    ):
        return certificate(
            groups[0],
            [bounds[::2].tolist(), bounds[1::2].tolist()],
            tuple(bool(v) for v in periodic),
            groups[1],
        )
    if planes is not None and bounds is not None and certificate is not None:
        return certificate(
            planes,
            [bounds[::2].tolist(), bounds[1::2].tolist()],
            tuple(bool(v) for v in periodic),
        )
    return planes is not None and all(
        not enabled or all(plane[axis] == 0 for plane in planes)
        for axis, enabled in enumerate(periodic)
    )


def _affine_planes(
    field: Scalar,
    scalar_mode: Sequence[bool],
    *,
    canonical: bool = False,
) -> dict[str, Scalar] | None:
    if not isinstance(field, _meshers.CompiledField) or any(scalar_mode):
        return None
    halfspaces = getattr(field, "_halfspaces", None)
    if halfspaces is None:
        return None
    if canonical:
        halfspaces = sorted(halfspaces)
    planes: dict[str, Scalar] = {}
    for index, coefficients in enumerate(halfspaces):

        def plane(
            x: Points,
            y: Points,
            z: Points,
            coefficients: Sequence[float] = coefficients,
        ) -> Points:
            return coefficients[0] * x + coefficients[1] * y + coefficients[2] * z + coefficients[3]

        planes[f"plane{index}"] = plane
    return planes


def _retained_constraints(
    field: Scalar,
    bounds: np.ndarray | None = None,
) -> dict[str, Scalar] | None:
    if not isinstance(field, _meshers.CompiledField) or not getattr(
        _meshers,
        "RETAINED_INTERSECTIONS",
        False,
    ):
        return None
    # This is the internal Python/native interface, gated by its capability above.
    fields = field._intersection_fields()  # ruff: ignore[private-member-access]
    if fields is None:
        return None
    if getattr(field, "_sdf_intersection", False):
        if bounds is None:
            return None
        # Check floating-point eligibility at the same 64 preflight sites as
        # the scalar generator. Refuse this optimization on an observed
        # zero-sublevel mismatch from squared-norm underflow.
        u = np.array([0.0, 0.23, 0.61, 1.0])
        grid = np.array(np.meshgrid(u, u, u, indexing="ij")).reshape(3, -1).T
        points = bounds[::2] + grid * (bounds[1::2] - bounds[::2])
        parent = field.values(points)
        branches = np.max(np.stack([branch.values(points) for branch in fields]), axis=0)
        if not np.array_equal(parent <= 0, branches <= 0):
            return None
    constraints = {}
    for index, branch in enumerate(fields):
        planes = _affine_planes(branch, (), canonical=True)
        if planes is None:
            constraints[f"branch{index}"] = branch
        else:
            constraints.update((f"branch{index}_{name}", plane) for name, plane in planes.items())
    return constraints if len(constraints) <= _meshers.MAX_CONSTRAINTS else None


def _check_cancel(cancel: _meshers.CancellationToken | None) -> None:
    if cancel is not None and cancel.cancelled:
        msg = "generation cancelled"
        raise _meshers.CancelledError(msg)


def _prepare(field: Scalar, *, compile_expression: bool) -> Scalar:
    if callable(field) and not isinstance(field, _meshers.CompiledField) and compile_expression:
        try:
            return compile_field(field)
        except TraceError as error:
            warnings.warn(f"Using NumPy callbacks: {error}", RuntimeWarning, stacklevel=3)
    return field


def _values(field: Scalar, points: Points) -> Points:
    if isinstance(field, _meshers.CompiledField):
        return field.values(points)
    if isinstance(field, str):
        x, y, z = (2 * np.pi * points).T
        if field == "gyroid":
            return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)
        if field == "schwarz_p":
            return np.cos(x) + np.cos(y) + np.cos(z)
        return (
            np.sin(x) * np.sin(y) * np.sin(z)
            + np.sin(x) * np.cos(y) * np.cos(z)
            + np.cos(x) * np.sin(y) * np.cos(z)
            + np.cos(x) * np.cos(y) * np.sin(z)
        )
    return np.broadcast_to(np.asarray(field(*points.T), dtype=np.float64), (len(points),))


def _mapped(points: Points, mapping: list[Scalar] | None) -> Points:
    return points if mapping is None else np.column_stack([_values(f, points) for f in mapping])


def _jacobian(points: Points, mapping: list[Scalar] | None, step: Points) -> Points:
    if mapping is None:
        return np.broadcast_to(np.eye(3), (len(points), 3, 3))
    return np.stack(
        [
            (
                _mapped(points + np.eye(3)[a] * step[a], mapping)
                - _mapped(points - np.eye(3)[a] * step[a], mapping)
            )
            / (2 * step[a])
            for a in range(3)
        ],
        axis=2,
    )


def _periodic_extension(field: Scalar, bounds: Points, closed_axes: Sequence[bool]) -> Scalar:
    """Return a periodic extension for sampled or analytic closed-chart fields."""
    if not any(closed_axes):
        return field

    def wrapped(x: Points, y: Points, z: Points) -> Points:
        points = np.column_stack((x, y, z))
        for axis, closed in enumerate(closed_axes):
            if closed:
                lo, hi = bounds[2 * axis : 2 * axis + 2]
                points[:, axis] = lo + np.mod(points[:, axis] - lo, hi - lo)
        return _values(field, points)

    return wrapped


def _gradient(field: Scalar, points: Points, step: Points) -> Points:
    """Return central differences in the field's parameter coordinates."""
    return np.column_stack([
        (
            _values(field, points + np.eye(3)[a] * step[a])
            - _values(field, points - np.eye(3)[a] * step[a])
        )
        / (2 * step[a])
        for a in range(3)
    ])


def _project_parameters(
    proposed: Points,
    origin: Points,
    fields: list[Scalar],
    *,
    locked: npt.NDArray[np.bool_],
    step: Points,
) -> Points | None:
    """Return a joint projection preserving cap locks, or None if unresolved."""
    point = proposed.copy()
    point[locked] = origin[locked]
    precision = float(np.max(step)) * 1e-7
    for _ in range(24):
        value = np.array([_values(f, point[None])[0] for f in fields])
        gradient = np.array([_gradient(f, point[None], step)[0] for f in fields])
        magnitude = np.linalg.norm(gradient, axis=1)
        if np.any(magnitude <= 0) or not np.all(np.isfinite(magnitude)):
            return None
        gradient /= magnitude[:, None]
        value /= magnitude
        if np.max(abs(value)) <= precision:
            return point
        gradient[:, locked] = 0
        delta = np.linalg.lstsq(gradient, value, rcond=None)[0]
        if not np.all(np.isfinite(delta)) or np.linalg.norm(gradient @ delta - value) > precision:
            return None
        point -= delta
        point[locked] = origin[locked]
    return None


_FACE_SAMPLES = np.array([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
    [0.5, 0.5, 0],
    [0.5, 0, 0.5],
    [0, 0.5, 0.5],
    [1 / 3, 1 / 3, 1 / 3],
])


@dataclass(frozen=True, slots=True)
class _MapInversion:
    mapping: list[Scalar] | None
    step: Points
    extent: float
    inverse_tolerance: float

    def distances(self, field: Scalar, p: Points, physical: Points) -> Points:
        # Resolve the map in physical units, accounting for representable
        # coordinate precision when the mesh is far from the origin.
        roundoff = (
            8 * np.finfo(np.float64).eps * np.maximum(np.max(np.abs(physical), axis=1), self.extent)
        )
        inverse_limit = np.maximum(self.inverse_tolerance, roundoff)
        for _ in range(8):
            residual = _mapped(p, self.mapping) - physical
            if np.all(np.linalg.norm(residual, axis=1) <= inverse_limit):
                break
            p -= np.linalg.solve(_jacobian(p, self.mapping, self.step), residual[..., None])[..., 0]
        inverse_error = np.linalg.norm(_mapped(p, self.mapping) - physical, axis=1)
        if not np.all(np.isfinite(inverse_error)) or np.any(inverse_error > inverse_limit):
            msg = "coordinate map could not be inverted for surface checks"
            raise _meshers.MeshingError(msg)
        gradient = _gradient(field, p, self.step)
        physical_gradient = (
            gradient
            if self.mapping is None
            else np.linalg.solve(
                np.swapaxes(_jacobian(p, self.mapping, self.step), 1, 2),
                gradient[..., None],
            )[..., 0]
        )
        magnitude = np.linalg.norm(physical_gradient, axis=1)
        value = np.abs(_values(field, p))
        distance = np.divide(value, magnitude, out=np.full_like(value, np.inf), where=magnitude > 0)
        # Do not hide unresolved physical displacement in the distance gate.
        distance += inverse_error
        if not np.all(np.isfinite(distance)):
            msg = "nonfinite or singular sampled surface distance"
            raise _meshers.MeshingError(msg)
        return distance


def _sample_chart(
    parameters: Points,
    chunk: Indices,
    sample_weights: Points,
    bounds: Points,
    closed_axes: Sequence[bool],
) -> Points:
    chart = parameters[chunk]
    for axis, closed in enumerate(closed_axes):
        if closed:
            period = bounds[2 * axis + 1] - bounds[2 * axis]
            delta = chart[:, :, axis] - chart[:, :1, axis]
            chart[:, :, axis] = chart[:, :1, axis] + delta - np.rint(delta / period) * period
    return np.einsum("sj,sjk->sk", sample_weights, chart)


def _wall_faces(result: _meshers._IntersectionResult, field_count: int) -> Iterator[Indices]:
    masks = result["constraint_masks"][result["surface"]]
    labels = masks[:, 0] & masks[:, 1] & masks[:, 2]
    for index in range(field_count):
        yield np.flatnonzero(labels & (1 << index))


def _sample_batches(selection: Indices, batch_size: int) -> Iterator[tuple[Indices, Points]]:
    count = len(selection) * len(_FACE_SAMPLES)
    limit = min(batch_size, 512 * len(_FACE_SAMPLES))
    for start in range(0, count, limit):
        samples = np.arange(start, min(start + limit, count))
        yield selection[samples // len(_FACE_SAMPLES)], _FACE_SAMPLES[samples % len(_FACE_SAMPLES)]


@dataclass(frozen=True, slots=True, kw_only=True)
class _SurfaceAudit:
    mapping: list[Scalar] | None
    bounds: Points
    closed_axes: Sequence[bool]
    geometry_tolerance: float
    batch_size: int
    cancel: _meshers.CancellationToken | None

    def errors(
        self,
        result: _meshers._IntersectionResult,
        fields: list[Scalar],
        *,
        errors_by_face: np.ndarray | None = None,
        face_selections: list[Indices] | None = None,
    ) -> list[float]:
        """Sample physical facets against analytic constraints, including map inversion.

        Returns:
            Largest sampled distance estimate per constraint.

        """
        step = (self.bounds[1::2] - self.bounds[::2]) * 1e-5
        points, parameters, faces = result["points"], result["parameters"], result["surface"]
        extent = float(np.max(np.ptp(points, axis=0))) if len(points) else 0.0
        inverse_tolerance = min(self.geometry_tolerance * 0.01, extent * 1e-12)
        inverse = _MapInversion(self.mapping, step, extent, inverse_tolerance)
        groups = _wall_faces(result, len(fields)) if face_selections is None else face_selections
        errors = []
        for field, selection in zip(fields, groups, strict=True):
            maximum = 0.0
            for indices, sample_weights in _sample_batches(selection, self.batch_size):
                _check_cancel(self.cancel)
                chunk = faces[indices]
                physical = np.einsum("sj,sjk->sk", sample_weights, points[chunk])
                p = _sample_chart(parameters, chunk, sample_weights, self.bounds, self.closed_axes)
                distance = inverse.distances(field, p, physical)
                if errors_by_face is not None:
                    np.maximum.at(errors_by_face, indices, distance)
                maximum = max(maximum, float(np.max(distance, initial=0)))
            errors.append(maximum)
        _check_cancel(self.cancel)
        return errors

    def check_walls(
        self,
        result: _meshers._IntersectionResult,
        fields: Mapping[str, Scalar],
        cells: Sequence[int],
        minimum_quality: float,
        *,
        can_fit: bool,
    ) -> dict[str, float]:
        diagnostics = result["diagnostics"]
        wall_fields = list(fields.values())
        try:
            errors = self.errors(result, wall_fields)
        except np.linalg.LinAlgError as error:
            msg = "singular coordinate map during surface checks"
            raise _meshers.MeshingError(msg) from error
        if max(errors) > self.geometry_tolerance and can_fit:
            fitting_start = time.perf_counter()
            moved = _BoundaryFit(result, wall_fields, self, minimum_quality).fit(cells)
            diagnostics["geometry_vertex_moves"] = moved
            diagnostics["surface_fitting_seconds"] = time.perf_counter() - fitting_start
            errors = self.errors(result, wall_fields)
            _refresh_linear_metrics(result, fitted=bool(moved))
        if max(errors) > self.geometry_tolerance:
            msg = f"sampled surface error {max(errors):g} exceeds {self.geometry_tolerance:g}"
            raise _meshers.MeshingError(msg)
        return dict(zip(fields, errors, strict=True))

    def caps(self, result: _meshers._IntersectionResult, tags: np.ndarray) -> dict[str, float]:
        """Return physical distance estimates for the chart caps present in the mesh.

        Raises:
            MeshingError: A measured cap exceeds tolerance or its map is singular.
        """
        if self.mapping is None:
            return {}
        cap_fields: list[Scalar] = []
        selections = []
        cap_names = []
        for side, bound in enumerate(self.bounds):
            axis = side // 2
            indices = np.flatnonzero(tags == side + 1)
            if self.closed_axes[axis] or not len(indices):
                continue

            def cap(
                x: Points,
                y: Points,
                z: Points,
                axis: int = axis,
                bound: float = bound,
            ) -> Points:
                return (x, y, z)[axis] - bound

            cap_fields.append(cap)
            selections.append(indices)
            cap_names.append(("x-", "x+", "y-", "y+", "z-", "z+")[side])
        try:
            cap_values = self.errors(result, cap_fields, face_selections=selections)
        except np.linalg.LinAlgError as error:
            msg = "singular coordinate map during cap checks"
            raise _meshers.MeshingError(msg) from error
        for name, value in zip(cap_names, cap_values, strict=True):
            if value > self.geometry_tolerance:
                msg = (
                    f"sampled mapped cap {name} error {value:g} exceeds {self.geometry_tolerance:g}"
                )
                raise _meshers.MeshingError(msg)
        return dict(zip(cap_names, cap_values, strict=True))


def _signed_quality(vertices: Points) -> Points:
    """Return signed MMG quality for each linear tetrahedron."""
    determinant = np.linalg.det(vertices[:, 1:] - vertices[:, :1])
    edges = sum(
        np.sum((vertices[:, a] - vertices[:, b]) ** 2, axis=1)
        for a in range(4)
        for b in range(a + 1, 4)
    )
    return np.sign(determinant) * np.sqrt(432 * determinant**2 / edges**3)


def _refresh_linear_metrics(result: _meshers._IntersectionResult, *, fitted: bool) -> None:
    """Update volume and any native shape statistics invalidated by fitting."""
    diagnostics = result["diagnostics"]
    vertices = result["points"][result["tetrahedra"]]
    if fitted:
        final_quality = _signed_quality(vertices)
        diagnostics["minimum_mmg_quality"] = float(np.min(final_quality))
        diagnostics["elements_below_quality_01"] = int(
            np.count_nonzero(final_quality < _LOW_QUALITY),
        )
    diagnostics["volume"] = float(np.linalg.det(vertices[:, 1:] - vertices[:, :1]).sum() / 6)


def _incident_cells(connectivity: Indices, node_count: int) -> list[list[int]]:
    incident: list[list[int]] = [[] for _ in range(node_count)]
    for index, cell in enumerate(connectivity):
        for node in cell:
            incident[node].append(index)
    return incident


def _trial_positions(origin: Points, cell: Points) -> Iterator[Points]:
    yield origin.copy()
    for axis in range(SPATIAL_DIMENSIONS):
        for factor in (-0.5, -0.2, -0.05, 0.05, 0.2, 0.5):
            candidate = origin.copy()
            candidate[axis] += factor * cell[axis]
            yield candidate


@dataclass
class _BoundaryFit:
    result: _meshers._IntersectionResult
    fields: list[Scalar]
    audit: _SurfaceAudit
    minimum_quality: float

    def measure(self, faces: Indices) -> Points:
        local = self.result.copy()
        local["surface"] = faces
        errors = np.zeros(len(faces))
        self.audit.errors(local, self.fields, errors_by_face=errors)
        return errors

    def feasible(self, candidate: Points, step: Points) -> bool:
        if any(
            not closed
            and not self.audit.bounds[2 * axis]
            <= candidate[axis]
            <= self.audit.bounds[2 * axis + 1]
            for axis, closed in enumerate(self.audit.closed_axes)
        ):
            return False
        return not any(
            _values(f, candidate[None])[0]
            > 1e-10 * np.linalg.norm(_gradient(f, candidate[None], step)[0])
            for f in self.fields
        )

    def best_position(
        self,
        node: int,
        tets: Indices,
        faces: Indices,
        cell: Points,
        step: Points,
    ) -> Points | None:
        points = self.result["points"]
        origin = self.result["parameters"][node].copy()
        locked = np.any(origin[:, None] == self.audit.bounds.reshape(3, 2), axis=1) & ~np.asarray(
            self.audit.closed_axes,
        )
        active = [
            f for j, f in enumerate(self.fields) if self.result["constraint_masks"][node] & (1 << j)
        ]
        if not active:
            return None
        score = float(np.max(self.measure(faces)))
        old = points[node].copy()
        best = None
        try:
            for trial in _trial_positions(origin, cell):
                candidate = _project_parameters(trial, origin, active, locked=locked, step=step)
                if candidate is None or not self.feasible(candidate, step):
                    continue
                points[node] = _mapped(candidate[None], self.audit.mapping)[0]
                quality = _signed_quality(points[tets])
                if not np.all(np.isfinite(quality)) or np.min(quality) < max(
                    self.minimum_quality,
                    1e-14,
                ):
                    continue
                new = float(np.max(self.measure(faces)))
                if new < score - 1e-10:
                    score, best = new, candidate
        finally:
            points[node] = old
        return best

    def fit(self, cells: Sequence[int]) -> int:
        faces, tets = self.result["surface"], self.result["tetrahedra"]
        node_count = len(self.result["points"])
        incident = _incident_cells(tets, node_count)
        boundary = _incident_cells(faces, node_count)
        span = self.audit.bounds[1::2] - self.audit.bounds[::2]
        step, cell = span * 1e-5, span / np.asarray(cells)
        moved = 0
        for _ in range(8):
            _check_cancel(self.audit.cancel)
            errors = self.measure(faces)
            if not np.any(errors > self.audit.geometry_tolerance):
                break
            ranked = np.argsort(errors)[::-1]
            vertices = dict.fromkeys(
                faces[ranked[errors[ranked] > self.audit.geometry_tolerance]].ravel(),
            )
            accepted = 0
            for node in list(vertices)[:128]:
                _check_cancel(self.audit.cancel)
                best = self.best_position(
                    int(node),
                    tets[incident[node]],
                    faces[boundary[node]],
                    cell,
                    step,
                )
                if best is not None:
                    for axis, closed in enumerate(self.audit.closed_axes):
                        if closed:
                            low = self.audit.bounds[2 * axis]
                            best[axis] = low + (best[axis] - low) % span[axis]
                    self.result["parameters"][node] = best
                    self.result["points"][node] = _mapped(best[None], self.audit.mapping)[0]
                    accepted += 1
            moved += accepted
            if not accepted:
                break
        return moved


def _check_rigid_pairs(
    points: Points,
    low: Indices,
    high: Indices,
    transform: Points,
    scale: float,
) -> None:
    mapped = points[low] @ transform[:3, :3].T + transform[:3, 3]
    # Bound rounding in the rigid transform, including cancellation between
    # large rotated coordinates and translations. Keep the geometry limit.
    roundoff = (
        8
        * np.finfo(np.float64).eps
        * (
            np.abs(points[low]) @ np.abs(transform[:3, :3]).T
            + np.abs(transform[:3, 3])
            + np.abs(points[high])
        )
    )
    if not np.all(np.isfinite(mapped)) or not np.all(
        np.abs(mapped - points[high]) <= scale * 1e-10 + roundoff,
    ):
        msg = "physical periodic nodes do not match the supplied transform"
        raise _meshers.MeshingError(msg)


def _periodic_metadata(
    result: _meshers._IntersectionResult,
    bounds: Points,
    *,
    periodic: Sequence[bool],
    transforms: Points,
) -> tuple[npt.NDArray[np.uint8], tuple[Indices, ...]]:
    points, parameters, faces = result["points"], result["parameters"], result["surface"]
    tags = np.zeros(len(faces), dtype=np.uint8)
    pairs = []
    for axis in range(3):
        low = parameters[:, axis] == bounds[2 * axis]
        high = parameters[:, axis] == bounds[2 * axis + 1]
        tags[np.all(low[faces], axis=1)] = 2 * axis + 1
        tags[np.all(high[faces], axis=1)] = 2 * axis + 2
        if not periodic[axis]:
            pairs.append(np.empty((0, 2), dtype=np.int64))
            continue
        other = [a for a in range(3) if a != axis]
        a, b = np.flatnonzero(low), np.flatnonzero(high)
        a = a[np.lexsort(parameters[a][:, other].T)]
        b = b[np.lexsort(parameters[b][:, other].T)]
        scale = max(float(np.max(bounds[1::2] - bounds[::2])), 1.0)
        if len(a) != len(b) or not np.allclose(
            parameters[a][:, other],
            parameters[b][:, other],
            rtol=0,
            atol=scale * 1e-11,
        ):
            msg = "periodic parameter nodes do not match"
            raise _meshers.MeshingError(msg)
        _check_rigid_pairs(points, a, b, transforms[axis], scale)
        correspondence = np.arange(len(points))
        correspondence[a] = b
        left = np.sort(correspondence[faces[np.all(low[faces], axis=1)]], axis=1)
        right = np.sort(faces[np.all(high[faces], axis=1)], axis=1)
        if {tuple(face) for face in left} != {tuple(face) for face in right}:
            msg = "periodic end triangles do not match"
            raise _meshers.MeshingError(msg)
        pairs.append(np.column_stack((a, b)))
    return tags, tuple(pairs)


def generate_intersection(  # ruff: ignore[too-many-arguments] - published keyword options.
    constraints: Mapping[str, Callable | _meshers.CompiledField | str],
    *,
    bounds: Sequence[float] = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
    cells: int | Sequence[int] = 24,
    coordinate_map: Callable | None = None,
    periodic: Sequence[bool] = (False, False, False),
    periodic_transforms: Mapping[int, Sequence[Sequence[float]]] | None = None,
    closed_axes: Sequence[bool] = (False, False, False),
    adapt: bool = False,
    inverse_map: Callable | None = None,
    physical_bounds: Sequence[float] | None = None,
    seed_tolerance: float = 0.01,
    min_length: float | None = None,
    max_length: float | None = None,
    gradation: float = 0.3,
    adapt_iterations: int = 4,
    quality_guard: float = 0.3,
    adapt_tolerance: float | None = None,
    geometry_tolerance: float | None = None,
    minimum_quality: float = 0.0,
    optimize_passes: int = 4,
    snap: float = 0.2,
    max_tetrahedra: int = 2_000_000,
    compile: bool = True,  # ruff: ignore[builtin-argument-shadowing] - public keyword.
    batch_size: int = 4096,
    cancel: _meshers.CancellationToken | None = None,
) -> Mesh:
    """Mesh the intersection of named negative sublevel sets, preserving their curves.

    Constraints are evaluated in parameter coordinates. A coordinate map returns
    a tuple ``(X(x,y,z), Y(x,y,z), Z(x,y,z))`` and acts on the background before
    cutting; it must preserve orientation and be locally invertible. Straight-sided
    tetrahedra approximate the mapped geometry. This path currently uses one worker.

    Args:
        constraints: Ordered mapping of 2-32 unique nonempty names to fields
            (8 without wide-constraints). Keep
            every field <= 0. Use two constraints for the sides of a TPMS sheet.
        bounds: Parameter box ``(xmin,xmax,ymin,ymax,zmin,zmax)``.
        cells: Background intervals, 4-128 per parameter axis.
        coordinate_map: Optional parameter-to-physical map; default identity.
        periodic: Matching low/high faces along the three parameter axes.
        periodic_transforms: Rigid homogeneous 4x4 matrices keyed by periodic axis
            (0, 1, 2), mapping low to high physical points. Required for each mapped
            periodic axis. Identity maps default to box translations. Node and
            triangle correspondence are checked after optimization.
        closed_axes: Identify coincident periodic caps before final repair. Requires
            identity periodic transforms; closed seams have no external node pairs.
        adapt: Adapt a validated mapped seed in physical space. Requires both maps
            and physical bounds; external periodic faces are not supported.
        inverse_map: Physical-to-parameter coordinates for mapped adaptation.
        physical_bounds: Physical bounding box used for adaptation.
        seed_tolerance: Separate sampled distance limit for the starting mesh.
        min_length: Smallest target length during adaptation.
        max_length: Largest target length during adaptation.
        gradation: Maximum target-length growth per unit distance.
        adapt_iterations: Number of size-driven adaptation rounds.
        quality_guard: Local quality floor for adaptive operations.
        adapt_tolerance: Adaptive distance target; defaults to the final tolerance.
        geometry_tolerance: Optional maximum sampled physical surface-distance
            estimate. None skips this relatively expensive check. This is an
            acceptance gate, not adaptive refinement or a Hausdorff bound.
        minimum_quality: Required final minimum MMG quality, 0-1; zero disables gate.
        optimize_passes: Integrated constrained optimization rounds, 0-20.
        snap: Joint background snap fraction, 0-0.2.
        max_tetrahedra: Output element budget.
        compile: Automatically compile supported field and map expressions.
        batch_size: Maximum callback batch size, 1-65536.
        cancel: Cooperative cancellation token.

    Returns:
        Mesh with physical points, parameter coordinates, constraint masks/names,
        feature edges, validated periodic pairs/transforms and quality diagnostics.

    Raises:
        ValueError: An option, name, transform or callback output is invalid.
        MeshingError: Construction, geometry, topology, quality or periodic checks fail.
        CancelledError: Generation is cancelled.

    Note:
        Constraints remain separate throughout construction. Optimization can move
        vertices on their assigned analytic surfaces and change the facet approximation.
        Acute input angles can make a requested minimum quality unattainable.
    """  # ruff: ignore[docstring-extraneous-exception] - includes native errors.
    if adapt:
        return _adapt_mapped_intersection(
            constraints,
            _MappedDomain(
                bounds=bounds,
                cells=cells,
                mapping=coordinate_map,
                inverse_map=inverse_map,
                physical_bounds=physical_bounds,
                periodic=periodic,
                transforms=periodic_transforms,
                closed_axes=closed_axes,
            ),
            _MeshRequirements(
                seed_tolerance=seed_tolerance,
                final_tolerance=geometry_tolerance,
                minimum_quality=minimum_quality,
                snap=snap,
            ),
            _IntersectionControls(
                optimize_passes=optimize_passes,
                max_tetrahedra=max_tetrahedra,
                compile=compile,
                batch_size=batch_size,
                cancel=cancel,
            ),
            _AdaptationSizing(
                min_length=min_length,
                max_length=max_length,
                gradation=gradation,
                adapt_iterations=adapt_iterations,
                quality_guard=quality_guard,
                adapt_tolerance=adapt_tolerance,
            ),
        )
    _check_constraint_names(constraints)
    b = _checked_bounds(bounds)
    cells = _checked_cells(cells)
    _check_quality_options(optimize_passes, snap, minimum_quality, geometry_tolerance)
    transforms = _periodic_transforms(b, periodic, coordinate_map, periodic_transforms)
    _check_closed_axes(closed_axes, periodic, transforms)
    start = time.perf_counter()
    fields = [_prepare(f, compile_expression=compile) for f in constraints.values()]
    mapping = None
    if coordinate_map is not None:
        mapping = [
            _prepare(lambda x, y, z, a=a: coordinate_map(x, y, z)[a], compile_expression=compile)
            for a in range(3)
        ]
    preparation = time.perf_counter() - start
    result = _meshers.generate_intersection(
        fields,
        mapping,
        [b[::2].tolist(), b[1::2].tolist()],
        cells,
        tuple(periodic),
        tuple(closed_axes),
        minimum_quality,
        max_tetrahedra,
        optimize_passes,
        snap,
        batch_size,
        cancel,
    )
    _check_cancel(cancel)
    tags, pairs = _periodic_metadata(
        result,
        b,
        periodic=tuple(p and not c for p, c in zip(periodic, closed_axes, strict=True)),
        transforms=transforms,
    )
    diagnostics = result["diagnostics"]
    diagnostics["preparation_seconds"] = preparation
    diagnostics["evaluator"] = (
        "compiled" if all(isinstance(f, _meshers.CompiledField) for f in fields) else "mixed"
    )
    diagnostics["empty"] = len(result["tetrahedra"]) == 0
    _refresh_linear_metrics(result, fitted=False)
    diagnostics["sampled_surface_error"] = None
    if geometry_tolerance is not None:
        wall_fields = {
            name: _periodic_extension(f, b, closed_axes)
            for name, f in zip(constraints, fields, strict=True)
        }
        if mapping is not None:
            mapping = [_periodic_extension(f, b, closed_axes) for f in mapping]
        audit = _SurfaceAudit(
            mapping=mapping,
            bounds=b,
            closed_axes=closed_axes,
            geometry_tolerance=geometry_tolerance,
            batch_size=batch_size,
            cancel=cancel,
        )
        errors = audit.check_walls(
            result,
            wall_fields,
            cells,
            minimum_quality,
            can_fit=not any(p and not c for p, c in zip(periodic, closed_axes, strict=True)),
        )
        caps = audit.caps(result, tags)
        diagnostics["surface_errors"] = errors
        diagnostics["mapped_cap_errors"] = caps
        diagnostics["sampled_surface_error"] = max([*errors.values(), *caps.values()])
    diagnostics["workflow_seconds"] = time.perf_counter() - start
    _check_cancel(cancel)
    return Mesh(
        **result,
        boundary_tags=tags,
        periodic_pairs=pairs,
        periodic_transforms=transforms,
        constraint_names=tuple(constraints),
    )


def _check_constraint_names(constraints: Mapping[str, Scalar]) -> None:
    maximum = _meshers.MAX_CONSTRAINTS
    if not isinstance(constraints, Mapping) or not _MIN_CONSTRAINTS <= len(constraints) <= maximum:
        msg = f"constraints must be a mapping of 2-{maximum} named fields"
        raise ValueError(msg)
    if any(not isinstance(name, str) or not name.strip() for name in constraints):
        msg = "constraint names must be nonempty strings"
        raise ValueError(msg)


def _checked_bounds(bounds: Sequence[float]) -> Points:
    b = np.asarray(bounds, dtype=np.float64)
    if b.shape != (6,) or not np.all(np.isfinite(b)) or np.any(b[1::2] <= b[::2]):
        msg = "bounds must contain three finite increasing intervals"
        raise ValueError(msg)
    return b


def _checked_cells(cells: int | Sequence[int]) -> Sequence[int]:
    cells = (int(cells),) * 3 if isinstance(cells, (int, np.integer)) else cells
    if len(cells) != SPATIAL_DIMENSIONS or any(
        not isinstance(v, (int, np.integer)) or not _MIN_GRID_CELLS <= v <= _MAX_GRID_CELLS
        for v in cells
    ):
        msg = "cells must contain three integers in 4..128"
        raise ValueError(msg)
    return cells


def _check_quality_options(
    optimize_passes: int,
    snap: float,
    minimum_quality: float,
    geometry_tolerance: float | None,
) -> None:
    if not isinstance(optimize_passes, int) or not 0 <= optimize_passes <= MAX_INTERSECTION_PASSES:
        msg = "optimize_passes must be an integer in 0..20"
        raise ValueError(msg)
    if not np.isfinite(snap) or not 0 <= snap <= _MAX_SNAP:
        msg = "snap must be finite and in 0..0.2"
        raise ValueError(msg)
    if not np.isfinite(minimum_quality) or not 0 <= minimum_quality <= 1:
        msg = "minimum_quality must be finite and in 0..1"
        raise ValueError(msg)
    if geometry_tolerance is not None and (
        not np.isfinite(geometry_tolerance) or geometry_tolerance <= 0
    ):
        msg = "geometry_tolerance must be positive and finite, or None"
        raise ValueError(msg)


def _check_closed_axes(
    closed_axes: Sequence[bool],
    periodic: Sequence[bool],
    transforms: Points,
) -> None:
    if len(closed_axes) != SPATIAL_DIMENSIONS or any(
        not isinstance(v, (bool, np.bool_)) for v in closed_axes
    ):
        msg = "closed_axes must contain three booleans"
        raise ValueError(msg)
    if any(
        closed_axes[a]
        and (not periodic[a] or not np.allclose(transforms[a], np.eye(4), rtol=0, atol=1e-12))
        for a in range(3)
    ):
        msg = "closed axes require periodic identity transforms"
        raise ValueError(msg)


def _check_periodic_axes(periodic: Sequence[bool]) -> None:
    if len(periodic) != SPATIAL_DIMENSIONS or any(
        not isinstance(v, (bool, np.bool_)) for v in periodic
    ):
        msg = "periodic must contain three booleans"
        raise ValueError(msg)


def _periodic_transforms(
    bounds: Points,
    periodic: Sequence[bool],
    coordinate_map: Callable | None,
    periodic_transforms: Mapping[int, Sequence[Sequence[float]]] | None,
) -> Points:
    _check_periodic_axes(periodic)
    transforms = np.tile(np.eye(4), (3, 1, 1))
    given = {} if periodic_transforms is None else dict(periodic_transforms)
    if any(a not in {0, 1, 2} or not periodic[a] for a in given):
        msg = "periodic_transforms keys must be enabled axes 0, 1, 2"
        raise ValueError(msg)
    for axis in range(3):
        if periodic[axis]:
            if coordinate_map is not None and axis not in given:
                msg = "mapped periodic axes require explicit periodic_transforms"
                raise ValueError(msg)
            transforms[axis, axis, 3] = bounds[2 * axis + 1] - bounds[2 * axis]
        if axis in given:
            t = np.asarray(given[axis], dtype=np.float64)
            if (
                t.shape != (4, 4)
                or not np.all(np.isfinite(t))
                or not np.allclose(t[3], [0, 0, 0, 1], atol=1e-12, rtol=0)
                or not np.allclose(t[:3, :3].T @ t[:3, :3], np.eye(3), atol=1e-12, rtol=0)
                or not np.isclose(np.linalg.det(t[:3, :3]), 1, atol=1e-12, rtol=0)
            ):
                msg = "periodic transforms must be proper rigid homogeneous matrices"
                raise ValueError(msg)
            transforms[axis] = t
    return transforms


def adapt_intersection(  # ruff: ignore[too-many-arguments] - published keyword options.
    seed: Mesh,
    constraints: Mapping[str, Scalar],
    *,
    bounds: Sequence[float],
    constraint_masks: npt.ArrayLike | None = None,
    seed_tolerance: float = 0.01,
    geometry_tolerance: float = 0.01,
    minimum_quality: float = 0.1,
    min_length: float | None = None,
    max_length: float | None = None,
    gradation: float = 0.3,
    adapt_iterations: int = 4,
    quality_guard: float = 0.3,
    adapt_tolerance: float | None = None,
    optimize_passes: int = 4,
    max_tetrahedra: int = 250_000,
    compile: bool = True,  # ruff: ignore[builtin-argument-shadowing] - public keyword.
    batch_size: int = 4096,
    cancel: _meshers.CancellationToken | None = None,
) -> Mesh:
    """Adapt a nonperiodic Tet4 seed against named physical constraints.

    The seed must pass `seed_tolerance`. The returned mesh must independently
    pass `geometry_tolerance` and `minimum_quality`. Constraint-mask bit j refers
    to the jth supplied field. Explicit masks allow parameter-domain cap labels
    to be replaced by their physical counterparts. No input arrays are modified.

    Returns:
        A physical Mesh with retained constraint names and new vertex masks.

    Raises:
        ValueError: Input arrays, labels, bounds or options are invalid.
        MeshingError: A seed or final geometry, topology or quality gate fails.
        CancelledError: The cancellation token is set.
    """  # ruff: ignore[docstring-extraneous-exception] - includes native errors.
    _check_constraint_names(constraints)
    b = _checked_bounds(bounds)
    _check_quality_options(optimize_passes, 0, minimum_quality, geometry_tolerance)
    if seed.tetrahedra.shape[1:] != (_TET4_NODES,):
        msg = "seed must contain Tet4 cells"
        raise ValueError(msg)
    if any(len(pairs) for pairs in seed.periodic_pairs):
        msg = "seed adaptation requires nonperiodic physical geometry"
        raise ValueError(msg)
    if constraint_masks is None:
        if seed.constraint_masks is None or seed.constraint_names != tuple(constraints):
            msg = "seed constraint names must match the supplied physical fields"
            raise ValueError(msg)
        constraint_masks = seed.constraint_masks
    masks = np.asarray(constraint_masks)
    if (
        masks.shape != (len(seed.points),)
        or masks.dtype.kind not in "iu"
        or np.any(masks < 0)
        or np.any(masks >= 2**_meshers.MAX_CONSTRAINTS)
    ):
        msg = "seed constraint masks must be bounded integer bits, one per vertex"
        raise ValueError(msg)
    cells = np.asarray(seed.tetrahedra)
    if cells.dtype.kind not in "iu" or np.any(cells < 0) or np.any(cells >= len(seed.points)):
        msg = "seed connectivity must contain valid integer node indices"
        raise ValueError(msg)
    start = time.perf_counter()
    fields = [_prepare(f, compile_expression=compile) for f in constraints.values()]
    preparation = time.perf_counter() - start
    _check_cancel(cancel)
    result = _meshers.adapt_intersection(
        fields,
        np.asarray(seed.points, dtype=np.float64),
        np.asarray(cells, dtype=np.int64),
        np.asarray(masks, dtype=np.uint32),
        [b[::2].tolist(), b[1::2].tolist()],
        geometry_tolerance,
        minimum_quality,
        max_tetrahedra,
        optimize_passes,
        (
            geometry_tolerance if min_length is None else float(min_length),
            float(np.min(b[1::2] - b[::2])) / 4 if max_length is None else float(max_length),
            float(gradation),
            0.0,
            int(adapt_iterations),
            float(quality_guard),
            geometry_tolerance if adapt_tolerance is None else float(adapt_tolerance),
        ),
        seed_tolerance,
        batch_size,
        cancel,
    )
    result["diagnostics"].update(
        preparation_seconds=preparation,
        workflow_seconds=time.perf_counter() - start,
        evaluator="compiled"
        if all(isinstance(f, _meshers.CompiledField) for f in fields)
        else "mixed",
        threads=1,
    )
    return Mesh(
        points=result["points"],
        tetrahedra=result["tetrahedra"],
        surface=result["surface"],
        boundary_tags=result["boundary_tags"],
        diagnostics=result["diagnostics"],
        periodic_pairs=tuple(result["periodic_pairs"]),
        constraint_masks=result["constraint_masks"],
        constraint_names=tuple(constraints),
    )


def _mapped_starting_mesh(
    build: Callable[[Sequence[int], float], Mesh],
    requested_cells: tuple[int, ...],
    closed_axes: Sequence[bool],
    minimum_quality: float,
) -> tuple[Mesh, tuple[int, ...], float]:
    """Retry only bounded closed-axis folds and starting-mesh quality failures.

    Returns:
        The accepted seed, its resolved grid and its starting quality requirement.

    Raises:
        MeshingError: The mapped starting mesh still fails its admission checks.
    """
    resolved_cells = requested_cells
    seed_quality = minimum_quality
    for attempt in range(_MAX_SEED_ATTEMPTS):
        try:
            seed = build(resolved_cells, seed_quality)
        except _meshers.MeshingError as error:
            if attempt + 1 == _MAX_SEED_ATTEMPTS:
                raise
            refined = tuple(
                min(_MAX_GRID_CELLS, 2 * n) if closed else n
                for n, closed in zip(requested_cells, closed_axes, strict=True)
            )
            if (
                str(error) == "coordinate map inverted a background tetrahedron"
                and resolved_cells == requested_cells
                and refined != requested_cells
            ):
                resolved_cells = refined
            elif str(error).startswith("minimum MMG quality ") and seed_quality > 0:
                seed_quality = 0.0
            else:
                raise
        else:
            return seed, resolved_cells, seed_quality
    msg = "mapped starting-mesh recovery exhausted"
    raise _meshers.MeshingError(msg)


@dataclass(frozen=True, slots=True, kw_only=True)
class _MappedDomain:
    bounds: Sequence[float]
    cells: int | Sequence[int]
    mapping: Callable | None
    inverse_map: Callable | None
    physical_bounds: Sequence[float] | None
    periodic: Sequence[bool]
    transforms: Mapping[int, Sequence[Sequence[float]]] | None
    closed_axes: Sequence[bool]


class _IntersectionControls(TypedDict):
    optimize_passes: int
    max_tetrahedra: int
    compile: bool
    batch_size: int
    cancel: _meshers.CancellationToken | None


@dataclass(frozen=True, slots=True, kw_only=True)
class _MeshRequirements:
    seed_tolerance: float
    final_tolerance: float | None
    minimum_quality: float
    snap: float


class _AdaptationSizing(TypedDict):
    min_length: float | None
    max_length: float | None
    gradation: float
    adapt_iterations: int
    quality_guard: float
    adapt_tolerance: float | None


def _physical_cap(inverse_map: Callable, axis: int, side: int, value: float) -> Callable:
    def cap(x: Points, y: Points, z: Points) -> Points:
        return (1 if side else -1) * (inverse_map(x, y, z)[axis] - value)

    return cap


def _physical_constraints(
    seed: Mesh,
    constraints: Mapping[str, Scalar],
    domain: _MappedDomain,
    inverse_map: Callable,
) -> tuple[dict[str, Scalar], npt.NDArray[np.uint32]]:
    """Return physical fields and masks for the walls and surviving chart caps.

    Raises:
        MeshingError: The mapped seed lost its labels or chart coordinates.
    """
    if seed.constraint_masks is None or seed.parameters is None:
        msg = "mapped seed lost its constraint labels or parameters"
        raise _meshers.MeshingError(msg)
    masks = seed.constraint_masks.astype(np.uint32, copy=True)
    fields = {}
    for name, field in constraints.items():

        def physical(x: Points, y: Points, z: Points, *, field: Scalar = field) -> Points:
            p = inverse_map(x, y, z)
            return (
                cast("Points", field(*p))
                if callable(field) and not isinstance(field, _meshers.CompiledField)
                else _values(field, np.column_stack(p))
            )

        fields[name] = physical
    for axis, closed in enumerate(domain.closed_axes):
        if closed:
            continue
        for side in (0, 1):
            value = domain.bounds[2 * axis + side]
            masks[seed.parameters[:, axis] == value] |= 1 << len(fields)

            name = f"domain_{axis}_{side}"
            while name in fields:
                name = "_" + name
            fields[name] = _physical_cap(inverse_map, axis, side, value)
    return fields, masks


def _adapt_mapped_intersection(
    constraints: Mapping[str, Scalar],
    domain: _MappedDomain,
    requirements: _MeshRequirements,
    controls: _IntersectionControls,
    sizing: _AdaptationSizing,
) -> Mesh:
    """Return physical adaptation of the named mapped solid.

    Raises:
        ValueError: The inverse map, physical bounds or final tolerance is absent.
    """
    inverse_map = domain.inverse_map
    tolerance = requirements.final_tolerance
    if (
        domain.mapping is None
        or inverse_map is None
        or domain.physical_bounds is None
        or tolerance is None
    ):
        msg = "mapped adaptation requires both maps, physical bounds and a final tolerance"
        raise ValueError(msg)
    if any(p and not c for p, c in zip(domain.periodic, domain.closed_axes, strict=True)):
        msg = "mapped adaptation requires nonperiodic physical geometry"
        raise ValueError(msg)
    if len(constraints) + 2 * sum(not c for c in domain.closed_axes) > _meshers.MAX_CONSTRAINTS:
        msg = "mapped adaptation exceeds the named-constraint limit with domain caps"
        raise ValueError(msg)
    start = time.perf_counter()
    cells = domain.cells
    requested_cells = (int(cells),) * 3 if isinstance(cells, (int, np.integer)) else tuple(cells)

    def starting_mesh(seed_cells: Sequence[int], seed_quality: float) -> Mesh:
        return generate_intersection(
            constraints,
            bounds=domain.bounds,
            cells=seed_cells,
            coordinate_map=domain.mapping,
            periodic=domain.periodic,
            periodic_transforms=domain.transforms,
            closed_axes=domain.closed_axes,
            geometry_tolerance=requirements.seed_tolerance,
            minimum_quality=seed_quality,
            snap=requirements.snap,
            **controls,
        )

    seed, resolved_cells, seed_quality = _mapped_starting_mesh(
        starting_mesh,
        requested_cells,
        domain.closed_axes,
        requirements.minimum_quality,
    )
    seed_seconds = time.perf_counter() - start
    fields, masks = _physical_constraints(seed, constraints, domain, inverse_map)
    result = adapt_intersection(
        seed,
        fields,
        bounds=domain.physical_bounds,
        constraint_masks=masks,
        seed_tolerance=requirements.seed_tolerance,
        geometry_tolerance=tolerance,
        minimum_quality=requirements.minimum_quality,
        **controls,
        **sizing,
    )
    result.diagnostics.update(
        seed_seconds=seed_seconds,
        seed_tetrahedra=len(seed.tetrahedra),
        seed_nodes=len(seed.points),
        seed_geometry_tolerance=requirements.seed_tolerance,
        seed_minimum_mmg_quality=seed.diagnostics["minimum_mmg_quality"],
        requested_seed_cells=tuple(int(n) for n in requested_cells),
        seed_cells=tuple(int(n) for n in resolved_cells),
        seed_quality_requirement=float(seed_quality),
        full_workflow_seconds=time.perf_counter() - start,
    )
    return result
