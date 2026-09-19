"""Named implicit intersections with retained curves and mapped periodic ends."""

from __future__ import annotations

import time
import warnings
from collections.abc import Callable, Mapping, Sequence
from itertools import starmap
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from . import _meshers
from .expression import TraceError, compile_field

if TYPE_CHECKING:
    from . import Mesh

Scalar = Callable | _meshers.CompiledField | str
Points = npt.NDArray[np.float64]
Indices = npt.NDArray[np.int64]


def _check_cancel(cancel: _meshers.CancellationToken | None) -> None:
    if cancel is not None and cancel.cancelled:
        raise _meshers.CancelledError("generation cancelled")


def _prepare(field: Scalar, *, compile: bool) -> Scalar:
    if callable(field) and not isinstance(field, _meshers.CompiledField) and compile:
        try:
            return compile_field(field)
        except TraceError as error:
            warnings.warn(f"Using NumPy callbacks: {error}", RuntimeWarning, stacklevel=3)
    return field


def _values(field: Scalar, points: Points) -> Points:
    if isinstance(field, _meshers.CompiledField):
        return np.fromiter(starmap(field, points), dtype=np.float64, count=len(points))
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


def _surface_errors(  # ruff: ignore[too-many-locals] - explicit physical inversion and gradient calculation.
    result: dict,
    fields: list[Scalar],
    mapping: list[Scalar] | None,
    bounds: Points,
    cancel: _meshers.CancellationToken | None,
    *,
    batch_size: int,
    geometry_tolerance: float,
) -> list[float]:
    """Sample physical facets against analytic constraints, including map inversion.

    Returns:
        Largest sampled distance estimate per constraint.

    Raises:
        MeshingError: The map or distance estimate is singular or invalid.
    """
    step = (bounds[1::2] - bounds[::2]) * 1e-5
    points, parameters, faces = result["points"], result["parameters"], result["surface"]
    extent = float(np.max(np.ptp(points, axis=0))) if len(points) else 0.0
    inverse_tolerance = min(geometry_tolerance * 0.01, extent * 1e-12)
    masks = result["constraint_masks"][faces]
    labels = masks[:, 0] & masks[:, 1] & masks[:, 2]
    weights = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [0.5, 0.5, 0],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5],
        [1 / 3, 1 / 3, 1 / 3],
    ])
    errors = []
    for j, field in enumerate(fields):
        selected = faces[(labels & (1 << j)) != 0]
        maximum = 0.0
        # Bound sample arrays, not just triangle counts: each face has seven samples.
        count = len(selected) * len(weights)
        limit = min(batch_size, 512 * len(weights))
        for start in range(0, count, limit):
            _check_cancel(cancel)
            samples = np.arange(start, min(start + limit, count))
            chunk = selected[samples // len(weights)]
            sample_weights = weights[samples % len(weights)]
            physical = np.einsum("sj,sjk->sk", sample_weights, points[chunk])
            p = np.einsum("sj,sjk->sk", sample_weights, parameters[chunk])
            # Resolve the map in physical units, accounting for representable
            # coordinate precision when the mesh is far from the origin.
            roundoff = (
                8 * np.finfo(np.float64).eps * np.maximum(np.max(np.abs(physical), axis=1), extent)
            )
            inverse_limit = np.maximum(inverse_tolerance, roundoff)
            for _ in range(8):
                residual = _mapped(p, mapping) - physical
                if np.all(np.linalg.norm(residual, axis=1) <= inverse_limit):
                    break
                p -= np.linalg.solve(_jacobian(p, mapping, step), residual[..., None])[..., 0]
            inverse_error = np.linalg.norm(_mapped(p, mapping) - physical, axis=1)
            if not np.all(np.isfinite(inverse_error)) or np.any(inverse_error > inverse_limit):
                raise _meshers.MeshingError(
                    "coordinate map could not be inverted for surface checks"
                )
            gradient = np.column_stack([
                (
                    _values(field, p + np.eye(3)[a] * step[a])
                    - _values(field, p - np.eye(3)[a] * step[a])
                )
                / (2 * step[a])
                for a in range(3)
            ])
            physical_gradient = np.linalg.solve(
                np.swapaxes(_jacobian(p, mapping, step), 1, 2), gradient[..., None]
            )[..., 0]
            magnitude = np.linalg.norm(physical_gradient, axis=1)
            value = np.abs(_values(field, p))
            distance = np.divide(
                value, magnitude, out=np.full_like(value, np.inf), where=magnitude > 0
            )
            # Do not hide unresolved physical displacement in the distance gate.
            distance += inverse_error
            if not np.all(np.isfinite(distance)):
                raise _meshers.MeshingError("nonfinite or singular sampled surface distance")
            maximum = max(maximum, float(np.max(distance, initial=0)))
        errors.append(maximum)
    _check_cancel(cancel)
    return errors


def _periodic_metadata(
    points: Points,
    parameters: Points,
    faces: Indices,
    bounds: Points,
    *,
    periodic: Sequence[bool],
    transforms: Points,
) -> tuple[npt.NDArray[np.uint8], tuple[Indices, ...]]:
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
            parameters[a][:, other], parameters[b][:, other], rtol=0, atol=scale * 1e-11
        ):
            raise _meshers.MeshingError("periodic parameter nodes do not match")
        transform = transforms[axis]
        mapped = points[a] @ transform[:3, :3].T + transform[:3, 3]
        # Bound rounding in the rigid transform, including cancellation between
        # large rotated coordinates and translations. Keep the geometry limit.
        roundoff = (
            8
            * np.finfo(np.float64).eps
            * (
                np.abs(points[a]) @ np.abs(transform[:3, :3]).T
                + np.abs(transform[:3, 3])
                + np.abs(points[b])
            )
        )
        if not np.all(np.isfinite(mapped)) or not np.all(
            np.abs(mapped - points[b]) <= scale * 1e-10 + roundoff
        ):
            raise _meshers.MeshingError(
                "physical periodic nodes do not match the supplied transform"
            )
        correspondence = np.arange(len(points))
        correspondence[a] = b
        left = np.sort(correspondence[faces[np.all(low[faces], axis=1)]], axis=1)
        right = np.sort(faces[np.all(high[faces], axis=1)], axis=1)
        if set(map(tuple, left)) != set(map(tuple, right)):
            raise _meshers.MeshingError("periodic end triangles do not match")
        pairs.append(np.column_stack((a, b)))
    return tags, tuple(pairs)


def generate_intersection(
    constraints: Mapping[str, Callable | _meshers.CompiledField | str],
    *,
    bounds: Sequence[float] = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
    cells: int | Sequence[int] = 24,
    coordinate_map: Callable | None = None,
    periodic: Sequence[bool] = (False, False, False),
    periodic_transforms: Mapping[int, Sequence[Sequence[float]]] | None = None,
    geometry_tolerance: float | None = None,
    minimum_quality: float = 0.0,
    optimize_passes: int = 4,
    snap: float = 0.2,
    max_tetrahedra: int = 2_000_000,
    compile: bool = True,
    batch_size: int = 4096,
    cancel: _meshers.CancellationToken | None = None,
) -> Mesh:
    """Mesh the intersection of named negative sublevel sets, preserving their curves.

    Constraints are evaluated in parameter coordinates. A coordinate map returns
    a tuple ``(X(x,y,z), Y(x,y,z), Z(x,y,z))`` and acts on the background before
    cutting; it must preserve orientation and be locally invertible. Straight-sided
    tetrahedra approximate the mapped geometry. This path currently uses one worker.

    Args:
        constraints: Ordered mapping of 2-8 unique nonempty names to fields. Keep
            every field <= 0. Use two constraints for the sides of a TPMS sheet.
        bounds: Parameter box ``(xmin,xmax,ymin,ymax,zmin,zmax)``.
        cells: Background intervals, 4-128 per parameter axis.
        coordinate_map: Optional parameter-to-physical map; default identity.
        periodic: Matching low/high faces along the three parameter axes.
        periodic_transforms: Rigid homogeneous 4x4 matrices keyed by periodic axis
            (0, 1, 2), mapping low to high physical points. Required for each mapped
            periodic axis. Identity maps default to box translations. Node and
            triangle correspondence are checked after optimization.
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
    """
    from . import Mesh

    b, cells, transforms = _validate_options(
        constraints,
        bounds=bounds,
        cells=cells,
        periodic=periodic,
        coordinate_map=coordinate_map,
        periodic_transforms=periodic_transforms,
        optimize_passes=optimize_passes,
        snap=snap,
        minimum_quality=minimum_quality,
        geometry_tolerance=geometry_tolerance,
    )
    start = time.perf_counter()
    fields = [_prepare(f, compile=compile) for f in constraints.values()]
    mapping = None
    if coordinate_map is not None:
        mapping = [
            _prepare(lambda x, y, z, a=a: coordinate_map(x, y, z)[a], compile=compile)
            for a in range(3)
        ]
    preparation = time.perf_counter() - start
    result = _meshers.generate_intersection(
        fields,
        mapping,
        [b[::2].tolist(), b[1::2].tolist()],
        cells,
        tuple(periodic),
        minimum_quality,
        max_tetrahedra,
        optimize_passes,
        snap,
        batch_size,
        cancel,
    )
    _check_cancel(cancel)
    tags, pairs = _periodic_metadata(
        result["points"],
        result["parameters"],
        result["surface"],
        b,
        periodic=periodic,
        transforms=transforms,
    )
    diagnostics = result["diagnostics"]
    diagnostics["preparation_seconds"] = preparation
    diagnostics["evaluator"] = (
        "compiled" if all(isinstance(f, _meshers.CompiledField) for f in fields) else "mixed"
    )
    diagnostics["empty"] = len(result["tetrahedra"]) == 0
    vertices = result["points"][result["tetrahedra"]]
    diagnostics["volume"] = float(np.linalg.det(vertices[:, 1:] - vertices[:, :1]).sum() / 6)
    diagnostics["sampled_surface_error"] = None
    if geometry_tolerance is not None:
        try:
            errors = _surface_errors(
                result,
                fields,
                mapping,
                b,
                cancel,
                batch_size=batch_size,
                geometry_tolerance=geometry_tolerance,
            )
        except np.linalg.LinAlgError as error:
            raise _meshers.MeshingError("singular coordinate map during surface checks") from error
        diagnostics["surface_errors"] = dict(zip(constraints, errors, strict=True))
        diagnostics["sampled_surface_error"] = max(errors)
        if max(errors) > geometry_tolerance:
            raise _meshers.MeshingError(
                f"sampled surface error {max(errors):g} exceeds {geometry_tolerance:g}"
            )
    _check_cancel(cancel)
    return Mesh(
        **result,
        boundary_tags=tags,
        periodic_pairs=pairs,
        periodic_transforms=transforms,
        constraint_names=tuple(constraints),
    )


def _validate_options(  # ruff: ignore[complex-structure, too-many-branches] - explicit independent public option checks.
    constraints: Mapping[str, Scalar],
    *,
    bounds: Sequence[float],
    cells: int | Sequence[int],
    periodic: Sequence[bool],
    coordinate_map: Callable | None,
    periodic_transforms: Mapping[int, Sequence[Sequence[float]]] | None,
    optimize_passes: int,
    snap: float,
    minimum_quality: float,
    geometry_tolerance: float | None,
) -> tuple[Points, Sequence[int], Points]:
    if not isinstance(constraints, Mapping) or not 2 <= len(constraints) <= 8:
        raise ValueError("constraints must be a mapping of 2-8 named fields")
    if any(not isinstance(name, str) or not name.strip() for name in constraints):
        raise ValueError("constraint names must be nonempty strings")
    b = np.asarray(bounds, dtype=np.float64)
    if b.shape != (6,) or not np.all(np.isfinite(b)) or np.any(b[1::2] <= b[::2]):
        raise ValueError("bounds must contain three finite increasing intervals")
    if isinstance(cells, (int, np.integer)):
        cells = (int(cells),) * 3
    if len(cells) != 3 or any(
        not isinstance(v, (int, np.integer)) or not 4 <= v <= 128 for v in cells
    ):
        raise ValueError("cells must contain three integers in 4..128")
    if len(periodic) != 3 or any(not isinstance(v, (bool, np.bool_)) for v in periodic):
        raise ValueError("periodic must contain three booleans")
    if not isinstance(optimize_passes, int) or not 0 <= optimize_passes <= 20:
        raise ValueError("optimize_passes must be an integer in 0..20")
    if not np.isfinite(snap) or not 0 <= snap <= 0.2:
        raise ValueError("snap must be finite and in 0..0.2")
    if not np.isfinite(minimum_quality) or not 0 <= minimum_quality <= 1:
        raise ValueError("minimum_quality must be finite and in 0..1")
    if geometry_tolerance is not None and (
        not np.isfinite(geometry_tolerance) or geometry_tolerance <= 0
    ):
        raise ValueError("geometry_tolerance must be positive and finite, or None")
    transforms = np.tile(np.eye(4), (3, 1, 1))
    given = {} if periodic_transforms is None else dict(periodic_transforms)
    if any(a not in {0, 1, 2} or not periodic[a] for a in given):
        raise ValueError("periodic_transforms keys must be enabled axes 0, 1, 2")
    for axis in range(3):
        if periodic[axis]:
            if coordinate_map is not None and axis not in given:
                raise ValueError("mapped periodic axes require explicit periodic_transforms")
            transforms[axis, axis, 3] = b[2 * axis + 1] - b[2 * axis]
        if axis in given:
            t = np.asarray(given[axis], dtype=np.float64)
            if (
                t.shape != (4, 4)
                or not np.all(np.isfinite(t))
                or not np.allclose(t[3], [0, 0, 0, 1], atol=1e-12, rtol=0)
                or not np.allclose(t[:3, :3].T @ t[:3, :3], np.eye(3), atol=1e-12, rtol=0)
                or not np.isclose(np.linalg.det(t[:3, :3]), 1, atol=1e-12, rtol=0)
            ):
                raise ValueError("periodic transforms must be proper rigid homogeneous matrices")
            transforms[axis] = t
    return b, cells, transforms
