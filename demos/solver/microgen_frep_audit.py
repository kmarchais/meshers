# Copyright (c) 2026 Kevin Marchais

"""Independent checks on Tet4 cells and field data read back by VTK."""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np
import pyvista as pv
from microgen import Cylinder
from periodic_bench import (
    Coordinates,
    Indices,
    gradient,
    tet_shapes,
    wall_reference,
    wall_to_mesh,
)
from periodic_geometry import _face_connectivity  # ruff: ignore[import-private-name] - same oriented Tet4 audit.
from scipy.ndimage import label
from scipy.spatial import ConvexHull

if TYPE_CHECKING:
    from microgen_frep_cases import Case

TOLERANCE = 0.01
MANIFOLD_INCIDENCE = 2
CAP_EPSILON = 1e-12
QUALITY = 0.10
VOLUME_LIMIT = 0.10
FACE_SAMPLES = np.array(
    [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]],
    dtype=float,
)
FACE_SAMPLES /= FACE_SAMPLES.sum(axis=1, keepdims=True)


def _triangles(faces: Indices) -> set[tuple]:
    return set(map(tuple, np.sort(faces, axis=1).tolist()))


def tet_topology(tets: Indices, faces: Indices, signed: Coordinates) -> dict:
    """Return face-connected volumes and the recovered surface incidence checks."""
    boundary, volumes, conforms = _face_connectivity(tets, 6 * signed)
    edges = np.sort(boundary[:, [(0, 1), (1, 2), (0, 2)]].reshape(-1, 2), axis=1)
    return {
        "boundary_matches_cells": len(faces) == len(boundary)
        and _triangles(boundary) == _triangles(faces),
        "oriented_face_incidence": conforms,
        "boundary_bad_edge_incidence": int(
            np.sum(np.unique(edges, axis=0, return_counts=True)[1] != MANIFOLD_INCIDENCE),
        ),
        "component_volumes": sorted(volumes.tolist(), reverse=True),
    }


def _periodic(grid: pv.UnstructuredGrid, tags: Indices, faces: Indices) -> bool:
    for axis, name in enumerate("xyz"):
        low, high = faces[tags == 1 + 2 * axis], faces[tags == 2 + 2 * axis]
        if not len(low) and not len(high):
            continue
        pairs = grid.field_data.get(f"periodic_pairs_{name}")
        shift = grid.field_data.get(f"periodic_shift_{name}")
        if pairs is None or shift is None or not np.array_equal(shift.reshape(3), np.eye(3)[axis]):
            return False
        image = dict(pairs.tolist())
        if (
            len(image) != len(pairs)
            or set(image) != set(low.ravel())
            or set(image.values()) != set(high.ravel())
            or len(set(image.values())) != len(pairs)
            or not np.allclose(
                grid.points[pairs[:, 1]] - grid.points[pairs[:, 0]],
                shift,
                rtol=0,
                atol=1e-13,
            )
        ):
            return False
        if _triangles(np.vectorize(image.__getitem__)(low)) != _triangles(high):
            return False
    return True


def _bisect(case: Case, lo: Coordinates, hi: Coordinates, f0: Coordinates) -> Coordinates:
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        fm = case.field(mid)
        same = ((fm < 0) == (f0 < 0)) & (fm != 0)
        lo[same], hi[~same] = mid[same], mid[~same]
    return hi


def _wall_distance(
    case: Case,
    samples: Coordinates,
    normal: Coordinates,
    tolerance: float,
    lock: int | None = None,
) -> Coordinates:
    """Bracket a root inside the domain for each sample.

    Returns:
        Upper distance bounds, or infinity where the search found no bracket.
    """
    field = case.field
    start = field(samples)
    distance = np.full(len(start), np.inf)
    distance[start == 0] = 0
    step = gradient(field, samples)
    if lock is not None:
        step[:, lock] = 0
    for direction in (step, normal):
        length = np.linalg.norm(direction, axis=1)
        unit = direction / np.maximum(length[:, None], 1e-300)
        for radius in (tolerance, 4 * tolerance):
            for sign in (-1, 1):
                test = (distance > tolerance) & (length > 0)
                if not np.any(test):
                    continue
                indices = np.flatnonzero(test)
                lo = samples[test]
                hi = np.clip(lo + radius * sign * unit[test], -0.5, 0.5)
                f0, f1 = start[test], field(hi)
                hit = ((f0 < 0) != (f1 < 0)) | (f1 == 0)
                if not np.any(hit):
                    continue
                indices, lo, hi, f0 = indices[hit], lo[hit], hi[hit], f0[hit]
                hi = _bisect(case, lo, hi, f0)
                distance[indices] = np.minimum(
                    distance[indices],
                    np.linalg.norm(hi - samples[indices], axis=1),
                )
    return distance


def planar_surface(
    vertices: Coordinates,
    facets: Indices,
    points: Coordinates,
    wall: Indices,
) -> dict:
    """Return two-way sampled distances for a prescribed triangulated polyhedron."""
    # Include each exact vertex and 16 subdivisions of all triangle edges/interiors.
    weights = np.array([(a, b, 16 - a - b) for a in range(17) for b in range(17 - a)]) / 16
    ref = np.unique(np.einsum("sk,nkd->nsd", weights, vertices[facets]).reshape(-1, 3), axis=0)
    mesh_samples = np.einsum("sk,nkd->nsd", FACE_SAMPLES, points[wall]).reshape(-1, 3)
    forward = wall_to_mesh(mesh_samples, vertices, facets)
    reverse = wall_to_mesh(ref, points, wall)
    return {
        "sampled_wall_forward_max": float(forward.max()),
        "wall_unbracketed": 0,
        "wall_samples": len(forward),
        "sampled_wall_reverse_max": float(reverse.max()),
        "wall_reference_samples": len(reverse),
        "max_vertex_to_mesh_distance": float(wall_to_mesh(vertices, points, wall).max()),
    }


def sphere_directions() -> Coordinates:
    """Return a Fibonacci sphere plus the six axis directions."""
    golden = (1 + np.sqrt(5)) / 2
    indices = np.arange(4096)
    z = 1 - 2 * (indices + 0.5) / len(indices)
    radial, angle = np.sqrt(1 - z * z), (2 * np.pi / golden) * indices
    return np.vstack((
        np.eye(3),
        -np.eye(3),
        np.column_stack((radial * np.cos(angle), radial * np.sin(angle), z)),
    ))


def _csg_spheres(case: Case, points: Coordinates, wall: Indices) -> dict:
    # Both source radii are .23, with centres on x at +/-.12. Their intersection
    # is a circle at x=0. A clipped sphere patch's nearest point is either its
    # unrestricted radial projection or its circular boundary.
    radius, dx = 0.23, 0.12
    circle_radius = np.sqrt(radius**2 - dx**2)
    samples = np.einsum("sk,nkd->nsd", FACE_SAMPLES, points[wall]).reshape(-1, 3)
    forward = np.hypot(samples[:, 0], np.linalg.norm(samples[:, 1:], axis=1) - circle_radius)
    phi = 2 * np.pi * np.arange(256) / 256
    ref = [
        np.column_stack((
            np.zeros_like(phi),
            circle_radius * np.cos(phi),
            circle_radius * np.sin(phi),
        )),
    ]
    for side in (-1, 1):
        center = np.array([side * dx, 0, 0])
        relative = samples - center
        norm = np.linalg.norm(relative, axis=1)
        projection = center + radius * relative / np.maximum(norm[:, None], 1e-300)
        shell = center + radius * sphere_directions()
        # The lens retains the left sphere's right half and the right sphere's
        # left half. Difference retains the left half of each intersecting sphere.
        keep = -side if case.two_sphere_operation == "lens" else -1
        valid = keep * projection[:, 0] >= 0
        forward[valid] = np.minimum(forward[valid], np.abs(norm[valid] - radius))
        ref.append(shell[keep * shell[:, 0] >= 0])
    reverse = wall_to_mesh(np.vstack(ref), points, wall)
    return {
        "wall_reference": "analytic clipped spheres and intersection circle",
        "sampled_wall_forward_max": float(forward.max()),
        "wall_unbracketed": 0,
        "wall_samples": len(forward),
        "sampled_wall_reverse_max": float(reverse.max()),
        "wall_reference_samples": len(reverse),
    }


def _cylinder_surface(shape: Cylinder, points: Coordinates, wall: Indices) -> dict:
    half, radius = shape.height / 2, shape.radius
    samples = np.einsum("sk,nkd->nsd", FACE_SAMPLES, points[wall]).reshape(-1, 3)
    local = shape.orientation.inv().apply(samples - shape.center)
    axial = abs(local[:, 0]) - half
    radial = np.linalg.norm(local[:, 1:], axis=1) - radius
    side = np.hypot(radial, np.maximum(0, axial))
    end = np.hypot(np.maximum(0, radial), axial)
    forward = np.minimum(side, end)
    angle = 2 * np.pi * np.arange(256) / 256
    unit = np.array((np.cos(angle), np.sin(angle))).T

    def ring(x: float, r: float) -> Coordinates:
        return (
            shape.orientation.apply(np.column_stack((np.full(len(unit), x), r * unit)))
            + shape.center
        )

    rims = np.vstack((ring(-half, radius), ring(half, radius)))
    ref = np.vstack([
        rims,
        *(ring(x, radius) for x in np.linspace(-half, half, 33)[1:-1]),
        *(ring(x, r) for x in (-half, half) for r in np.linspace(0, radius, 17)[:-1]),
    ])
    reverse = wall_to_mesh(ref, points, wall)
    return {
        "wall_reference": "analytic cylinder side, ends and both rims",
        "sampled_wall_forward_max": float(forward.max()),
        "wall_unbracketed": 0,
        "wall_samples": len(forward),
        "sampled_wall_reverse_max": float(reverse.max()),
        "wall_reference_samples": len(reverse),
        "max_rim_to_mesh_distance": float(reverse[: len(rims)].max()),
    }


def _surface(
    case: Case,
    points: Coordinates,
    faces: Indices,
    tags: Indices,
    tolerance: float,
) -> dict:
    wall = faces[tags == 0]
    if isinstance(case.shape, Cylinder) and len(wall):
        return _cylinder_surface(case.shape, points, wall)
    if case.two_sphere_operation is not None and len(wall):
        return _csg_spheres(case, points, wall)
    if case.hull_vertices is not None and len(wall):
        facets = ConvexHull(case.hull_vertices).simplices
        return {
            "wall_reference": "analytic convex hull",
            **planar_surface(case.hull_vertices, facets, points, wall),
        }
    ref = wall_reference(case.field, grid=40)
    if not len(wall) or not len(ref):
        msg = "No interior wall or no independently sampled sign changes"
        raise ValueError(msg)
    forward = []
    for chunk in np.array_split(wall, max(1, len(wall) // 256)):
        p = points[chunk]
        normal = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        samples = np.einsum("sk,nkd->nsd", FACE_SAMPLES, p).reshape(-1, 3)
        forward.extend(
            _wall_distance(case, samples, np.repeat(normal, len(FACE_SAMPLES), axis=0), tolerance),
        )
    measured = np.asarray(forward)
    return {
        "sampled_wall_forward_max": float(measured.max()) if np.isfinite(measured).all() else None,
        "wall_unbracketed": int((~np.isfinite(measured)).sum()),
        "wall_samples": len(measured),
        "sampled_wall_reverse_max": float(wall_to_mesh(ref, points, wall).max()),
        "wall_reference_samples": len(ref),
    }


def _caps(case: Case, points: Coordinates, faces: Indices, tags: Indices, tolerance: float) -> dict:
    coord = (np.arange(40) + 0.5) / 40 - 0.5
    reverse_max, overshoot_max, mislabeled = 0.0, 0.0, 0
    for axis, side in itertools.product(range(3), range(2)):
        tag = 1 + 2 * axis + side
        cap = faces[tags == tag]
        other = [i for i in range(3) if i != axis]
        probes = np.zeros((len(coord) ** 2, 3))
        probes[:, other] = np.reshape(np.meshgrid(coord, coord), (2, -1)).T
        probes[:, axis] = side - 0.5
        # Some Microgen solids include the cell box in the field itself.
        inside = probes.copy()
        inside[:, axis] *= 1 - 1e-9
        ref = probes[case.field(inside) < 0]
        if len(ref):
            gap = wall_to_mesh(ref, points, cap).max() if len(cap) else np.inf
            reverse_max = max(reverse_max, float(gap))
        if len(cap):
            p = points[cap]
            mislabeled += int(
                np.sum(np.max(np.abs(p[:, :, axis] - (side - 0.5)), axis=1) > CAP_EPSILON),
            )
            samples = np.einsum("sk,nkd->nsd", FACE_SAMPLES, p).reshape(-1, 3)
            samples[:, axis] *= 1 - 1e-9
            outside = samples[case.field(samples) > 0]
            if len(outside):
                distance = _wall_distance(
                    case,
                    outside,
                    np.zeros_like(outside),
                    tolerance,
                    lock=axis,
                )
                overshoot_max = max(overshoot_max, float(distance.max()))
    return {
        "sampled_cap_reverse_max": reverse_max if np.isfinite(reverse_max) else None,
        "sampled_cap_overshoot_max": overshoot_max if np.isfinite(overshoot_max) else None,
        "mislabeled_cap_faces": mislabeled,
    }


def _voxel_components(case: Case, n: int) -> list[float]:
    axis = (np.arange(n) + 0.5) / n - 0.5
    points = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    inside = np.concatenate([case.field(chunk) < 0 for chunk in np.array_split(points, n)])
    labels, _ = label(inside.reshape(n, n, n))
    return sorted((np.bincount(labels.ravel())[1:] / n**3).tolist(), reverse=True)


def _volume_reference(case: Case) -> dict:
    coarse = _voxel_components(case, 64)
    fine = _voxel_components(case, 96)
    if case.expected_volume is not None:
        volume = case.expected_volume
        change = 0.0
    else:
        volume = sum(fine)
        change = abs(sum(coarse) / volume - 1)
    count = case.expected_components
    if count is None and len(coarse) == len(fine):
        count = len(fine)
    return {
        "reference_kind": "analytic" if case.expected_volume is not None else "midpoint64/96",
        "reference_volume": volume,
        "volume_grid_relative_change": change,
        "voxel_component_volumes_64": coarse,
        "voxel_component_volumes_96": fine,
        "expected_components": count,
        "component_reference_kind": "analytic"
        if case.expected_components is not None
        else "sampled",
    }


def audit(
    case: Case,
    grid: pv.UnstructuredGrid,
    *,
    tolerance: float = TOLERANCE,
    volume_limit: float = VOLUME_LIMIT,
) -> tuple[dict, list[str]]:
    """Measure orientation, connectivity, bidirectional wall/cap distance, and volume.

    Returns:
        Measurements on the exported mesh and names of missed limits.
    """
    tets = grid.cells_dict[pv.CellType.TETRA]
    faces = grid.field_data["meshers_boundary_faces"]
    tags = grid.field_data["meshers_boundary_tags"]
    signed, quality, angles = tet_shapes(grid.points, tets)
    reference = _volume_reference(case)
    topo = tet_topology(tets, faces, signed)
    measured = {
        "nodes": grid.n_points,
        "tetrahedra": len(tets),
        "inverted": int(np.sum(signed <= 0)),
        "minimum_quality": float(quality.min()),
        "quality_p05_median": np.percentile(quality, [5, 50]).tolist(),
        "dihedral_range": [float(angles.min()), float(angles.max())],
        "volume": float(signed.sum()),
        **reference,
        "relative_volume_error": float(signed.sum() / reference["reference_volume"] - 1),
        "tags": dict(zip(*[v.tolist() for v in np.unique(tags, return_counts=True)], strict=True)),
        "periodic_metadata_valid": _periodic(grid, tags, faces) if case.periodic else None,
        **topo,
        **_surface(case, grid.points, faces, tags, tolerance),
        **_caps(case, grid.points, faces, tags, tolerance),
    }
    limits = {
        "volume": abs(measured["relative_volume_error"]) + measured["volume_grid_relative_change"]
        <= volume_limit,
        "positive": measured["inverted"] == 0,
        "quality": measured["minimum_quality"] >= QUALITY,
        "component_count": reference["expected_components"]
        in {None, len(topo["component_volumes"])},
        "conformity": topo["oriented_face_incidence"] and not topo["boundary_bad_edge_incidence"],
        "labels": topo["boundary_matches_cells"]
        and not measured["mislabeled_cap_faces"]
        and set(tags).issubset(range(7)),
        "periodicity": measured["periodic_metadata_valid"] is not False,
    }
    if case.component_targets is not None:
        actual = np.asarray(measured["component_volumes"])
        error = (
            float(np.max(abs(actual / case.component_targets - 1)))
            if len(actual) == len(case.component_targets)
            else None
        )
        measured["analytic_component_volumes"] = case.component_targets.tolist()
        measured["max_relative_component_volume_error"] = error
        limits["component_volumes"] = error is not None and error <= volume_limit
    for key in (
        "sampled_wall_forward_max",
        "sampled_wall_reverse_max",
        "sampled_cap_reverse_max",
        "sampled_cap_overshoot_max",
    ):
        limits[key] = measured[key] is not None and measured[key] <= tolerance
    return measured, [name for name, passed in limits.items() if not passed]
