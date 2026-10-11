"""Printing and visualization diagnostics on complete triangle surfaces."""

import numpy as np

from .geometry import geometry_metrics, surface_topology
from .surface_intersections import check


def validate_surface(points, faces, case, max_pairs=2_000_000, target_points=None):
    points, faces = np.asarray(points, float), np.asarray(faces, np.int64)
    if points.ndim != 2 or points.shape[1] != 3 or faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("Expected points (N,3) and triangles (M,3)")
    if (
        not len(faces)
        or not np.isfinite(points).all()
        or faces.min() < 0
        or faces.max() >= len(points)
    ):
        raise ValueError("Empty or malformed triangle surface")
    x = points / case.scale
    if np.any(np.abs(x) > 0.5 + 1e-8):
        raise ValueError("Surface leaves the declared specimen box")
    surface_topology(faces)
    v = x[faces]
    normals = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    twice_area = np.linalg.norm(normals, axis=1)
    if np.any(twice_area <= 1e-15):
        raise ValueError("Degenerate surface triangle")
    intersections = check(x, faces, max_pairs)
    if intersections["intersecting_pairs"]:
        raise ValueError(f"Nonadjacent surface triangles intersect: {intersections['first_pairs']}")
    signed_volume = np.sum(np.einsum("ij,ij->i", v[:, 0], np.cross(v[:, 1], v[:, 2]))) / 6
    if signed_volume <= 0:
        raise ValueError("Nonpositive signed surface volume; expected outward winding")
    edges_squared = sum(
        np.sum((v[:, a] - v[:, b]) ** 2, axis=1) for a, b in [(0, 1), (1, 2), (2, 0)]
    )
    shape = 2 * np.sqrt(3) * twice_area / edges_squared
    return {
        "nodes": len(np.unique(faces)),
        "triangles": len(faces),
        "volume": float(signed_volume * case.scale**3),
        "solid_fraction": float(signed_volume),
        "surface_shape_min": float(shape.min()),
        "surface_shape_p05": float(np.percentile(shape, 5)),
        "geometry_status": "checked-at-declared-tolerances",
        "surface_intersection_scope": "Floating-point nonadjacent triangle test; shared-vertex pairs excluded. Not a global intersection certificate.",
        **intersections,
        **geometry_metrics(points, faces, case, target_points),
    }
