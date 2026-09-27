"""Join identical periodic unit-cell meshes without repeating optimization."""

from collections.abc import Sequence
from itertools import product
import time

import numpy as np


def tile_periodic(mesh, *, bounds: Sequence[float], repeats: Sequence[int]):
    """Repeat a certified periodic surface or volume cell on an integer lattice.

    The source field and band must be identical in every cell. `bounds` gives
    the source cell's six coordinates; `repeats` gives x, y, z copy counts.
    """
    from . import Mesh, SurfaceMesh

    start = time.perf_counter()

    if not isinstance(mesh, (Mesh, SurfaceMesh)):
        raise TypeError("mesh must be a Mesh or SurfaceMesh")
    if isinstance(mesh, Mesh) and (
        mesh.parameters is not None
        or mesh.constraint_masks is not None
        or mesh.periodic_transforms is not None
    ):
        raise ValueError("mapped and intersection meshes cannot be tiled")
    b = np.asarray(bounds, dtype=np.float64)
    if b.shape != (6,) or not np.all(np.isfinite(b)) or np.any(b[1::2] <= b[::2]):
        raise ValueError("bounds must contain three finite increasing intervals")
    if len(repeats) != 3 or any(
        not isinstance(v, (int, np.integer)) or v < 1 for v in repeats
    ):
        raise ValueError("repeats must contain three positive integers")
    repeats = tuple(int(v) for v in repeats)
    lo = b[::2]
    hi = b[1::2]
    width = hi - lo
    if np.any(mesh.points < lo - 1e-10) or np.any(mesh.points > hi + 1e-10):
        raise ValueError("mesh points lie outside source cell bounds")
    source_keys = np.rint((mesh.points - lo) / width * 1e11).astype(np.int64)
    if len(np.unique(source_keys, axis=0)) != len(mesh.points):
        raise ValueError("distinct source nodes are too close to tile safely")
    labels = mesh.boundary_tags if isinstance(mesh, Mesh) else mesh.labels
    faces = mesh.surface if isinstance(mesh, Mesh) else mesh.triangles
    cap_offset = 1 if isinstance(mesh, Mesh) else 2
    for axis in range(3):
        if repeats[axis] == 1:
            continue
        caps = []
        for side in range(2):
            cap_points = mesh.points[
                faces[labels == cap_offset + 2 * axis + side]
            ].copy()
            cap_points[:, :, axis] = lo[axis]
            keys = np.rint((cap_points - lo) / width * 1e11).astype(np.int64)
            caps.append({tuple(sorted(map(tuple, face))) for face in keys})
        if caps[0] != caps[1]:
            raise ValueError(f"source cell caps do not match on axis {axis}")
        if isinstance(mesh, Mesh) and caps[0] and len(mesh.periodic_pairs[axis]) == 0:
            raise ValueError(f"source cell has no periodic pairs on axis {axis}")
    all_points = []
    all_faces = []
    all_labels = []
    all_tets = []
    for index in product(*(range(count) for count in repeats)):
        offset = len(all_points) * len(mesh.points)
        shift = np.asarray(index) * width
        all_points.append(mesh.points + shift)
        keep = np.ones(len(faces), dtype=bool)
        for axis in range(3):
            if index[axis] > 0:
                keep &= labels != cap_offset + 2 * axis
            if index[axis] + 1 < repeats[axis]:
                keep &= labels != cap_offset + 2 * axis + 1
        all_faces.append(faces[keep] + offset)
        all_labels.append(labels[keep])
        if isinstance(mesh, Mesh):
            all_tets.append(mesh.tetrahedra + offset)
    copied_points = np.concatenate(all_points)
    # Integer local coordinates merge only translated copies of the same node.
    keys = np.rint((copied_points - lo) / width * 1e11).astype(np.int64)
    _, first, mapping = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    points = copied_points[first]
    surface = mapping[np.concatenate(all_faces)]
    boundary_tags = np.concatenate(all_labels)
    if isinstance(mesh, SurfaceMesh):
        used = np.unique(surface)
        compact = np.full(len(points), -1, dtype=np.int64)
        compact[used] = np.arange(len(used))
        return SurfaceMesh(
            points[used],
            compact[surface],
            boundary_tags,
            {
                **mesh.diagnostics,
                "repeated_cells": int(np.prod(repeats)),
                "tiling_seconds": time.perf_counter() - start,
                "seconds": mesh.diagnostics["seconds"] + time.perf_counter() - start,
            },
        )
    tetrahedra = mapping[np.concatenate(all_tets)]
    pairs = []
    global_hi = lo + np.asarray(repeats) * width
    point_keys = np.rint((points - lo) / width * 1e11).astype(np.int64)
    for axis in range(3):
        if len(mesh.periodic_pairs[axis]) == 0:
            pairs.append(np.empty((0, 2), dtype=np.int64))
            continue
        low = np.flatnonzero(np.isclose(points[:, axis], lo[axis], rtol=0, atol=1e-10))
        high = np.flatnonzero(
            np.isclose(points[:, axis], global_hi[axis], rtol=0, atol=1e-10)
        )
        high_map = {tuple(np.delete(point_keys[i], axis)): int(i) for i in high}
        matched = [
            (int(i), high_map[tuple(np.delete(point_keys[i], axis))]) for i in low
        ]
        if len(matched) != len(high):
            raise ValueError(f"tiled periodic nodes do not match on axis {axis}")
        pairs.append(np.asarray(matched, dtype=np.int64).reshape(-1, 2))
    diagnostics = {
        **mesh.diagnostics,
        "tiling_seconds": time.perf_counter() - start,
        "seconds": mesh.diagnostics["seconds"] + time.perf_counter() - start,
        "volume": mesh.diagnostics["volume"] * int(np.prod(repeats)),
        "repeated_cells": int(np.prod(repeats)),
        "elements_below_quality_01": mesh.diagnostics["elements_below_quality_01"]
        * int(np.prod(repeats)),
    }
    return Mesh(
        points,
        tetrahedra,
        surface,
        boundary_tags,
        tuple(pairs),
        diagnostics,
    )
