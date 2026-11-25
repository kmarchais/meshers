# Copyright (c) 2026 Kevin Marchais

"""Check the Tet4 readback against the raw-level solid in its unit cube."""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np
import pyvista as pv
from microgen_spatial_cases import BORE_APOTHEM, BORE_CENTER
from periodic_bench import gradient, mesh_to_wall, tet_shapes, wall_reference, wall_to_mesh
from periodic_geometry import _face_connectivity  # ruff: ignore[import-private-name] - reuse the oriented-face audit.
from scipy.ndimage import label
from scipy.stats import qmc
from vtkmodules.vtkCommonCore import vtkIdList
from vtkmodules.vtkCommonDataModel import vtkStaticCellLocator

if TYPE_CHECKING:
    from collections.abc import Callable

    from microgen_spatial_cases import Specimen
    from numpy.typing import NDArray

from microgen_frep_audit import FACE_SAMPLES, QUALITY
from microgen_frep_audit import TOLERANCE as TOL
from microgen_frep_audit import VOLUME_LIMIT as VOLUME_REL_LIMIT

EPS = 1e-12
CAP_PLANE_ATOL = 1e-10
EDGE_INCIDENCE = 2
BARYCENTRIC_ATOL = 1e-9
HALF = 0.5


def region(case: Specimen, p: NDArray) -> NDArray:
    """Return the x or normalized depth bin at each point."""
    if case.name == "distance-infill":
        return np.minimum(3, (4 * (HALF - np.abs(p).max(axis=1)) / case.depth_max).astype(int))
    return np.minimum(3, (4 * (p[:, 0] + HALF)).astype(int))


def fractions(mask: NDArray, bins: NDArray) -> NDArray:
    """Return the sampled material fraction in each of the four subregions."""
    return np.bincount(bins, weights=mask, minlength=4) / np.bincount(bins, minlength=4)


def reference(case: Specimen, p: NDArray) -> NDArray:
    """Return the phase using the upstream raw world gyroid."""
    xyz = case.world(p)
    f = case.upstream.raw_field(*xyz.T)
    if case.name == "distance-infill":
        # The normalized field cancels between the public upper and lower parts.
        upper = case.upstream.as_upper_skeletal().require_field()(*xyz.T)
        lower = case.upstream.as_lower_skeletal().require_field()(*xyz.T)
        h = (upper + lower) / 2
    else:
        h = case.offset(*p.T) / 2
    phase = np.abs(f) - h if case.part == "sheet" else h - f
    if case.name == "bored-infill":
        # Classification uses the actual upstream polyhedron, independent of the halfspaces.
        # The outer cube is enforced by the meshing bounds; cap samples belong to that cube.
        hole = pv.PolyData(xyz).compute_implicit_distance(case.upstream.obj)["implicit_distance"]
        return np.where(hole > CAP_PLANE_ATOL, 1.0, phase)
    return phase


def material(case: Specimen, grid: pv.UnstructuredGrid) -> dict:
    """Return reference and VTK occupancy at the same reproducible points."""
    result = []
    for seed in (1733, 5077):
        q = qmc.Sobol(3, scramble=True, seed=seed).random_base2(16) - HALF
        bins = region(case, q)
        inside = reference(case, q) < 0
        retained = grid.find_containing_cell(case.world(q)) >= 0
        result.append({
            "seed": seed,
            "samples": len(q),
            "bin_counts": np.bincount(bins).tolist(),
            "reference_volume": float(inside.mean()),
            "mesh_volume_sampled": float(retained.mean()),
            "reference_fractions": fractions(inside, bins).tolist(),
            "mesh_fractions": fractions(retained, bins).tolist(),
            "false_material_fractions": fractions(retained & ~inside, bins).tolist(),
            "missing_material_fractions": fractions(inside & ~retained, bins).tolist(),
        })
    return {"occupancy": result}


def features(case: Specimen, grid: pv.UnstructuredGrid) -> dict:
    """Return missed material and void chord midpoints on 24 off-grid lines."""
    net = qmc.Sobol(2, scramble=True, seed=897).random_base2(3) - HALF
    probes = []
    expected = []
    for axis, ends in itertools.product(range(3), net):
        q = np.empty((513, 3))
        q[:, axis] = np.linspace(-HALF + 1e-6, HALF - 1e-6, len(q))
        q[:, np.arange(3) != axis] = ends
        state = reference(case, q) < 0
        runs = np.split(np.arange(len(q)), np.flatnonzero(np.diff(state)) + 1)
        for run in runs:
            if np.ptp(q[run, axis]) > 2 * TOL:
                i = run[len(run) // 2]
                probes.append(q[i])
                expected.append(state[i])
    q = np.asarray(probes)
    material_here = np.asarray(expected)
    found = grid.find_containing_cell(case.world(q)) >= 0
    return {
        "solid_chords": int(material_here.sum()),
        "void_chords": int((~material_here).sum()),
        "missing_solid_by_bin": np.bincount(
            region(case, q[material_here & ~found]),
            minlength=4,
        ).tolist(),
        "filled_void_by_bin": np.bincount(
            region(case, q[~material_here & found]),
            minlength=4,
        ).tolist(),
    }


def bracket_faces(field: Callable, tri: NDArray) -> NDArray:
    """Return per-triangle upper bounds from explicit brackets on the closed solid."""
    weights = FACE_SAMPLES
    p = np.einsum("sk,nkd->nsd", weights, tri).reshape(-1, 3)
    value = field(p)
    normal = np.repeat(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), len(weights), axis=0)
    directions = [
        normal,
        gradient(field, p),
        *(np.broadcast_to(axis, p.shape) for axis in np.eye(3)),
    ]
    distances = np.where(value == 0, 0.0, np.inf)
    for direction in directions:
        unit = direction / np.maximum(np.linalg.norm(direction, axis=1)[:, None], 1e-300)
        for reach in (TOL / 2, -TOL / 2, TOL, -TOL):
            # A root at an endpoint or between opposite signs is a physical witness.
            ids = np.flatnonzero(distances > TOL)
            if not len(ids):
                break
            lo = p[ids]
            hi = lo + reach * unit[ids]
            end = field(hi)
            crossed = ((value[ids] < 0) != (end < 0)) | (end == 0)
            ids, lo, hi = ids[crossed], lo[crossed], hi[crossed]
            if not len(ids):
                continue
            for _ in range(22):
                mid = (lo + hi) / 2
                fm = field(mid)
                same = ((fm < 0) == (value[ids] < 0)) & (fm != 0)
                lo[same], hi[~same] = mid[same], mid[~same]
            distances[ids] = np.linalg.norm(hi - p[ids], axis=1)
    return distances.reshape(len(tri), len(weights)).max(axis=1)


def fibonacci_directions(count: int) -> NDArray:
    """Return search rays distributed across the sphere, including Cartesian directions."""
    i = np.arange(count)
    z = 1 - (2 * i + 1) / count
    r = np.sqrt(1 - z * z)
    angle = np.pi * (3 - np.sqrt(5)) * i
    return np.vstack((
        np.eye(3),
        -np.eye(3),
        np.column_stack((r * np.cos(angle), r * np.sin(angle), z)),
    ))


def bracket(field: Callable, points: NDArray, radius: float, ray_count: int = 256) -> NDArray:
    """Return upper bounds from brackets strictly inside each physical ball."""
    reach = radius - 4 * EPS
    rays = fibonacci_directions(ray_count)
    delta = np.concatenate([np.zeros((1, 3)), rays * (reach / 2), rays * reach])
    result = np.full(len(points), np.inf)
    for start in range(0, len(points), 32):
        p = points[start : start + 32]
        locations = p[:, None, :] + delta
        values = field(locations.reshape(-1, 3)).reshape(len(p), -1)
        low = values.argmin(axis=1)
        high = values.argmax(axis=1)
        rows = np.arange(len(p))
        hit = (values[rows, low] < -EPS) & (values[rows, high] > EPS)
        if not hit.any():
            continue
        a, b = locations[rows, low][hit], locations[rows, high][hit]
        # Retain strict signs at both endpoints, even near clipped intersections.
        for _ in range(25):
            mid = (a + b) / 2
            val = field(mid)
            a[val < -EPS] = mid[val < -EPS]
            b[val > EPS] = mid[val > EPS]
        result[start + np.flatnonzero(hit)] = EPS + np.maximum(
            np.linalg.norm(a - p[hit], axis=1),
            np.linalg.norm(b - p[hit], axis=1),
        )
    return result


def forward_distances(closed: Callable, points: NDArray, boundary: NDArray) -> dict:
    """Return distances from mesh samples to the closed solid in the rigid cube frame."""
    forward = mesh_to_wall(closed, points, boundary)
    doubtful = (forward > TOL) | ~np.isfinite(forward)
    newton_unresolved = int(np.count_nonzero(~np.isfinite(forward)))
    if np.any(doubtful):
        forward[doubtful] = np.fmin(
            forward[doubtful],
            bracket_faces(closed, points[boundary[doubtful]]),
        )
    # At cap/wall intersections the gradient changes discontinuously. Search
    # for opposite signs inside the physical ball around each doubtful probe.
    suspect = (forward > TOL) | ~np.isfinite(forward)
    if suspect.any():
        triangles = points[boundary[suspect]]
        probes = np.einsum("sk,nkd->nsd", FACE_SAMPLES, triangles)
        unique, inverse = np.unique(probes.reshape(-1, 3), axis=0, return_inverse=True)
        upper = bracket(closed, unique, TOL)
        remaining = ~np.isfinite(upper)
        upper[remaining] = bracket(closed, unique[remaining], TOL, 4096)
        forward[suspect] = upper[inverse].reshape(len(triangles), -1).max(axis=1)
    return {
        "ball_bracket_rechecked_faces": int(suspect.sum()),
        "exterior_control_clear": not np.isfinite(
            bracket(closed, np.array([[0.56, 0.11, 0.21]]), TOL)[0],
        ),
        "newton_unresolved_faces": newton_unresolved,
        "faces_rechecked_by_bracket": int(doubtful.sum()),
        "forward_max": float(forward.max()) if np.isfinite(forward).all() else None,
        "forward_unresolved_faces": int(np.count_nonzero(~np.isfinite(forward))),
    }


def boundary_distances(
    case: Specimen,
    points: NDArray,
    boundary: NDArray,
    faces: NDArray,
    tags: NDArray,
) -> dict:
    """Return the two sampled physical distances and the off-plane cap count."""

    # Check both the requested phase and the envelope when measuring surface distance.
    def phase(p: NDArray) -> NDArray:
        return case.level(*p.T)

    def closed(p: NDArray) -> NDArray:
        return np.maximum(phase(p), np.abs(p).max(axis=1) - HALF)

    walls = wall_reference(phase, grid=40)
    net = qmc.Sobol(2, scramble=True, seed=1749).random_base2(12) - HALF
    reverse = {}
    mislabeled = 0
    for tag in range(7):
        hit = faces[tags == tag]
        samples = walls
        if tag:
            axis = (tag - 1) // 2
            side = -HALF if tag % 2 else HALF
            cap = np.empty((len(net), 3))
            cap[:, axis] = side
            cap[:, np.arange(3) != axis] = net
            samples = cap[phase(cap) < 0]
            mislabeled += int(
                np.count_nonzero(np.max(abs(points[hit, axis] - side), axis=1) > CAP_PLANE_ATOL),
            )
        distance = wall_to_mesh(samples, points, hit) if len(hit) and len(samples) else []
        reverse[tag] = {
            "triangles": len(hit),
            "probes": len(samples),
            "max": float(max(distance)) if len(distance) else None,
        }
    return {
        "bore_footprint": bore_probes(case, points, faces, tags)
        if case.name == "bored-infill"
        else None,
        "cap_off_plane": mislabeled,
        **forward_distances(closed, points, boundary),
        "reverse_by_tag": reverse,
    }


def bore_probes(case: Specimen, points: NDArray, faces: NDArray, tags: NDArray) -> dict:
    """Return reverse distances to the material footprint on each bore face."""
    net = qmc.Sobol(2, scramble=True, seed=4207).random_base2(10) - HALF
    distance = []
    for angle in np.pi / 8 + np.arange(8) * np.pi / 4:
        normal = np.array([np.cos(angle), np.sin(angle)])
        tangent = np.array([-normal[1], normal[0]])
        width = 2 * BORE_APOTHEM * np.tan(np.pi / 8)
        xy = BORE_CENTER + BORE_APOTHEM * normal + (width * net[:, :1]) * tangent
        xyz = np.column_stack((xy, net[:, 1]))
        solid = case.upstream.raw_field(*xyz.T) > case.offset(*xyz.T) / 2
        if solid.any():
            distance.extend(wall_to_mesh(xyz[solid], points, faces[tags == 0]))
    return {"probes": len(distance), "max": max(distance) if distance else None}


def sampled_overlap(grid: pv.UnstructuredGrid, tet: NDArray, signed: NDArray) -> dict:
    """Return counts of sampled tetrahedron centroids strictly inside another Tet4."""
    locator = vtkStaticCellLocator()
    locator.SetDataSet(grid)
    locator.BuildLocator()
    ids = vtkIdList()
    indices = np.linspace(0, len(tet) - 1, min(1024, len(tet)), dtype=int)
    overlaps = []
    for i in indices:
        center = grid.points[tet[i]].mean(axis=0)
        ids.Reset()
        locator.FindCellsWithinBounds(np.repeat(center, 2), ids)
        candidates = np.array([ids.GetId(k) for k in range(ids.GetNumberOfIds())], dtype=int)
        candidates = candidates[(candidates != i) & (signed[candidates] != 0)]
        if candidates.size:
            vertices = grid.points[tet[candidates]]
            systems = (vertices[:, 1:] - vertices[:, :1]).transpose(0, 2, 1)
            w = np.linalg.solve(systems, (center - vertices[:, 0])[:, :, None])[:, :, 0]
            hits = candidates[
                (w.min(axis=1) > BARYCENTRIC_ATOL) & (w.sum(axis=1) < 1 - BARYCENTRIC_ATOL)
            ]
            overlaps.extend((int(i), int(hit)) for hit in hits)
    return {
        "probes": len(indices),
        "pairs": len(overlaps),
        "first_pair": overlaps[0] if overlaps else None,
    }


def voxel_components(case: Specimen, n: int) -> list[float]:
    """Return face-connected material volumes on an independent midpoint grid."""
    axis = (np.arange(n) + 0.5) / n - HALF
    x, y = np.meshgrid(axis, axis, indexing="ij")
    occupied = np.stack([case.level(x, y, z) < 0 for z in axis], axis=-1)
    labels, _ = label(occupied)
    return sorted((np.bincount(labels.ravel())[1:] / n**3).tolist(), reverse=True)


def audit(case: Specimen, grid: pv.UnstructuredGrid) -> dict:
    """Return boundary, signed cell volume, material and void measurements."""
    tet = grid.cells_dict[pv.CellType.TETRA]
    signed, quality, angles = tet_shapes(grid.points, tet)
    boundary, component_volumes, conforms = _face_connectivity(tet, 6 * signed)
    faces = np.asarray(grid.field_data["meshers_boundary_faces"])
    tags = np.asarray(grid.field_data["meshers_boundary_tags"])
    edges = np.sort(boundary[:, [(0, 1), (1, 2), (2, 0)]].reshape(-1, 2), axis=1)
    incidence = np.unique(edges, axis=0, return_counts=True)[1]
    occupancy = material(case, grid)
    return {
        "nodes": grid.n_points,
        "tetrahedra": len(tet),
        "signed_volume": float(signed.sum()),
        "signed_component_volumes": sorted(component_volumes.tolist(), reverse=True),
        "voxel_component_volumes": {str(n): voxel_components(case, n) for n in (64, 128)},
        "sampled_overlap": sampled_overlap(grid, tet, signed),
        "nonpositive_tets": int(np.count_nonzero(signed <= 0)),
        "quality_min_p05_median": np.percentile(quality, [0, 5, 50]).tolist(),
        "dihedral_min_max": [float(angles.min()), float(angles.max())],
        "oriented_face_incidence": conforms,
        "boundary_bad_edges": int(np.count_nonzero(incidence != EDGE_INCIDENCE)),
        "boundary_matches_exposed_faces": len(boundary) == len(faces)
        and set(map(tuple, np.sort(faces, axis=1))) == set(map(tuple, np.sort(boundary, axis=1))),
        "unknown_labels": int(np.count_nonzero(~np.isin(tags, np.arange(7)))),
        **boundary_distances(case, case.local(grid.points), boundary, faces, tags),
        **occupancy,
        **features(case, grid),
    }


def missed_checks(metrics: dict) -> list[str]:
    """Return named misses of the declared geometry screen on a returned mesh."""
    reverse = metrics["reverse_by_tag"]
    signed = metrics["signed_volume"]
    refs = [row["reference_volume"] for row in metrics["occupancy"]]
    spread = abs(refs[0] - refs[1])
    bore = metrics["bore_footprint"]
    counts = [len(volumes) for volumes in metrics["voxel_component_volumes"].values()]
    checks = {
        "sampled_component_count": all(
            n == len(metrics["signed_component_volumes"]) for n in counts
        ),
        "positive_tet4": metrics["nonpositive_tets"] == 0,
        "sampled_overlap": metrics["sampled_overlap"]["pairs"] == 0,
        "minimum_quality": metrics["quality_min_p05_median"][0] >= QUALITY,
        "conformity": metrics["oriented_face_incidence"] and metrics["boundary_bad_edges"] == 0,
        "boundary_labels": metrics["boundary_matches_exposed_faces"]
        and metrics["unknown_labels"] == 0
        and metrics["cap_off_plane"] == 0,
        "forward_distance": metrics["exterior_control_clear"]
        and metrics["forward_max"] is not None
        and metrics["forward_max"] <= TOL,
        "reverse_distance": all(
            row["max"] is not None and row["max"] <= TOL for row in reverse.values()
        ),
        "bore_distance": bore is None or (bore["max"] is not None and bore["max"] <= TOL),
        "volume": max(abs(signed - ref) + spread for ref in refs) <= VOLUME_REL_LIMIT * min(refs),
        "material_and_void_chords": not any(
            metrics["missing_solid_by_bin"] + metrics["filled_void_by_bin"],
        ),
    }
    return [name for name, passed in checks.items() if not passed]
