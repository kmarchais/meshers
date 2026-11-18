"""Fixed-boundary Schwarz control from research commit 9d156a29."""

import itertools
import time

import numpy as np
import pyvista as pv
from scipy.optimize import minimize

FACES = np.array(((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)))


def quality(points, cells):
    """Return the normalized MMG quality on each linear tetrahedron."""
    coordinates = points[cells]
    volume = abs(np.linalg.det(coordinates[:, 1:] - coordinates[:, :1])) / 6
    squared = sum(
        ((coordinates[:, first] - coordinates[:, second]) ** 2).sum(axis=1)
        for first, second in itertools.combinations(range(4), 2)
    )
    return (12 * (3 * volume) ** (2 / 3) / squared) ** 1.5


def bounds(triangles):
    """Return the quality ceiling imposed by each fixed opposite face."""
    area = (
        np.linalg.norm(
            np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]), axis=1
        )
        / 2
    )
    squared_edges = sum(
        np.sum((triangles[:, first] - triangles[:, second]) ** 2, axis=1)
        for first, second in ((0, 1), (1, 2), (2, 0))
    )
    return 4 * np.sqrt(3) * area / squared_edges, squared_edges


def relocate_star(points, cells, node, target=0.1):  # ruff: ignore[too-many-locals] - preserve the archived numerical control.
    """Return the candidate vertex and diagnostics for its incident tetrahedra.

    Raises:
        ValueError: The vertex star is empty or degenerate, or the target is invalid.
    """
    selected = np.flatnonzero((cells == node).any(axis=1))
    if not len(selected) or not target > 0:
        raise ValueError("Relocation needs an incident tetrahedron and a positive quality target")
    star = cells[selected]
    positions = np.argmax(star == node, axis=1)
    fixed = np.array([
        cell[np.arange(4) != position] for cell, position in zip(star, positions, strict=True)
    ])
    ceiling, _ = bounds(points[fixed])
    initial = quality(points, star)
    constraint_target = target + 1e-9
    row = {
        "node": int(node),
        "incident_tets": len(star),
        "quality_target": target,
        "initial_minimum_quality": float(initial.min()),
        "fixed_opposite_face_upper_bound": float(ceiling.min()),
        "target_attained": False,
        "status": "fixed-face-obstructed",
        "numerical_quality_guard": 1e-9,
    }
    if ceiling.min() < target:
        return points[node].copy(), row
    origin = points[node].copy()
    scale = float(np.linalg.norm(points[np.unique(star)] - origin, axis=1).max())
    if not scale > 0:
        raise ValueError("Degenerate vertex star")
    coordinates = (points[star] - origin) / scale
    initial_squared = sum(
        np.sum((coordinates[:, first] - coordinates[:, second]) ** 2, axis=1)
        for first in range(4)
        for second in range(first)
    )

    def constraints(displacement):
        moved = coordinates.copy()
        moved[np.arange(len(star)), positions] = displacement
        determinants = np.linalg.det(moved[:, 1:] - moved[:, :1])
        squared = sum(
            np.sum((moved[:, first] - moved[:, second]) ** 2, axis=1)
            for first in range(4)
            for second in range(first)
        )
        return (
            determinants / initial_squared**1.5
            - constraint_target / (12 * np.sqrt(3)) * (squared / initial_squared) ** 1.5
        )

    started = time.monotonic()
    result = minimize(
        lambda displacement: np.dot(displacement, displacement),
        np.zeros(3),
        method="SLSQP",
        constraints={"type": "ineq", "fun": constraints},
        options={"maxiter": 300, "ftol": 1e-13},
    )
    candidate = origin + scale * result.x
    moved = points.copy()
    moved[node] = candidate
    final = quality(moved, star)
    determinants = np.linalg.det(moved[star[:, 1:]] - moved[star[:, :1]])
    attained = bool(determinants.min() > 0 and final.min() >= target)
    row.update(
        status="feasible-numerically" if attained else "not-found-not-impossibility-proof",
        target_attained=attained,
        solver_success=bool(result.success),
        solver_message=str(result.message),
        iterations=int(result.nit),
        seconds=time.monotonic() - started,
        final_minimum_quality=float(final.min()),
        maximum_move=float(np.linalg.norm(candidate - origin)),
        minimum_normalized_constraint=float(constraints(result.x).min()),
        candidate_coordinates=candidate.tolist(),
        minimum_signed_determinant=float(determinants.min()),
    )
    return candidate, row


def audit(points, cells, surface):  # ruff: ignore[too-many-locals] - preserve the archived numerical control.
    """Return the facet incidence and orientation checks for the input mesh.

    Raises:
        ValueError: A tetrahedron has repeated or out-of-range point indices.
    """
    cells = cells[:, :4]
    surface = surface[:, :3]
    if cells.min() < 0 or cells.max() >= len(points):
        raise ValueError("Tetrahedron point index out of range")
    if (np.diff(np.sort(cells, axis=1), axis=1) == 0).any():
        raise ValueError("Tetrahedron contains repeated vertices")
    oriented = cells[:, FACES].reshape(-1, 3)
    inversions = (oriented[:, 0] > oriented[:, 1]).astype(int)
    inversions += oriented[:, 0] > oriented[:, 2]
    inversions += oriented[:, 1] > oriented[:, 2]
    signs = np.where(inversions % 2, -1, 1)
    sorted_faces = np.sort(oriented, axis=1)
    ordering = np.lexsort(sorted_faces[:, ::-1].T)
    ordered = sorted_faces[ordering]
    starts = np.r_[0, np.flatnonzero((ordered[1:] != ordered[:-1]).any(axis=1)) + 1]
    counts = np.diff(np.r_[starts, len(ordered)])
    sums = np.add.reduceat(signs[ordering], starts)
    unique = ordered[starts]
    conflicts = np.flatnonzero((counts == 2) & (sums != 0))
    extracted = {tuple(map(int, face)) for face in unique[counts == 1]}
    provided = {tuple(sorted(map(int, face))) for face in surface}
    examples = []
    for group in conflicts[:10]:
        rows = ordering[starts[group] : starts[group] + counts[group]]
        indices = rows // 4
        face = unique[group]
        triangle = points[face]
        normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
        opposite = [int(next(iter(set(cells[index]) - set(face)))) for index in indices]
        distances = [float(np.dot(points[node] - triangle[0], normal)) for node in opposite]
        examples.append({
            "face": face.tolist(),
            "incident_elements": indices.tolist(),
            "orientation_sum": int(sums[group]),
            "opposite_vertices": opposite,
            "signed_plane_products": distances,
            "same_physical_side": bool(distances[0] * distances[1] > 0),
        })
    determinants = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])
    duplicate_cells = len(cells) - len(np.unique(np.sort(cells, axis=1), axis=0))
    duplicate_surface = len(surface) - len(provided)
    matched = extracted == provided
    coherent = (
        not len(conflicts)
        and not (counts > 2).any()
        and not duplicate_cells
        and not duplicate_surface
        and matched
    )
    return {
        "nodes": len(points),
        "elements": len(cells),
        "unique_faces": len(unique),
        "interior_faces": int((counts == 2).sum()),
        "extracted_boundary_faces": len(extracted),
        "same_orientation_shared_faces": len(conflicts),
        "nonmanifold_faces": int((counts > 2).sum()),
        "duplicate_cells": duplicate_cells,
        "duplicate_surface_faces": duplicate_surface,
        "provided_boundary_matches_extraction": matched,
        "missing_provided_faces": len(extracted - provided),
        "extra_provided_faces": len(provided - extracted),
        "nonpositive_corner_determinants": int((determinants <= 0).sum()),
        "minimum_corner_determinant": float(determinants.min()),
        "oriented_facet_consistency_passed": bool(coherent),
        "conflict_examples": examples,
        "scope": "Oriented corner-face incidence and same-side neighbor checks",
    }


def check_repair(points, moved, cells, boundary):
    """Check that relocation preserves the boundary and material volume.

    Returns:
        Signed volume change and smallest final signed determinant.

    Raises:
        ValueError: A repair reversed an element or changed the volume.
    """
    np.testing.assert_array_equal(moved[boundary], points[boundary])
    old_determinants = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])
    new_determinants = np.linalg.det(moved[cells[:, 1:]] - moved[cells[:, :1]])
    if not np.all(new_determinants > 0):
        raise ValueError("Repair produced a nonpositive or nonfinite determinant")
    volume_difference = float((new_determinants.sum() - old_determinants.sum()) / 6)
    if not abs(volume_difference) < 1e-12:
        raise ValueError(f"Repair changed signed volume by {volume_difference}")
    return volume_difference, float(new_determinants.min())


def repair(points, cells, surface, *, maximum_attempts=400):
    """Return the repaired points and the diagnostics from the two-sweep repair.

    Raises:
        ValueError: The mesh fails a topology, repair progress, orientation or volume check.
    """
    if not audit(points, cells, surface)["oriented_facet_consistency_passed"]:
        raise ValueError("Source fails oriented facet consistency; relocation is unsupported")
    moved = points.copy()
    boundary = np.zeros(len(points), dtype=bool)
    boundary[np.unique(surface)] = True
    attempts, accepted = [], []
    initial = quality(moved, cells)
    measured = initial.copy()
    for sweep in range(2):
        accepted_before = len(accepted)
        for element in np.argsort(measured):
            if measured[element] >= 0.1 or len(attempts) >= maximum_attempts:
                continue
            free = cells[element][~boundary[cells[element]]]
            if not len(free):
                attempts.append({
                    "element": int(element),
                    "sweep": sweep,
                    "status": "all-boundary",
                    "target_attained": False,
                })
            for node in free[: maximum_attempts - len(attempts)]:
                candidate, record = relocate_star(moved, cells, int(node), target=0.100001)
                record.update(element=int(element), sweep=sweep)
                attempts.append(record)
                if record["target_attained"]:
                    bad_before = int((measured < 0.1).sum())
                    moved[node] = candidate
                    measured = quality(moved, cells)
                    if int((measured < 0.1).sum()) >= bad_before:
                        raise ValueError(
                            "Accepted relocation did not reduce the number of poor elements"
                        )
                    accepted.append(int(node))
                    break
        if (
            not (measured < 0.1).any()
            or len(accepted) == accepted_before
            or len(attempts) >= maximum_attempts
        ):
            break
    volume_difference, minimum_determinant = check_repair(points, moved, cells, boundary)
    return moved, {
        "initial_bad_elements": int((initial < 0.1).sum()),
        "remaining_bad_elements": int((measured < 0.1).sum()),
        "initial_minimum_quality": float(initial.min()),
        "final_minimum_quality": float(measured.min()),
        "accepted_moves": len(accepted),
        "moved_vertices": len(set(accepted)),
        "attempts": attempts,
        "maximum_attempts": maximum_attempts,
        "maximum_sweeps": 2,
        "boundary_coordinates_unchanged": True,
        "signed_volume_difference": volume_difference,
        "minimum_signed_determinant": minimum_determinant,
    }


def repaired_grid(grid, row):
    """Return the reoriented volume after the frozen interior relocation control."""
    skin = grid.extract_surface(pass_pointid=True, algorithm="dataset_surface")
    faces = skin.point_data["vtkOriginalPointIds"][skin.regular_faces]
    cells = grid.cells_dict[pv.CellType.TETRA].copy()
    points = np.asarray(grid.points, dtype=float)
    det = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])
    reverse = det < 0
    cells[reverse, :2] = cells[reverse, :2][:, ::-1]
    moved, row["repair"] = repair(points, cells, faces)
    return pv.UnstructuredGrid({pv.CellType.TETRA: cells}, moved)
