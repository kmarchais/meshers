"""Small binary signed-face solve with independent integer and volume replay."""

from itertools import combinations
import warnings

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from bistellar_probe import oriented_boundary


def solve_cavity(points, original, target, *, maximum_vertices=7, minimum_relative_determinant=0):
    vertices = np.unique(original)
    if maximum_vertices not in (7, 9):
        raise ValueError("Only seven- or nine-vertex stages are supported")
    if len(vertices) > maximum_vertices:
        raise ValueError("Expanded cavity exceeds seven-vertex limit")
    candidates = np.array(list(combinations(map(int, vertices), 4)))
    old_volume = float(np.linalg.det(points[original[:, 1:]] - points[original[:, :1]]).sum())
    assert old_volume > 0
    determinants = np.linalg.det(points[candidates[:, 1:]] - points[candidates[:, :1]])
    retained = abs(determinants) > minimum_relative_determinant * old_volume
    candidates = candidates[retained]
    determinants = determinants[retained]
    negative = determinants < 0
    candidates[negative, 1], candidates[negative, 2] = candidates[negative, 2].copy(), candidates[negative, 1].copy()
    assert len(candidates) <= 126
    boundaries = [oriented_boundary(cell[None]) for cell in candidates]
    keys = sorted(set(target).union(*(set(boundary) for boundary in boundaries)))
    matrix = np.array([[boundary.get(key, 0) for boundary in boundaries] for key in keys], dtype=np.int16)
    desired = np.array([target.get(key, 0) for key in keys])
    limits = np.array([1 if key in target else 2 for key in keys])
    with warnings.catch_warnings(record=True) as caught:
        result = milp(np.ones(len(candidates)), integrality=np.ones(len(candidates)),
            bounds=Bounds(0, 1), constraints=[LinearConstraint(matrix, desired, desired),
                LinearConstraint(abs(matrix), 0, limits)],
            options={"time_limit": 2, "node_limit": 1000, "threads": 1, "mip_rel_gap": 0})
    report = {"vertices": vertices.tolist(), "candidate_tetrahedra": len(candidates),
        "solver_status": int(result.status), "solver_message": result.message,
        "solver_warnings": [str(item.message) for item in caught], "passed": False}
    if result.x is None:
        return None, report
    rounded = np.rint(result.x).astype(np.int16)
    if np.max(abs(result.x - rounded)) > 1e-8 or not np.array_equal(matrix @ rounded, desired) or np.any(abs(matrix) @ rounded > limits):
        return None, dict(report, replay_rejected="integer-face-equations")
    selected = candidates[rounded.astype(bool)]
    new_determinants = np.linalg.det(points[selected[:, 1:]] - points[selected[:, :1]])
    assert old_volume > 0
    difference = float(abs(new_determinants.sum() - old_volume) / old_volume)
    if difference >= 1e-10 or not (new_determinants > 0).all() or oriented_boundary(selected) != target:
        return None, dict(report, replay_rejected="volume-or-oriented-boundary", relative_volume_difference=difference)
    report.update(passed=True, alternative_tetrahedra=len(selected),
        relative_volume_difference=difference, integer_face_equations_exact=True,
        minimum_corner_determinant=float(new_determinants.min()))
    return selected, report
