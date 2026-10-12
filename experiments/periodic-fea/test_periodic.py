"""Checks for failure modes that would invalidate a periodic comparison."""

import numpy as np
import pytest
from periodic import periodic_pairs, tensor_difference


def cube():
    from microgen_benchmark.geometry import boundary

    points = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ],
        float,
    )
    tets = np.array(
        [[0, 1, 3, 7], [0, 3, 2, 7], [0, 2, 6, 7], [0, 6, 4, 7], [0, 4, 5, 7], [0, 5, 1, 7]]
    )
    return points, boundary(points, tets)


def test_periodic_corners_form_one_class():
    points, surface = cube()
    pairs, classes, diagnostics = periodic_pairs(points, surface)
    assert [len(pair) for pair in pairs] == [4, 4, 4]
    assert len(np.unique(classes)) == 1
    assert max(diagnostics["periodic_coordinate_errors"]) == 0


def test_changed_opposite_diagonal_is_rejected():
    points, surface = cube()
    high = np.all(points[surface, 0] == 1, axis=1)
    changed = np.vstack([surface[~high], [[1, 3, 5], [3, 7, 5]]])
    with pytest.raises(ValueError, match="triangulations differ"):
        periodic_pairs(points, changed)


def test_unmatched_face_node_is_rejected():
    points, surface = cube()
    points[7, 1] -= 1e-5
    with pytest.raises(ValueError):
        periodic_pairs(points, surface)


def test_tensor_score_uses_coupled_strain_combinations():
    reference = np.diag([3, 3, 3, 1, 1, 1])
    candidate = reference.copy().astype(float)
    candidate[0, 3] = candidate[3, 0] = 0.1
    result = tensor_difference(candidate, reference)
    assert result["worst_relative_strain_energy_difference"] == pytest.approx(0.1 / np.sqrt(3))
    assert tensor_difference(reference, reference)["worst_relative_strain_energy_difference"] == 0


def test_tensor_score_rejects_singular_reference():
    with pytest.raises(ValueError, match="positive definite"):
        tensor_difference(np.eye(6), np.zeros((6, 6)))
