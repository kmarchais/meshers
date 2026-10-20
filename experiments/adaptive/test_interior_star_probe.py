"""Fixed-neighbor constraints and volume checks for the one-vertex diagnostic."""

import itertools

import numpy as np

from interior_star_probe import relocate_star
from mmg_quality_localization import quality


def star():
    points = np.vstack(([0., 0., .97], np.eye(3), -np.eye(3)))
    cells = np.array([[0, 1 + choices[0] * 3, 2 + choices[1] * 3, 3 + choices[2] * 3]
        for choices in itertools.product((0, 1), repeat=3)])
    negative = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]]) < 0
    cells[negative, 1], cells[negative, 2] = cells[negative, 2].copy(), cells[negative, 1].copy()
    return points, cells


def test_complete_star_reaches_target_without_boundary_motion():
    points, cells = star()
    assert quality(points, cells).min() < .1
    candidate, row = relocate_star(points, cells, 0, target=.2)
    moved = points.copy()
    moved[0] = candidate
    assert row["target_attained"]
    assert quality(moved, cells).min() >= .2 - 1e-10
    np.testing.assert_array_equal(moved[1:], points[1:])
    old_volume = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]]).sum()
    new_volume = np.linalg.det(moved[cells[:, 1:]] - moved[cells[:, :1]]).sum()
    assert abs(new_volume - old_volume) < 1e-12


def test_fixed_opposite_face_obstruction_stops_before_optimization():
    points = np.array([[0., 0., 1.], [0., 0., 0.], [1., 0., 0.], [.5, .001, 0.]])
    cells = np.array([[0, 1, 2, 3]])
    if np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])[0] < 0:
        cells[:, [1, 2]] = cells[:, [2, 1]]
    candidate, row = relocate_star(points, cells, 0)
    assert row["status"] == "fixed-face-obstructed" and not row["target_attained"]
    np.testing.assert_array_equal(candidate, points[0])
