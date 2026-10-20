"""Whole-mesh repair cannot move boundary vertices or accept new bad tets."""

import itertools

import numpy as np

from interior_relocation_trial import repair
from test_interior_star_probe import star


def test_complete_interior_star_repairs_with_fixed_boundary():
    points, cells = star()
    moved, report = repair(points, cells, cells[:, 1:], maximum_attempts=4)
    assert report["initial_bad_elements"] > 0 and report["remaining_bad_elements"] == 0
    assert report["accepted_moves"] == 1
    np.testing.assert_array_equal(moved[1:], points[1:])


def test_all_boundary_bad_tet_stays_rejected():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [.5, .5, .001]])
    cells = np.arange(4)[None]
    surface = np.array(list(itertools.combinations(range(4), 3)))
    moved, report = repair(points, cells, surface)
    assert report["remaining_bad_elements"] == 1 and report["accepted_moves"] == 0
    assert report["attempts"][0]["status"] == "all-boundary"
    np.testing.assert_array_equal(moved, points)
