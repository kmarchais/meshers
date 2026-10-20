"""Enumeration checks every fixed-vertex simplex subset in one tiny cavity."""

import numpy as np

from bistellar_probe import oriented_boundary
from cavity_triangulation import enumerate_cavity, replacement_boundary


def saved_nonconical_star():
    points = np.array([[-.13333333333333336, -.5, -.4], [-.13333333333333336, -.4666666666666667, -.4],
        [-.13333333333333336, -.5, -.3666666666666667], [-.16334780789045183, -.5, -.4],
        [-.15617647183609332, -.48950980516942666, -.4], [-.164405327992552, -.5, -.39773866132588537]])
    cells = np.array([[2, 4, 3, 5], [1, 3, 2, 4], [0, 2, 1, 3]])
    return points, cells


def test_original_cavity_is_reproduced_without_new_points():
    points, cells = saved_nonconical_star()
    target = oriented_boundary(cells)
    alternative, row = enumerate_cavity(points, cells, target)
    assert row["passed"] and oriented_boundary(alternative) == target
    assert set(alternative.ravel()) <= set(cells.ravel())


def test_nonconical_cap_replacement_is_checked_not_assumed():
    points, cells = saved_nonconical_star()
    target = replacement_boundary(points, cells, np.array([[2, 3, 5], [0, 2, 3]]), np.array([[0, 2, 5], [0, 5, 3]]))
    alternative, row = enumerate_cavity(points, cells, target)
    assert row["subsets_checked"] <= 32767
    assert alternative is None and not row["passed"]
    assert row["oriented_boundary_matches"] == 0


def test_changed_boundary_has_a_positive_convex_control():
    points = np.array([[0, 0, -.5], [.4, 0, -.5], [.4, .4, -.5], [0, .4, -.5], [.2, .2, 0]])
    cells = np.array([[0, 1, 2, 4], [0, 2, 3, 4]])
    target = replacement_boundary(points, cells, np.array([[0, 1, 2], [0, 2, 3]]), np.array([[0, 1, 3], [1, 2, 3]]))
    alternative, row = enumerate_cavity(points, cells, target)
    assert row["passed"] and oriented_boundary(alternative) == target


def test_unclosed_boundary_is_rejected():
    points, cells = saved_nonconical_star()
    target = oriented_boundary(cells)
    target.pop(next(iter(target)))
    alternative, row = enumerate_cavity(points, cells, target)
    assert alternative is None and not row["passed"]
