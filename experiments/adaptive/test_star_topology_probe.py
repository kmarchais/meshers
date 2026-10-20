"""Analytic apex and containment checks do not change the cavity domain."""

import numpy as np

from bistellar_probe import oriented_boundary
from cavity_milp import solve_cavity
from mmg_quality_localization import quality
from star_topology_probe import interior_margin, optimal_apex


def test_optimal_apex_uses_the_original_material_side():
    points = np.array([[0, 0, 0], [1, 0, 0], [.5, np.sqrt(3) / 2, 0], [.5, .2, 1.]])
    apex = optimal_apex(points, [0, 1, 2], 3)
    augmented = np.vstack((points, apex))
    assert apex[2] > 0
    np.testing.assert_allclose(quality(augmented, np.array([[0, 1, 2, 4]])), 1, atol=1e-13)
    points[3, 2] = -1
    assert optimal_apex(points, [0, 1, 2], 3)[2] < 0


def test_containment_and_extra_vertex_preserve_original_volume():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1.], [.1, .1, .1]])
    cells = np.arange(4)[None]
    assert interior_margin(points, cells, points[4]) > 0
    assert interior_margin(points, cells, np.array([2., 2., 2.])) < 0
    target = oriented_boundary(cells)
    alternative, row = solve_cavity(points, cells, target, maximum_vertices=16, extra_vertices=[4])
    assert row["passed"] and 4 in row["vertices"]
    assert row["relative_volume_difference"] < 1e-12
    assert oriented_boundary(alternative) == target
