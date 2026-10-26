"""Directionality and periodicity of the geometry-only refinement metric."""

import numpy as np

from periodic_crease_grading import directional_metric
from crease_directional_refinement import edge_metric_squared


def test_curve_metric_refines_transverse_edges_and_keeps_periodic_copies_equal():
    # A curve on x=-.5 whose wall normal is y has tangent z.
    cloud = np.array([[0., .5, .6]])
    def wall(x, y, z):
        return y
    def scalar(x, y, z):
        return np.full_like(x, .01)
    metric = directional_metric(scalar, cloud, wall, maximum=.1, ratio=5.)
    values = metric(np.array([[-.5, 0, .1], [.5, 0, .1]]))
    np.testing.assert_array_equal(values[0], values[1])
    assert np.linalg.eigvalsh(values).min() > 0
    radial, tangent = np.array([.03, 0, 0]), np.array([0, 0, .03])
    assert radial@values[0]@radial > (4/3)**2
    assert tangent@values[0]@tangent < (4/3)**2


def test_fine_curve_endpoint_is_not_hidden_by_a_coarse_edge_midpoint():
    cloud = np.array([[0., .5, .6]])
    def wall(x, y, z):
        return y
    def scalar(x, y, z):
        return np.maximum(.002, .1*np.minimum(abs(x+.5)/(1/26), 1.)**.75)
    metric = directional_metric(scalar, cloud, wall, maximum=.1, ratio=5.)
    points = np.array([[-.5, 0, .1], [-.48, 0, .1], [-.49, 0, .1]])
    delta = points[1]-points[0]
    midpoint_score = delta@metric(points[2:])[0]@delta
    assert midpoint_score < (4/3)**2
    assert edge_metric_squared(points, np.array([[0, 1]]), np.array([2]), metric)[0] > (4/3)**2
