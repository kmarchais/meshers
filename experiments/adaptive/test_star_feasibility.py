"""Analytic derivatives, concavity, and incompatible fixed-neighbor controls."""

import numpy as np

from cap_quality_bound import bounds
from mmg_quality_localization import quality
from star_feasibility import StarConstraints, checked_upper_bound, diagnose


def paired_tets():
    triangle = np.array([[-.4, -.8 * np.sqrt(3) / 6], [.4, -.8 * np.sqrt(3) / 6], [0, .8 * np.sqrt(3) / 3]])
    points = np.vstack((np.zeros(3), np.column_stack((triangle, np.full(3, -.4))),
        np.column_stack((triangle, np.full(3, .4)))))
    cells = np.array([[0, 1, 2, 3], [0, 4, 5, 6]])
    negative = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]]) < 0
    cells[negative, 1], cells[negative, 2] = cells[negative, 2].copy(), cells[negative, 1].copy()
    return points, cells


def test_analytic_gradient_and_global_concave_tangent():
    points, cells = paired_tets()
    model = StarConstraints(points, cells, 0, .99)
    position = np.array([.03, -.02, .015])
    values, gradients = model.evaluate(position)
    for axis in range(3):
        step = np.eye(3)[axis] * 1e-6
        numeric = (model.evaluate(position + step)[0] - model.evaluate(position - step)[0]) / 2e-6
        np.testing.assert_allclose(gradients[:, axis], numeric, atol=1e-9, rtol=1e-8)
    for query in np.random.default_rng(20261020).uniform(model.lower, model.upper, (256, 3)):
        tangent = values + gradients @ (query - position)
        assert np.all(model.evaluate(query)[0] <= tangent + 1e-12)


def test_joint_star_can_be_obstructed_when_each_face_is_good():
    points, cells = paired_tets()
    assert bounds(points[cells[:, 1:]])[0].min() > .99
    result = diagnose(points, cells, 0, .99)
    assert result["fixed_neighbor_obstruction"] and not result["positive_candidate_at_target"]
    model = StarConstraints(points, cells, 0, .99)
    bound = checked_upper_bound(model, np.asarray(result["normalized_expansion_point"]), result["weights"])
    assert bound < -1e-5


def test_feasible_star_has_no_false_obstruction():
    points, cells = paired_tets()
    result = diagnose(points, cells, 0, .5)
    assert result["positive_candidate_at_target"] and not result["fixed_neighbor_obstruction"]
    assert result["checked_tangent_upper_bound"] >= 0


def test_concave_constraints_match_independent_mmg_quality():
    points, cells = paired_tets()
    target = .7
    model = StarConstraints(points, cells, 0, target)
    for position in np.random.default_rng(20261020).uniform(-.1, .1, (32, 3)):
        moved = points.copy()
        moved[0] = model.origin + model.scale * position
        determinants = np.linalg.det(moved[cells[:, 1:]] - moved[cells[:, :1]]) / model.scale**3
        signed_quality = quality(moved, cells) * np.sign(determinants)
        expected = determinants / model.squared**1.5 * (1 - target / signed_quality)
        np.testing.assert_allclose(model.evaluate(position)[0], expected, atol=1e-13, rtol=1e-12)
