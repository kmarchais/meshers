"""Planar cap constraints retain periodic translations and convex upper bounds."""

import numpy as np

import microgen_workflow as workflow
from cap_plane_star import PlaneStar, diagnose
from cap_quality_bound import caps


def cube():
    _, _, meshers, *_ = workflow.imports()
    mesh = meshers.generate(lambda x, *_: x * 0 - 1, bounds=(-.5, .5) * 3, cells=2, threads=1,
        periodic=(True, True, True))
    nodes = np.flatnonzero((abs(mesh.points[:, 0]) == .5) & (abs(mesh.points[:, 1:]).max(axis=1) < 1e-12))
    assert len(nodes) == 2
    faces, _ = caps({"points": mesh.points, "surface": mesh.surface})
    return mesh.points, mesh.tetrahedra, faces, nodes


def test_combined_gradients_and_concave_tangent_are_correct():
    points, cells, faces, nodes = cube()
    model = PlaneStar(points, cells, faces, nodes, 0, .1)
    position = np.array([.01, -.02])
    values, derivative = model.evaluate(position)
    for axis in range(2):
        step = np.eye(2)[axis] * 1e-6
        numeric = (model.evaluate(position + step)[0] - model.evaluate(position - step)[0]) / 2e-6
        np.testing.assert_allclose(derivative[:, axis], numeric, atol=1e-9, rtol=1e-7)
    for query in np.random.default_rng(6).uniform(model.lower, model.upper, (32, 2)):
        assert np.all(model.evaluate(query)[0] <= values + derivative @ (query - position) + 1e-12)


def test_feasible_and_obstructed_periodic_planar_stars():
    points, cells, faces, nodes = cube()
    passed = diagnose(points, cells, faces, nodes, 0, .1)
    assert passed["positive_candidate_at_cap_target"] and not passed["fixed_neighbor_obstruction"]
    moved = np.array(passed["candidate_coordinates"])
    np.testing.assert_allclose(moved[1] - moved[0], points[nodes[1]] - points[nodes[0]], atol=1e-15)
    rejected = diagnose(points, cells, faces, nodes, 0, .99)
    assert rejected["fixed_neighbor_obstruction"] and not rejected["positive_candidate_at_cap_target"]
