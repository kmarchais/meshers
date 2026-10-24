"""Joint derivatives and a known repair must work before the spinodoid search."""

import numpy as np

import microgen_workflow as workflow
from cap_coupled_patch import CoupledPatch, search
from cap_quality_bound import caps


def cube(distort=False):
    _, _, meshers, *_ = workflow.imports()
    mesh = meshers.generate(lambda x, *_: x * 0 - 1, bounds=(-.5, .5) * 3, cells=4,
        periodic=(True, True, True), threads=1)
    points = mesh.points.copy()
    orbits = [np.flatnonzero((abs(points[:, 0]) == .5) & (abs(points[:, 1] - y) < 1e-12)
        & (abs(points[:, 2]) < 1e-12)) for y in (0., .25)]
    assert all(len(orbit) == 2 for orbit in orbits)
    if distort:
        points[orbits[0], 1] += .23
    faces, _ = caps({"points": points, "surface": mesh.surface})
    return points, mesh.tetrahedra, faces, orbits


def test_joint_derivatives_include_cross_terms_and_periodic_copies():
    points, cells, faces, orbits = cube()
    model = CoupledPatch(points, cells, faces, orbits, 0)
    query = np.array([.01, -.02, -.015, .005])
    values, gradient = model.evaluate(query)
    for axis in range(4):
        step = np.eye(4)[axis] * 1e-6
        numeric = (model.evaluate(query + step)[0] - model.evaluate(query - step)[0]) / 2e-6
        np.testing.assert_allclose(gradient[:, axis], numeric, atol=1e-9, rtol=1e-7)
    other = model.evaluate(query + np.array([.01, 0., 0., -.01]))[1]
    assert np.linalg.norm(other - gradient) > 1e-3
    moved = model.moved(query)
    for orbit in orbits:
        np.testing.assert_allclose(moved[orbit[1]] - moved[orbit[0]], points[orbit[1]] - points[orbit[0]], atol=1e-15)


def test_known_distorted_periodic_patch_can_be_repaired_without_source_mutation():
    points, cells, faces, orbits = cube(distort=True)
    original = points.copy()
    model = CoupledPatch(points, cells, faces, orbits, 0, target=.3)
    assert not model.check(np.zeros(4))["feasible"]
    chosen, attempts = search(model)
    assert chosen is not None and attempts[-1]["feasible"]
    assert model.check(chosen)["minimum_cap_quality"] >= .3
    np.testing.assert_array_equal(points, original)
