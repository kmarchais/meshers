"""Independent polynomial, weak-form, continuity and numbering controls."""

import numpy as np
import pytest

import microgen_workflow as w

w.imports()
from hierarchical_residual_trial import basis, build_modes, integrate, topology, trace_check, trial


def test_cubic_bubbles_vanish_at_old_nodes_and_have_the_reported_derivatives():
    from fedoo.lib_elements.tetrahedron import Tet10

    v, _ = basis(Tet10().xi_nd)
    assert abs(v).max() < 1e-15
    points = np.random.default_rng(17).dirichlet(np.ones(4), size=17)[:, :3]
    _, derivative = basis(points)
    for a in range(3):
        step = np.eye(3)[a]*1e-7
        measured = (basis(points+step)[0]-basis(points-step)[0])/2e-7
        np.testing.assert_allclose(measured, derivative[:, a], atol=3e-10, rtol=1e-7)


@pytest.fixture(scope='module')
def cube():
    _, _, meshers, *_ = w.imports()
    from homogenization import strain_tensor
    from periodic import periodic_pairs
    from quadrature_audit import duffy_rule

    mesh = meshers.generate(lambda x, *_: x*0-1, bounds=(-.5,.5)*3, cells=4,
        periodic=(True,True,True), order=2, threads=1)
    fields = np.array([mesh.points@strain_tensor(load).T for load in np.eye(6)])
    _, classes, _ = periodic_pairs(mesh.points, mesh.surface)
    return mesh, fields, classes, duffy_rule(4)


def test_periodic_affine_patch_has_zero_hierarchical_residual(cube):
    from homogenization import MATERIAL

    mesh, fields, _, rule = cube
    space = topology(mesh.points, mesh.tetrahedra, mesh.surface)
    coefficients, residual, _ = build_modes(mesh.points, mesh.tetrahedra, fields, space, rule)
    assert abs(residual).max() < 1e-12
    energy, _, checked = integrate(mesh.points, mesh.tetrahedra, fields, space, coefficients, rule)
    np.testing.assert_allclose(energy[:6, :6], MATERIAL, atol=1e-12, rtol=1e-12)
    assert abs(checked['volume']-1) < 1e-12
    assert trace_check(space, coefficients)['maximum_trace_jump'] < 1e-12


def test_scattered_residual_matches_direct_cross_energy_and_renumbering(cube):
    mesh, affine, classes, rule = cube
    noise = np.random.default_rng(31).normal(size=(6, len(mesh.points), 3))
    fields = affine+.01*noise[:, classes]
    space = topology(mesh.points, mesh.tetrahedra, mesh.surface)
    coefficients, residual, _ = build_modes(mesh.points, mesh.tetrahedra, fields, space, rule)
    energy, _, _ = integrate(mesh.points, mesh.tetrahedra, fields, space, coefficients, rule)
    cross = np.einsum('fcl,fcm->lm', residual, coefficients)
    np.testing.assert_allclose(energy[:6, 6:], cross, atol=1e-12, rtol=1e-10)
    assert trace_check(space, coefficients)['maximum_trace_jump'] < 1e-12
    optimized = trial(energy)
    assert optimized['maximum_relative_gain'] > 0
    assert optimized['direct_energy_identity_relative_error'] < 1e-10
    order = np.random.default_rng(87).permutation(len(mesh.points))
    inverse = np.empty(len(order), dtype=int); inverse[order] = np.arange(len(order))
    points, cells, surface = mesh.points[order], inverse[mesh.tetrahedra], inverse[mesh.surface]
    other = topology(points, cells, surface)
    c, _, _ = build_modes(points, cells, fields[:, order], other, rule)
    result, _, _ = integrate(points, cells, fields[:, order], other, c, rule)
    np.testing.assert_allclose(result, energy, atol=1e-12, rtol=1e-10)
    assert trace_check(other, c)['maximum_trace_jump'] < 1e-12
