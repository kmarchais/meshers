"""Cross derivatives and coupled trial minimization have independent controls."""

import numpy as np

import microgen_workflow as workflow
from hp_cross_energy import combine, integrate


def test_opposite_periodic_errors_cancel_and_mixed_terms_are_retained():
    rng = np.random.default_rng(14)
    base = rng.normal(size=(6, 6))
    base = base.T @ base + np.eye(6)
    error = rng.normal(size=(6, 6))
    error = error.T @ error + np.eye(6)
    result = combine(np.block([[base + error, base - error], [base - error, base + error]]))
    np.testing.assert_allclose(result["trial_energy_tensor"], base, atol=1e-12)
    np.testing.assert_allclose(result["coefficients"], np.eye(6) / 2, atol=1e-12)


def test_nonidentity_child_maps_reproduce_exact_polynomial_cross_integrals():
    workflow.imports()
    from cubic_reference import register_cubic
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4
    from quadrature_audit import duffy_rule

    p2, p3, linear = Tet10(), register_cubic()(), Tet4()
    maps = np.array([np.eye(4), np.eye(4)])
    maps[0, 1] = [.5, .5, 0., 0.]
    maps[1, 0] = [.5, .5, 0., 0.]
    points = (linear.shape_function(p2.xi_nd)[None] @ maps @ linear.xi_nd).reshape(-1, 3)
    cells = np.arange(20).reshape(2, 10)
    vectors = np.random.default_rng(2).normal(size=(6, 3))
    h = vectors[:, None] * points[None, :, :1]**2
    p = vectors[:, None] * p3.xi_nd[None, :, :1]**3
    matrix, volume, minimum, coupling = integrate(points, cells, h, np.arange(20)[None], p, np.zeros(2, dtype=int), maps,
        duffy_rule(4), return_coupling=True)
    metric = vectors @ np.diag([3., 1., 1.]) @ vectors.T
    expected = np.block([[metric / 15, metric / 20], [metric / 20, metric * 3 / 70]])
    np.testing.assert_allclose(matrix, expected, atol=1e-13, rtol=1e-12)
    np.testing.assert_allclose(volume, 1 / 6, atol=1e-14)
    assert minimum > 0
    np.testing.assert_allclose(h.transpose(0, 2, 1).reshape(6, -1) @ coupling,
        matrix[:6, 6:] - matrix[:6, :6], atol=1e-13)
    np.testing.assert_allclose(coupling.reshape(3, len(points), 6).sum(axis=1), 0., atol=1e-13)
    extended, _, _, mean = integrate(points, cells, h, np.arange(20)[None], p, np.zeros(2, dtype=int), maps,
        duffy_rule(4), additional_fields=2 * h, return_stress=True)
    basis = np.hstack((np.eye(12), np.vstack((2 * np.eye(6), np.zeros((6, 6))))))
    np.testing.assert_allclose(extended, basis.T @ expected @ basis, atol=1e-13)
    packed = np.array([3 * vectors[:, 0], vectors[:, 0], vectors[:, 0], vectors[:, 1], vectors[:, 2], np.zeros(6)])
    np.testing.assert_allclose(mean, np.hstack((packed / 12, packed / 20, packed / 6)), atol=1e-13)
