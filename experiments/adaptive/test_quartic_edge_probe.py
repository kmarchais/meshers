"""Algebraic controls for the experimental quartic edge modes."""

import numpy as np

from quartic_edge_probe import EDGES, edge_basis, strain_basis


def test_quartic_derivatives_match_centered_differences():
    barycentric = np.random.default_rng(20261020).dirichlet(np.ones(4), 32)
    _, derivatives = edge_basis(barycentric)
    for axis in range(4):
        step = np.eye(4)[axis] * 1e-6
        numeric = (edge_basis(barycentric + step)[0] - edge_basis(barycentric - step)[0]) / 2e-6
        np.testing.assert_allclose(derivatives[:, :, axis], numeric, atol=1e-9, rtol=1e-8)


def test_edge_values_vanish_on_other_faces_and_at_vertices():
    np.testing.assert_array_equal(edge_basis(np.eye(4))[0], np.zeros((4, 6)))
    barycentric = np.random.default_rng(20261020).dirichlet(np.ones(4), 32)
    for edge, (first, second) in enumerate(EDGES):
        for endpoint in (first, second):
            face = barycentric.copy()
            face[:, endpoint] = 0
            np.testing.assert_array_equal(edge_basis(face)[0][:, edge], np.zeros(len(face)))
    swapped = barycentric[:, [1, 0, 2, 3]]
    np.testing.assert_allclose(edge_basis(barycentric)[0][:, 0], edge_basis(swapped)[0][:, 0], atol=1e-14)


def test_engineering_strain_basis_matches_symmetric_gradient():
    gradients = np.random.default_rng(20261020).normal(size=(4, 6, 3))
    coefficients = np.random.default_rng(20261021).normal(size=(4, 6, 3))
    actual = np.einsum("teai,tei->ta", strain_basis(gradients), coefficients)
    gradient = np.einsum("tei,tej->tij", coefficients, gradients)
    expected = np.stack([gradient[:, 0, 0], gradient[:, 1, 1], gradient[:, 2, 2],
        gradient[:, 0, 1] + gradient[:, 1, 0], gradient[:, 0, 2] + gradient[:, 2, 0],
        gradient[:, 1, 2] + gradient[:, 2, 1]], axis=1)
    np.testing.assert_allclose(actual, expected, atol=1e-14)
