"""A nonzero nodal lift residual must enter the full coupled solution."""

import numpy as np

from hp_ranked_solve import lifted_condensation


def test_general_lift_matches_full_augmented_system_and_ignoring_residual_fails():
    rng = np.random.default_rng(53)
    a = rng.normal(size=(9, 9)); a = a.T @ a + np.eye(9)
    coupling = rng.normal(size=(9, 6)) / 10
    nodal_residual = rng.normal(size=(9, 6)) / 10
    mode_residual = rng.normal(size=(6, 6)) / 10
    gram = coupling.T @ np.linalg.solve(a, coupling) + np.eye(6)
    h = 4 * np.eye(6)
    cross = h + mode_residual.T
    p = gram - h + cross + cross.T
    energy = np.block([[h, cross], [cross.T, p]])
    u, z, tensor, _, _ = lifted_condensation(energy, coupling, nodal_residual,
        np.linalg.solve(a, coupling), np.linalg.solve(a, nodal_residual))
    matrix = np.block([[a, coupling], [coupling.T, gram]])
    rhs = np.vstack((nodal_residual, mode_residual))
    direct = -np.linalg.solve(matrix, rhs)
    np.testing.assert_allclose(np.vstack((u, z)), direct, atol=1e-13)
    np.testing.assert_allclose(tensor, h + rhs.T @ direct, atol=1e-13)
    wrong = lifted_condensation(energy, coupling, np.zeros_like(nodal_residual),
        np.linalg.solve(a, coupling), np.zeros_like(nodal_residual))[2]
    assert np.linalg.norm(wrong - tensor) > .01
