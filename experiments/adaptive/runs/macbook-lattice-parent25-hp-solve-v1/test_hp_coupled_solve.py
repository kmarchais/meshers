"""Small Schur solve must equal the complete augmented system."""

import numpy as np
import pytest

from hp_coupled_solve import condense


def test_schur_correction_matches_complete_dense_augmented_system():
    rng = np.random.default_rng(25)
    a = rng.normal(size=(9, 9)); a = a.T @ a + np.eye(9)
    coupling = rng.normal(size=(9, 6)) / 10
    responses = np.linalg.solve(a, coupling)
    gram = coupling.T @ responses + np.eye(6)
    residual = rng.normal(size=(6, 6)) / 10
    h = np.eye(6) * 3
    cross = residual.T + h
    p = gram - h + cross + cross.T
    energy = np.block([[h, cross], [cross.T, p]])
    coefficients, tensor, _, _ = condense(energy, coupling, responses)
    augmented = np.block([[a, coupling], [coupling.T, gram]])
    rhs = np.vstack((np.zeros((9, 6)), residual))
    direct = -np.linalg.solve(augmented, rhs)
    np.testing.assert_allclose(coefficients, direct[9:], atol=1e-13)
    np.testing.assert_allclose(-responses @ coefficients, direct[:9], atol=1e-13)
    np.testing.assert_allclose(tensor, h + rhs.T @ direct, atol=1e-13)
    with pytest.raises(AssertionError):
        condense(energy, coupling, 10000 * responses)
