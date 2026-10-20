"""Synthetic controls for three-level full-tensor refinement diagnostics."""

import numpy as np
import pytest

from cubic_refinement import fit_refinement


@pytest.mark.parametrize("order", [.5, 2., 4.])
def test_known_common_matrix_order_and_limit_are_recovered(order):
    lengths = np.array([1.2, 1.1, 1.])
    limit = np.diag([3., 4., 5., 1., 1.5, 2.])
    coefficient = np.ones((6, 6)) * .002 + np.eye(6) * .02
    tensors = np.array([limit + coefficient * length**order for length in lengths])
    row = fit_refinement(tensors, lengths, 1e-9)
    assert row["compatible"]
    assert abs(row["observed_order"] - order) < 1e-9
    np.testing.assert_allclose(row["extrapolated_stiffness"], limit, atol=1e-12, rtol=0)


def test_reversed_corrections_are_not_positive_order_convergence():
    limit = np.eye(6) * 2
    tensors = [limit, limit + np.eye(6) * .1, limit]
    row = fit_refinement(tensors, [1.2, 1.1, 1.], 1e-9)
    assert not row["compatible"] and row["best_correction_ratio"] < 0
    assert "Correction directions reverse" in row["reasons"]


def test_nonproportional_matrix_corrections_are_rejected():
    limit = np.eye(6) * 2
    middle = limit + np.diag([0., .1, 0., 0., 0., 0.])
    coarse = middle + np.diag([.1, 0., 0., 0., 0., 0.])
    row = fit_refinement([coarse, middle, limit], [1.2, 1.1, 1.], 1e-9)
    assert not row["compatible"] and row["matrix_model_relative_residual"] > .1


def test_quadrature_unresolved_corrections_do_not_promote_a_limit():
    lengths = np.array([1.2, 1.1, 1.])
    tensors = np.array([np.eye(6) * (2 + 1e-7 * length**2) for length in lengths])
    row = fit_refinement(tensors, lengths, 1e-5)
    assert not row["compatible"] and not row["signal_above_ten_times_quadrature"]
    assert row["extrapolated_stiffness"] is None


def test_nonrefining_lengths_are_rejected():
    with pytest.raises(AssertionError):
        fit_refinement([np.eye(6)] * 3, [1., 1.1, 1.2], 1e-9)
