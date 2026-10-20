"""The pair bound cannot replace full-matrix checking of a common witness."""

import numpy as np
import pytest

from reference_envelope import assess


@pytest.mark.parametrize("ratio,contradiction", [(1.01, False), (1.04, True)])
def test_scalar_scaled_models_have_exact_minimax_bound(ratio, contradiction):
    report = assess([np.eye(6), ratio * np.eye(6)])
    expected = (ratio - 1) / (ratio + 1)
    assert abs(report["pairwise_lower_bound"] - expected) < 1e-13
    assert report["contradiction_passed"] is contradiction
    assert abs(report["harmonic_trial_maximum_disagreement"] - expected) < 1e-13


def test_mixed_load_direction_and_congruence_invariance():
    basis = np.eye(6)
    basis[:2, :2] = [[.8, -.6], [.6, .8]]
    high = basis @ np.diag([1.04, 1, 1, 1, 1, 1]) @ basis.T
    report = assess([np.eye(6), high])
    assert report["contradiction_passed"]
    direction = np.asarray(report["mixed_strain_direction"])
    assert abs(direction[0]) > .1 and abs(direction[1]) > .1
    transform = np.diag([1, 2, 3, 4, 5, 6])
    changed = assess([transform.T @ matrix @ transform for matrix in (np.eye(6), high)])
    assert abs(changed["pairwise_lower_bound"] - report["pairwise_lower_bound"]) < 1e-13


def test_identical_models_and_invalid_stiffness():
    report = assess([np.eye(6), np.eye(6)])
    assert report["pairwise_lower_bound"] == 0 and report["harmonic_trial_passes_all_models"]
    with pytest.raises(np.linalg.LinAlgError):
        assess([np.eye(6), -np.eye(6)])
