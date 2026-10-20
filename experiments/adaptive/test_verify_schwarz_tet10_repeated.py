"""The independent verifier must detect mixed-load errors, not just diagonals."""

import numpy as np

from verify_schwarz_tet10_repeated import model_errors


def test_mixed_load_off_diagonal_error_is_detected():
    reference = np.eye(6)
    candidate = reference.copy()
    candidate[0, 1] = candidate[1, 0] = .002
    np.testing.assert_array_equal(candidate.diagonal(), reference.diagonal())
    errors = model_errors(candidate, reference, reference, [.8, .9])
    np.testing.assert_allclose(errors, [.002] * 8, atol=1e-14, rtol=0)


def test_all_order_and_length_models_are_used():
    errors = model_errors(np.eye(6), np.eye(6) * 1.01, np.eye(6), [.8, .9])
    assert len(errors) == 8 and len(set(errors)) == 8
