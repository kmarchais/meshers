"""The independent verifier must detect mixed-load errors, not just diagonals."""

import numpy as np
import pytest

import microgen_workflow as workflow
from verify_schwarz_tet10_repeated import model_errors, verify_hash


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


def test_measured_driver_can_be_verified_after_source_changes(tmp_path):
    source = tmp_path / "schwarz_tet10_live.py"
    archive = tmp_path / "warmup/driver.py"
    archive.parent.mkdir()
    source.write_text("measured source")
    expected = workflow.digest(source)
    assert not verify_hash(source, expected, tmp_path)
    archive.write_text(source.read_text())
    source.write_text("later source")
    assert verify_hash(source, expected, tmp_path)
    archive.write_text("wrong archive")
    with pytest.raises(AssertionError):
        verify_hash(source, expected, tmp_path)


def test_changed_data_never_uses_source_archive_fallback(tmp_path):
    data = tmp_path / "fields.npz"
    data.write_bytes(b"original")
    expected = workflow.digest(data)
    data.write_bytes(b"changed")
    with pytest.raises(AssertionError):
        verify_hash(data, expected, tmp_path)
