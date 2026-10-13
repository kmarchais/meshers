"""Analytic checks for local complementary-energy comparisons and point location."""

import numpy as np
import pytest
from accuracy_metrics import field_score, probe_cells


def test_energy_norm_scales_and_includes_mixed_loads():
    reference = np.eye(6)[None].repeat(3, axis=0)
    assert field_score(reference, reference)["worst_relative_stress_energy_norm_difference"] == 0
    assert field_score(reference * 1.1, reference)["worst_relative_stress_energy_norm_difference"] == pytest.approx(.1)
    candidate = reference.copy()
    # Two identical error columns produce a larger mixed-load error than either unit load.
    candidate[:, 3, 0] += .1
    candidate[:, 3, 1] += .1
    score = field_score(candidate, reference)
    assert score["worst_relative_stress_energy_norm_difference"] > max(score["per_unit_load_relative_stress_energy_norm_difference"])


def test_probe_reports_cell_membership_and_void():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    ids, valid = probe_cells(points, np.array([[0, 1, 2, 3]]), np.array([[.1, .1, .1], [2, 2, 2]]))
    assert valid.tolist() == [True, False]
    assert ids[0] == 0


def test_stress_score_is_invariant_to_load_basis_and_units():
    rng = np.random.default_rng(20261013)
    reference = rng.normal(size=(12, 6, 6))
    candidate = reference + .03 * rng.normal(size=reference.shape)
    basis = rng.normal(size=(6, 6)) + 4 * np.eye(6)
    expected = field_score(candidate, reference)["worst_relative_stress_energy_norm_difference"]
    actual = field_score(candidate @ basis, reference @ basis)["worst_relative_stress_energy_norm_difference"]
    assert actual == pytest.approx(expected, rel=1e-12)
