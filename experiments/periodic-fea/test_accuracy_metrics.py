"""Analytic checks for local complementary-energy comparisons and point location."""

import numpy as np
import pytest
from accuracy_metrics import field_score, probe_cells
from reference_sensitivity import read


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


@pytest.mark.parametrize("failed_metric", ["sampled_geometry_distance", "sampled_reverse_distance",
                                         "sampled_cap_overshoot", "relative_volume_error_with_uncertainty"])
def test_common_gate_checks_both_directions_caps_and_volume(tmp_path, failed_metric):
    import json

    row = {"status": "ok", "homogenization_audited": True, "quality_eligible": True, "periodic_eligible": True,
           "sampled_geometry_distance": .0001, "sampled_reverse_distance": .0001,
           "sampled_cap_overshoot": .0001, "sampled_outward_orientation": True,
           "relative_volume_error_with_uncertainty": .001}
    row[failed_metric] = .02
    path = tmp_path / "result.json"
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError, match="Unadmitted"):
        read(path, geometry_gate=.002)


def test_failed_post_solve_check_cannot_supply_an_accuracy_candidate(tmp_path):
    import json

    row = {"status": "failed", "geometry_eligible": True, "homogenization_audited": True,
           "quality_eligible": True, "periodic_eligible": True}
    path = tmp_path / "result.json"
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError, match="Unadmitted"):
        read(path)


def test_automatic_solver_bypass_cannot_supply_an_accuracy_candidate(tmp_path):
    import json

    row = {"status": "ok", "geometry_eligible": True, "homogenization_audited": True,
           "quality_eligible": True, "periodic_eligible": True,
           "homogenization": {"strain_controlled_loads": [{}] * 6, "automatic_stress_controlled_solves": []}}
    path = tmp_path / "result.json"
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError, match="twelve audited"):
        read(path)
