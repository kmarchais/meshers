"""Keep failed solves and weakened quadrature checks out of tensor comparisons."""

import json
import sys

import pytest

import microgen_workflow as workflow
import stiffness_pair


@pytest.mark.parametrize("updates", [None, {"status": "failed"},
    {"quadrature_sensitivity_passed": False}, {"quadrature_sensitivity_limit": 0.01},
    {"qualified": False}])
def test_four_point_comparison_admission(tmp_path, monkeypatch, updates):
    control = workflow.HERE / "runs/g3-lattice-uniform34-quadrature4-v1/result.json"
    reference = workflow.HERE / "runs/g3-lattice-adapt24-tol0028-v1/result.json"
    row = json.loads(control.read_text())
    assert row["order"] == 2 and row["solver_quadrature_points"] == 4
    assert row["qualified"] and row["quadrature_sensitivity_passed"]
    if updates:
        row.update(updates)
    first = tmp_path / "first.json"
    workflow.save(first, row)
    output = tmp_path / "comparison.json"
    monkeypatch.setattr(sys, "argv", ["stiffness_pair.py", "--first", str(first),
        "--second", str(reference), "--output", str(output)])
    if updates is None:
        stiffness_pair.main()
        assert json.loads(output.read_text())["difference"]["worst_relative_strain_energy_difference"] < 0.01
    else:
        with pytest.raises(AssertionError):
            stiffness_pair.main()
        assert not output.exists()
