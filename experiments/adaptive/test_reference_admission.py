"""Keep failed solves and weakened quadrature checks out of tensor comparisons."""

import json
import sys

import pytest

import microgen_workflow as workflow
import stiffness_pair
from sizing_pilot import enforce_mesh_cap

pytest.importorskip("fedoo")
pytest.importorskip("microgen")
pytest.importorskip("pyvista")


def test_post_adaptation_mesh_cap():
    enforce_mesh_cap(83_574, workflow.CAP)
    with pytest.raises(ValueError, match="mesh element cap"):
        enforce_mesh_cap(83_574, 424_754)
    with pytest.raises(ValueError, match="mesh element cap"):
        enforce_mesh_cap(111_087, 567_821)


@pytest.mark.parametrize("updates", [None, {"status": "failed"},
    {"quadrature_sensitivity_passed": False}, {"quadrature_sensitivity_limit": 0.01},
    {"qualified": False}, {"elements": 424_754}])
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


@pytest.mark.parametrize(("source", "expected"), [
    ("g3-schwarz-microgen-straight-tet10-v1", True), ("g3-schwarz-curved14-v1", False)])
def test_supplementary_stricter_target(tmp_path, monkeypatch, source, expected):
    output = tmp_path / "comparison.json"
    monkeypatch.setattr(sys, "argv", ["stiffness_pair.py", "--first",
        str(workflow.HERE / "runs" / source / "result.json"), "--reference-summary",
        str(workflow.HERE / "runs/g3-microgen-defaults-v4/summary.json"),
        "--orders", "1", "2", "3", "4", "--target", "0.001", "--output", str(output)])
    stiffness_pair.main()
    report = json.loads(output.read_text())
    assert report["stiffness_target"] == 0.001
    assert report["passes_target_under_all_declared_models"] is expected
    assert report["passes_one_percent_under_all_declared_models"]
