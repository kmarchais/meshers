"""Controls for the matched-cost screen and the held-out solver trial."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
import schwarz_trials


def test_only_admitted_settings_enter_rotating_repeats(tmp_path, monkeypatch):
    calls = []

    def record(_out, label, arm):
        calls.append((label, arm))
        status = schwarz_trials.TARGET if arm in {"h050", "h030", "native"} else "geometry-rejected"
        if label == "r1-h030":
            status = "outside-reference-target"
        return {"status": status, "wall_seconds": 3}

    monkeypatch.setattr(schwarz_trials, "measure", record)
    schwarz_trials.replay("sizing", tmp_path / "trial")
    result = json.loads((tmp_path / "trial" / "done.json").read_bytes())
    assert result["median_wall_seconds_at_target"]["h030"] is None
    assert result["fastest_measured_microgen"] == "h050"
    assert [a for _, a in calls[4:]] == [
        "native",
        "h050",
        "h030",
        "h050",
        "h030",
        "native",
        "h030",
        "native",
        "h050",
    ]


@pytest.mark.parametrize("admitted", [False, True])
def test_native_bracket_requires_an_admitted_candidate(tmp_path, monkeypatch, admitted):
    calls = []

    def record(_out, _label, arm):
        calls.append(arm)
        return {
            "status": schwarz_trials.TARGET if admitted else "geometry-rejected",
            "wall_seconds": 1,
        }

    monkeypatch.setattr(schwarz_trials, "measure", record)
    schwarz_trials.replay("bracket", tmp_path / "trial")
    expected = ["near-target", "h050", "h050", "near-target", "near-target", "h050"]
    assert calls == ["near-target", *(expected if admitted else [])]


@pytest.mark.parametrize("change", ["none", "mesh", "hardware", "tensor"])
def test_phase_speedup_requires_matching_inputs_and_stiffness(tmp_path, monkeypatch, change):
    def record(_out, label, _arm):
        altered = label == "r2-shifted-rcm"
        return {
            "status": schwarz_trials.TARGET,
            "wall_seconds": 2 if label.endswith("rcm") else 3,
            "solver_mesh_sha256": "other" if altered and change == "mesh" else "original",
            "hardware": {"cpu_count": 4 if altered and change == "hardware" else 8},
            "fea": {
                "stiffness_engineering_voigt": [[1.01 if altered and change == "tensor" else 1.0]]
            },
        }

    monkeypatch.setattr(schwarz_trials, "measure", record)
    schwarz_trials.replay("phase", tmp_path / "trial")
    result = json.loads((tmp_path / "trial" / "done.json").read_bytes())
    assert result["identical_solver_mesh"] is (change != "mesh")
    assert result["identical_hardware"] is (change != "hardware")
    if change == "none":
        assert result["wall_change_percent"] == pytest.approx(-100 / 3)
    else:
        assert result["wall_change_percent"] is None
