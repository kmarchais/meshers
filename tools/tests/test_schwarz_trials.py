"""Controls for the matched-cost screen and the held-out solver trial."""

import json
import sys
from pathlib import Path

import numpy as np
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


def test_shifted_wall_uses_the_translated_level_set():
    pytest.importorskip("pyvista")
    pytest.importorskip("scipy")
    from periodic_geometry import field

    points = np.array([[0.125, 0.25, -0.25], [-0.47, 0.093, 0.29]])
    np.testing.assert_allclose(
        field(points, "schwarz-shifted"),
        field(points + np.array([0.125, 0, 0]), "schwarz"),
        atol=1e-14,
    )


def test_reordered_cg_returns_original_degrees_of_freedom():
    pytest.importorskip("fedoo")
    pytest.importorskip("meshers")
    from fedoo_homogenization import ReorderedCG
    from scipy.sparse import csr_matrix

    matrix = csr_matrix([
        [4.0, 0.0, -1.0, 0.0],
        [0.0, 2.0, 0.0, -1.0],
        [-1.0, 0.0, 5.0, -1.0],
        [0.0, -1.0, -1.0, 3.0],
    ])
    rhs = np.array([2.0, -3.0, 4.0, 1.0])
    solver = ReorderedCG()
    np.testing.assert_allclose(matrix @ solver(matrix, rhs), rhs, atol=1e-12)
    assert solver.records[-1]["original_relative_residual"] < 1e-12
