import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("fedoo")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
import schwarz_workflow as workflow


@pytest.mark.parametrize(
    ("status", "expected"),
    [("geometry-rejected", 1), ("outside-reference-target", 1), ("within-reference-target", 0)],
)
def test_candidate_exit_code_and_record(tmp_path, monkeypatch, status, expected):
    output = tmp_path / "candidate"
    monkeypatch.setattr(sys, "argv", ["study", "native-candidate", str(output)])
    monkeypatch.setattr(workflow, "hardware_snapshot", lambda: {"platform": "test"})

    def run(route, cells, folder, row, *, native_settings):
        assert route == "native-candidate"
        assert cells == 24
        assert "microgen" not in row["versions"]
        row["status"] = status

    monkeypatch.setattr(workflow, "run", run)
    assert workflow.main() == expected
    assert json.loads((output / "result.json").read_text())["status"] == status
