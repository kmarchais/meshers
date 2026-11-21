# Copyright (c) 2026 Kevin Marchais

"""A failed repeat disqualifies the level from the matched-target timing table."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
import schwarz_series
from schwarz_series import TARGET, summarize


def test_times_require_every_repeat_to_meet_target():
    success = {"route": "native", "cells": 16, "status": TARGET, "workflow_seconds": 12.0}
    rejected = {**success, "status": "geometry-rejected", "workflow_seconds": 1.0}

    def level(rows):
        return next(r for r in summarize(rows, 2) if r["route"] == "native" and r["cells"] == 16)

    assert level([success])["median_workflow_seconds_at_target"] is None
    assert level([rejected, success])["median_workflow_seconds_at_target"] is None
    assert (
        level([success, {**success, "workflow_seconds": 20}])["median_workflow_seconds_at_target"]
        == 16
    )


def test_incomplete_and_contradictory_workers_do_not_abort_series(tmp_path, monkeypatch):
    calls = 0

    def observe(cmd, _env, _log, _timeout):
        nonlocal calls
        folder = Path(cmd[3])
        folder.mkdir()
        path = folder / "result.json"
        calls += 1
        if calls == 1:
            path.write_text("{", encoding="utf8")
            return "timeout", 42
        record = {
            "route": cmd[2],
            "cells": int(cmd[5]) if len(cmd) > 4 else None,
            "status": TARGET,
            "workflow_seconds": 2,
        }
        path.write_text(json.dumps(record), encoding="utf8")
        return int(calls == 2), 42

    monkeypatch.setattr(schwarz_series, "observe_worker", observe)
    output = tmp_path / "series"
    schwarz_series.run_series(output, 1, 10)
    summary = json.loads((output / "summary.json").read_bytes())
    assert summary["completed_workers"] == 9
    rejected = [r for r in summary["levels"] if r["statuses"] == ["failed"]]
    assert len(rejected) == 2
    assert all(r["median_workflow_seconds_at_target"] is None for r in rejected)
    first = json.loads((output / "r0-native-12" / "result.json").read_bytes())
    assert "timeout" in first["error"]
    assert first["sampled_peak_tree_rss_bytes"] == 42
