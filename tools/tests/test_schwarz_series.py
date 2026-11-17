"""A failed repeat disqualifies the level from the matched-target timing table."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
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
