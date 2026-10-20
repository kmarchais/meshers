"""Reproduce the late audit finding at actual standard-Job transitions."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path.home() / "bin"))
import pueue_progress

from job_parts import worker_part


@pytest.mark.parametrize("fail", [False, True])
def test_worker_part_retains_values_at_completion(tmp_path, monkeypatch, fail):
    records = []
    def capture(fraction, message="", *, eta_s=None, **values):
        records.append({"message": message, "values": values})
    monkeypatch.setenv("PUEUE_PROGRESS_FILE", str(tmp_path / "isolated-progress.json"))
    monkeypatch.setattr(pueue_progress, "report", capture)
    with pueue_progress.Job(parts=["worker"], every_s=0) as job:
        try:
            with worker_part(job, "worker") as part:
                part.update(.5, "working", rss_gib=1., peak_gib=1.)
                if fail:
                    raise ValueError("controlled failure")
        except ValueError:
            assert fail
    assert records[-2]["values"]
    assert records[-1]["values"]["memory_ceiling_gib"] == 8
    assert records[-1]["values"]["elapsed_s"] >= 0
    assert ("failed" in records[-1]["message"]) == fail
