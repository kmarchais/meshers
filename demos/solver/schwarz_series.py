# Copyright (c) 2026 Kevin Marchais

"""Replay alternating Schwarz P refinement levels in fresh, single-threaded workers."""

from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import subprocess  # ruff: ignore[suspicious-subprocess-import] - isolated benchmark workers.
import sys
import time
from contextlib import suppress
from pathlib import Path
from typing import TYPE_CHECKING, Literal, NotRequired, TypedDict, cast

import psutil

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping, Sequence
    from typing import BinaryIO


class WorkerRecord(TypedDict):
    """Fields shared by successful, rejected and failed refinement attempts."""

    route: str
    cells: int | None
    status: str
    workflow_seconds: float
    error: NotRequired[str]
    sampled_peak_tree_rss_bytes: NotRequired[int]
    repeat: NotRequired[int]


WORKFLOW = Path(__file__).with_name("schwarz_workflow.py")
LEVELS = {
    "native": (12, 16, 24),
    "microgen-default": (None, 24, 30),
    "microgen-tuned": (18, 24, 30),
}
TARGET = "within-reference-target"


def summarize(rows: Sequence[WorkerRecord], repeats: int) -> list[dict[str, object]]:
    """Return observed target costs only for levels passing every repetition."""
    result: list[dict[str, object]] = []
    for route, levels in LEVELS.items():
        for cells in levels:
            group = [r for r in rows if (r["route"], r["cells"]) == (route, cells)]
            passed = [r for r in group if r["status"] == TARGET]
            complete = len(group) == repeats and len(passed) == repeats
            result.append({
                "route": route,
                "cells": cells,
                "observations": len(group),
                "within_target": len(passed),
                "statuses": [r["status"] for r in group],
                "median_workflow_seconds_at_target": statistics.median(
                    r["workflow_seconds"] for r in passed
                )
                if complete
                else None,
            })
    return result


def observe_worker(
    cmd: Sequence[str],
    env: Mapping[str, str],
    log: BinaryIO,
    timeout: float,
) -> tuple[int | Literal["timeout"], int]:
    """Return worker exit status and peak sampled RSS including MMG children."""
    deadline = time.monotonic() + timeout
    peak = 0
    with subprocess.Popen(cmd, env=env, stdout=log, stderr=subprocess.STDOUT) as worker:  # ruff: ignore[subprocess-without-shell-equals-true] - fixed local worker.
        parent = psutil.Process(worker.pid)
        while worker.poll() is None:
            processes = [parent]
            with suppress(psutil.NoSuchProcess):
                processes.extend(parent.children(recursive=True))
            rss = 0
            for process in processes:
                with suppress(psutil.NoSuchProcess):
                    rss += process.memory_info().rss
            peak = max(peak, rss)
            if time.monotonic() >= deadline:
                for process in reversed(processes):
                    with suppress(psutil.NoSuchProcess):
                        process.kill()
                worker.wait()
                return "timeout", peak
            with suppress(subprocess.TimeoutExpired):
                worker.wait(timeout=0.25)
        return worker.returncode, peak


def worker_record(
    path: Path,
    route: str,
    cells: int | None,
    exit_status: int | Literal["timeout"],
    elapsed: float,
) -> WorkerRecord:
    """Read the worker output.

    Returns:
        A failed observation for incomplete or contradictory output, or the worker result.
    """
    try:
        record = json.loads(path.read_bytes())
    except (FileNotFoundError, ValueError):
        record = None
    if not isinstance(record, dict):
        record = {}
    duration = record.get("workflow_seconds")
    valid_time = isinstance(duration, (int, float)) and math.isfinite(duration) and duration >= 0
    valid_case = record.get("route") == route and "cells" in record and record["cells"] == cells
    outcomes = (TARGET, "outside-reference-target", "geometry-rejected", "failed")
    if not valid_case or not valid_time or record.get("status") not in outcomes:
        return {
            "route": route,
            "cells": cells,
            "status": "failed",
            "error": f"Worker exit {exit_status} without a complete result",
            "workflow_seconds": elapsed,
        }
    if exit_status and record["status"] != "failed":
        record.update(status="failed", error=f"Worker exited with {exit_status}")
    return cast("WorkerRecord", record)


def run_series(
    output: Path,
    repeats: int,
    timeout: float,
    on_step: Callable[[int, str], None] | None = None,
) -> None:
    """Run a fixed refinement series in new directories and retain each worker outcome.

    Raises:
        ValueError: Repetitions and timeout must be positive.
    """
    if repeats < 1 or timeout <= 0:
        msg = "Expected positive repetitions and timeout"
        raise ValueError(msg)
    output.mkdir()
    env = {
        **os.environ,
        **dict.fromkeys(
            (
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "VTK_SMP_MAX_THREADS",
                "RAYON_NUM_THREADS",
            ),
            "1",
        ),
    }
    rows: list[WorkerRecord] = []
    for repeat in range(repeats):
        for level in range(3):
            routes = list(LEVELS)
            offset = (repeat + level) % len(routes)
            for route in routes[offset:] + routes[:offset]:
                cells = LEVELS[route][level]
                name = f"r{repeat}-{route}-{cells or 'default'}"
                folder = output / name
                cmd = [sys.executable, str(WORKFLOW), route, str(folder)]
                if cells is not None:
                    cmd.extend(("--cells", str(cells)))
                if on_step:
                    on_step(len(rows), name)
                start = time.perf_counter()
                with (output / (name + ".log")).open("wb") as log:
                    status, peak = observe_worker(cmd, env, log, timeout)
                path = folder / "result.json"
                record = worker_record(path, route, cells, status, time.perf_counter() - start)
                folder.mkdir(exist_ok=True)
                record["sampled_peak_tree_rss_bytes"] = peak
                path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf8")
                rows.append({**record, "repeat": repeat})
                summary = {
                    "repetitions": repeats,
                    "timeout_seconds": timeout,
                    "completed_workers": len(rows),
                    "levels": summarize(rows, repeats),
                }
                (output / "summary.json").write_text(
                    json.dumps(summary, indent=2) + "\n",
                    encoding="utf8",
                )


def main() -> None:
    """Run the same timed series from a clean output directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    run_series(args.output, args.repeats, args.timeout)


if __name__ == "__main__":
    main()
