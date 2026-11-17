"""Replay alternating Schwarz P refinement levels in fresh, single-threaded workers."""

import argparse
import json
import os
import statistics
import subprocess  # ruff: ignore[suspicious-subprocess-import] - isolated benchmark workers.
import sys
import time
from contextlib import suppress
from pathlib import Path

import psutil

LEVELS = {
    "native": (12, 16, 24),
    "microgen-default": (None, 24, 30),
    "microgen-tuned": (18, 24, 30),
}
TARGET = "within-reference-target"


def summarize(rows, repeats):
    """Return observed target costs only for levels passing every repetition."""
    result = []
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


def observe_worker(cmd, env, log, timeout):
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


def run_series(output, repeats, timeout, on_step=None):
    """Run a fixed refinement series in new directories and retain each worker outcome.

    Raises:
        ValueError: Repetitions and timeout must be positive.
        RuntimeError: The worker exit code contradicts its recorded result.
    """
    if repeats < 1 or timeout <= 0:
        raise ValueError("Expected positive repetitions and timeout")
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
    script = Path(__file__).with_name("schwarz_workflow.py")
    rows = []
    for repeat in range(repeats):
        for level in range(3):
            routes = list(LEVELS)
            offset = (repeat + level) % len(routes)
            for route in routes[offset:] + routes[:offset]:
                cells = LEVELS[route][level]
                name = f"r{repeat}-{route}-{cells or 'default'}"
                folder = output / name
                cmd = [sys.executable, str(script), route, str(folder)]
                if cells is not None:
                    cmd.extend(("--cells", str(cells)))
                if on_step:
                    on_step(len(rows), name)
                start = time.perf_counter()
                with (output / (name + ".log")).open("wb") as log:
                    status, peak = observe_worker(cmd, env, log, timeout)
                path = folder / "result.json"
                if path.exists():
                    record = json.loads(path.read_text(encoding="utf8"))
                else:
                    folder.mkdir(exist_ok=True)
                    record = {
                        "route": route,
                        "cells": cells,
                        "status": "failed",
                        "error": f"Worker exit {status}",
                        "workflow_seconds": time.perf_counter() - start,
                    }
                record["sampled_peak_tree_rss_bytes"] = peak
                path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf8")
                if status and record["status"] != "failed":
                    raise RuntimeError(f"Worker {name} failed after writing an inconsistent result")
                rows.append({**record, "repeat": repeat})
                summary = {
                    "repetitions": repeats,
                    "timeout_seconds": timeout,
                    "completed_workers": len(rows),
                    "levels": summarize(rows, repeats),
                }
                (output / "summary.json").write_text(
                    json.dumps(summary, indent=2) + "\n", encoding="utf8"
                )


def main():
    """Run the same timed series from a clean output directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    run_series(args.output, args.repeats, args.timeout)


if __name__ == "__main__":
    main()
