"""Replay the resolution-30 sizing screen or the held-out phase experiment on one configuration."""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from statistics import median

from schwarz_series import TARGET, observe_worker

REPAIR_SIZES = {"h050": 0.05, "h040": 0.04, "h035": 0.035, "h030": 0.03}
THREAD_VARIABLES = (
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
    "RAYON_NUM_THREADS",
    "VTK_SMP_MAX_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)


def measure(output, label, arm):
    """Return one fresh worker's outcome, including process startup and peak tree RSS."""
    path = output / f"{label}.json"
    if arm in REPAIR_SIZES:
        args = ["repaired-solve", str(path), "--repair-hmax", str(REPAIR_SIZES[arm])]
    else:
        args = ["finest" if arm == "native" else arm, str(path)]
    command = [sys.executable, str(Path(__file__).with_name("schwarz_diagnostics.py")), *args]
    start = time.perf_counter()
    with (output / f"{label}.log").open("xb") as log:
        code, peak = observe_worker(
            command, {**os.environ, **dict.fromkeys(THREAD_VARIABLES, "1")}, log, 750
        )
    try:
        row = json.loads(path.read_bytes())
    except (FileNotFoundError, ValueError):
        row = {"status": "failed", "error": f"Worker exit {code} without a complete record"}
    if code and row["status"] != "failed":
        row.update(status="failed", error=f"Worker exited with {code}")
    row.update(wall_seconds=time.perf_counter() - start, sampled_peak_tree_rss_bytes=peak)
    path.write_text(json.dumps(row, indent=2, allow_nan=False) + "\n", encoding="utf8")
    return row


def summarize(repeats):
    """Return median wall times for arms meeting the target in all three repeats."""
    return {
        arm: median(r["wall_seconds"] for r in rows)
        if len(rows) == 3 and all(r["status"] == TARGET for r in rows)
        else None
        for arm, rows in repeats.items()
    }


def replay(study, output):
    """Run the four-setting screen and qualified repeats, or six frozen phase trials."""
    output.mkdir()
    if study == "sizing":
        qualified = [
            arm for arm in REPAIR_SIZES if measure(output, f"screen-{arm}", arm)["status"] == TARGET
        ]
        arms = ["native", *qualified] if qualified else []
        result = {"screened": REPAIR_SIZES, "qualified": qualified}
    else:
        arms = ["shifted-baseline", "shifted-rcm"]
        result = {"phase_x": 0.125, "order": "AB BA AB"}
    repeats = {arm: [] for arm in arms}
    for repeat in range(3):
        offset = repeat % len(arms) if arms else 0
        for arm in arms[offset:] + arms[:offset]:
            repeats[arm].append(measure(output, f"r{repeat}-{arm}", arm))
    times = summarize(repeats)
    result["median_wall_seconds_at_target"] = times
    if study == "sizing":
        eligible = {arm: t for arm, t in times.items() if arm in REPAIR_SIZES and t is not None}
        result["fastest_measured_microgen"] = min(eligible, key=eligible.get) if eligible else None
    else:
        hashes = [r.get("solver_mesh_sha256") for rows in repeats.values() for r in rows]
        result["identical_solver_mesh"] = None not in hashes and len(set(hashes)) == 1
    (output / "done.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", choices=("sizing", "phase"))
    parser.add_argument("output", type=Path)
    opts = parser.parse_args()
    replay(opts.study, opts.output)
