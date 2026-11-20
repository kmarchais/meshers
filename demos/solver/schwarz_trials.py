# Copyright (c) 2026 Kevin Marchais

"""Compare Schwarz P mesh costs at a fixed accuracy and solver ordering on a held-out phase."""

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
            command,
            {**os.environ, **dict.fromkeys(THREAD_VARIABLES, "1")},
            log,
            750,
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


def phase_comparison(repeats, times):
    """Return the change in wall time for admitted solves on the same mesh and hardware."""
    rows = [row for arm in repeats.values() for row in arm]
    hashes = [row.get("solver_mesh_sha256") for row in rows]
    identical = None not in hashes and len(set(hashes)) == 1
    configs = [row.get("hardware") for row in rows]
    same_config = configs[0] is not None and all(h == configs[0] for h in configs)
    baseline = times["shifted-baseline"]
    reordered = times["shifted-rcm"]
    result = {
        "identical_solver_mesh": identical,
        "identical_hardware": same_config,
        "largest_tensor_entry_difference": None,
        "wall_change_percent": None,
    }
    if identical and same_config and baseline and reordered:
        tensors = [
            [x for line in row["fea"]["stiffness_engineering_voigt"] for x in line] for row in rows
        ]
        difference = max(max(col) - min(col) for col in zip(*tensors, strict=True))
        result["largest_tensor_entry_difference"] = difference
        if difference <= 1e-9 * max(abs(x) for x in tensors[0]):
            result["wall_change_percent"] = 100 * (reordered / baseline - 1)
    return result


def replay(study, output):
    """Run fresh workers for the selected study and retain every admission decision."""
    output.mkdir()
    if study == "sizing":
        qualified = [
            arm for arm in REPAIR_SIZES if measure(output, f"screen-{arm}", arm)["status"] == TARGET
        ]
        arms = ["native", *qualified] if qualified else []
        result = {"screened": REPAIR_SIZES, "qualified": qualified}
    elif study == "bracket":
        candidate = measure(output, "screen-near-target", "near-target")
        arms = ["near-target", "h050"] if candidate["status"] == TARGET else []
        result = {"settings": ["near-target", "h050"], "screen_status": candidate["status"]}
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
    elif study == "phase":
        result.update(phase_comparison(repeats, times))
    (output / "done.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", choices=("sizing", "phase", "bracket"))
    parser.add_argument("output", type=Path)
    opts = parser.parse_args()
    replay(opts.study, opts.output)
