# Copyright (c) 2026 Kevin Marchais

"""Replay the dated cap study on one installed wheel and a saved baseline worker."""

import argparse
import datetime as dt
import hashlib
import json
import runpy
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

from tools import performance_worker as current
from tools.benchmark_hardware import hardware_snapshot
from tools.performance_gate import RELEASE_CASES
from tools.performance_worker import _measure_generation
from tools.tests.test_periodic_trace import (
    dropped_pair,
    flipped_diagonal,
    rotational_fixture,
    swapped_midpoints,
    unit_shift,
)


def verdict(sample, case, generated):
    with patch.dict(sample.__globals__, _measure_generation=lambda _: generated):
        try:
            return "pass", sample(case)
        except RuntimeError as error:
            return str(error), None


def measure(case, before):
    mesh, seconds, peak = _measure_generation(case)
    variants = {"original": mesh}
    if case in {"compiled_gyroid", "gyroid_tet10"}:
        variants["drop_pair"] = dropped_pair(mesh)
        if case == "compiled_gyroid":
            variants["flip_volume_and_surface_diagonal"] = flipped_diagonal(mesh)
        else:
            variants["swap_surface_midpoints"] = swapped_midpoints(mesh)
    results = {}
    for name, variant in variants.items():
        generated = variant, seconds, peak
        old, old_record = verdict(before, case, generated)
        new, new_record = verdict(current.sample, case, generated)
        results[name] = {"before": old, "after": new}
        if name == "original":
            results[name]["records_identical"] = new_record == old_record and new == "pass"
    trace_ms = []
    if case not in {"graded", "mapped_intersection"}:
        for _ in range(6):
            start = perf_counter()
            for axis in range(3):
                current.check_periodic_axis(mesh, axis, unit_shift(axis))
            trace_ms.append(1000 * (perf_counter() - start))
    return {
        "nodes": len(mesh.points),
        "tetrahedra": len(mesh.tetrahedra),
        "periodic_pairs": [len(p) for p in mesh.periodic_pairs],
        "generation_seconds": seconds,
        "three_axis_check_ms": trace_ms[1:],
        "controls": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_worker", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    baseline = runpy.run_path(str(args.baseline_worker))
    result = {
        "recorded_at": dt.datetime.now(dt.UTC).isoformat(),
        **hardware_snapshot(),
        **current.runtime_identity(),
        "baseline_worker_sha256_lf": hashlib.sha256(
            args.baseline_worker.read_text(encoding="utf-8").encode(),
        ).hexdigest(),
        "worker_sha256_lf": hashlib.sha256(
            Path(current.__file__).read_text(encoding="utf-8").encode(),
        ).hexdigest(),
        "cases": {case: measure(case, baseline["sample"]) for case in RELEASE_CASES},
    }
    mesh, transform = rotational_fixture()
    current.check_periodic_axis(mesh, 2, transform)
    result["mapped_rotation"] = {"pairs": len(mesh.periodic_pairs[2]), "radians": 1 / 3}
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for case, row in result["cases"].items():
        for name, control in row["controls"].items():
            assert control["before"] == "pass", (case, name, control)
            if name == "original":
                assert control["records_identical"], (case, control)
            else:
                assert control["after"] != "pass", (case, name, control)


if __name__ == "__main__":
    main()
