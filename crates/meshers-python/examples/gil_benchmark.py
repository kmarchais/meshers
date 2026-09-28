"""Measure concurrent independent meshes and native worker scaling separately."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import statistics
import time
from pathlib import Path

import numpy as np
import meshers


def field(x, y, z):
    x, y, z = x * (2 * np.pi), y * (2 * np.pi), z * (2 * np.pi)
    return (np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)) / 0.3


def graded(x, y, z):
    return field(x, y, z) / (1 + 0.3 * x)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepared = meshers.compile_field(field)
    graded_prepared = meshers.compile_field(graded)
    rows = []

    def run(evaluator, workers=1, repeat=1, small=False):
        f = field if evaluator == "callback" else graded_prepared if repeat == 2 else prepared
        result = meshers.generate(
            f, compile=False, cells=10 if small else 15 if repeat == 1 else 23,
            bounds=(-repeat / 2, repeat / 2) * 3, band=(-1, 1),
            periodic=(repeat == 1, True, True), optimize_passes=1 if small else 4,
            threads=workers, minimum_quality=0.0,
        )
        fingerprint = hashlib.sha256(result.points.tobytes() + result.tetrahedra.tobytes()).hexdigest()
        return {"fingerprint": fingerprint, "tets":len(result.tetrahedra), **result.diagnostics}

    # Same jobs and native worker count; only Python concurrency changes.
    for evaluator, small in [("compiled", False), ("callback", True)]:
        reference = run(evaluator, small=small)
        for concurrency in [1, 2, 4]:
            times = []
            for trial in range(3):
                start = time.perf_counter()
                with ThreadPoolExecutor(max_workers=concurrency) as executor:
                    outputs = list(executor.map(lambda _: run(evaluator, small=small), range(4)))
                times.append(time.perf_counter() - start)
                assert all(r["fingerprint"] == reference["fingerprint"] for r in outputs)
            row = {"kind":"four_independent_meshes", "evaluator":evaluator,
                   "small_callback_probe":small,"python_threads":concurrency,"native_workers":1,
                   "seconds":statistics.median(times),"times":times,"result":reference}
            rows.append(row)
            print(json.dumps(row), flush=True)
            args.output.write_text(json.dumps(rows, indent=2))

    # One larger graded job. Colored ordering is the same for 1, 2 and 4 workers.
    reference = run("compiled", repeat=2)
    for workers in [1, 2, 4]:
        times = []
        for trial in range(3):
            start = time.perf_counter()
            result = run("compiled", workers=workers, repeat=2)
            times.append(time.perf_counter() - start)
            assert result["fingerprint"] == reference["fingerprint"]
        row = {"kind":"graded_single_mesh", "evaluator":"compiled", "native_workers":workers,
               "seconds":statistics.median(times),"times":times,"result":result}
        rows.append(row)
        print(json.dumps(row), flush=True)
        args.output.write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
