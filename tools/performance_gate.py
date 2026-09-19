"""Compare installed base/head wheels on one Linux runner in alternating fresh processes."""

import argparse
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

CASES = ("compiled_gyroid", "callback_gyroid", "graded", "mapped_intersection")


def sample(case):
    """Measure a fresh process, then independently validate its mesh.

    Returns:
        Timing, peak resident memory and independently computed mesh metrics.

    Raises:
        RuntimeError: The mesh fails the benchmark's validity or quality checks.
    """
    import resource

    import meshers
    import numpy as np

    def gyroid(x, y, z):
        a, b, c = 2 * np.pi * x, 2 * np.pi * y, 2 * np.pi * z
        return np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)

    options = {"cells": 24, "optimize_passes": 2, "geometry_tolerance": 0.1}
    start = time.perf_counter()
    if case == "mapped_intersection":

        def mapping(x, y, z):
            return (3 + x) * np.cos(z / 3), (3 + x) * np.sin(z / 3), -y

        mesh = meshers.generate_intersection(
            {"right": lambda x, y, z: x - 0.501, "back": lambda x, y, z: y - 0.751},
            coordinate_map=mapping,
            **options,
        )
    else:
        field = gyroid
        if case == "graded":

            def field(x, y, z):
                return gyroid(x, y, z) / (0.35 + 0.15 * x)

        mesh = meshers.generate(
            field,
            band=(-1, 1) if case == "graded" else (-0.5, 0.5),
            periodic=(False,) * 3 if case == "graded" else (True,) * 3,
            compile=case != "callback_gyroid",
            threads=1,
            batch_size=4096,
            **options,
        )
    elapsed = time.perf_counter() - start
    # Capture OS high-water RSS before allocating arrays for independent checks.
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    points = mesh.points[mesh.tetrahedra]
    determinants = np.linalg.det(points[:, 1:] - points[:, :1])
    if not len(points) or not np.all(np.isfinite(mesh.points)) or not np.all(determinants > 0):
        raise RuntimeError("Nonempty, finite, positively oriented mesh required")
    edges = sum(
        np.sum((points[:, i] - points[:, j]) ** 2, axis=1)
        for i in range(4)
        for j in range(i + 1, 4)
    )
    quality = np.sqrt(432 * determinants**2 / edges**3)
    boundary_edges = np.sort(
        np.concatenate([mesh.surface[:, [0, 1]], mesh.surface[:, [1, 2]], mesh.surface[:, [2, 0]]]),
        axis=1,
    )
    _, incidence = np.unique(boundary_edges, axis=0, return_counts=True)
    if not np.all(incidence == 2):
        raise RuntimeError("Closed manifold boundary required")
    if case in {"compiled_gyroid", "callback_gyroid"}:
        for axis, pairs in enumerate(mesh.periodic_pairs):
            expected = np.zeros(3)
            expected[axis] = 1
            if not len(pairs) or not np.allclose(
                mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]], expected, atol=1e-10, rtol=0
            ):
                raise RuntimeError("Periodic correspondence lost")
    record = {
        "seconds": elapsed,
        "peak_rss_bytes": peak,
        "tetrahedra": len(points),
        "volume": float(determinants.sum() / 6),
        "quality_min": float(quality.min()),
        "quality_p05": float(np.quantile(quality, 0.05)),
        "quality_median": float(np.median(quality)),
        "surface_error": mesh.diagnostics["sampled_surface_error"],
    }
    if record["quality_min"] < 0.01 or record["surface_error"] > 0.1:
        raise RuntimeError("Benchmark geometry or quality floor failed")
    return record


def regressions(base, head):
    """Return failures; resource noise allowances never relax geometry validity."""
    failures = []
    for key in base:
        if (
            not math.isfinite(base[key])
            or not math.isfinite(head[key])
            or base[key] < 0
            or head[key] < 0
        ):
            failures.append(f"invalid {key}")
    if failures:
        return failures
    for key, ratio, allowance in (("seconds", 1.25, 0.02), ("peak_rss_bytes", 1.15, 8 * 1024**2)):
        if head[key] > base[key] * ratio + allowance:
            failures.append(f"{key} increased: {base[key]:.6g} -> {head[key]:.6g}")
    for key in ("quality_min", "quality_p05", "quality_median"):
        if head[key] < base[key] * 0.99 - 1e-12:
            failures.append(f"{key} decreased: {base[key]:.6g} -> {head[key]:.6g}")
    for key, tolerance in (("volume", 0.005), ("tetrahedra", 0.01)):
        if abs(head[key] - base[key]) > abs(base[key]) * tolerance + 1e-12:
            failures.append(f"{key} changed beyond {tolerance:.1%}")
    if head["surface_error"] > base["surface_error"] * 1.05 + 1e-10:
        failures.append("surface error increased by more than 5%")
    return failures


def main():
    """Collect repeated measurements and fail when a regression exceeds its allowance.

    Raises:
        SystemExit: At least one comparison exceeds its regression allowance.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", choices=CASES)
    parser.add_argument("--base-python")
    parser.add_argument("--head-python")
    parser.add_argument("--base-ref", default="local")
    parser.add_argument("--head-ref", default="local")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path, default=Path("performance.json"))
    args = parser.parse_args()
    if sys.platform != "linux":
        parser.error("Run in Linux or WSL; peak RSS units and CI thresholds target Linux")
    if args.sample:
        print(json.dumps(sample(args.sample), allow_nan=False))
        return
    if not args.base_python or not args.head_python or args.repeats < 3:
        parser.error("Two installed environments and at least three repetitions are required")
    result = {
        "platform": platform.platform(),
        "python": platform.python_version(),
        "base_ref": args.base_ref,
        "head_ref": args.head_ref,
        "cases": {},
    }
    failures = []
    script = str(Path(__file__).resolve())
    for case in CASES:
        runs = {"base": [], "head": []}
        for repetition in range(args.repeats + 1):
            order = ("base", "head") if repetition % 2 == 0 else ("head", "base")
            for label in order:
                python = args.base_python if label == "base" else args.head_python
                raw = subprocess.check_output(
                    [python, script, "--sample", case], text=True, timeout=120
                )
                if repetition:  # Discard the first pair to warm the OS page cache.
                    runs[label].append(json.loads(raw))
        medians = {
            label: {key: statistics.median(r[key] for r in rows) for key in rows[0]}
            for label, rows in runs.items()
        }
        errors = regressions(medians["base"], medians["head"])
        failures.extend(f"{case}: {error}" for error in errors)
        result["cases"][case] = {"runs": runs, "medians": medians, "failures": errors}
        print(case, "FAIL" if errors else "PASS", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    lines = [
        "| Case | Time ratio | RSS ratio | Minimum quality, base → head | Result |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for name, data in result["cases"].items():
        b, h = data["medians"]["base"], data["medians"]["head"]
        lines.append(
            f"| {name} | {h['seconds'] / b['seconds']:.3f} | "
            f"{h['peak_rss_bytes'] / b['peak_rss_bytes']:.3f} | "
            f"{b['quality_min']:.5f} → {h['quality_min']:.5f} | "
            f"{'FAIL' if data['failures'] else 'PASS'} |"
        )
    args.output.with_suffix(".md").write_text(
        "\n".join([*lines, "", *failures]) + "\n", encoding="utf-8"
    )
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    main()
