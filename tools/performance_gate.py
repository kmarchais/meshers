"""Compare installed base/head wheels on one Linux runner in alternating fresh processes."""

import argparse
import datetime
import hashlib
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

CASES = (
    "compiled_gyroid",
    "callback_gyroid",
    "graded",
    "mapped_intersection",
    "schwarz_p",
    "spinodoid",
    "smooth_lattice",
)
RELEASE_CASES = (
    *CASES,
    "adaptive_gyroid",
    "adaptive_spinodoid",
    "adaptive_lattice",
    "gyroid_tet10",
)
QUALITY = ("min", "p05", "p10", "median", "mean", "p90", "max")
CAP = 250_000


def cell_geometry(vertices):
    """Integrate Tet4/Tet10 volume and sample Jacobians independently of Meshers.

    Returns:
        Integrated physical volumes per element.

    Raises:
        RuntimeError: A sampled Jacobian is nonpositive.
    """
    import numpy as np

    corner = np.linalg.det(vertices[:, 1:4] - vertices[:, :1])
    if vertices.shape[1] == 4:
        return corner / 6
    edges = np.array([(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)])
    bary = np.vstack([
        np.full((1, 4), 0.25),
        np.full((4, 4), 1 / 6) + np.eye(4) / 3,
        np.eye(4),
        (np.eye(4)[edges[:, 0]] + np.eye(4)[edges[:, 1]]) / 2,
    ])
    dl = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    determinants = []
    for location in bary:
        derivative = np.vstack([
            (4 * location - 1)[:, None] * dl,
            4
            * (
                location[edges[:, 0], None] * dl[edges[:, 1]]
                + location[edges[:, 1], None] * dl[edges[:, 0]]
            ),
        ])
        determinant = np.linalg.det(np.einsum("tni,nj->tij", vertices, derivative))
        if not np.isfinite(determinant).all() or np.any(determinant <= 0):
            raise RuntimeError("Nonpositive sampled curved Jacobian")
        determinants.append(determinant)
    return np.array([-2 / 15, *([3 / 40] * 4)]) @ np.array(determinants[:5])


def sample(case):  # ruff: ignore[too-many-locals] - keep the timing boundary and independent checks together.
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

    def spinodoid(x, y, z):
        # Fixed periodic Fourier fixture; not a stochastic field generator.
        return sum(
            a * np.sin(2 * np.pi * (i * x + j * y + k * z) + phase)
            for i, j, k, a, phase in [
                (1, 1, 0, 0.6, 0.2),
                (0, 1, 1, 0.5, 0.7),
                (1, 0, -1, 0.45, 1.1),
                (2, 1, 0, 0.25, 0.4),
                (0, 2, -1, 0.2, 0.9),
                (1, -1, 2, 0.15, 1.4),
            ]
        )

    def lattice(x, y, z):
        a, b, c = (np.sin(np.pi * v) ** 2 for v in (x, y, z))

        def union(f, g):
            return np.minimum(f, g) - np.maximum(0.08 - np.abs(f - g), 0) ** 2 / 0.32

        return union(union(a + b - 0.24, b + c - 0.24), a + c - 0.24)

    options = {"cells": 24, "optimize_passes": 2, "geometry_tolerance": 0.1, "max_tetrahedra": CAP}
    if case.startswith("adaptive_"):
        options.update(
            adapt=True,
            min_length=0.005,
            max_length=0.1,
            adapt_tolerance=0.002,
            geometry_tolerance=0.003,
            gradation=0.6,
        )
    if case == "gyroid_tet10":
        options["order"] = 2
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
        field = (
            spinodoid
            if "spinodoid" in case
            else lattice
            if "lattice" in case
            else (
                lambda x, y, z: (
                    np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y) + np.cos(2 * np.pi * z)
                )
            )
            if case == "schwarz_p"
            else gyroid
        )
        if case == "graded":

            def field(x, y, z):
                return gyroid(x, y, z) / (0.35 + 0.15 * x)

        mesh = meshers.generate(
            field,
            band=(
                (-1, 1)
                if case == "graded"
                else None
                if case
                in {
                    "schwarz_p",
                    "spinodoid",
                    "smooth_lattice",
                    "adaptive_spinodoid",
                    "adaptive_lattice",
                }
                else (-0.5, 0.5)
            ),
            periodic=(False,) * 3 if case == "graded" else (True,) * 3,
            compile=case != "callback_gyroid",
            threads=1,
            batch_size=4096,
            **options,
        )
    elapsed = time.perf_counter() - start
    # Capture OS high-water RSS before allocating arrays for independent checks.
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (
        1 if sys.platform == "darwin" else 1024
    )
    if max(len(mesh.points), len(mesh.tetrahedra)) > CAP:
        raise RuntimeError("Benchmark node/tetrahedron cap exceeded")
    points = mesh.points[mesh.tetrahedra]
    determinants = np.linalg.det(points[:, 1:4] - points[:, :1])
    volumes = cell_geometry(points)
    if not len(points) or not np.all(np.isfinite(mesh.points)) or not np.all(determinants > 0):
        raise RuntimeError("Nonempty, finite, positively oriented mesh required")
    edges = sum(
        np.sum((points[:, i] - points[:, j]) ** 2, axis=1)
        for i in range(4)
        for j in range(i + 1, 4)
    )
    quality = np.sqrt(432 * determinants**2 / edges**3)
    surface = mesh.surface[:, :3]
    boundary_edges = np.sort(
        np.concatenate([surface[:, [0, 1]], surface[:, [1, 2]], surface[:, [2, 0]]]),
        axis=1,
    )
    _, incidence = np.unique(boundary_edges, axis=0, return_counts=True)
    if not np.all(incidence == 2):
        raise RuntimeError("Closed manifold boundary required")
    if case not in {"graded", "mapped_intersection"}:
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
        "nodes": len(mesh.points),
        "volume": float(volumes.sum()),
        **{
            f"quality_{name}": float(value)
            for name, value in zip(
                QUALITY,
                [
                    quality.min(),
                    np.quantile(quality, 0.05),
                    np.quantile(quality, 0.1),
                    np.median(quality),
                    quality.mean(),
                    np.quantile(quality, 0.9),
                    quality.max(),
                ],
                strict=True,
            )
        },
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
    for key in (f"quality_{name}" for name in QUALITY):
        if head[key] < base[key] * 0.99 - 1e-12:
            failures.append(f"{key} decreased: {base[key]:.6g} -> {head[key]:.6g}")
    for key, tolerance in (("volume", 0.005), ("tetrahedra", 0.01), ("nodes", 0.01)):
        if abs(head[key] - base[key]) > abs(base[key]) * tolerance + 1e-12:
            failures.append(f"{key} changed beyond {tolerance:.1%}")
    if head["surface_error"] > base["surface_error"] * 1.05 + 1e-10:
        failures.append("surface error increased by more than 5%")
    return failures


def runtime_identity():
    """Return versions and the actual loaded native binary fingerprint."""
    import importlib
    import importlib.metadata

    native_module = importlib.import_module("meshers._meshers")
    import numpy as np

    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "meshers": importlib.metadata.version("meshers"),
        "native_sha256": hashlib.sha256(Path(native_module.__file__).read_bytes()).hexdigest(),
    }


def case_regressions(base, head, head_runs):
    """Return resource-median and every-repetition quality regression failures."""
    failures = regressions(base, head)
    for index, row in enumerate(head_runs):
        for name in QUALITY:
            key = f"quality_{name}"
            if row[key] < base[key] * 0.99 - 1e-12:
                failures.append(f"sample {index + 1}: {key} decreased by more than 1%")
    return failures


def write_summary(result, output, failures):
    """Write the concise CI comparison beside the complete JSON artifact."""
    lines = [
        "| Case | Time ratio | RSS ratio | Minimum quality, base → head | Result |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for name, data in result["cases"].items():
        if data["medians"] is None:
            lines.append(f"| {name} | — | — | — | FAIL: sample failed |")
            continue
        b, h = data["medians"]["base"], data["medians"]["head"]
        lines.append(
            f"| {name} | {h['seconds'] / b['seconds']:.3f} | "
            f"{h['peak_rss_bytes'] / b['peak_rss_bytes']:.3f} | "
            f"{b['quality_min']:.5f} → {h['quality_min']:.5f} | "
            f"{'FAIL' if data['failures'] else 'PASS'} |"
        )
    lines += [
        "",
        "Quality statistics, base → candidate; linear corners for Tet10.",
        "",
        "| Case | Min | P05 | P10 | Median | Mean | P90 | Max |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for name, data in result["cases"].items():
        if data["medians"]:
            b, h = data["medians"]["base"], data["medians"]["head"]
            values = [f"{b[f'quality_{q}']:.4f} → {h[f'quality_{q}']:.4f}" for q in QUALITY]
            lines.append(f"| {name} | " + " | ".join(values) + " |")
    output.with_suffix(".md").write_text(
        "\n".join([*lines, "", *failures]) + "\n", encoding="utf-8"
    )


def configure_worker():
    """Select one numerical thread before importing NumPy in each fresh worker.

    Raises:
        RuntimeError: Single-thread Accelerate cannot be selected.
    """
    os.environ.update(
        dict.fromkeys(
            (
                "OPENBLAS_NUM_THREADS",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "RAYON_NUM_THREADS",
            ),
            "1",
        )
    )
    if sys.platform == "darwin":
        import ctypes

        blas = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
        blas.BLASSetThreading.argtypes = [ctypes.c_uint]
        if blas.BLASSetThreading(1) != 0:
            raise RuntimeError("Cannot select single-thread Accelerate")


def arguments():
    """Return the benchmark parser and its parsed command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", choices=RELEASE_CASES)
    parser.add_argument("--identity", action="store_true")
    parser.add_argument("--suite", choices=["pr", "release"], default="pr")
    parser.add_argument("--base-python")
    parser.add_argument("--head-python")
    parser.add_argument("--base-ref", default="local")
    parser.add_argument("--head-ref", default="local")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path, default=Path("performance.json"))
    return parser, parser.parse_args()


def main():
    """Collect repeated measurements and fail when a regression exceeds its allowance.

    Raises:
        SystemExit: At least one comparison exceeds its regression allowance.
    """
    parser, args = arguments()
    if sys.platform not in {"linux", "darwin"}:
        parser.error("OS peak-RSS reporting supports Linux and macOS")
    configure_worker()
    if args.identity:
        print(json.dumps(runtime_identity()))
        return
    if args.sample:
        print(json.dumps(sample(args.sample), allow_nan=False))
        return
    if not args.base_python or not args.head_python or args.repeats < 3:
        parser.error("Two installed environments and at least three repetitions are required")
    result = {
        "schema_version": 2,
        "fixture_version": 2,
        "suite": args.suite,
        "recorded_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "runner": os.environ.get("RUNNER_NAME", platform.node()),
        "machine": platform.machine(),
        "cpu_model": (
            subprocess.check_output(
                ["/usr/sbin/sysctl", "-n", "machdep.cpu.brand_string"], text=True
            ).strip()
            if sys.platform == "darwin"
            else platform.processor()
        ),
        "cpu_count": os.cpu_count(),
        "platform": platform.platform(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "threads": 1,
        "warmup_pairs": 1,
        "measured_pairs": args.repeats,
        "quality_basis": (
            "linear corners, including Tet10; volume integrates Tet10 mapping; "
            "surface_error is the linear parent diagnostic"
        ),
        "scope": (
            "Regression fixtures, not the matched-accuracy Microgen study or FEA certification"
        ),
        "python": platform.python_version(),
        "base_ref": args.base_ref,
        "head_ref": args.head_ref,
        "cases": {},
    }
    failures = []
    script = str(Path(__file__).resolve())
    result["runtime"] = {
        label: json.loads(
            subprocess.check_output([python, script, "--identity"], text=True, timeout=30)
        )
        for label, python in (("base", args.base_python), ("head", args.head_python))
    }
    if any(
        result["runtime"]["base"][key] != result["runtime"]["head"][key]
        for key in ("python", "numpy")
    ):
        parser.error("Both environments must use the same Python and NumPy versions")
    result["same_native_binary"] = (
        result["runtime"]["base"]["native_sha256"] == result["runtime"]["head"]["native_sha256"]
    )

    for case in RELEASE_CASES if args.suite == "release" else CASES:
        runs = {"base": [], "head": []}
        try:
            for repetition in range(args.repeats + 1):
                order = ("base", "head") if repetition % 2 == 0 else ("head", "base")
                for label in order:
                    python = args.base_python if label == "base" else args.head_python
                    raw = subprocess.check_output(
                        [python, script, "--sample", case], text=True, timeout=120
                    )
                    if repetition:  # Discard the first pair to warm the OS page cache.
                        runs[label].append(json.loads(raw))
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            message = f"{case}: sample failed: {error}"
            failures.append(message)
            result["cases"][case] = {"runs": runs, "failures": [message], "medians": None}
            print(message, flush=True)
            continue
        medians = {
            label: {key: statistics.median(r[key] for r in rows) for key in rows[0]}
            for label, rows in runs.items()
        }
        errors = case_regressions(medians["base"], medians["head"], runs["head"])
        failures.extend(f"{case}: {error}" for error in errors)
        result["cases"][case] = {"runs": runs, "medians": medians, "failures": errors}
        print(case, "FAIL" if errors else "PASS", flush=True)
    result["complete"] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_summary(result, args.output, failures)
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    main()
