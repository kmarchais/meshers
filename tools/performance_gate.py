# Copyright (c) 2026 Kevin Marchais

"""Compare installed base/head wheels on one configuration in alternating fresh processes."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import logging
import math
import os
import platform
import statistics
import sys
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

if sys.platform == "darwin":
    import ctypes

if not __package__ and __name__ == "__main__":
    from benchmark_hardware import hardware_snapshot
    from commands import CalledProcessError, TimeoutExpired, run_command
else:
    from tools.benchmark_hardware import hardware_snapshot
    from tools.commands import CalledProcessError, TimeoutExpired, run_command

_LOG = logging.getLogger(__name__)


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
SCHEMA_VERSION = 2
FIXTURE_VERSION = 5
# Versions 4 and 5 sample identical meshes. Version 5 changes the resource verdict.
MAPPED_CAP_VERSION = 4
PAIRED_RESOURCE_VERSION = 5
ResourceLimits = tuple[tuple[str, float, float], ...]
PRE_CALIBRATION_LIMITS: ResourceLimits = (
    ("seconds", 1.25, 0.02),
    ("peak_rss_bytes", 1.15, 8 * 1024**2),
)
RESOURCE_LIMITS_BY_VERSION: dict[int, ResourceLimits] = {
    2: PRE_CALIBRATION_LIMITS,
    3: PRE_CALIBRATION_LIMITS,
    4: PRE_CALIBRATION_LIMITS,
    PAIRED_RESOURCE_VERSION: (("seconds", 1.15, 0.01), ("peak_rss_bytes", 1.05, 4 * 1024**2)),
}
MIN_MEASURED_PAIRS = 3


class _CaseResult(TypedDict):
    runs: dict[str, list[dict[str, float]]]
    medians: dict[str, dict[str, float]] | None
    failures: list[str]


class _Summary(TypedDict):
    cases: dict[str, _CaseResult]


def _resource_failures(
    medians: Mapping[str, Mapping[str, float]],
    runs: Mapping[str, Sequence[Mapping[str, float]]],
    version: int,
) -> list[str]:
    limits = RESOURCE_LIMITS_BY_VERSION[version]
    base, head = medians["base"], medians["head"]
    if version < PAIRED_RESOURCE_VERSION:
        return [
            f"{key} increased: {base[key]:.6g} -> {head[key]:.6g}"
            for key, ratio, allowance in limits
            if head[key] > base[key] * ratio + allowance
        ]
    pairs = list(zip(runs["base"], runs["head"], strict=True))
    failures = []
    for key, ratio, allowance in limits:
        if any(not math.isfinite(row[key]) or row[key] < 0 for pair in pairs for row in pair):
            failures.append(f"invalid {key} sample")
            continue
        excess = statistics.median(h[key] - (b[key] * ratio + allowance) for b, h in pairs)
        if excess > 0:
            failures.append(f"{key} increased: median paired excess over limit {excess:.6g}")
    return failures


def _mesh_failures(base: Mapping[str, float], head: Mapping[str, float]) -> list[str]:
    failures = [
        f"invalid {key}"
        for key in base
        if not math.isfinite(base[key])
        or not math.isfinite(head[key])
        or base[key] < 0
        or head[key] < 0
    ]
    if failures:
        return failures
    failures.extend(
        f"{key} decreased: {base[key]:.6g} -> {head[key]:.6g}"
        for key in (f"quality_{name}" for name in QUALITY)
        if head[key] < base[key] * 0.99 - 1e-12
    )
    for key, tolerance in (("volume", 0.005), ("tetrahedra", 0.01), ("nodes", 0.01)):
        if abs(head[key] - base[key]) > abs(base[key]) * tolerance + 1e-12:
            failures.append(f"{key} changed beyond {tolerance:.1%}")
    if head["surface_error"] > base["surface_error"] * 1.05 + 1e-10:
        failures.append("surface error increased by more than 5%")
    if (
        "mapped_cap_error" in base
        and head["mapped_cap_error"] > base["mapped_cap_error"] * 1.05 + 1e-10
    ):
        failures.append("mapped cap error increased by more than 5%")
    return failures


def case_regressions(
    medians: Mapping[str, Mapping[str, float]],
    runs: Mapping[str, Sequence[Mapping[str, float]]],
    *,
    fixture_version: int = FIXTURE_VERSION,
) -> list[str]:
    """Return resource, geometry and per-sample quality failures for this fixture version."""
    base = medians["base"]
    failures = _resource_failures(medians, runs, fixture_version)
    failures.extend(_mesh_failures(base, medians["head"]))
    for index, row in enumerate(runs["head"]):
        for name in QUALITY:
            key = f"quality_{name}"
            if row[key] < base[key] * 0.99 - 1e-12:
                failures.append(f"sample {index + 1}: {key} decreased by more than 1%")
    return failures


def write_summary(result: _Summary, output: Path, failures: Sequence[str]) -> None:
    """Write the concise CI comparison beside the complete JSON artifact."""
    lines = [
        "Ratios summarize group medians. Resource limits use the median excess across paired runs.",
        "",
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
            f"{'FAIL' if data['failures'] else 'PASS'} |",
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
        "\n".join([*lines, "", *failures]) + "\n",
        encoding="utf-8",
    )


def configure_worker() -> None:
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
        ),
    )
    if sys.platform == "darwin":
        blas = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
        blas.BLASSetThreading.argtypes = [ctypes.c_uint]
        if blas.BLASSetThreading(1) != 0:
            msg = "Cannot select single-thread Accelerate"
            raise RuntimeError(msg)


def arguments() -> tuple[argparse.ArgumentParser, argparse.Namespace]:
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


def _worker_output(python: str, script: str, *options: str, timeout: int) -> str:
    """Return the JSON emitted by the selected wheel interpreter."""
    return run_command(python, script, *options, stdout_encoding="utf-8", timeout=timeout).stdout


def _collect_case(
    case: str,
    script: str,
    interpreters: Mapping[str, str],
    repeats: int,
) -> _CaseResult:
    runs: dict[str, list[dict[str, float]]] = {"base": [], "head": []}
    try:
        for repetition in range(repeats + 1):
            order = ("base", "head") if repetition % 2 == 0 else ("head", "base")
            for label in order:
                python = interpreters[label]
                raw = _worker_output(python, script, "--sample", case, timeout=120)
                if repetition:  # Discard the first pair to warm the OS page cache.
                    runs[label].append(json.loads(raw))
    except (CalledProcessError, TimeoutExpired) as error:
        reason = "timeout" if isinstance(error, TimeoutExpired) else f"exit {error.returncode}"
        message = f"{case}: sample failed ({reason})"
        return {"runs": runs, "failures": [message], "medians": None}
    medians = {
        label: {key: statistics.median(r[key] for r in rows) for key in rows[0]}
        for label, rows in runs.items()
    }
    errors = case_regressions(medians, runs)
    return {"runs": runs, "medians": medians, "failures": errors}


def compare(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Collect repeated measurements and fail when a regression exceeds its allowance.

    Raises:
        SystemExit: At least one comparison exceeds its regression allowance.
    """
    if not args.base_python or not args.head_python or args.repeats < MIN_MEASURED_PAIRS:
        parser.error("Two installed environments and at least three repetitions are required")
    summary: _Summary = {"cases": {}}
    result = {
        "schema_version": SCHEMA_VERSION,
        "fixture_version": FIXTURE_VERSION,
        "suite": args.suite,
        "recorded_at": dt.datetime.now(dt.UTC).isoformat(),
        **hardware_snapshot(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "worker_sha256": hashlib.sha256(
            Path(__file__).with_name("performance_worker.py").read_bytes(),
        ).hexdigest(),
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
        **summary,
    }
    failures = []
    script = str(Path(__file__).resolve())
    result["runtime"] = {
        label: json.loads(
            _worker_output(python, script, "--identity", timeout=30),
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
        row = _collect_case(
            case,
            script,
            {"base": args.base_python, "head": args.head_python},
            args.repeats,
        )
        summary["cases"][case] = row
        if row["medians"] is None:
            failures.extend(row["failures"])
            _LOG.error("%s", row["failures"][0])
        else:
            failures.extend(f"{case}: {error}" for error in row["failures"])
            _LOG.info("%s %s", case, "FAIL" if row["failures"] else "PASS")
    result["complete"] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_summary(summary, args.output, failures)
    if failures:
        raise SystemExit("\n".join(failures))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    parser, args = arguments()
    if sys.platform not in {"linux", "darwin", "win32"}:
        parser.error("OS peak-RSS reporting supports Linux, macOS and Windows")
    configure_worker()
    if args.identity or args.sample:
        if __package__:
            from tools.performance_worker import runtime_identity, sample
        else:
            from performance_worker import runtime_identity, sample

        json.dump(
            runtime_identity() if args.identity else sample(args.sample),
            sys.stdout,
            allow_nan=False,
        )
        sys.stdout.write("\n")
    else:
        compare(parser, args)
