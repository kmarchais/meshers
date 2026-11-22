# Copyright (c) 2026 Kevin Marchais

"""Archive immutable benchmark JSON and generate documentation pages and plots."""

import argparse
import hashlib
import html
import json
import math
import re
import statistics
from pathlib import Path

from tools.performance_gate import (
    CASES,
    FIXTURE_VERSION,
    MIN_MEASURED_PAIRS,
    QUALITY,
    RELEASE_CASES,
    SCHEMA_VERSION,
    case_regressions,
)

SUPPORTED_FIXTURE_VERSIONS = {2, 3, FIXTURE_VERSION}


def validate(data: dict) -> None:
    """Reject incomplete measurements rather than publishing an apparent success.

    Raises:
        ValueError: The record is incomplete or contains invalid measurements.
    """
    if data.get("schema_version") != SCHEMA_VERSION or not data.get("complete"):
        msg = f"Expected a complete schema-{SCHEMA_VERSION} benchmark"
        raise ValueError(msg)
    if data.get("fixture_version") not in SUPPORTED_FIXTURE_VERSIONS:
        msg = "Unsupported benchmark fixture version"
        raise ValueError(msg)
    expected = {"pr": CASES, "release": RELEASE_CASES}.get(data.get("suite"), ())
    if (
        set(data.get("cases", {})) != set(expected)
        or data.get("measured_pairs", 0) < MIN_MEASURED_PAIRS
    ):
        msg = "Incomplete fixture suite or insufficient repetitions"
        raise ValueError(msg)
    for ref in ("base_ref", "head_ref", "script_sha256"):
        if not re.fullmatch(r"[a-f0-9]{40,64}", data.get(ref, "")):
            msg = f"Missing immutable source identity: {ref}"
            raise ValueError(msg)
    worker_digest = data.get("worker_sha256", "")
    if "worker_sha256" in data and not (
        isinstance(worker_digest, str) and re.fullmatch(r"[a-f0-9]{64}", worker_digest)
    ):
        msg = "Invalid benchmark worker SHA-256"
        raise ValueError(msg)
    for name, case in data["cases"].items():
        needs_cap = name == "mapped_intersection" and data["fixture_version"] == FIXTURE_VERSION
        validate_case(case, data["measured_pairs"], needs_cap=needs_cap)


def validate_case(case: dict, measured_pairs: int, *, needs_cap: bool) -> None:
    """Check raw samples, summary medians and the recorded regression verdict.

    Raises:
        ValueError: Measurements or their verdict are inconsistent.
    """
    if case["medians"] is None:
        if not case["failures"]:
            msg = "Missing samples cannot pass"
            raise ValueError(msg)
        return
    for label in ("base", "head"):
        rows = case["runs"][label]
        if len(rows) != measured_pairs:
            msg = "Missing measured samples"
            raise ValueError(msg)
        for row in [*rows, case["medians"][label]]:
            required = (
                "seconds",
                "peak_rss_bytes",
                "tetrahedra",
                "nodes",
                "volume",
                "surface_error",
                *(("mapped_cap_error",) if needs_cap else ()),
                *(f"quality_{q}" for q in QUALITY),
            )
            if any(k not in row or not math.isfinite(row[k]) or row[k] < 0 for k in required):
                msg = "Invalid numeric measurement"
                raise ValueError(msg)
        if any(
            case["medians"][label][key] != statistics.median(r[key] for r in rows)
            for key in required
        ):
            msg = "Summary does not match the raw measurements"
            raise ValueError(msg)
    if case["failures"] != case_regressions(
        case["medians"]["base"],
        case["medians"]["head"],
        case["runs"]["head"],
    ):
        msg = "Regression verdict does not match the measurements"
        raise ValueError(msg)


def timing_plot(data: dict) -> str:
    """Return an accessible SVG of within-run base and candidate generation times."""
    cases = [(name, c["medians"]) for name, c in data["cases"].items() if c["medians"]]
    svg = [
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 {55 + len(cases) * 58}" '
            'role="img"><title>Within-run generation times in seconds</title>'
        ),
        (
            "<style>:root{color:#24413b;--base:#ad633b;--head:#0b8072}"
            "@media(prefers-color-scheme:dark){:root{color:#dcebe6;--base:#e5a679;--head:#76d6be}}"
            "</style>"
        ),
        '<g font-family="system-ui,sans-serif" font-size="13" fill="currentColor">',
        '<text x="20" y="25">Median seconds · upper bar: base · lower bar: candidate</text>',
    ]
    scale = max((m[k]["seconds"] for _, m in cases for k in ("base", "head")), default=0)
    for index, (name, medians) in enumerate(cases):
        y = 55 + index * 58
        svg.append(f'<text x="20" y="{y + 13}">{html.escape(name)}</text>')
        for offset, label, color in [(0, "base", "var(--base)"), (21, "head", "var(--head)")]:
            seconds = medians[label]["seconds"]
            width = 400 * seconds / scale if scale else 0
            svg.append(
                f'<rect x="220" y="{y + offset}" width="{width:.2f}" height="15" '
                f'fill="{color}"/><text x="{230 + width:.2f}" '
                f'y="{y + offset + 12}">{seconds:.4f} s</text>',
            )
    return "".join([*svg, "</g></svg>"])


def render(directory: Path) -> None:
    """Regenerate navigable Markdown and SVGs from the immutable JSON records."""
    records = []
    for path in sorted(directory.glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        data, label = record["data"], record["label"]
        validate(data)
        failed = any(c["failures"] for c in data["cases"].values())
        status = (
            "FAIL" if failed else "REFERENCE" if data["base_ref"] == data["head_ref"] else "PASS"
        )
        records.append((label, data, status))
        svg = path.with_suffix(".svg")
        svg.write_text(timing_plot(data), encoding="utf-8")
        lines = [
            f"# {label}",
            "",
            f"Result: **{status}**. Suite: `{data['suite']}`.",
            "",
            f"Source: `{data['base_ref']}` → `{data['head_ref']}`.",
            "",
            (
                f"Host: `{data['platform']}`, `{data['machine']}`, "
                f"{data['cpu_count']} CPUs. One numerical thread."
            ),
            "",
            (
                f"Recorded: {data['recorded_at']}. {data['measured_pairs']} measured fresh-process "
                "pairs after one discarded warmup pair."
            ),
            "",
            (
                f"[Raw samples, versions and native hashes]({path.name}) · "
                "[All records](../benchmark-history.md)"
            ),
            "",
        ]
        for name, case in data["cases"].items():
            lines.extend(f"- **{name}** — {error}" for error in case["failures"])
        lines += [
            "",
            f"![Median generation times in seconds on one common scale]({svg.name})",
            "",
            (
                "For Tet10, quality and sampled surface error describe the linear parent. "
                "Volume integrates the curved element."
            ),
            "",
        ]
        if data["same_native_binary"]:
            lines += [
                (
                    "**Both sides loaded the same native binary.** Timings still include Python "
                    "preparation. Inspect the source identities before interpreting a speedup "
                    "as a native-code improvement."
                ),
                "",
            ]
        lines += ['<details markdown="1"><summary>Median metrics by case</summary>', ""]
        for name, case in data["cases"].items():
            if case["medians"]:
                lines += [f"### {name}", ""]
                b, h = (case["medians"][k] for k in ("base", "head"))
                lines += ["| Metric | Base | Candidate |", "| --- | ---: | ---: |"]
                for key in (
                    "seconds",
                    "peak_rss_bytes",
                    "tetrahedra",
                    "nodes",
                    "volume",
                    "surface_error",
                    *(("mapped_cap_error",) if "mapped_cap_error" in b else ()),
                    *(f"quality_{q}" for q in QUALITY),
                ):
                    lines.append(f"| {key} | {b[key]:.6g} | {h[key]:.6g} |")
                lines.append("")
        lines.append("</details>")
        path.with_suffix(".md").write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    records.sort(key=lambda r: r[1]["recorded_at"])
    latest_plot = (
        [f"![Latest within-run timing comparison](benchmarks/{records[-1][0]}.svg)", ""]
        if records
        else []
    )
    index = [
        "---",
        "icon: lucide/history",
        "---",
        "",
        "# Benchmark history",
        "",
        "Open a record for its timing plot, quality gates and raw samples.",
        "",
        "Compare times within a run. A REFERENCE compares one source revision to itself.",
        "",
        *latest_plot,
        "| Record | Suite | Result | Platform |",
        "| --- | --- | --- | --- |",
    ]
    index += [
        f"| [{label}](benchmarks/{label}.md) | {d['suite']} | {status} | "
        f"{d['platform']} / {d['machine']} |"
        for label, d, status in records
    ]
    index += [
        "",
        (
            "See [regression protocol and release instructions](performance.md) "
            "for thresholds and how to add a release record."
        ),
    ]
    (directory.parent / "benchmark-history.md").write_text(
        "\n".join(index) + "\n",
        encoding="utf-8",
    )


def archive(source: Path, label: str, directory: Path) -> None:
    """Archive a measured run without overwriting any previous release evidence.

    Raises:
        ValueError: The label or input is invalid.
    """
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}", label):
        msg = "Use a plain version or release-candidate label"
        raise ValueError(msg)
    payload = source.read_bytes()
    data = json.loads(payload)
    data.pop("runner", None)  # Early records included the CI or local hostname.
    validate(data)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{label}.json"
    with path.open("x", encoding="utf-8") as file:
        json.dump(
            {"label": label, "input_sha256": hashlib.sha256(payload).hexdigest(), "data": data},
            file,
            indent=2,
            allow_nan=False,
        )
        file.write("\n")
    render(directory)


def main() -> None:
    """Archive a release measurement and rebuild its documentation history."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--directory", type=Path, default=Path("docs/site/benchmarks"))
    args = parser.parse_args()
    archive(args.input, args.label, args.directory)


if __name__ == "__main__":
    main()
