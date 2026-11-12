"""Archive immutable benchmark JSON and generate documentation pages and plots."""

import argparse
import hashlib
import html
import json
import math
import re
import statistics
from pathlib import Path

from tools.performance_gate import CASES, QUALITY, RELEASE_CASES, case_regressions


def validate(data):
    """Reject incomplete measurements rather than publishing an apparent success.

    Raises:
        ValueError: The record is incomplete or contains invalid measurements.
    """
    if data.get("schema_version") != 2 or not data.get("complete"):
        raise ValueError("Expected a complete schema-2 benchmark")
    expected = {"pr": CASES, "release": RELEASE_CASES}.get(data.get("suite"), ())
    if set(data.get("cases", {})) != set(expected) or data.get("measured_pairs", 0) < 3:
        raise ValueError("Incomplete fixture suite or insufficient repetitions")
    for ref in ("base_ref", "head_ref", "script_sha256"):
        if not re.fullmatch(r"[a-f0-9]{40,64}", data.get(ref, "")):
            raise ValueError(f"Missing immutable source identity: {ref}")
    for case in data["cases"].values():
        validate_case(case, data["measured_pairs"])


def validate_case(case, measured_pairs):
    """Check raw samples, summary medians and the recorded regression verdict.

    Raises:
        ValueError: Measurements or their verdict are inconsistent.
    """
    if case["medians"] is None:
        if not case["failures"]:
            raise ValueError("Missing samples cannot pass")
        return
    for label in ("base", "head"):
        rows = case["runs"][label]
        if len(rows) != measured_pairs:
            raise ValueError("Missing measured samples")
        for row in [*rows, case["medians"][label]]:
            required = (
                "seconds",
                "peak_rss_bytes",
                "tetrahedra",
                "nodes",
                "volume",
                "surface_error",
                *(f"quality_{q}" for q in QUALITY),
            )
            if any(not math.isfinite(row[k]) or row[k] < 0 for k in required):
                raise ValueError("Invalid numeric measurement")
        if any(
            case["medians"][label][key] != statistics.median(r[key] for r in rows)
            for key in required
        ):
            raise ValueError("Summary does not match the raw measurements")
    if case["failures"] != case_regressions(
        case["medians"]["base"], case["medians"]["head"], case["runs"]["head"]
    ):
        raise ValueError("Regression verdict does not match the measurements")


def timing_plot(data):
    """Return an accessible SVG of within-run base and candidate generation times."""
    cases = [(name, c["medians"]) for name, c in data["cases"].items() if c["medians"]]
    svg = [
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 {55 + len(cases) * 58}" '
            'role="img"><title>Within-run generation times in seconds</title>'
        ),
        '<rect width="100%" height="100%" fill="#0b202b"/>',
        '<g font-family="sans-serif" font-size="13" fill="#eef5f7">',
        '<text x="20" y="25">Base (blue) and candidate (teal), same host and settings</text>',
    ]
    for index, (name, medians) in enumerate(cases):
        scale = max(medians[k]["seconds"] for k in ("base", "head"))
        y = 55 + index * 58
        svg.append(f'<text x="20" y="{y + 13}">{html.escape(name)}</text>')
        for offset, label, color in [(0, "base", "#8eafff"), (21, "head", "#52d5b8")]:
            seconds = medians[label]["seconds"]
            width = 400 * seconds / scale if scale else 0
            svg.append(
                f'<rect x="220" y="{y + offset}" width="{width:.2f}" height="15" '
                f'fill="{color}"/><text x="{230 + width:.2f}" '
                f'y="{y + offset + 12}">{seconds:.4f} s</text>'
            )
    return "".join([*svg, "</g></svg>"])


def render(directory):
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
            "Absolute times from different hosts are not a release-to-release speed ranking.",
            "",
            f"![Within-run generation times]({svg.name})",
            "",
            (
                "Quality and sampled surface error describe linear corners/parent surfaces, "
                "including Tet10. Tet10 volumes integrate the curved mapping and its Jacobians "
                "are sampled for positivity. These are regression checks, not geometry or "
                "FEA certification."
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
        for name, case in data["cases"].items():
            lines += [f"## {name}", ""]
            if case["medians"]:
                b, h = (case["medians"][k] for k in ("base", "head"))
                lines += ["| Metric | Base | Candidate |", "| --- | ---: | ---: |"]
                for key in (
                    "seconds",
                    "peak_rss_bytes",
                    "tetrahedra",
                    "nodes",
                    "volume",
                    "surface_error",
                    *(f"quality_{q}" for q in QUALITY),
                ):
                    lines.append(f"| {key} | {b[key]:.6g} | {h[key]:.6g} |")
                lines.append("")
            lines += [f"- {error}" for error in case["failures"]]
        path.with_suffix(".md").write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    records.sort(key=lambda r: r[1]["recorded_at"])
    index = [
        "# Benchmark history",
        "",
        (
            "Versioned measurements from fixed regression fixtures. Failed records remain visible. "
            "Release candidates and same-binary references are explicitly labelled."
        ),
        "",
        (
            "Compare base/candidate measurements within each run. Different machines, suite "
            "versions or settings cannot establish an absolute performance trend. Raw JSON "
            "retains every measured sample and source/native identity."
        ),
        "",
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
        "\n".join(index) + "\n", encoding="utf-8"
    )


def archive(source, label, directory):
    """Archive a measured run without overwriting any previous release evidence.

    Raises:
        ValueError: The label or input is invalid.
    """
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}", label):
        raise ValueError("Use a plain version or release-candidate label")
    payload = source.read_bytes()
    data = json.loads(payload)
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


def main():
    """Archive a release measurement and rebuild its documentation history."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--directory", type=Path, default=Path("docs/site/benchmarks"))
    args = parser.parse_args()
    archive(args.input, args.label, args.directory)


if __name__ == "__main__":
    main()
