# Copyright (c) 2026 Kevin Marchais

"""Archive real-shaped records without hiding failures or replacing history."""

import copy
import json
import xml.etree.ElementTree as ET  # ruff: ignore[suspicious-xml-etree-import] - own generated SVG.

import pytest

from tools import benchmark_history as history
from tools import performance_gate as gate
from tools.benchmark_hardware import hardware_snapshot


def measurement():
    row = dict(
        seconds=1.0,
        peak_rss_bytes=1000000.0,
        tetrahedra=1000,
        nodes=300,
        volume=0.5,
        surface_error=0.01,
        **{f"quality_{q}": 0.7 for q in gate.QUALITY},
    )
    case = {
        "runs": {k: [copy.deepcopy(row) for _ in range(5)] for k in ("base", "head")},
        "medians": {k: copy.deepcopy(row) for k in ("base", "head")},
        "failures": [],
    }
    return {
        "schema_version": 2,
        "fixture_version": 2,
        "complete": True,
        "suite": "release",
        "measured_pairs": 5,
        "base_ref": "a" * 40,
        "head_ref": "a" * 40,
        "script_sha256": "c" * 64,
        "cases": {k: copy.deepcopy(case) for k in gate.RELEASE_CASES},
        "same_native_binary": True,
        "platform": "test host",
        "machine": "arm64",
        "cpu_count": 10,
        "recorded_at": "2026-11-12T00:00:00Z",
    }


def test_published_configuration_excludes_computer_name(monkeypatch):
    monkeypatch.setenv("RUNNER_NAME", "private-identifier")
    monkeypatch.setenv("COMPUTERNAME", "private-identifier")
    specs = hardware_snapshot()
    assert set(specs) == {"cpu_model", "cpu_count", "ram_bytes", "machine", "platform"}
    assert "private-identifier" not in json.dumps(specs)
    assert specs["cpu_count"] > 0
    assert specs["ram_bytes"] > 0


def test_archive_plots_and_immutable_history(tmp_path):
    source = tmp_path / "input.json"
    data = measurement()
    slow = data["cases"]["compiled_gyroid"]
    for label, runs in slow["runs"].items():
        for row in [*runs, slow["medians"][label]]:
            row["seconds"] *= 2
    source.write_text(json.dumps({**data, "runner": "private-identifier"}))
    directory = tmp_path / "benchmarks"
    history.archive(source, "v0.2.0-rc-reference", directory)
    page = (directory / "v0.2.0-rc-reference.md").read_text()
    assert "same native binary" in page
    assert "quality_mean" in page
    assert "quality_p90" in page
    assert "REFERENCE" in (tmp_path / "benchmark-history.md").read_text()
    chart = ET.parse(directory / "v0.2.0-rc-reference.svg")  # ruff: ignore[suspicious-xml-element-tree-usage] - own generated SVG.
    bars = chart.findall(".//{*}rect[@height='15']")
    widths = {float(bar.attrib["width"]) for bar in bars}
    assert len(widths) == 2
    assert max(widths) == 2 * min(widths)
    frozen = (directory / "v0.2.0-rc-reference.json").read_bytes()
    assert b"private-identifier" not in frozen
    with pytest.raises(FileExistsError):
        history.archive(source, "v0.2.0-rc-reference", directory)
    assert (directory / "v0.2.0-rc-reference.json").read_bytes() == frozen


@pytest.mark.parametrize(
    ("corruption", "error"),
    [
        ("missing", "Missing measured samples"),
        ("nonfinite", "Invalid numeric measurement"),
        ("verdict", "Regression verdict does not match"),
        ("median", "Summary does not match"),
    ],
)
def test_refuses_incomplete_or_false_success(corruption, error):
    data = measurement()
    case = data["cases"]["compiled_gyroid"]
    if corruption == "missing":
        case["runs"]["head"].pop()
    elif corruption == "nonfinite":
        case["runs"]["head"][0]["seconds"] = float("nan")
    elif corruption == "verdict":
        case["failures"] = ["pretend failure"]
    else:
        case["medians"]["head"]["seconds"] = 99.0
    with pytest.raises(ValueError, match=error):
        history.validate(data)


def test_preserves_real_regression_and_failed_worker(tmp_path):
    data = measurement()
    data["same_native_binary"] = False
    case = data["cases"]["compiled_gyroid"]
    for row in [*case["runs"]["head"], case["medians"]["head"]]:
        row["quality_min"] = 0.4
    case["failures"] = gate.case_regressions(
        case["medians"]["base"],
        case["medians"]["head"],
        case["runs"]["head"],
    )
    data["cases"]["adaptive_gyroid"] = {
        "runs": {"base": [], "head": []},
        "medians": None,
        "failures": ["sample failed"],
    }
    source = tmp_path / "input.json"
    source.write_text(json.dumps(data))
    history.archive(source, "v0.2.0-failed-control", tmp_path / "benchmarks")
    assert "FAIL" in (tmp_path / "benchmark-history.md").read_text()
    page = (tmp_path / "benchmarks/v0.2.0-failed-control.md").read_text()
    assert "quality_min decreased" in page
    assert "sample failed" in page


def test_mapped_cap_measurements_are_required_in_new_fixture_records():
    data = measurement()
    data["fixture_version"] = gate.FIXTURE_VERSION
    with pytest.raises(ValueError, match="Invalid numeric measurement"):
        history.validate(data)
    case = data["cases"]["mapped_intersection"]
    for label, rows in case["runs"].items():
        for row in [*rows, case["medians"][label]]:
            row["mapped_cap_error"] = 0.001
    history.validate(data)
