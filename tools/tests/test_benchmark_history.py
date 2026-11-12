"""Archive real-shaped records without hiding failures or replacing history."""

import copy
import json
import xml.etree.ElementTree as ET  # ruff: ignore[suspicious-xml-etree-import] - own generated SVG.

import pytest

from tools import benchmark_history as history
from tools import performance_gate as gate


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


def test_archive_plots_and_immutable_history(tmp_path):
    source = tmp_path / "input.json"
    source.write_text(json.dumps(measurement()))
    directory = tmp_path / "benchmarks"
    history.archive(source, "v0.2.0-rc-reference", directory)
    page = (directory / "v0.2.0-rc-reference.md").read_text()
    assert "same native binary" in page
    assert "quality_mean" in page
    assert "quality_p90" in page
    assert "REFERENCE" in (tmp_path / "benchmark-history.md").read_text()
    ET.parse(directory / "v0.2.0-rc-reference.svg")  # ruff: ignore[suspicious-xml-element-tree-usage] - own generated SVG.
    frozen = (directory / "v0.2.0-rc-reference.json").read_bytes()
    with pytest.raises(FileExistsError):
        history.archive(source, "v0.2.0-rc-reference", directory)
    assert (directory / "v0.2.0-rc-reference.json").read_bytes() == frozen


@pytest.mark.parametrize("corruption", ["missing", "nonfinite", "verdict", "median"])
def test_refuses_incomplete_or_false_success(corruption):
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
    with pytest.raises(ValueError):
        history.validate(data)


def test_preserves_real_regression_and_failed_worker(tmp_path):
    data = measurement()
    data["same_native_binary"] = False
    case = data["cases"]["compiled_gyroid"]
    for row in [*case["runs"]["head"], case["medians"]["head"]]:
        row["quality_min"] = 0.4
    case["failures"] = gate.case_regressions(
        case["medians"]["base"], case["medians"]["head"], case["runs"]["head"]
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
