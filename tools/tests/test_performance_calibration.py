# Copyright (c) 2026 Kevin Marchais

"""Replay measured controls against the production resource gate."""

import json
from pathlib import Path
from statistics import median

import pytest

from tools import benchmark_history as history
from tools import performance_gate as gate

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / "benchmarks/gate-calibration-20261121"


@pytest.mark.parametrize(
    "path",
    [STUDY / config / "same-commit.json" for config in ("linux-epyc-7763", "linux-ryzen-7-7730u")],
    ids=lambda p: p.parent.name,
)
def test_identical_source_passes_in_both_directions(path):
    record = json.loads(path.read_bytes())
    history.validate(record)
    assert record["base_ref"] == record["head_ref"]
    assert record["measured_pairs"] == 5
    assert set(record["cases"]) == set(gate.RELEASE_CASES)
    for name, case in record["cases"].items():
        for base, head in (("base", "head"), ("head", "base")):
            errors = gate.case_regressions(
                {"base": case["medians"][base], "head": case["medians"][head]},
                {"base": case["runs"][base], "head": case["runs"][head]},
            )
            assert not errors, (name, base, head, errors)


@pytest.mark.parametrize(
    "path",
    [STUDY / "linux-ryzen-7-7730u" / name for name in ("injection.json", "schwarz-injection.json")],
    ids=lambda p: p.stem,
)
def test_measured_delay_and_allocation_trip_only_the_new_resource_limits(path):
    experiment = json.loads(path.read_bytes())
    reference = json.loads(path.with_name("same-commit.json").read_bytes())
    assert experiment["delay_seconds"] == pytest.approx(
        0.2 * reference["cases"][experiment["case"]]["medians"]["base"]["seconds"],
    )
    assert experiment["extra_bytes"] == 12 * 1024**2
    runs = experiment["runs"]
    assert len(runs["base"]) == len(runs["head"]) == 5
    medians = {
        label: {key: median(row[key] for row in rows) for key in rows[0]}
        for label, rows in runs.items()
    }
    assert not gate.case_regressions(
        medians,
        runs,
        fixture_version=4,
    )
    errors = gate.case_regressions(medians, runs)
    assert [error.split()[0] for error in errors] == ["seconds", "peak_rss_bytes"]
