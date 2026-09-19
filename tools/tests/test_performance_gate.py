"""Prove the gate rejects each measured regression and accepts bounded noise."""

import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "performance_gate", Path(__file__).parents[1] / "performance_gate.py"
)
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)

BASE = {
    "seconds": 1.0,
    "peak_rss_bytes": 100_000_000,
    "tetrahedra": 1000,
    "volume": 0.5,
    "quality_min": 0.1,
    "quality_p05": 0.3,
    "quality_median": 0.7,
    "surface_error": 0.01,
}


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("seconds", 1.5),
        ("peak_rss_bytes", 150_000_000),
        ("tetrahedra", 800),
        ("volume", 0.45),
        ("quality_min", 0.08),
        ("quality_p05", 0.2),
        ("quality_median", 0.6),
        ("surface_error", 0.02),
        ("seconds", float("nan")),
    ],
)
def test_rejects_regression(key, value):
    assert gate.regressions(BASE, {**BASE, key: value})


def test_accepts_identical_and_bounded_resource_noise():
    assert not gate.regressions(BASE, BASE)
    assert not gate.regressions(BASE, {**BASE, "seconds": 1.1, "peak_rss_bytes": 105_000_000})
