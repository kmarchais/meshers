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
    "nodes": 300,
    "volume": 0.5,
    "quality_min": 0.1,
    "quality_p05": 0.3,
    "quality_p10": 0.4,
    "quality_mean": 0.65,
    "quality_p90": 0.9,
    "quality_max": 1.0,
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
        ("quality_mean", 0.5),
        ("quality_p10", 0.2),
        ("quality_p90", 0.7),
        ("quality_max", 0.8),
        ("nodes", 100),
        ("surface_error", 0.02),
        ("seconds", float("nan")),
    ],
)
def test_rejects_regression(key, value):
    assert gate.regressions(BASE, {**BASE, key: value})


def test_accepts_identical_and_bounded_resource_noise():
    assert not gate.regressions(BASE, BASE)
    assert not gate.regressions(BASE, {**BASE, "seconds": 1.1, "peak_rss_bytes": 105_000_000})


def test_curved_geometry_integration_and_inversion():
    import numpy as np

    corners = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    edges = np.array([(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)])
    points = np.vstack([corners, corners[edges].mean(axis=1)])[None]
    np.testing.assert_allclose(gate.cell_geometry(points), [1 / 6])
    points[0, 4] = [2.0, 2.0, 2.0]
    with pytest.raises(RuntimeError, match="Jacobian"):
        gate.cell_geometry(points)


def test_a_single_bad_quality_sample_cannot_hide_in_passing_medians():
    rows = [dict(BASE) for _ in range(5)]
    rows[2]["quality_p10"] = 0.1
    assert not gate.regressions(BASE, BASE)
    assert "sample 3: quality_p10" in gate.case_regressions(BASE, BASE, rows)[0]
