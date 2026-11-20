# Copyright (c) 2026 Kevin Marchais

"""Prove the gate rejects each measured regression and accepts bounded noise."""

import os
import runpy
import sys
import venv
from pathlib import Path

import numpy as np
import pytest

from tools import performance_gate as gate
from tools.commands import run_command
from tools.performance_worker import cell_geometry

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


def test_driver_loads_without_worker_dependencies(monkeypatch):
    for dependency in (
        "meshers",
        "numpy",
        "psutil",
        "tools.performance_worker",
        "performance_worker",
    ):
        monkeypatch.setitem(sys.modules, dependency, None)
    runpy.run_path(gate.__file__)


def test_command_uses_selected_virtual_environment(tmp_path):
    environment = tmp_path / "wheel environment"
    venv.EnvBuilder(with_pip=False, symlinks=os.name != "nt").create(environment)
    python = environment / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    script = tmp_path / "identity.py"
    script.write_text("import sys\nprint(sys.prefix)\n", encoding="utf-8")
    assert (
        Path(
            run_command(
                python.name,
                str(script),
                env={**os.environ, "PATH": str(python.parent)},
                timeout=30,
                stdout_encoding="utf-8",
            ).stdout.strip(),
        )
        == environment
    )


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
    corners = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    edges = np.array([(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)])
    points = np.vstack([corners, corners[edges].mean(axis=1)])[None]
    np.testing.assert_allclose(cell_geometry(points), [1 / 6])
    points[0, 4] = [2.0, 2.0, 2.0]
    with pytest.raises(RuntimeError, match="Jacobian"):
        cell_geometry(points)


def test_a_single_bad_quality_sample_cannot_hide_in_passing_medians():
    rows = [dict(BASE) for _ in range(5)]
    rows[2]["quality_p10"] = 0.1
    assert not gate.regressions(BASE, BASE)
    assert "sample 3: quality_p10" in gate.case_regressions(BASE, BASE, rows)[0]


def test_fold_between_quadrature_points():
    corners = np.vstack([np.zeros(3), np.eye(3)])
    edges = np.array([(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)])
    r, s, t = np.vstack([corners, corners[edges].mean(axis=1)]).T
    folded = np.column_stack([r, (r - 0.38) * s, (r - 0.42) * t])
    with pytest.raises(RuntimeError, match="curved Jacobian"):
        cell_geometry(folded[None])
    # This has nearly the same corner maps but a positive interior minimum.
    valid = np.column_stack([r, (r - 0.4) * s - 0.025 * t, (r - 0.4) * t + 0.025 * s])
    assert cell_geometry(valid[None])[0] == pytest.approx(
        (0.1 - 0.8 / 4 + 0.16 + 0.025**2) / 6,
    )


def test_mapped_cap_error_is_gated_independently_of_larger_wall_error():
    base = {**BASE, "mapped_cap_error": 0.0001}
    head = {**base, "mapped_cap_error": 0.0002}
    assert gate.regressions(base, head) == ["mapped cap error increased by more than 5%"]


@pytest.mark.skipif(sys.platform != "win32", reason="Windows working-set measurement")
def test_identity_without_windows_memory_sampler(monkeypatch):
    monkeypatch.setitem(sys.modules, "psutil", None)
    worker = runpy.run_path(str(Path(gate.__file__).with_name("performance_worker.py")))
    assert worker["runtime_identity"]()["native_sha256"]
    with pytest.raises(ModuleNotFoundError, match="Windows benchmark memory sampling"):
        worker["peak_rss_bytes"]()
