"""Geometry and solver admission controls for the historical implicit solids."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
import implicit_coverage as coverage


@pytest.mark.parametrize(("geometry_ok", "quality_ok"), [(False, True), (True, False)])
def test_rejected_geometry_never_runs_fea(monkeypatch, geometry_ok, quality_ok):
    mesh = SimpleNamespace(points=np.eye(4, 3), tetrahedra=np.array([[0, 1, 2, 3]]))
    monkeypatch.setattr(coverage.meshers, "generate", lambda *_args, **_kw: mesh)
    case = SimpleNamespace(id="gyroid", band=(-1, 1), field=None)
    oracle = SimpleNamespace(
        reference=SimpleNamespace(boundary_samples=lambda _case: np.zeros((1, 3))),
        geometry=SimpleNamespace(validate=lambda *_args: ({"solid_fraction": 0.16}, None)),
        qualification=SimpleNamespace(
            admission=lambda *_args: {
                "geometry_eligible": geometry_ok,
                "quality_eligible": quality_ok,
                "admission_reasons": [],
                "relative_volume_error_with_uncertainty": 0.01,
            },
        ),
    )
    record = coverage.measure(case, {}, oracle, fea=True)
    assert record["status"] == "geometry-rejected"
    assert "fea" not in record


@pytest.mark.parametrize(("error", "status"), [(0.019, "qualified"), (0.021, "outside-target")])
def test_energy_target_is_separate_from_geometry(error, status):
    metrics = dict.fromkeys(("patch_max_error", "free_residual", "assembly_s", "solve_s"), 0.0)
    metrics.update(worst_energy_error=error, free_dofs=6)
    mesh = SimpleNamespace(points=None, tetrahedra=None)
    result = coverage.mechanics(lambda *_args: (metrics, {}), mesh, None, 1)
    assert result["status"] == status
    assert result["worst_energy_error"] == error
