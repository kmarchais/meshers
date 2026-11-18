"""Known-answer controls for the independent periodic-volume observer."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
pytest.importorskip("scipy")
pytest.importorskip("pyvista")
from periodic_geometry import BASELINE, audit, boundary_metrics, field, forward, roots


@pytest.mark.parametrize(
    ("specimen", "level"), [("gyroid", -0.55), ("gyroid", 0.55), ("schwarz", 0.0)]
)
def test_analytic_wall_samples(specimen, level):
    lines = np.random.default_rng(24).uniform(-0.5, 0.5, (200, 3))
    for axis in range(3):
        coords, valid = roots(lines, specimen, level, axis)
        assert valid.any()
        for column in coords.T:
            probes = lines[valid].copy()
            probes[:, axis] = column[valid]
            np.testing.assert_allclose(field(probes, specimen), level, atol=1e-12)
            assert forward(probes, specimen, level, range(3)).max() < 1e-12


def test_missing_patch_on_cap():
    points = np.array([[-0.5, -0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    net = np.random.default_rng(971).uniform(-0.5, 0.5, (256, 2))
    whole, _ = boundary_metrics(points, faces, "schwarz", net)
    damaged, _ = boundary_metrics(points, faces[:1], "schwarz", net)
    key = "z+0.5"
    assert whole[key]["reference_to_mesh"]["maximum"] < 1e-12
    assert damaged[key]["reference_to_mesh"]["maximum"] > 0.1


def test_detached_tetrahedron_fails_component_check():
    # A second positive tet with no shared face must fail even if all probes miss it.
    import pyvista as pv

    first = np.array([[-0.4, -0.4, -0.4]]) + 0.1 * np.vstack([np.zeros(3), np.eye(3)])
    extra = 0.2495 + 0.00004 * np.vstack([np.zeros(3), np.eye(3)])
    baseline = json.loads(BASELINE.read_bytes())
    for points, expected in [(first, 1), (np.vstack([first, extra]), 2)]:
        cells = np.arange(len(points)).reshape(-1, 4)
        grid = pv.UnstructuredGrid({pv.CellType.TETRA: cells}, points)
        result = audit(grid, "schwarz", baseline, power=4)
        assert result["face_connected_components"] == expected
        assert ("face_connected_components" in result["failed_checks"]) == (expected == 2)
        if expected == 2:
            assert not result["sampled_admission"]
