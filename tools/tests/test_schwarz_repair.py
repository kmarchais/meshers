"""Check that the historical repair normalizes orientation and keeps boundary nodes fixed."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
pv = pytest.importorskip("pyvista")
pytest.importorskip("scipy")


def test_negative_tet_is_reoriented_without_moving_its_boundary():
    from historical_schwarz_repair import repaired_grid

    points = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    grid = pv.UnstructuredGrid({pv.CellType.TETRA: np.array([[1, 0, 2, 3]])}, points)
    record = {}
    repaired = repaired_grid(grid, record)
    vertices = repaired.points[repaired.cells_dict[pv.CellType.TETRA][0]]
    assert np.linalg.det(vertices[1:] - vertices[0]) == pytest.approx(1)
    np.testing.assert_array_equal(repaired.points, points)
    assert record["repair"]["accepted_moves"] == 0
