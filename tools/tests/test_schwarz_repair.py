"""Check that the historical repair normalizes orientation and keeps boundary nodes fixed."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
pytest.importorskip("pyvista")
pytest.importorskip("scipy")
import pyvista as pv
from historical_schwarz_repair import audit, check_repair, repaired_grid


def test_negative_tet_is_reoriented_without_moving_its_boundary():
    points = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    grid = pv.UnstructuredGrid({pv.CellType.TETRA: np.array([[1, 0, 2, 3]])}, points)
    record = {}
    repaired = repaired_grid(grid, record)
    vertices = repaired.points[repaired.cells_dict[pv.CellType.TETRA][0]]
    assert np.linalg.det(vertices[1:] - vertices[0]) == pytest.approx(1)
    np.testing.assert_array_equal(repaired.points, points)
    assert record["repair"]["accepted_moves"] == 0


def test_repeated_tet_vertices_are_rejected():
    cells = np.array([[0, 1, 2, 2]])
    with pytest.raises(ValueError, match="repeated vertices"):
        audit(np.eye(3), cells, np.array([[0, 1, 2]]))


@pytest.mark.parametrize(("scale", "message"), [(-1, "determinant"), (1.1, "volume")])
def test_invalid_repair_fails_even_with_optimized_python(scale, message):
    points = np.vstack([np.zeros(3), np.eye(3)])
    moved = points.copy()
    moved[1, 0] *= scale
    with pytest.raises(ValueError, match=message):
        check_repair(points, moved, np.array([[0, 1, 2, 3]]), np.array([True, False, True, True]))
