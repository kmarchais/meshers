"""Explicit seed replay must preserve coordinates and reject obstructed input."""

import json

import numpy as np
import pytest

from microgen_tuned import load_prepared_seed


def test_prepared_seed_preserves_arrays_and_requires_cap_and_rank_checks(tmp_path):
    points = np.array([[0, 0, 0], [.2, 0, 0], [0, .2, 0], [0, 0, .2]])
    cells = np.array([[0, 1, 2, 3]])
    np.savez(tmp_path / "mesh.npz", points=points, tetrahedra=cells)
    row = {"case": "gyroid", "original_points_identical": True, "wall_triangles_identical": True,
        "cap_faces_below_quality_gate": 0, "minimum_replacement_relative_determinant": .001}
    (tmp_path / "result.json").write_text(json.dumps(row))
    grid = load_prepared_seed(tmp_path)
    np.testing.assert_array_equal(grid.points, points)
    np.testing.assert_array_equal(grid.cells_dict[10], cells)
    row["cap_faces_below_quality_gate"] = 1
    (tmp_path / "result.json").write_text(json.dumps(row))
    with pytest.raises(ValueError, match="obstructing caps"):
        load_prepared_seed(tmp_path)
