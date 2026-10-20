"""Small binary solve agrees with exhaustive six-vertex controls."""

import numpy as np
import pytest

from bistellar_probe import oriented_boundary
from cavity_milp import solve_cavity
from cavity_triangulation import replacement_boundary
from test_cavity_triangulation import saved_nonconical_star


def test_original_cavity_has_exact_integer_face_replay():
    points, cells = saved_nonconical_star()
    target = oriented_boundary(cells)
    alternative, row = solve_cavity(points, cells, target)
    assert row["passed"] and row["integer_face_equations_exact"]
    assert oriented_boundary(alternative) == target


def test_original_obstruction_agrees_with_exhaustive_result():
    points, cells = saved_nonconical_star()
    target = replacement_boundary(points, cells, np.array([[2, 3, 5], [0, 2, 3]]), np.array([[0, 2, 5], [0, 5, 3]]))
    alternative, row = solve_cavity(points, cells, target)
    assert alternative is None and not row["passed"] and row["solver_status"] == 2


def test_seven_vertex_limit_is_enforced():
    points = np.zeros((8, 3))
    with pytest.raises(ValueError, match="seven-vertex"):
        solve_cavity(points, np.arange(8).reshape(2, 4), {})
