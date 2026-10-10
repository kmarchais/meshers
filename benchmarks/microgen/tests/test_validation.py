import numpy as np
import pytest

from microgen_benchmark.cases import CASES
from microgen_benchmark.geometry import canonical_mesh, overlapping_pair, validate


def tet():
    return np.array([[0.0, 0, 0], [0.3, 0, 0], [0, 0.3, 0], [0, 0, 0.3]])


def test_overlap_and_contact_are_distinct():
    a = tet()
    assert overlapping_pair(np.array([a, a + 0.02]))[0] == [0, 1]
    assert overlapping_pair(np.array([a, a + 0.3]))[0] is None


def test_incomplete_validation_cannot_pass():
    a = tet()
    with pytest.raises(RuntimeError, match="incomplete"):
        overlapping_pair(np.array([a, a + 0.01, a + 0.02]), max_pairs=1)


def test_duplicate_cells_rejected():
    p, t, _, _ = canonical_mesh(tet(), [[0, 1, 2, 3], [0, 1, 2, 3]])
    with pytest.raises(ValueError, match="Duplicate"):
        validate(p, t, CASES[0])


def test_outside_specimen_rejected():
    p, t, _, _ = canonical_mesh(tet() + 1, [[0, 1, 2, 3]])
    with pytest.raises(ValueError, match="specimen"):
        validate(p, t, CASES[0])


def test_vertex_only_material_contact_rejected():
    a = tet()
    points = np.vstack([a, -a[1:]])
    p, t, _, _ = canonical_mesh(points, [[0, 1, 2, 3], [0, 4, 5, 6]])
    with pytest.raises(ValueError, match="touch only"):
        validate(p, t, CASES[0])
