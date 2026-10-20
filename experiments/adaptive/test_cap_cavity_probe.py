"""Volume alternatives must preserve the complete non-cap oriented boundary."""

import numpy as np

from cap_cavity_probe import cavity
from cap_flip_probe import flip_pair


def star():
    points = np.array([[0, 0, -.5], [.4, 0, -.5], [.4, .4, -.5], [0, .4, -.5], [.2, .2, 0]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    cells = np.column_stack((faces, [4, 4]))
    return points, cells, faces


def test_complete_boundary_star_replacement_preserves_volume():
    points, cells, faces = star()
    alternative_faces, _ = flip_pair(points, faces, [0, 1])
    alternative, row = cavity(points, cells, faces, alternative_faces, [0, 2])
    assert row["cavity_passed"] and row["relative_volume_difference"] < 1e-12
    assert alternative.shape == (2, 4)


def test_incomplete_edge_star_is_rejected():
    points, cells, faces = star()
    cells = np.vstack((cells, cells))
    alternative_faces, _ = flip_pair(points, faces, [0, 1])
    alternative, row = cavity(points, cells, faces, alternative_faces, [0, 2])
    assert alternative is None and row["status"] == "not-two-or-three-tet-edge-star"


def test_wrong_cap_patch_is_rejected():
    points, cells, faces = star()
    alternative, row = cavity(points, cells, faces, np.array([[0, 1, 4], [1, 2, 4]]), [0, 2])
    assert alternative is None and not row["cavity_passed"]


def test_complete_three_tet_star_replacement_preserves_volume():
    points, _, faces = star()
    points[4] = [.3, .2, -.2]
    points = np.vstack((points, [.2, .3, -.2]))
    cells = np.array([[0, 2, 1, 4], [0, 2, 4, 5], [0, 2, 5, 3]])
    negative = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]]) < 0
    cells[negative, 1], cells[negative, 2] = cells[negative, 2].copy(), cells[negative, 1].copy()
    alternative_faces, _ = flip_pair(points, faces, [0, 1])
    alternative, row = cavity(points, cells, faces, alternative_faces, [0, 2])
    assert row["cavity_passed"] and row["relative_volume_difference"] < 1e-12
    assert alternative.shape == (4, 4)


def test_roundoff_scale_cone_is_not_accepted_by_rank_guard():
    points = np.array([[.3666666666666667, -.5, -.09999999999999998], [.3666666666666667, -.4666666666666667, -.09999999999999998],
        [.3666666666666667, -.5, -.06666666666666665], [.33665219210954817, -.5, -.09999999999999998],
        [.34382352816390666, -.4895098051694267, -.09999999999999998], [.3406272003984264, -.5, -.09270613293490694]])
    cells = np.array([[2, 4, 3, 5], [1, 3, 2, 4], [0, 2, 1, 3]])
    old_faces, new_faces = np.array([[0, 2, 3], [2, 5, 3]]), np.array([[0, 2, 5], [0, 5, 3]])
    alternative, _ = cavity(points, cells, old_faces, new_faces, [2, 3])
    determinants = np.linalg.det(points[alternative[:, 1:]] - points[alternative[:, :1]])
    assert determinants.min() / determinants.sum() < 1e-12
    alternative, row = cavity(points, cells, old_faces, new_faces, [2, 3], minimum_relative_determinant=1e-12)
    assert alternative is None and not row["cavity_passed"]
