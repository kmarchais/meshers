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
    cells = np.vstack((cells, cells[:1]))
    alternative_faces, _ = flip_pair(points, faces, [0, 1])
    alternative, row = cavity(points, cells, faces, alternative_faces, [0, 2])
    assert alternative is None and row["status"] == "not-two-tet-edge-star"


def test_wrong_cap_patch_is_rejected():
    points, cells, faces = star()
    alternative, row = cavity(points, cells, faces, np.array([[0, 1, 4], [1, 2, 4]]), [0, 2])
    assert alternative is None and not row["cavity_passed"]
