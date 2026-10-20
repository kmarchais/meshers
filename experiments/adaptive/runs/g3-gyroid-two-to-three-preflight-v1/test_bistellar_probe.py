"""Five-vertex cavity geometry and topology admission controls."""

import numpy as np

from bistellar_probe import flip_candidate, oriented_boundary


def bipyramid():
    points = np.array([[-.4, -.2, 0], [.4, -.2, 0], [0, .4, 0], [0, 0, .3], [0, 0, -.3]])
    cells = np.array([[0, 1, 2, 3], [0, 2, 1, 4]])
    assert np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]]).min() > 0
    return points, cells


def test_convex_cavity_preserves_oriented_boundary_and_volume():
    points, cells = bipyramid()
    alternative, row = flip_candidate(points, cells, [0, 1])
    assert row["local_target_passed"] and alternative.shape == (3, 4)
    assert oriented_boundary(cells) == oriented_boundary(alternative)
    assert row["relative_cavity_volume_difference"] < 1e-12


def test_nonconvex_alternative_is_rejected_by_oriented_boundary():
    points, cells = bipyramid()
    points[4, :2] = [.49, .49]
    alternative, row = flip_candidate(points, cells, [0, 1])
    assert alternative is None and row["status"] == "oriented-cavity-boundary-mismatch"


def test_existing_new_edge_is_rejected():
    points, cells = bipyramid()
    cells = np.vstack((cells, [3, 4, 0, 1]))
    alternative, row = flip_candidate(points, cells, [0, 1])
    assert alternative is None and row["status"] == "new-edge-already-present"
