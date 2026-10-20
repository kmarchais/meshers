"""Cap alternatives retain the exact oriented outer polygon and area."""

import numpy as np

from cap_flip_probe import flip_pair, triangle_boundary


def square():
    points = np.array([[0, 0, -.5], [.4, 0, -.5], [.4, .4, -.5], [0, .4, -.5]])
    return points, np.array([[0, 1, 2], [0, 2, 3]])


def test_convex_flip_preserves_oriented_patch_boundary():
    points, triangles = square()
    alternative, row = flip_pair(points, triangles, [0, 1])
    assert row["local_target_passed"] and row["relative_area_difference"] < 1e-12
    assert triangle_boundary(alternative) == triangle_boundary(triangles)


def test_concave_patch_is_rejected():
    points, triangles = square()
    points[2, :2] = [.08, .08]
    alternative, row = flip_pair(points, triangles, [0, 1])
    assert alternative is None and row["status"] == "nonconvex-or-boundary-mismatch"


def test_existing_alternative_edge_is_rejected():
    points, triangles = square()
    triangles = np.vstack((triangles, [1, 3, 0]))
    alternative, row = flip_pair(points, triangles, [0, 1])
    assert alternative is None and row["status"] == "alternative-edge-already-present"


def test_non_coplanar_patch_is_rejected():
    points, triangles = square()
    points[3, 2] += .01
    alternative, row = flip_pair(points, triangles, [0, 1])
    assert alternative is None and row["status"] == "different-cap-planes"


def test_thin_convex_patch_is_not_promoted_by_area_preservation():
    points, triangles = square()
    points[:, 1] *= .01
    alternative, row = flip_pair(points, triangles, [0, 1])
    assert alternative is not None and not row["local_target_passed"]
    assert row["status"] == "cap-ceiling-rejected"
