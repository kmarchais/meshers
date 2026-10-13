"""Verify exact partition and boundary preservation independently of the marker."""

import numpy as np
import pytest
from body_refine import split_interior_edges, verify_preserved_boundary
from microgen_benchmark.geometry import boundary, canonical_mesh


def test_cube_body_diagonal_split_preserves_boundary_volume_and_affine_field():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
                       [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], float) - .5
    tets = np.array([[0, 1, 3, 7], [0, 3, 2, 7], [0, 2, 6, 7],
                     [0, 6, 4, 7], [0, 4, 5, 7], [0, 5, 1, 7]])
    points, tets, _, _ = canonical_mesh(points, tets)
    original_boundary = boundary(points, tets)
    edges = np.unique(np.sort(tets[:, [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]].reshape(-1, 2), axis=1), axis=0)
    refined_points, children, info = split_interior_edges(points, tets, original_boundary, np.ones(len(edges)), 1)
    assert info["inserted_body_nodes"] == 1
    assert len(children) == 2 * len(tets)
    np.testing.assert_array_equal(refined_points[:len(points)], points)
    assert {tuple(sorted(t)) for t in boundary(refined_points, children)} == {tuple(sorted(t)) for t in original_boundary}
    refined_boundary = boundary(refined_points, children)
    verify_preserved_boundary(points, original_boundary, refined_points, refined_boundary)
    # Triangle order and cyclic permutations preserve the same oriented surface.
    verify_preserved_boundary(points, original_boundary, refined_points, np.roll(refined_boundary[::-1], 1, axis=1))
    moved = refined_points.copy()
    moved[0, 0] += .001
    with pytest.raises(ValueError, match="original vertex"):
        verify_preserved_boundary(points, original_boundary, moved, refined_boundary)
    with pytest.raises(ValueError, match="boundary triangle"):
        verify_preserved_boundary(points, original_boundary, refined_points, refined_boundary[:, ::-1])
    volumes = np.linalg.det(refined_points[children[:, 1:]] - refined_points[children[:, :1]]) / 6
    assert volumes.min() > 0
    np.testing.assert_allclose(volumes.sum(), 1, atol=1e-14)
    # Midpoint interpolation retains every affine displacement exactly.
    matrix = np.array([[1, .2, .3], [.4, 2, .5], [.6, .7, 3]])
    displacement = points @ matrix.T
    np.testing.assert_allclose((displacement[0] + displacement[7]) / 2,
                               refined_points[-1] @ matrix.T, atol=1e-14)


def test_zero_indicator_does_not_refine():
    points = np.array([[0, 0, 0], [.3, 0, 0], [0, .3, 0], [0, 0, .3]], float)
    tets = np.array([[0, 1, 2, 3]])
    p, t, info = split_interior_edges(points, tets, boundary(points, tets), np.zeros(6), 2)
    assert info["inserted_body_nodes"] == 0
    np.testing.assert_array_equal(p, points)
    np.testing.assert_array_equal(t, tets)
