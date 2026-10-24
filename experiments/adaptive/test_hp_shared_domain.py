"""Affine map composition must retain child-to-parent ordering."""

import numpy as np

from hp_shared_domain import compose


def test_composition_reproduces_noncommuting_bisection_and_permutation():
    first = np.eye(4)
    first[1] = [.5, .5, 0., 0.]
    second = np.eye(4)[[2, 0, 1, 3]]
    parents, maps = compose(np.array([9, 4]), np.array([np.eye(4), first]), np.array([1]), second[None])
    barycentric = np.array([.1, .2, .3, .4])
    np.testing.assert_array_equal(parents, [4])
    np.testing.assert_allclose(barycentric @ maps[0], (barycentric @ second) @ first, atol=1e-15)
    assert not np.allclose(maps[0], first @ second)
