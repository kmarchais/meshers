"""Nodal reconstruction must evaluate the original cubic, not a P2 interpolant."""

import numpy as np

import microgen_workflow as workflow
from hp_coupled_audit import cubic_nodal_values


def test_cubic_values_at_child_nodes_reproduce_the_original_polynomial():
    workflow.imports()
    from cubic_reference import register_cubic
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4

    p1, p2, p3 = Tet4(), Tet10(), register_cubic()()
    maps = np.array([np.eye(4), np.eye(4)])
    maps[0, 1] = [.5, .5, 0., 0.]
    maps[1, 0] = [.5, .5, 0., 0.]
    points = (p1.shape_function(p2.xi_nd)[None] @ maps @ p1.xi_nd).reshape(-1, 3)
    vectors = np.arange(18).reshape(6, 3) + 1.
    scalar = lambda p: p[:, 0]**3 + 2 * p[:, 1] * p[:, 2]**2
    fields = vectors[:, None] * scalar(p3.xi_nd)[None, :, None]
    actual = cubic_nodal_values(points, np.arange(20).reshape(2, 10), np.arange(20)[None], fields,
        np.zeros(2, dtype=int), maps)
    expected = vectors[:, None] * scalar(points)[None, :, None]
    np.testing.assert_allclose(actual, expected, atol=1e-13)
