"""Affine elevation preserves geometry, shared faces and the preallocation cap."""

import numpy as np

from microgen_tet10_control import straight_elevation
from test_facet_incidence_audit import surface


def test_shared_edge_elevation_and_node_cap():
    points = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.], [0.,0.,-1.]])
    cells = np.array([[0,1,2,3], [0,2,1,4]])
    faces = surface(cells)
    result = straight_elevation(points, cells, faces)
    assert result["nodes"] == 14 and result["within_node_cap"]
    np.testing.assert_array_equal(result["points"][:5], points)
    for cell, old in zip(result["tetrahedra"], cells, strict=True):
        for offset, (a,b) in enumerate(((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))):
            np.testing.assert_array_equal(result["points"][cell[4+offset]], (points[old[a]]+points[old[b]])/2)
    assert set(result["tetrahedra"][0,4:]) & set(result["tetrahedra"][1,4:])
    rejected = straight_elevation(points, cells, faces, node_cap=13)
    assert rejected == {"nodes": 14, "within_node_cap": False}
