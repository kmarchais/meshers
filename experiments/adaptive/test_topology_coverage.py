"""Independent controls for analytic fixtures and boundary topology accounting."""

import numpy as np
import pytest

from topology_coverage import boundary_topology, fixture


@pytest.mark.parametrize('name', ('torus', 'shell', 'disjoint', 'octahedron'))
def test_reference_points_are_on_the_exact_surface(name):
    field, queries, volume, components, euler = fixture(name, np)
    assert np.max(np.abs(field(*queries.T))) < 3e-16
    assert volume > 0 and components > 0
    assert euler in (0, 2, 4)


def test_boundary_topology_counts_disconnected_components_and_torus_handle():
    tetrahedron = np.array(((0, 1, 2), (0, 3, 1), (1, 3, 2), (0, 2, 3)))
    assert boundary_topology(tetrahedron, np) == (1, 2)
    assert boundary_topology(np.concatenate((tetrahedron, tetrahedron+4)), np) == (2, 4)
    faces = []
    for i in range(4):
        for j in range(4):
            a, b, c, d = 4*i+j, 4*((i+1)%4)+j, 4*((i+1)%4)+(j+1)%4, 4*i+(j+1)%4
            faces.extend(((a, b, c), (a, c, d)))
    assert boundary_topology(np.array(faces), np) == (1, 0)
