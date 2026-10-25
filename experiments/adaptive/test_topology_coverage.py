"""Independent controls for analytic fixtures and boundary topology accounting."""

import numpy as np
import pytest

from topology_coverage import bored_prism_constraints, bored_prism_distance, boundary_topology, fixture, notched_prism


@pytest.mark.parametrize('name', ('torus', 'shell', 'disjoint', 'octahedron', 'notched-prism', 'bored-prism'))
def test_reference_points_are_on_the_exact_surface(name):
    field, queries, volume, components, euler = fixture(name, np)
    closed = np.maximum(field(*queries.T), np.max(np.abs(queries), axis=1)-.5)
    assert np.max(np.abs(closed)) < 3e-16
    assert volume > 0 and components > 0
    assert euler in (0, 2, 4)


def test_notched_prism_has_the_exact_oriented_volume_and_boundary():
    points, surface = notched_prism(np)
    triangles = points[surface]
    volume = np.einsum('ij,ij->i', triangles[:, 0],
        np.cross(triangles[:, 1], triangles[:, 2])).sum()/6
    assert abs(volume-.24) < 1e-15
    assert boundary_topology(surface, np) == (1, 2)
    edges = np.sort(surface[:, ((0, 1), (1, 2), (2, 0))].reshape(-1, 2), axis=1)
    assert np.all(np.unique(edges, axis=0, return_counts=True)[1] == 2)


def test_bored_prism_distance_handles_void_caps_and_exterior_corners():
    _, queries, volume, components, euler = fixture('bored-prism', np)
    assert np.max(bored_prism_distance(queries, np)) < 1e-15
    probes = np.array(((.03, -.02, 0.), (.03, -.02, .3), (.03, -.02, .4),
        (.45, .45, .4), (.3, .3, .4), (.3, .3, 0.)))
    expected = (.12, .12, np.hypot(.12, .1), np.sqrt(.03), .1, .05)
    np.testing.assert_allclose(bored_prism_distance(probes, np), expected, atol=1e-15)
    assert abs(volume-(.7*.7*.6-np.pi*.12**2*.6)) < 1e-15
    assert (components, euler) == (1, 0)


def test_named_bored_prism_is_the_identical_scalar_domain():
    field, _, _, _, _ = fixture('bored-prism', np)
    points = np.random.default_rng(42).uniform(-.5, .5, (1024, 3))
    named = bored_prism_constraints(np)
    assert len(named) == 7
    np.testing.assert_allclose(np.max([f(*points.T) for f in named.values()], axis=0), field(*points.T), atol=1e-16)


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


def test_periodic_column_preserves_exact_domain_and_cap_pairing():
    from topology_coverage import audit_column_periodicity, periodic_column, periodic_column_constraints

    points, faces = periodic_column(np)
    field, queries, volume, components, euler = fixture('periodic-column', np)
    closed = np.maximum(field(*queries.T), np.max(abs(queries), axis=1)-.5)
    assert abs(closed).max() < 2e-16
    triangles = points[faces]
    actual_volume = np.einsum('ij,ij->i', triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])).sum()/6
    assert abs(actual_volume-volume) < 1e-15
    assert boundary_topology(faces, np) == (components, euler) == (1, 2)
    assert audit_column_periodicity(points, faces, np)['pair_count'] == 4
    rng = np.random.default_rng(42)
    probes = rng.uniform(-.5, .5, (1024, 3))
    np.testing.assert_allclose(field(*probes.T), np.max([f(*probes.T) for f in periodic_column_constraints(np).values()], axis=0), atol=1e-16)
    shifted = probes.copy();shifted[:, 2] += 1
    np.testing.assert_array_equal(field(*probes.T), field(*shifted.T))
    with pytest.raises(AssertionError):
        audit_column_periodicity(points, faces[1:], np)
