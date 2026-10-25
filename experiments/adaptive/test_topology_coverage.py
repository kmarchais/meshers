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


def test_periodic_wedge_has_exact_volume_matching_caps_and_an_unchanged_normal_jump():
    from topology_coverage import audit_column_periodicity, mapped_wedge_input, periodic_wedge

    points, faces = periodic_wedge(np)
    field, queries, volume, components, euler = fixture('periodic-wedge', np)
    closed = np.maximum(field(*queries.T), np.max(abs(queries), axis=1)-.5)
    assert abs(closed).max() < 2e-16
    triangles = points[faces]
    actual_volume = np.einsum('ij,ij->i', triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])).sum()/6
    assert abs(actual_volume-volume) < 1e-15 and abs(volume-.22) < 1e-15
    assert boundary_topology(faces, np) == (components, euler) == (1, 2)
    assert audit_column_periodicity(points, faces, np)['pair_count'] == 4
    constraints, mapping, transforms = mapped_wedge_input(np)
    probes = np.random.default_rng(42).uniform(-.5, .5, (1024, 3))
    physical = np.column_stack(mapping(*probes.T))
    parameter_inside = np.max([f(*probes.T) for f in constraints.values()], axis=0) <= 0
    np.testing.assert_array_equal(field(*physical.T) <= 0, parameter_inside)
    low = probes.copy(); low[:, 2] = -.5
    high = low.copy(); high[:, 2] = .5
    np.testing.assert_array_equal(field(*low.T), field(*high.T))
    np.testing.assert_allclose(np.column_stack(mapping(*high.T))-np.column_stack(mapping(*low.T)),
        np.broadcast_to((0, 0, 1), high.shape), atol=1e-15)
    assert transforms[2][2][3] == 1
    limits = np.array(((.2, 0., -.5+1e-6), (.2, 0., .5-1e-6)))
    gradients = []
    for axis in range(3):
        step = np.eye(3)[axis]*1e-7
        gradients.append((field(*(limits+step).T)-field(*(limits-step).T))/2e-7)
    np.testing.assert_allclose(np.column_stack(gradients), ((1, 0, -.3), (1, 0, .3)), atol=1e-8)
    cosine = (1-.3**2)/(1+.3**2)
    assert abs(np.degrees(np.arccos(cosine))-2*np.degrees(np.arctan(.3))) < 1e-12
    assert np.degrees(np.arccos(cosine)) > 30


def test_rotated_periodic_wedge_preserves_volume_caps_and_physical_domain():
    from topology_coverage import audit_column_periodicity, periodic_wedge

    angle = .37
    c, s = np.cos(angle), np.sin(angle)
    rotation = np.array(((c, -s, 0.), (s, c, 0.), (0., 0., 1.)))
    translation = np.array((.02, -.03, 0.))
    original, faces = periodic_wedge(np)
    points = original@rotation.T+translation
    assert np.max(abs(points[:, :2])) < .5
    assert audit_column_periodicity(points, faces, np)['matching_caps']
    triangles = points[faces]
    volume = np.einsum('ij,ij->i', triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])).sum()/6
    assert abs(volume-.22) < 1e-15
    field, queries, _, _, _ = fixture('periodic-wedge', np)
    transformed = queries@rotation.T+translation
    recovered = (transformed-translation)@rotation
    closed = np.maximum(field(*recovered.T), np.max(abs(transformed), axis=1)-.5)
    assert abs(closed).max() < 5e-16
    rng = np.random.default_rng(42)
    probes = rng.uniform(-.5, .5, (1024, 3))
    probes[:, 2] = -.5
    opposite = probes.copy(); opposite[:, 2] = .5
    np.testing.assert_array_equal(field(*((probes-translation)@rotation).T),
        field(*((opposite-translation)@rotation).T))


def test_periodic_notched_wedge_has_exact_volume_caps_and_reentrant_boundary():
    from topology_coverage import audit_column_periodicity, periodic_notched_wedge

    p, f = periodic_notched_wedge(np)
    field, queries, volume, components, euler = fixture('periodic-notched-wedge', np)
    assert abs(np.maximum(field(*queries.T), np.max(abs(queries), axis=1)-.5)).max() < 3e-16
    triangles = p[f]
    actual = np.einsum('ij,ij->i', triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])).sum()/6
    assert abs(volume-(.22-.15*(.275-.05))) < 1e-15
    assert abs(actual-volume) < 1e-15
    assert boundary_topology(f, np) == (components, euler) == (1, 2)
    edges = np.sort(f[:, ((0,1),(1,2),(2,0))].reshape(-1,2), axis=1)
    assert np.all(np.unique(edges, axis=0, return_counts=True)[1] == 2)
    periodic = audit_column_periodicity(p, f, np)
    assert periodic['pair_count'] == 6 and periodic['cap_count'] == 4
    assert field(.1,.1,0.) > 0 and field(.0,.1,0.) < 0 and field(.1,.0,0.) < 0
    low = np.random.default_rng(42).uniform(-.5,.5,(1024,3)); low[:,2] = -.5
    high = low.copy(); high[:,2] = .5
    np.testing.assert_array_equal(field(*low.T), field(*high.T))
