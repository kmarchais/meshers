"""A requested short edge must wait for longer edges in its incident tets."""

import numpy as np
import microgen_workflow as w
from longest_edge_closure import closure_scores
from nested_refinement import EDGES, refine, interpolation_check


def test_short_edge_request_is_propagated_to_a_locally_longest_edge():
    w.imports()
    from fedoo.lib_elements.tetrahedron import Tet10
    from test_nested_refinement import FACES

    points = Tet10().xi_nd.copy()
    cells = np.arange(10)[None]
    mids = {tuple(sorted(pair)): 4+i for i, pair in enumerate(EDGES)}
    surface = np.array([list(face)+[mids[tuple(sorted((face[a], face[b])))]
        for a, b in ((0, 1), (1, 2), (2, 0))] for face in FACES])
    request = {(0, 2): 10.}
    scores, report = closure_scores(points, cells, surface, request, periodic=False)
    assert report['closure_orbits'] > report['requested_orbits'] == 1
    assert (0, 2) not in scores
    for edge, value in scores.items():
        assert np.linalg.norm(points[edge[1]]-points[edge[0]]) >= np.sqrt(2)-1e-12
        assert value == 10.
    fields = np.random.default_rng(26).normal(size=(6, 10, 3))
    result = refine(points, cells, surface, fields, periodic=False, edge_scores=scores)
    assert result['diagnostics']['accepted_edge_orbits'] == 1
    interpolation_check(points, cells, fields, result)


def test_periodic_closure_keeps_whole_orbits_and_selects_only_local_longest_edges():
    _, _, meshers, *_ = w.imports()
    from periodic import periodic_pairs
    mesh = meshers.generate(lambda x, y, z: x*0-1, bounds=(-.5, .5)*3,
        cells=4, periodic=(True,)*3, order=2, threads=1)
    points, cells, surface = mesh.points, mesh.tetrahedra, mesh.surface
    edges = np.unique(np.sort(cells[:, EDGES].reshape(-1, 2), axis=1), axis=0)
    request = {tuple(map(int, edge)): 1. for edge in edges[::11]}
    scores, report = closure_scores(points, cells, surface, request)
    assert report['locally_longest_leaf_orbits'] > 0
    _, classes, _ = periodic_pairs(points, surface)
    selected = {tuple(sorted(classes[list(edge)])) for edge in scores}
    for edge in edges:
        assert (tuple(map(int, edge)) in scores) == (tuple(sorted(classes[edge])) in selected)
    for cell in cells:
        lengths = [np.linalg.norm(points[cell[b]]-points[cell[a]]) for a, b in EDGES]
        for edge, length in zip(EDGES, lengths, strict=True):
            if tuple(sorted(cell[edge])) in scores:
                assert length >= max(lengths)*(1-1e-12)
