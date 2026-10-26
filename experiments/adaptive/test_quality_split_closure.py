"""Do not accept a fine-target split that leaves a child unable to refine."""

import json
from pathlib import Path
import numpy as np

import microgen_workflow as w
from nested_refinement import EDGES, FACES, positive_quality, refine
from quality_split_closure import select_with_quality_closure


def surface():
    mids = {tuple(sorted(pair)): 4+i for i, pair in enumerate(EDGES)}
    return np.array([list(face)+[mids[tuple(sorted((face[a], face[b])))]
        for a,b in ((0,1),(1,2),(2,0))] for face in FACES])


def fine_metric(points):
    return np.broadcast_to(np.eye(3)/.002**2, (len(points),3,3))


def legal_splits(points, cells):
    allowed = np.ones((len(cells),6), dtype=bool)
    for offset, endpoints in enumerate(EDGES):
        for endpoint in endpoints:
            xyz = points[cells[:,:4]].copy()
            xyz[:,endpoint] = points[cells[:,offset+4]]
            allowed[:,offset] &= positive_quality(xyz) >= .1+1e-12
    return allowed


def test_real_admitted_greedy_split_must_not_trap_a_child():
    w.imports()
    fixture = json.loads((Path(__file__).parent/'fixtures/quality-split-trap.json').read_text())
    points, cells = np.array(fixture['parent_points']), np.arange(10)[None]
    requested = tuple(sorted(fixture['requested_local_edge']))
    # This is the original failure: all immediate children pass q, but one
    # cannot be bisected along any of its six edges under the same q gate.
    result = refine(points, cells, surface(), np.zeros((6,10,3)), periodic=False,
        edge_scores={requested:1.})
    assert len(result['tetrahedra']) == 2
    assert positive_quality(result['points'][result['tetrahedra'][:,:4]]).min() >= .1
    assert (~legal_splits(result['points'],result['tetrahedra']).any(axis=1)).any()
    selected, _ = select_with_quality_closure(points,cells,surface(),{requested:1.},
        fine_metric,periodic=False)
    assert requested not in selected, 'Immediate child quality admitted a refinement dead end'


def test_regular_short_edge_does_not_need_longest_edge_propagation():
    w.imports()
    from fedoo.lib_elements.tetrahedron import Tet10
    points, cells = Tet10().xi_nd.copy(), np.arange(10)[None]
    requested = (0,2)
    selected, report = select_with_quality_closure(points,cells,surface(),{requested:1.},
        fine_metric,periodic=False)
    assert requested in selected
    assert report['direct_requested_orbits'] == 1


def test_finished_metric_does_not_require_an_unneeded_future_split():
    w.imports()
    fixture = json.loads((Path(__file__).parent/'fixtures/quality-split-trap.json').read_text())
    points,cells = np.array(fixture['parent_points']),np.arange(10)[None]
    requested = tuple(sorted(fixture['requested_local_edge']))
    def coarse(points):
        return np.broadcast_to(np.eye(3), (len(points),3,3))
    selected,_ = select_with_quality_closure(points,cells,surface(),{requested:1.},
        coarse,periodic=False)
    assert requested in selected


def test_periodic_split_selection_retains_exact_maps_and_macro_jumps():
    _,_,meshers,*_ = w.imports()
    from crease_directional_refinement import edge_metric_squared
    from hp_shared_domain import periodic_jumps
    from homogenization import strain_tensor
    from nested_refinement import interpolation_check
    mesh = meshers.generate(lambda x,y,z:x*0-1,bounds=(-.5,.5)*3,cells=4,
        periodic=(True,)*3,order=2,threads=1)
    points,cells,faces = mesh.points,mesh.tetrahedra,mesh.surface
    pairs = np.sort(cells[:,EDGES].reshape(-1,2),axis=1)
    edges,first = np.unique(pairs,axis=0,return_index=True)
    mids = cells[:,4:].ravel()[first]
    def metric(points):
        return np.broadcast_to(np.diag([400.,400.,16.]),(len(points),3,3))
    values = np.maximum(edge_metric_squared(points,edges,mids,metric)-(4/3)**2,0)
    scores = {tuple(map(int,e)):float(v) for e,v in zip(edges,values,strict=True) if v>0}
    selected,_ = select_with_quality_closure(points,cells,faces,scores,metric)
    fields = np.array([points@strain_tensor(load).T for load in np.eye(6)])
    result = refine(points,cells,faces,fields,edge_scores=selected)
    assert result['diagnostics']['accepted_edge_orbits'] > 0
    interpolation_check(points,cells,fields,result)
    periodic_jumps(result['points'],result['surface'],result['displacement'])
