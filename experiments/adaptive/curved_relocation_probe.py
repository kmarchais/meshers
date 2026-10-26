"""Local interior Tet10 relocation feasibility with the cavity boundary fixed."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np

import microgen_workflow as w
from nested_refinement import EDGES, FACES, positive_quality
from quadratic_edge_removal_probe import boundary, measure, oriented


POLICY = dict(lattice_vertices=256, other_case_vertices=64, proposals_per_vertex=1,
    displacement='Half the vector to the mean of adjacent corners, as in native relocate',
    edge_midpoint_displacement='Half the corner displacement',
    minimum_corner_quality_gain=1e-6, median_corner_quality_must_not_decrease=True,
    sampled_motion_fractions=[.25, .5, .75, 1.], boundary_coordinates_fixed=True,
    fields_used=False, global_mesh_changed=False)


def propose(points, cells, node, surface_nodes):
    if node in surface_nodes:
        return None
    assert (cells[:, :4] == node).any(axis=1).all()
    position = np.argmax(cells[:, :4] == node, axis=1)
    adjacent = np.unique(cells[:, :4])
    adjacent = adjacent[adjacent != node]
    delta = .5 * (points[adjacent].mean(axis=0)-points[node])
    weights = np.zeros(cells.shape)
    weights[np.arange(len(cells)), position] = 1.
    weights[:, 4:] = .5*np.any(EDGES[None] == position[:, None, None], axis=2)
    moved = set(cells[weights > 0])
    assert not moved.intersection(surface_nodes)
    assert not moved.intersection({n for face in boundary(cells) for n in face})
    before = points[cells]
    after = before + weights[:, :, None]*delta
    q0, q1 = positive_quality(before[:, :4]), positive_quality(after[:, :4])
    if q1.min() <= q0.min()+1e-6 or np.median(q1) < np.median(q0):
        return None
    if abs(after).max() > .5+1e-12:
        return None
    return dict(node=int(node), cells=cells, before=before, after=after, delta=delta,
        initial_minimum_quality=float(q0.min()), final_minimum_quality=float(q1.min()),
        initial_median_quality=float(np.median(q0)), final_median_quality=float(np.median(q1)),
        all_recorded_local_quantiles_no_worse=bool(np.all(np.quantile(q1,[0,.01,.1,.5,.9,.99,1]) >=
            np.quantile(q0,[0,.01,.1,.5,.9,.99,1])-1e-12)))


def controls():
    # A complete warped octahedron star around a displaced interior vertex.
    points = [np.array([.07,.03,.01])]
    points.extend(.3*np.vstack((np.eye(3),-np.eye(3))))
    corners = [oriented(np.asarray(points),[0,1+i,2+j,3+k])
        for i in (0,3) for j in (0,3) for k in (0,3)]
    mids, cells = {}, []
    for c in corners:
        row = list(c)
        for a,b in EDGES:
            edge = tuple(sorted((int(c[a]),int(c[b]))))
            if edge not in mids:
                mids[edge] = len(points)
                points.append((points[edge[0]]+points[edge[1]])/2)
            row.append(mids[edge])
        cells.append(row)
    points, cells = np.array(points), np.array(cells)
    points[:,2] += .03*points[:,0]**2
    outer = {n for f in boundary(cells) for n in f}
    assert propose(points,cells,1,outer) is None
    row = propose(points,cells,0,outer)
    assert row is not None
    _, old_volume = measure(row['before'])
    for fraction in POLICY['sampled_motion_fractions']:
        minimum, volume = measure(row['before']+fraction*(row['after']-row['before']))
        assert minimum.min() > 0 and abs(volume.sum()/old_volume.sum()-1) < 1e-10


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    w.imports(); controls()
    receipt_path = w.ROOT/'handoff/macbook/adaptive-interior-validation-20261025.json'
    receipt = json.loads(receipt_path.read_text())
    sources = {name:w.ROOT/receipt['cases'][name]['run'] for name in ('gyroid','schwarz','spinodoid')}
    sources['lattice'] = args.source
    paths = [Path(__file__), w.HERE/'quadratic_edge_removal_probe.py',
        w.HERE/'nested_refinement.py', w.HERE/'quadrature_audit.py', receipt_path]
    for i,p in enumerate(paths):
        if p.suffix == '.py':shutil.copyfile(p,args.run/f'input-{i}-{p.name}')
    results = {}
    for label,source in sources.items():
        source_row = json.loads((source/'result.json').read_text())
        assert source_row['qualified'] and source_row['geometry']['geometry_admitted']
        with np.load(source/'mesh.npz') as mesh:
            points, cells, surface = (mesh[k] for k in ('points','tetrahedra','surface'))
        assert cells.shape[1] == 10 and surface.shape[1] == 6
        q = positive_quality(points[cells[:, :4]])
        boundary_nodes = set(surface.ravel())
        rank = np.full(len(points),np.inf)
        np.minimum.at(rank,cells[:,:4].ravel(),np.repeat(q,4))
        eligible = np.array([u for u in np.argsort(rank,kind='stable') if
            np.isfinite(rank[u]) and u not in boundary_nodes])
        selected = eligible[:POLICY['lattice_vertices'] if label == 'lattice' else POLICY['other_case_vertices']]
        proposals = []
        for node in selected:
            ids = np.flatnonzero((cells[:,:4] == node).any(axis=1))
            row = propose(points,cells[ids],int(node),boundary_nodes)
            if row is not None:proposals.append(row)
        admitted = []
        for row in proposals:
            _,old_volume = measure(row['before'],positive_rule=False)
            minimum, maximum_volume_error = np.inf, 0.
            for fraction in POLICY['sampled_motion_fractions']:
                jacobian,volume = measure(row['before']+fraction*(row['after']-row['before']))
                minimum = min(minimum,float(jacobian.min()))
                maximum_volume_error = max(maximum_volume_error,abs(float(volume.sum()/old_volume.sum()-1)))
            if minimum <= 0 or maximum_volume_error >= 1e-10:continue
            admitted.append({k:v for k,v in row.items() if k not in ('cells','before','after','delta')} |
                dict(delta=row['delta'].tolist(), incident_tets=len(row['cells']),
                    minimum_sampled_path_jacobian=minimum, maximum_relative_volume_difference=maximum_volume_error))
        results[label] = dict(source=str(source), eligible_interior_corners=len(eligible),
            examined_vertices=len(selected), corner_quality_improving_proposals=len(proposals),
            curved_path_admitted_local_proposals=len(admitted),
            all_local_quantiles_no_worse=sum(row['all_recorded_local_quantiles_no_worse'] for row in admitted),
            proposals=admitted)
        w.save(args.run/'partial-results.json',results)
        paths.extend(source/name for name in ('mesh.npz','result.json'))
    w.save(args.run/'result.json',dict(passed=True, warped_star_and_boundary_controls_passed=True,
        policy=POLICY, cases=results, old_approximation_space_preserved=False,
        global_mesh_changed=False, global_admission_performed=False, new_fields=False,
        scope='Potentially overlapping local star moves only. Exact cavity Tri6 traces and sampled positive curved Jacobians along the tested path, not a global mesh, nested FE transfer, accuracy result or timing claim.',
        inputs={str(p):w.digest(p) for p in paths}))
    print('Locally feasible curved moves:', {k:v['curved_path_admitted_local_proposals'] for k,v in results.items()},flush=True)


if __name__ == '__main__':main()
