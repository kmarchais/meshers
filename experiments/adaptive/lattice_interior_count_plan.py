"""Freeze an interior-node request from admitted counts, without stiffness data."""

import argparse
import json
import math
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--seed', type=Path, required=True)
    parser.add_argument('--actual', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    import numpy as np

    paths = [Path(__file__), args.seed/'result.json', args.seed/'mesh.npz',
        args.actual/'result.json', args.actual/'mesh.npz']
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    seed, actual = [json.loads((p/'result.json').read_text()) for p in (args.seed, args.actual)]
    for row in (seed, actual):
        assert row['case'] == 'smooth-lattice' and row['qualified'] and row['order'] == 2
        assert row['geometry']['geometry_admitted'] and max(row['nodes'], row['elements']) <= w.CAP
    assert seed.get('interior_refine_nodes', 0) == 0
    mesh = np.load(args.seed/'mesh.npz')
    points, tets, surface = mesh['points'], mesh['tetrahedra'][:, :4], mesh['surface'][:, :3]
    corners = np.unique(tets)
    np.testing.assert_array_equal(corners, np.arange(len(corners)))
    v = len(corners)
    edges = np.sort(np.concatenate([tets[:, [i, j]] for i in range(4) for j in range(i+1, 4)]), axis=1)
    edges, incidence = np.unique(edges, axis=0, return_counts=True)
    assert len(points) == v + len(edges)
    boundary = np.sort(np.concatenate([surface[:, [i, j]] for i, j in ((0, 1), (0, 2), (1, 2))]), axis=1)
    midpoint = points[edges].mean(axis=1)
    eligible = (~np.isin(edges[:, 0]*v+edges[:, 1], boundary[:, 0]*v+boundary[:, 1])
        & np.all((midpoint > -.5+1e-9) & (midpoint < .5-1e-9), axis=1))
    m = incidence[eligible]
    assert len(m) and m.min() >= 3
    lower, upper = math.ceil(.98*actual['elements']), actual['elements']
    center = .99*actual['elements']
    desired = round((center-seed['elements'])/float(m.mean()))
    safe = min((upper-seed['elements'])//int(m.max()),
        (actual['nodes']-seed['nodes'])//(int(m.max())+2))
    request = min(desired, safe)
    bounds = [seed['elements']+request*int(n) for n in (m.min(), m.max())]
    assert request > 0 and lower <= bounds[0] <= bounds[1] <= upper
    result = dict(seed=str(args.seed), actual=str(args.actual), requested_interior_nodes=request,
        seed_nodes=seed['nodes'], seed_tets=seed['elements'], eligible_edges=len(m),
        edge_star_valence_min=int(m.min()), edge_star_valence_max=int(m.max()),
        edge_star_valence_mean=float(m.mean()), target_tet_interval=[lower, upper],
        target_tet_center=center, tet_bounds_if_request_filled=bounds,
        node_upper_bound_if_request_filled=seed['nodes']+request*(int(m.max())+2),
        selection_uses_stiffness=False, new_solve=False, qualified=False,
        scope='Count plan only. Requires filled request, native preflight, all gates, resource check, and fresh six-load accuracy.',
        physical_domain='Unchanged scalar specimen; verify curved boundary after native regeneration.',
        quadratic_space_nesting_claim=False, inputs=hashes)
    shutil.copyfile(__file__, args.run/'analysis-driver.py')
    w.save(args.run/'result.json', result)
    assert all(w.digest(Path(p)) == value for p, value in hashes.items())
    print(f"Frozen interior request {request}; conditional tetrahedra bounds {bounds}", flush=True)


if __name__ == '__main__':
    main()
