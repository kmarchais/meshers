"""Read-only allocation and quality evidence for choosing the next Meshers change."""
import argparse
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np

import microgen_workflow as w
from nested_refinement import EDGES, child_maps, positive_quality
from crease_directional_refinement import edge_metric_squared
from periodic_crease_grading import make_metric, directional_metric
from periodic_crease_audit import features


def cut_edges(cells, parents, mappings):
    counts = np.bincount(parents, minlength=len(cells))
    assert ((counts == 1) | (counts == 2)).all()
    refined = np.flatnonzero(counts == 2)
    _, first = np.unique(parents, return_index=True)
    maps = mappings[first[refined]]
    fractional = np.sum(maps == .5, axis=2) == 2
    assert (fractional.sum(axis=1) == 1).all()
    bary = maps[np.arange(len(maps)), fractional.argmax(axis=1)]
    local = np.sort(np.argsort(bary, axis=1)[:, -2:], axis=1)
    offset = np.argmax(np.all(local[:, None] == np.sort(EDGES, axis=1)[None], axis=2), axis=1)
    assert np.array_equal(np.sort(EDGES[offset], axis=1), local)
    edges = np.sort(cells[refined[:, None], local], axis=1)
    mids = cells[refined, 4 + offset]
    return refined, edges, mids


def controls():
    assert quantiles(np.array([])) is None
    cells = np.arange(10)[None]
    for edge in range(6):
        parents = np.zeros(2, dtype=int)
        selected, pairs, mids = cut_edges(cells, parents, np.asarray(child_maps(edge)))
        np.testing.assert_array_equal(selected, [0])
        np.testing.assert_array_equal(pairs, np.sort(EDGES[edge:edge+1], axis=1))
        np.testing.assert_array_equal(mids, [4+edge])


def quantiles(values):
    if len(values) == 0:
        return None
    return dict(zip(('min', 'p01', 'p10', 'median', 'p90', 'p99', 'max'),
        map(float, np.quantile(values, [0, .01, .1, .5, .9, .99, 1]))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--run', required=True, type=Path)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    w.imports()
    from periodic import periodic_pairs
    controls()
    paths = [Path(__file__), w.HERE/'nested_refinement.py',
        w.HERE/'crease_directional_refinement.py', w.HERE/'periodic_crease_grading.py',
        w.HERE/'periodic_crease_audit.py', w.HERE/'microgen_workflow.py']
    for i, path in enumerate(paths):
        shutil.copyfile(path, args.run/f'input-{i}-{path.name}')
    source = json.loads((args.source/'result.json').read_text())
    assert source['qualified'] and source['quadratic_quality_pass_requested']
    seed = args.source/'feature-seed'
    with np.load(seed/'mesh.npz') as mesh:
        old = {k: mesh[k] for k in mesh.files}
    case, _, _ = w.specimen('smooth-lattice')
    scalar, metadata, cloud = make_metric(SimpleNamespace(points=old['points'],
        tetrahedra=old['tetrahedra'], surface=old['surface']), case.function, .002, .1, 26)
    assert metadata['law'] == source['crease_bisection']['law']
    frozen_cloud = w.HERE/'runs/macbook-lattice-crease-distance-preflight-v1/crease-cloud.npz'
    with np.load(frozen_cloud) as mesh:
        np.testing.assert_array_equal(cloud, mesh['canonical_points'])
    paths.extend((frozen_cloud, seed/'mesh.npz', args.source/'result.json'))
    tensor = directional_metric(scalar, cloud, case.function, .1, 5.)
    history = json.loads((args.source/'directional-refinement/history.json').read_text())
    levels = []
    for record in history:
        directory = args.source/'directional-refinement'/f'level-{record["step"]}'
        with np.load(directory/'mesh.npz') as mesh:
            new = {k: mesh[k] for k in mesh.files}
        with np.load(directory/'parent-map.npz') as mapping:
            parents, maps = mapping['parents'], mapping['parent_barycentric_maps']
        np.testing.assert_array_equal(old['points'], new['points'][:len(old['points'])])
        refined, edges, mids = cut_edges(old['tetrahedra'], parents, maps)
        assert len(refined) == len(new['tetrahedra'])-len(old['tetrahedra'])
        ratio = np.sqrt(edge_metric_squared(old['points'], edges, mids, tensor))
        unique, indices = np.unique(edges, axis=0, return_index=True)
        parent_ratio = np.full(len(old['tetrahedra']), np.nan)
        parent_ratio[refined] = ratio
        ids = new['tetrahedra'].ravel()
        added = ids >= len(old['points'])
        owner = np.repeat(parent_ratio[parents], 10)[added]
        assert np.isfinite(owner).all()
        added_ids = ids[added]-len(old['points'])
        n = len(new['points'])-len(old['points'])
        low, high = np.full(n, np.inf), np.zeros(n)
        np.minimum.at(low, added_ids, owner)
        np.maximum.at(high, added_ids, owner)
        assert np.isfinite(low).all() and (high > 0).all()
        limit = 4/3*(1+1e-12)
        levels.append(dict(level=record['step'], nodes_added=n, tets_added=len(refined),
            physical_cut_edges=len(unique), cut_metric_ratio_quantiles=quantiles(ratio),
            cut_edges_already_below_split_threshold=int((ratio[indices] <= limit).sum()),
            tets_added_on_already_below_threshold_cuts=int((ratio <= limit).sum()),
            nodes_added_only_for_already_below_threshold_cuts=int((high <= limit).sum()),
            nodes_added_with_mixed_cut_attribution=int(((low <= limit) & (high > limit)).sum())))
        paths.extend(directory/name for name in ('mesh.npz', 'parent-map.npz', 'result.json'))
        w.save(args.run/'closure-levels.json', levels)
        old = new
    # This is the saved pre-pass mesh; changed topology is never called nested.
    with np.load(args.source/'quality-source/mesh.npz') as mesh:
        for key in old: np.testing.assert_array_equal(old[key], mesh[key])
    quality_routes = {}
    assessment = json.loads((w.HERE/'runs/macbook-lattice-uniform37-assessment-v1/result.json').read_text())
    routes = dict(compact=seed, uniform37=w.ROOT/assessment['source'],
        longest_edge=args.source/'quality-source', curved_quality=args.source,
        actual=w.HERE/'runs/macbook-lattice-actual-live-solve-v1')
    for label, directory in routes.items():
        with np.load(directory/'mesh.npz') as mesh:
            points, cells, surface = (mesh[k] for k in ('points', 'tetrahedra', 'surface'))
        feature = features(dict(points=points, tetrahedra=cells, surface=surface))
        groups = feature['groups']
        # Explicitly disjoint partition, unlike the overlapping feature groups.
        curves = groups['curve_vertex_incident']
        bulk = groups['neither_wall_nor_cap_vertex']
        masks = dict(curve_incident=curves, other_boundary=~curves & ~bulk, interior=bulk)
        assert np.all(sum(mask.astype(int) for mask in masks.values()) == 1)
        q = positive_quality(points[cells[:, :4]])
        row = dict(nodes=len(points), tets=len(cells), overall_quality=quantiles(q), regions={})
        for name, mask in masks.items():
            row['regions'][name] = dict(tets=int(mask.sum()), quality=quantiles(q[mask]),
                fraction_of_all_tets_below_05=float(((q < .5) & mask).sum()/max(1, (q < .5).sum())))
        if label == 'curved_quality':
            pairs = np.sort(cells[:, EDGES].reshape(-1, 2), axis=1)
            edges, first = np.unique(pairs, axis=0, return_index=True)
            mids = cells[:, 4:].ravel()[first]
            ratio = np.sqrt(edge_metric_squared(points, edges, mids, tensor))
            interior_node = np.ones(len(points), dtype=bool)
            interior_node[surface.ravel()] = False
            interior_edges = interior_node[edges].all(axis=1)
            short = (ratio < .8) & interior_edges
            row['metric'] = dict(edges=len(edges), above_split_threshold=int((ratio > 4/3).sum()),
                maximum_ratio=float(ratio.max()), short_interior_edges=int(short.sum()),
                short_all_edges=int((ratio < .8).sum()),
                scope='Length-only coarsening candidates, not link/quality/new-edge/curved-Jacobian-admitted collapses')
            _, _, periodic = periodic_pairs(points, surface)
            row['periodic_pair_counts'] = periodic['periodic_pair_counts']
        paths.append(directory/'mesh.npz')
        quality_routes[label] = row
    sums = {key: sum(row[key] for row in levels) for key in (
        'nodes_added', 'tets_added', 'physical_cut_edges', 'cut_edges_already_below_split_threshold',
        'tets_added_on_already_below_threshold_cuts', 'nodes_added_only_for_already_below_threshold_cuts',
        'nodes_added_with_mixed_cut_attribution')}
    paths.extend((args.source/'directional-refinement/history.json',
        w.HERE/'runs/macbook-lattice-uniform37-assessment-v1/result.json'))
    result = dict(passed=True, source=str(args.source), controls_passed=6, levels=levels,
        closure_totals=sums, routes=quality_routes, identical_frozen_metric_cloud=True,
        scope='Saved-mesh geometry-only allocation and corner-quality diagnostics. Extra closure cuts may be necessary for admissible refinement; they are not assumed removable. No new mesh, FE solution, timing win or reference promotion.',
        inputs={str(p):w.digest(p) for p in paths})
    w.save(args.run/'result.json', result)
    print('Closure totals:', sums, flush=True)
    print('Remaining metric:', quality_routes['curved_quality']['metric'], flush=True)


if __name__ == '__main__':
    main()
