"""Audit saved exact-map refinement and identify quality-blocked metric edges."""

import argparse
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np

import microgen_workflow as w


def prolong(points, cells, root_cells, fields, parents, maps):
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4

    element, linear = Tet10(), Tet4()
    nodes, first = np.unique(cells.ravel(), return_index=True)
    np.testing.assert_array_equal(nodes, np.arange(len(points)))
    children, local = np.divmod(first, 10)
    barycentric = linear.shape_function(element.xi_nd)
    values = np.empty((6, len(points), 3))
    for start in range(0, len(points), 2048):
        take = slice(start, start+2048)
        root = np.einsum('na,nab->nb', barycentric[local[take]], maps[children[take]])
        shape = element.shape_function(root @ linear.xi_nd)
        values[:, take] = np.einsum('tn,ltnd->ltd', shape,
            fields[:, root_cells[parents[children[take]]]])
    return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--field-audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    w.imports()
    from hp_shared_domain import compose, periodic_jumps
    from nested_refinement import interpolation_check, EDGES, positive_quality
    from quadrature_audit import duffy_rule
    from tet10_integration import integrated_tensors
    from periodic import tensor_difference, periodic_pairs
    from periodic_crease_grading import make_metric, directional_metric
    from crease_directional_refinement import edge_metric_squared
    from gradient_tie_probe import symmetric_lattice

    inputs = [args.source/n for n in ('result.json', 'mesh.npz', 'fields.npz')]
    inputs += [args.parent/n for n in ('result.json', 'mesh.npz', 'fields.npz')]
    inputs += [args.source/'feature-seed/mesh.npz',
        args.field_audit/(args.source.name+'.json')]
    inputs += [w.HERE/n for n in ('directional_nesting_audit.py', 'nested_refinement.py',
        'hp_shared_domain.py', 'quadrature_audit.py', 'tet10_integration.py',
        'periodic_crease_grading.py', 'crease_directional_refinement.py')]
    row = json.loads((args.source/'result.json').read_text())
    root = np.load(args.parent/'mesh.npz')
    seed = np.load(args.source/'feature-seed/mesh.npz')
    assert all(np.array_equal(root[k], seed[k]) for k in root.files)
    mesh = np.load(args.source/'mesh.npz')
    points, cells, surface = (mesh[k] for k in ('points', 'tetrahedra', 'surface'))
    parents, maps = np.arange(len(root['tetrahedra'])), np.repeat(np.eye(4)[None], len(root['tetrahedra']), axis=0)
    for level in row['crease_bisection']['history']:
        path = args.source/'directional-refinement'/f'level-{level["step"]}'/'parent-map.npz'
        assert w.digest(path) == level['parent_map_sha256']
        inputs.append(path)
        mapping = np.load(path)
        parents, maps = compose(parents, maps, mapping['parents'], mapping['parent_barycentric_maps'])
    hashes = {str(p): w.digest(p) for p in inputs}
    w.save(args.output.parent/'analysis-inputs.json', hashes)
    shutil.copyfile(__file__, args.output.parent/'analysis-driver.py')
    old = np.load(args.parent/'fields.npz')['displacement']
    lifted = prolong(points, cells, root['tetrahedra'], old, parents, maps)
    check = interpolation_check(root['points'], root['tetrahedra'], old,
        dict(points=points, tetrahedra=cells, displacement=lifted,
            parents=parents, parent_barycentric_maps=maps))
    periodic = periodic_jumps(points, surface, lifted)
    solved = np.load(args.source/'fields.npz')['displacement']
    # Constant gauges have zero energy; keep the original fields for the identity.
    _, parent_energy, _ = integrated_tensors(points, cells, lifted, duffy_rule(5))
    _, correction_energy, _ = integrated_tensors(points, cells, lifted-solved, duffy_rule(5))
    audited = json.loads((args.field_audit/(args.source.name+'.json')).read_text())
    fine_energy = np.asarray(audited['rules'][-1]['fixed_field_energy_tensor'])
    balance = tensor_difference(parent_energy, fine_energy+correction_energy)
    assert balance['worst_relative_strain_energy_difference'] <= 1e-4
    assert np.linalg.eigvalsh(correction_energy).min() >= -1e-12
    parent_row = json.loads((args.parent/'result.json').read_text())
    inherited = tensor_difference(parent_energy, np.asarray(parent_row['stiffness']))
    assert inherited['worst_relative_strain_energy_difference'] <= 1e-4
    def field(x, y, z):
        return np.minimum(symmetric_lattice(x, y, z), .02)
    scalar, _, cloud = make_metric(SimpleNamespace(**dict(seed)), field, .002, .1, 26)
    tensor = directional_metric(scalar, cloud, field, .1, 5.)
    pairs = np.sort(cells[:, EDGES].reshape(-1, 2), axis=1)
    edges, first, inverse = np.unique(pairs, axis=0, return_index=True, return_inverse=True)
    mids = cells[:, 4:].ravel()[first]
    metric = edge_metric_squared(points, edges, mids, tensor)
    wanted = metric > (4/3)**2
    assert int(wanted.sum()) == row['crease_bisection']['remaining_edges_above_metric_split_threshold']
    allowed = np.ones((len(cells), 6), dtype=bool)
    for offset, endpoints in enumerate(EDGES):
        for endpoint in endpoints:
            coordinates = points[cells[:, :4]].copy()
            coordinates[:, endpoint] = points[cells[:, offset+4]]
            allowed[:, offset] &= positive_quality(coordinates) >= .1+1e-12
    incidence = inverse.reshape(-1, 6)
    blocked = wanted[incidence] & ~allowed
    blocking_cells = blocked.any(axis=1)
    _, classes, _ = periodic_pairs(points, surface)
    groups = {}
    for i, edge in enumerate(edges):
        groups.setdefault(tuple(sorted(classes[edge])), []).append(i)
    bad_edges = np.bincount(inverse, weights=(~allowed).ravel(), minlength=len(edges)) > 0
    wanted_groups = [g for g in groups.values() if wanted[g].any()]
    obstruction = dict(requested_edges=int(wanted.sum()),
        requested_orbits=len(wanted_groups),
        quality_blocked_requested_orbits=sum(bool(bad_edges[g].any()) for g in wanted_groups),
        blocking_cells=int(blocking_cells.sum()),
        blocking_cells_with_some_legal_alternative=int((blocking_cells & allowed.any(axis=1)).sum()),
        blocking_cells_with_no_legal_edge=int((blocking_cells & ~allowed.any(axis=1)).sum()),
        scope='Legal alternatives are per-element only; no conforming periodic closure is inferred')
    w.save(args.output, dict(source=str(args.source), parent=str(args.parent),
        inputs=hashes, exact_seed_arrays=True, mapping_and_field_reproduction=check,
        periodic=periodic, galerkin_energy_balance=balance, inherited_energy_difference=inherited,
        correction_energy_tensor=correction_energy.tolist(), quality_obstruction=obstruction,
        passed=True, reference_promoted=False, mesh_changed=False,
        scope='Read-only fixed-domain nesting and obstruction audit; no new solve or sizing choice'))
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    print(f'Nested energy balance {balance["worst_relative_strain_energy_difference"]:.3g}; '
        f'{obstruction["quality_blocked_requested_orbits"]}/{len(wanted_groups)} requested orbits quality-blocked', flush=True)


if __name__ == '__main__':
    main()
