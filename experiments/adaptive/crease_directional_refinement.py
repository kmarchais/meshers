"""Geometry-only directional refinement with exact parent Tet10 map restriction."""

import numpy as np
import shutil

import microgen_workflow as w
from nested_refinement import EDGES, refine, interpolation_check
from periodic_crease_grading import make_metric, directional_metric


def edge_metric_squared(points, edges, mids, tensor):
    delta = points[edges[:, 1]]-points[edges[:, 0]]
    matrices = tensor(points)
    eigenvalues = np.linalg.eigvalsh(matrices)
    assert np.isfinite(matrices).all() and eigenvalues.min() > 0
    assert (eigenvalues[:, -1]/eigenvalues[:, 0]).max() <= 25*(1+1e-12)
    return np.maximum.reduce([np.einsum('ni,nij,nj->n', delta, matrices[ids], delta)
        for ids in (edges[:, 0], mids, edges[:, 1])])


def generate(seed, function, args, seed_volume):
    from types import SimpleNamespace
    from homogenization import strain_tensor
    from hp_shared_domain import periodic_jumps
    from periodic import periodic_pairs

    scalar, metadata, cloud = make_metric(seed, function, args.min_length, args.max_length, args.cells)
    tensor = directional_metric(scalar, cloud, function, args.max_length, 5.)
    points, cells, surface = seed.points, seed.tetrahedra, seed.surface
    # These nonlinear periodic controls test prolongation; they are not FE solutions.
    controls = np.array([points@strain_tensor(load).T for load in np.eye(6)])
    for load in range(6):
        axes = np.roll(np.arange(3), load % 3)
        modulation = .01*np.sin(2*np.pi*points[:, axes[0]])*np.cos(2*np.pi*points[:, axes[1]])
        controls[load, :, load % 3] += modulation
    # Guard the fixed-point loop; disjoint stars may need separate directional passes.
    maximum_passes = 3*int(np.ceil(np.log2((1/args.cells)/args.min_length)))
    if getattr(args, 'crease_longest_closure', False) or getattr(args, 'crease_quality_closure', False):
        # Closure must first split longer neighbors. The 15-pass pilot still
        # accepted 507 orbits at its guard, without a quality/budget rejection.
        maximum_passes *= 4
    assert maximum_passes > 0
    history = []
    directory = args.run/'directional-refinement'
    directory.mkdir(exist_ok=False)
    for step in range(maximum_passes):
        pairs = np.sort(cells[:, EDGES].reshape(-1, 2), axis=1)
        edges, first = np.unique(pairs, axis=0, return_index=True)
        mids = cells[:, 4:].ravel()[first]
        squared = edge_metric_squared(points, edges, mids, tensor)
        values = np.maximum(squared-(4/3)**2, 0)
        scores = {tuple(map(int, edge)): float(value) for edge, value in zip(edges, values, strict=True) if value > 0}
        if not scores:
            history.append(dict(step=step+1, status='metric-satisfied', nodes=len(points), tets=len(cells)))
            break
        closure = None
        if getattr(args, 'crease_quality_closure', False):
            from quality_split_closure import select_with_quality_closure
            scores, closure = select_with_quality_closure(points,cells,surface,scores,tensor)
            closure_key = 'quality_split_closure'
        elif getattr(args, 'crease_longest_closure', False):
            from longest_edge_closure import closure_scores
            scores, closure = closure_scores(points, cells, surface, scores)
            closure_key = 'longest_edge_closure'
        result = refine(points, cells, surface, controls, maximum_nodes=w.CAP,
            maximum_tets=w.CAP, minimum_quality=.1, edge_scores=scores)
        check = interpolation_check(points, cells, controls, result)
        diagnostics = result['diagnostics']
        diagnostics['selection'] = 'Geometry-only directional edge-metric excess per added-node upper bound'
        if closure is not None:
            diagnostics[closure_key] = closure
            diagnostics['selection'] = closure['selection']+'; inherited metric excess per node bound'
        diagnostics['fe_residual_scores_used'] = False
        diagnostics['score_units'] = 'Squared normalized metric-length excess, dimensionless'
        points, cells, surface, controls = (result[k] for k in ('points', 'tetrahedra', 'surface', 'displacement'))
        assert max(len(points), len(cells)) <= w.CAP
        part = directory/f'level-{step+1}'
        part.mkdir(exist_ok=False)
        for suffix in ('-samples.npz', '-specimen.json'):
            shutil.copyfile(args.run/(args.case+suffix), part/(args.case+suffix))
        np.savez_compressed(part/'mesh.npz', points=points, tetrahedra=cells, surface=surface)
        np.savez_compressed(part/'parent-map.npz', parents=result['parents'],
            parent_barycentric_maps=result['parent_barycentric_maps'])
        geometry = w.geometry(part, args.case, points, cells, surface)
        record = dict(step=step+1, nodes=len(points), tets=len(cells), diagnostics=diagnostics,
            mapping_and_field_reproduction=check, geometry=geometry,
            relative_volume_change=geometry['volume']/seed_volume-1,
            parent_map_sha256=w.digest(part/'parent-map.npz'),
            mesh_sha256=w.digest(part/'mesh.npz'))
        w.save(part/'result.json', record)
        periodic = periodic_jumps(points, surface, controls)
        record['control_periodicity'] = periodic
        history.append(record)
        w.save(directory/'history.json', history)
        assert geometry['geometry_admitted'] and abs(record['relative_volume_change']) < 1e-10
        if diagnostics['accepted_edge_orbits'] == 0:
            if getattr(args, 'crease_quality_closure', False):
                history[-1]['status'] = ('node-or-tet-budget-blocked' if diagnostics['budget_rejected_orbits']
                    else 'quality-continuation-blocked' if not scores else 'no-admissible-periodic-orbit')
            else:
                history[-1]['status'] = 'no-admissible-orbit-within-remaining-budget'
            break
    w.save(directory/'history.json', history)
    np.savez_compressed(directory/'prolonged-controls.npz', displacement=controls)
    pairs = np.sort(cells[:, EDGES].reshape(-1, 2), axis=1)
    edges, first = np.unique(pairs, axis=0, return_index=True)
    mids = cells[:, 4:].ravel()[first]
    remaining_ratio = np.sqrt(edge_metric_squared(points, edges, mids, tensor))
    _, _, final_periodic = periodic_pairs(points, surface)
    quality_policy = None
    if getattr(args, 'crease_quality_closure', False):
        from quality_split_closure import QUALITY_POLICY
        quality_policy = QUALITY_POLICY
    metadata.update(maximum_passes=maximum_passes, directional_ratio=5., history=history,
        longest_edge_closure=getattr(args, 'crease_longest_closure', False),
        quality_split_closure=getattr(args, 'crease_quality_closure', False),
        quality_split_policy=quality_policy,
        edge_metric_sampling='Maximum at both endpoints and quadratic edge midpoint',
        remaining_edges_above_metric_split_threshold=int((remaining_ratio > 4/3).sum()),
        maximum_remaining_metric_edge_ratio=float(remaining_ratio.max()),
        metric_completely_satisfied=bool((remaining_ratio <= 4/3).all()),
        control_fields='Six affine macrostrains plus deterministic nonlinear periodic nodal perturbations; no equilibrium fields',
        saved_fe_fields_replayed=False, geometry_map_restricted_exactly=True,
        approximation_spaces_nested=True, scalar_gradation_applies_to_seed_only=True,
        final_periodic=final_periodic, native_only=False, python_exact_map_refinement=True)
    return SimpleNamespace(points=points, tetrahedra=cells, surface=surface,
        diagnostics=dict(source='Fresh native seed plus experimental exact-map directional refinement',
            seed=seed.diagnostics, directional_refinement=metadata)), metadata
