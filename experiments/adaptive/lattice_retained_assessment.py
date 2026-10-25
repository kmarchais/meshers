"""Assess one fully field-audited lattice tensor against all 24 frozen models."""

import argparse
import hashlib
import json
from pathlib import Path

import microgen_workflow as w
from cubic_refinement import MODEL_REPORTS


def canonical_ordered_mesh(mesh, np):
    points = mesh['points']
    assert len(np.unique(points, axis=0)) == len(points)
    order = np.lexsort(points[:, ::-1].T)
    ranks = np.empty(len(points), dtype=np.int64)
    ranks[order] = np.arange(len(points))
    result = dict(points=points[order])
    for key in ('tetrahedra', 'surface'):
        # Preserve all local corner/mid-edge associations; sort only whole rows.
        cells = ranks[mesh[key]]
        result[key] = cells[np.lexsort(cells[:, ::-1].T)]
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--preflight', type=Path, required=True)
    parser.add_argument('--field-audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    np, *_ = w.imports()
    from periodic import tensor_difference
    source = args.source
    row = json.loads((source/'result.json').read_text())
    assert row['case'] == 'smooth-lattice' and row['order'] == 2
    assert row['qualified'] and row['geometry']['geometry_admitted']
    assert row['status'] == 'solved-needs-reference' and row['completed_loads'] == 6
    assert len(row['solver_records']) == 6
    assert max(r['relative_residual'] for r in row['solver_records']) < 2e-9
    assert row['nodes'] <= w.CAP and row['elements'] <= w.CAP
    manifest = json.loads((source/'manifest.json').read_text())
    assert manifest['architecture'] == 'arm64' and manifest['threads'] == 1
    assert manifest['native_sha256'] == manifest['native_sha256_after']
    worker = json.loads((source/'worker.json').read_text())
    assert worker['returncode'] == 0 and worker['limit_failure'] is None
    assert worker['memory_bytes'] == w.MEMORY and worker['memory_enforcement'] == 'process-tree-rss-watchdog'
    audit_manifest = json.loads((args.field_audit/'manifest.json').read_text())
    audited = {Path(p).resolve(): h for p, h in audit_manifest['inputs'].items()}
    hashes = {str(Path(__file__).resolve()): w.digest(Path(__file__))}
    for path, expected in audit_manifest['source_sha256'].items():
        assert w.digest(Path(path)) == expected
        hashes[path] = expected
    for name in ('result.json', 'mesh.npz', 'fields.npz'):
        path = source/name
        assert w.digest(path) == audited[path.resolve()]
        hashes[str(path)] = w.digest(path)
    own_path = args.field_audit/(source.name+'-consistency.json')
    own = json.loads(own_path.read_text())
    assert max(own[k] for k in ('reaction_relative_error', 'hill_mandel_relative_error',
        'integrated_stiffness_skew_relative')) < 1e-7
    assert own['stress_reaction_mixed_load_difference']['worst_relative_strain_energy_difference'] < 1e-7
    assert max(own['periodic_jump_max'], own['periodic_mid_edge_jump_max']) < 1e-9
    quad_path = args.field_audit/(source.name+'.json')
    quad = json.loads(quad_path.read_text())
    quad_worker = args.field_audit/(source.name+'-worker.json')
    assert json.loads(quad_worker.read_text())['returncode'] == 0
    assert [r['points'] for r in quad['rules']] == [64, 125]
    sensitivity = max(r['energy_vs_solved_tensor']['worst_relative_strain_energy_difference'] for r in quad['rules'])
    assert sensitivity <= 1e-4
    assert quad['rules'][-1]['energy_vs_previous_rule']['worst_relative_strain_energy_difference'] <= 1e-4
    previous = w.HERE/'runs/g3-lattice-adapt26-tol0028-amg-solve-v1'
    old = json.loads((previous/'result.json').read_text())
    for key in ('case', 'cells', 'min_length', 'order', 'native_tolerance', 'optimization_passes',
        'adaptive', 'adapt_tolerance', 'solve_node_cap', 'solver_quadrature_points'):
        assert old[key] == row[key]
    assert old['void_clamp']['positive_ceiling'] == row['void_clamp']['positive_ceiling'] == .02
    assert old['symmetric_union']['identity'] == row['symmetric_union']['identity']
    for suffix in ('-samples.npz', '-specimen.json'):
        name = row['case']+suffix
        assert w.digest(source/name) == w.digest(previous/name) == w.digest(args.preflight/name)
        hashes[str(source/name)] = w.digest(source/name)
        hashes[str(previous/name)] = w.digest(previous/name)
    meshes = [np.load(run/'mesh.npz') for run in (source, args.preflight)]
    canonical = [canonical_ordered_mesh(mesh, np) for mesh in meshes]
    assert all(np.array_equal(canonical[0][key], canonical[1][key]) for key in canonical[0])
    fields = np.load(source/'fields.npz')['displacement']
    assert fields.shape == (6, row['nodes'], 3) and np.isfinite(fields).all()
    tensor = np.asarray(row['stiffness'])
    assert tensor.shape == (6, 6) and np.isfinite(tensor).all() and np.linalg.eigvalsh(tensor).min() > 0
    models, model_bytes = [], []
    for filename in MODEL_REPORTS:
        path = w.HERE/'runs'/filename
        report = json.loads(path.read_text())
        assert report['case'] == row['case'] and len(report['models']) == 8 and report['stiffness_target'] == .01
        hashes[str(path)] = w.digest(path)
        for index, model in enumerate(report['models']):
            matrix = np.asarray(model['stiffness'], dtype=np.float64)
            assert matrix.shape == (6, 6) and np.isfinite(matrix).all() and np.linalg.eigvalsh(matrix).min() > 0
            model_bytes.append(hashlib.sha256(matrix.tobytes()).hexdigest())
            models.append(dict(source=str(path), index=index, length_proxy=model['length_proxy'],
                assumed_order=model['assumed_order'], difference=tensor_difference(tensor, matrix),
                previous_difference=tensor_difference(np.asarray(old['stiffness']), matrix)))
    assert len(models) == 24
    errors = [m['difference']['worst_relative_strain_energy_difference'] for m in models]
    old_errors = [m['previous_difference']['worst_relative_strain_energy_difference'] for m in models]
    for path in (source/'manifest.json', source/'worker.json', args.preflight/'result.json',
        args.preflight/'mesh.npz', previous/'result.json', previous/'mesh.npz', previous/'fields.npz',
        own_path, quad_path, quad_worker, args.field_audit/'manifest.json'):
        hashes[str(path)] = w.digest(path)
    passed = max(errors) <= .01
    w.save(args.output, dict(source=str(source), field_audit_passed=True, complete_six_load_fields=True,
        fixed_recipe_and_specimen=True, preflight_mesh_identical_up_to_numbering=True,
        full_tet10_and_tri6_local_associations_preserved=True,
        retained_model_count=24, retained_model_tensor_hashes=model_bytes, retained_models_unchanged=True,
        models=models, minimum_model_disagreement=min(errors), worst_model_disagreement=max(errors),
        previous_worst_model_disagreement=max(old_errors), passes_one_percent_all_models=passed,
        quadrature_sensitivity=sensitivity, quadrature_sensitivity_limit=1e-4,
        workflow_seconds=row['total_workflow_s'], peak_rss_gib=worker['peak_rss_bytes']/1024**3,
        geometry=row['geometry'], previous_geometry=old['geometry'],
        tensor_change_from_previous=tensor_difference(tensor, np.asarray(old['stiffness'])),
        reference_promoted=False, comparative_timing_claim=False, full_objective_complete=False,
        inputs=hashes, scope='Fresh same-host lattice workflow and complete retained-model assessment; no continuum, local stress or four-case timing victory'))
    print(f'All 24 models retained: worst disagreement {max(errors):.9%}; passes 1% {passed}; no reference promoted', flush=True)


if __name__ == '__main__':
    main()
