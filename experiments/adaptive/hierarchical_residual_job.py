"""One immutable-field cubic residual trial; no meshing or equilibrium solve."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--assessment', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, *_ = w.imports()
    from fedoo.lib_elements.tetrahedron import Tet10
    from periodic import tensor_difference
    from scipy.linalg import eigh
    import hierarchical_residual_trial as ht
    from periodic_crease_audit import features
    from quadrature_audit import duffy_rule
    checked = json.loads(args.assessment.read_text())
    assert checked['field_audit_passed'] and checked['complete_six_load_fields']
    assert checked['retained_model_count'] == 24 and checked['retained_models_unchanged']
    assert Path(checked['source']).resolve() == args.source.resolve()
    paths = [Path(__file__), Path(ht.__file__), w.HERE/'test_hierarchical_residual_trial.py',
        w.HERE/'hp_cross_energy.py', w.HERE/'quadrature_audit.py', w.HERE/'periodic_crease_audit.py',
        Path(w.__file__), args.assessment,
        w.ROOT.parent/'periodic-microgen/experiments/periodic-fea/periodic.py']
    paths.extend(args.source/name for name in ('result.json', 'mesh.npz', 'fields.npz',
        'manifest.json', 'smooth-lattice-samples.npz', 'smooth-lattice-specimen.json'))
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    for path, expected in checked['inputs'].items():
        assert w.digest(Path(path)) == expected, path
        hashes[str(Path(path).resolve())] = expected
    w.save(args.run/'analysis-inputs.json', hashes)
    for i, path in enumerate(paths):
        if path.suffix == '.py':
            shutil.copyfile(path, args.run/f'input-{i}-{path.name}')
    row = json.loads((args.source/'result.json').read_text())
    mesh = np.load(args.source/'mesh.npz')
    points, cells, surface = (mesh[k] for k in ('points', 'tetrahedra', 'surface'))
    fields = np.load(args.source/'fields.npz')['displacement']
    assert fields.shape == (6, len(points), 3) and len(points) <= 250000 and len(cells) <= 250000
    assert abs(ht.basis(Tet10().xi_nd)[0]).max() < 1e-15
    space = ht.topology(points, cells, surface)
    w.save(args.run/'preflight.json', dict(nodes=len(points), elements=len(cells),
        geometry_unchanged=True, parent_fields_unchanged=True, new_mesh_created=False,
        cubic_edge_features=space['edges'], cubic_face_features=space['faces'],
        feature_coefficient_entries=space['features']*3*6,
        extra_scalar_trial_unknowns=6, full_cubic_space_solved=False,
        build_rule=64, positive_audit_rules=[64,125], gain_precision_limit=.01,
        gain_signal_screen=.0005, stiffness_target=.01,
        scope='Six fixed global trial functions, not independent higher-order nodal unknowns or a native mesh candidate'))
    print('Build six same-domain periodic residual modes with positive 64-point integration', flush=True)
    coefficients, residual, construction = ht.build_modes(points, cells, fields, space, duffy_rule(4))
    direct_cross = np.einsum('fcl,fcm->lm', residual, coefficients)
    del residual
    traces = ht.trace_check(space, coefficients)
    np.savez_compressed(args.run/'modes.npz', feature_ids=space['ids'], signs=space['signs'],
        coefficients=coefficients, corner_classes=space['corner_classes'], face_ids=space['face_ids'])
    rules, local_grams = [], []
    previous_gain = None
    for order in (4,5):
        print(f'Integrate fixed parent/mode fields with {order**3} positive points', flush=True)
        energy, local, geometry = ht.integrate(points, cells, fields, space, coefficients, duffy_rule(order))
        if order == 4:
            np.testing.assert_allclose(energy[:6,6:], direct_cross, atol=1e-12, rtol=1e-8)
        parent_difference = tensor_difference(energy[:6,:6], np.asarray(row['stiffness']))
        assert parent_difference['worst_relative_strain_energy_difference'] <= 1e-4
        assert abs(geometry['volume']/row['geometry']['volume']-1) <= 1e-8
        optimized = ht.trial(energy)
        gain = np.asarray(optimized['gain_tensor'])
        change = None if previous_gain is None else float(np.linalg.norm(gain-previous_gain)/np.linalg.norm(gain))
        rules.append(dict(points=order**3, geometry=geometry, parent_vs_reactions=parent_difference,
            trial=optimized, gain_vs_previous_rule=change, parent_energy_tensor=energy[:6,:6].tolist()))
        local_grams.append(local)
        previous_gain = gain
    assert rules[-1]['gain_vs_previous_rule'] <= .01
    alpha = np.asarray(rules[-1]['trial']['coefficients'])
    local_energy = [np.einsum('mi,tmn,nj->tij', alpha, g, alpha, optimize=True) for g in local_grams]
    assert np.linalg.norm(local_energy[1]-local_energy[0])/np.linalg.norm(local_energy[1]) <= .01
    gain = np.asarray(rules[-1]['trial']['gain_tensor'])
    np.testing.assert_allclose(local_energy[-1].sum(axis=0), gain, atol=1e-12, rtol=1e-8)
    parent = np.asarray(rules[-1]['parent_energy_tensor'])
    _, directions = eigh(gain, parent)
    witness = directions[:, -1]
    scalars = dict(six_load_normalized_trace=np.einsum('ij,tji->t', np.linalg.inv(parent), local_energy[-1]),
        largest_gain_direction=np.einsum('i,tij,j->t', witness, local_energy[-1], witness))
    groups = features(mesh)['groups']
    localization = {}
    for name, scalar in scalars.items():
        assert scalar.min() > -1e-12 and scalar.sum() > 0
        localization[name] = {label: dict(elements=int(mask.sum()),
            energy_fraction=float(scalar[mask].sum()/scalar.sum())) for label, mask in groups.items()}
    for order, local in zip((4,5), local_energy, strict=True):
        np.savez_compressed(args.run/f'local-{order**3}.npz', correction_energy=local)
    trial_tensor = np.asarray(rules[-1]['trial']['trial_energy_tensor'])
    models = []
    for model in checked['models']:
        path = Path(model['source'])
        reference = np.asarray(json.loads(path.read_text())['models'][model['index']]['stiffness'])
        models.append(dict(source=str(path), index=model['index'], difference=tensor_difference(trial_tensor, reference)))
    assert len(models) == 24
    worst = max(m['difference']['worst_relative_strain_energy_difference'] for m in models)
    gain_max = rules[-1]['trial']['maximum_relative_gain']
    for name in ('modes.npz', 'local-64.npz', 'local-125.npz'):
        hashes[str((args.run/name).resolve())] = w.digest(args.run/name)
    w.save(args.run/'result.json', dict(passed=True, source=str(args.source), native_candidate=False,
        geometry_unchanged=True, old_approximation_space_preserved=True, original_nodes=len(points),
        original_elements=len(cells), extra_scalar_trial_unknowns=6,
        potential_cubic_features=space['features'], mode_construction=construction,
        complete_cubic_trace_check=traces, rules=rules, gain_signal_passed=gain_max >= .0005,
        gain_signal_screen=.0005, positive_rule_precision_passed=True,
        largest_gain_witness=witness.tolist(), localization=localization,
        retained_model_count=24, models=models, trial_worst_model_disagreement=worst,
        trial_passes_one_percent_all_models=worst <= .01, equilibrium_solve_performed=False,
        reference_promoted=False, comparative_timing_claim=False, inputs=hashes,
        scope='Six-mode admissible fixed-domain energy trial, not a newly equilibrated tensor, native workflow or continuum/local stress certificate'))
    assert all(w.digest(Path(p)) == h for p,h in hashes.items())
    print(f'Fixed-domain trial gain {gain_max:.9%}; all-model disagreement {worst:.9%}; no reference promoted', flush=True)


if __name__ == '__main__':
    main()
