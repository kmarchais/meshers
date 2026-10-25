"""Screen residual modes of the complete audited field on its unchanged domain."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--audit', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, *_ = w.imports()
    import hierarchical_residual_trial as ht
    from hp_coupled_preflight import inherited_inputs
    from periodic import tensor_difference
    from quadrature_audit import duffy_rule
    audited = json.loads((args.audit/'result.json').read_text())
    assert audited['passed'] and audited['equilibrium']['fresh_operator']
    solved_path = Path(audited['solve'])
    solved = json.loads((solved_path/'result.json').read_text())
    source, original_trial = Path(solved['source']), Path(solved['trial'])
    paths = [Path(__file__),Path(ht.__file__),w.HERE/'test_hierarchical_residual_trial.py',
             args.audit/'result.json',solved_path/'coupled-fields.npz',original_trial/'modes.npz']
    hashes = {**inherited_inputs(args.audit), **{str(p.resolve()):w.digest(p) for p in paths}}
    w.save(args.run/'analysis-inputs.json',hashes)
    for i,path in enumerate(paths):
        if path.suffix == '.py':
            shutil.copyfile(path,args.run/f'input-{i}-{path.name}')
    mesh = np.load(source/'mesh.npz')
    points,cells,surface = (mesh[k] for k in ('points','tetrahedra','surface'))
    fields = np.load(source/'fields.npz')['displacement']
    corrected = np.load(solved_path/'coupled-fields.npz')
    alpha = corrected['mode_coefficients']
    nodal = fields-np.einsum('mni,ml->lni',corrected['response_displacement'],alpha)
    original = np.load(original_trial/'modes.npz')
    parent_coefficients = original['coefficients']@alpha
    space = ht.topology(points,cells,surface)
    np.testing.assert_array_equal(space['ids'],original['feature_ids'])
    np.testing.assert_array_equal(space['signs'],original['signs'])
    print('Construct residual modes of the complete enriched field using positive64',flush=True)
    coefficients,residual,construction = ht.build_modes(points,cells,nodal,space,duffy_rule(4),
        parent_coefficients=parent_coefficients)
    direct = np.einsum('fcl,fcm->lm',residual,coefficients)
    del residual
    traces = ht.trace_check(space,coefficients)
    np.savez_compressed(args.run/'modes.npz',coefficients=coefficients,feature_ids=space['ids'],signs=space['signs'])
    rules = []
    previous = None
    for order in (4,5):
        print(f'Integrate complete parent and new residual modes with positive{order**3}',flush=True)
        energy,local,geometry = ht.integrate(points,cells,nodal,space,coefficients,duffy_rule(order),
            parent_coefficients=parent_coefficients)
        if order == 4:
            np.testing.assert_allclose(energy[:6,6:],direct,atol=1e-12,rtol=1e-8)
        parent_difference = tensor_difference(energy[:6,:6],np.asarray(solved['predicted_energy_tensor']))
        assert parent_difference['worst_relative_strain_energy_difference'] < 1e-4
        trial = ht.trial(energy)
        gain = np.asarray(trial['gain_tensor'])
        change = None if previous is None else float(np.linalg.norm(gain-previous)/np.linalg.norm(gain))
        if change is not None:
            assert change <= .01
        a = np.asarray(trial['coefficients'])
        local_energy = np.einsum('mi,tmn,nj->tij',a,local,a,optimize=True)
        np.testing.assert_allclose(local_energy.sum(axis=0),gain,atol=1e-12,rtol=1e-8)
        np.savez_compressed(args.run/f'local-{order**3}.npz',correction_energy=local_energy)
        rules.append(dict(points=order**3,geometry=geometry,parent_vs_audited=parent_difference,
            trial=trial,gain_vs_previous=change))
        previous = gain
    models = []
    for model in audited['models']:
        path = Path(model['source'])
        reference = np.asarray(json.loads(path.read_text())['models'][model['index']]['stiffness'])
        models.append(dict(source=str(path),index=model['index'],difference=
            tensor_difference(np.asarray(rules[-1]['trial']['trial_energy_tensor']),reference)))
    assert len(models) == 24
    worst = max(m['difference']['worst_relative_strain_energy_difference'] for m in models)
    gain_max = rules[-1]['trial']['maximum_relative_gain']
    for name in ('modes.npz','local-64.npz','local-125.npz'):
        path = args.run/name
        hashes[str(path.resolve())] = w.digest(path)
    w.save(args.run/'result.json',dict(passed=True,audit=str(args.audit),source=str(source),
        physical_nodes=len(points),physical_tets=len(cells),geometry_unchanged=True,
        parent_field_includes_all_bubbles=True,old_approximation_space_preserved=True,
        new_trial_scalar_unknowns=6,polynomial_degree_unchanged=True,construction=construction,
        traces=traces,rules=rules,maximum_additional_trial_gain=gain_max,
        gain_signal_screen=.0005,gain_signal_passed=gain_max >= .0005,retained_model_count=24,
        models=models,trial_worst_model_disagreement=worst,trial_passes_one_percent_all_models=worst <= .01,
        equilibrium_solve_performed=False,reference_promoted=False,native_candidate_promoted=False,
        inputs=hashes,scope='Residual screen of complete equilibrated enriched field; a fixed-parent trial, not remaining-error bound or mesh workflow'))
    assert all(w.digest(Path(p)) == h for p,h in hashes.items())
    print(f'Additional fixed-parent trial gain {gain_max:.9%}; all-model disagreement {worst:.9%}',flush=True)


if __name__ == '__main__':
    main()
