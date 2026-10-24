"""Positive-rule saved-field trial on the verified parent25 physical domain."""

import argparse
import json
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.linalg import eigvalsh
import microgen_workflow as workflow
from hp_parent25_probe import PREPARED, H_SOLVE, P_MESH, P_SOLVE


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--benchmark', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from hp_batched_energy import integrate_batched
    from hp_cross_energy import combine
    from hp_coupled_preflight import inherited_inputs
    from quadrature_audit import duffy_rule
    from periodic import tensor_difference
    from periodic_crease_audit import features
    from cubic_refinement import MODEL_REPORTS

    benchmark = json.loads((args.benchmark/'result.json').read_text())
    assert benchmark['passed'] and benchmark['conservative_positive_rule_time_estimate_s'] < 600
    inputs = inherited_inputs(args.benchmark)
    paths = [Path(__file__), args.benchmark/'result.json', workflow.HERE/'hp_batched_energy.py',
        workflow.HERE/'periodic_crease_audit.py', workflow.ROOT/'crates/meshers-python/examples/periodic.py']
    paths.extend(workflow.HERE/'runs'/name for name in MODEL_REPORTS)
    inputs.update({str(path): workflow.digest(path) for path in paths})
    workflow.save(args.output.parent/'analysis-inputs.json', inputs)
    shutil.copyfile(__file__, args.output.parent/'analysis-driver.py')
    shutil.copyfile(workflow.HERE/'hp_batched_energy.py', args.output.parent/'hp_batched_energy.py')
    mesh, cubic = (dict(np.load(directory/'mesh.npz')) for directory in (H_SOLVE, P_MESH))
    h, p = (np.load(directory/'fields.npz')['displacement'] for directory in (H_SOLVE, P_SOLVE))
    maps = dict(np.load(PREPARED/'prolongation.npz'))
    groups = features(mesh)['groups']
    records = []
    for order in (5, 6):
        begin = time.perf_counter()
        energy, volume, minimum, local, volumes = integrate_batched(mesh['points'], mesh['tetrahedra'], h,
            cubic['tetrahedra'], p, maps['parents'], maps['parent_barycentric_maps'], duffy_rule(order))
        combined = combine(energy)
        coefficients = np.asarray(combined['coefficients'])
        local_gain = np.einsum('ij,tik,kl->tjl', coefficients, local, coefficients, optimize=True)
        np.testing.assert_allclose(local_gain.sum(axis=0), combined['gain_tensor'], atol=1e-12, rtol=1e-9)
        normalized_gain = np.einsum('ij,tji->t', np.linalg.inv(energy[:6,:6]), local_gain)
        assert normalized_gain.min() > -1e-12
        localization = {name: dict(tets=int(mask.sum()), volume_fraction=float(volumes[mask].sum()/volume),
            normalized_gain_fraction=float(normalized_gain[mask].sum()/normalized_gain.sum())) for name, mask in groups.items()}
        record = dict(positive_rule_points=order**3, seconds=time.perf_counter()-begin, energy_matrix=energy.tolist(),
            volume=volume, minimum_sampled_jacobian=minimum, localization=localization, **combined)
        for name, block, source in (('nested', energy[:6,:6], H_SOLVE), ('cubic', energy[6:,6:], P_SOLVE)):
            difference = tensor_difference(block, np.asarray(json.loads((source/'result.json').read_text())['stiffness']))
            assert difference['worst_relative_strain_energy_difference'] < 1e-4
            record[name+'_own_energy_difference'] = difference
        models = []
        for name in MODEL_REPORTS:
            path = workflow.HERE/'runs'/name
            for index, model in enumerate(json.loads(path.read_text())['models']):
                tensor = np.asarray(model['stiffness'])
                values = eigvalsh(np.asarray(combined['trial_energy_tensor'])-tensor, tensor)
                models.append(dict(source=str(path), index=index, signed_relative_energy_eigenvalues=values.tolist(),
                    worst_relative_difference=float(abs(values).max())))
        assert len(models) == 24
        record.update(models=models, retained_model_count=24,
            worst_retained_model_disagreement=max(row['worst_relative_difference'] for row in models))
        np.savez_compressed(args.output.parent/f'local-rule-{order**3}.npz', mode_gram=local, volumes=volumes,
            normalized_gain=normalized_gain)
        record['local_archive_sha256'] = workflow.digest(args.output.parent/f'local-rule-{order**3}.npz')
        workflow.save(args.output.parent/f'rule-{order**3}.json', record)
        records.append(record)
        print(f'Rule {order**3}: gain {combined["maximum_relative_gain"]:.9g}, 24-model error {record["worst_retained_model_disagreement"]:.9g}', flush=True)
    gain_precision = float(np.linalg.norm(np.asarray(records[1]['gain_tensor'])-records[0]['gain_tensor'])/np.linalg.norm(records[1]['gain_tensor']))
    energy_precision = tensor_difference(np.asarray(records[1]['trial_energy_tensor']), np.asarray(records[0]['trial_energy_tensor']))
    assert gain_precision <= .01 and energy_precision['worst_relative_strain_energy_difference'] < 1e-4
    workflow.save(args.output, dict(passed=True, benchmark=str(args.benchmark), rules=records,
        gain_rule_relative_change=gain_precision, positive_rule_energy_difference=energy_precision,
        precision_gate=.01, new_equilibrium_solve=False, reference_promoted=False, native_candidate_promoted=False,
        scope='Fixed-domain saved-field trial in six periodic p-minus-h modes. Localized gain is the reduction from this subspace, not remaining error or continuum certification. No live workflow timing claim.'))
    assert all(workflow.digest(Path(name)) == expected for name, expected in inputs.items())


if __name__ == '__main__':
    main()
