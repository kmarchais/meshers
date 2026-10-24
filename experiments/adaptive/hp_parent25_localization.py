"""Localize the audited same-domain coupled correction, without another solve."""

import argparse
import json
from pathlib import Path
import shutil
import time

import numpy as np
import microgen_workflow as workflow


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from hp_coupled_preflight import inherited_inputs
    from hp_batched_energy import integrate_batched
    from hp_cross_energy import integrate
    from periodic_crease_audit import features
    from quadrature_audit import duffy_rule

    audit = json.loads((args.audit/'result.json').read_text())
    assert audit['passed'] and audit['parent25'] and audit['retained_model_count'] == 24
    source = Path(audit['source'])
    solved = json.loads((source/'result.json').read_text())
    assert workflow.digest(source/'hp-fields.npz') == solved['field_archive_sha256']
    consistent = json.loads((Path(solved['consistent_source'])/'result.json').read_text())
    trial = json.loads((Path(consistent['source'])/'result.json').read_text())
    inputs = inherited_inputs(args.audit)
    paths = [Path(__file__), args.audit/'result.json']
    paths.extend(workflow.HERE/name for name in ('hp_batched_energy.py','hp_cross_energy.py','periodic_crease_audit.py','quadrature_audit.py'))
    inputs.update({str(path):workflow.digest(path) for path in paths})
    workflow.save(args.output.parent/'analysis-inputs.json', inputs)
    shutil.copyfile(__file__, args.output.parent/'analysis-driver.py')
    mesh, cubic = np.load(Path(solved['nested_source'])/'mesh.npz'), np.load(solved['cubic_mesh'])
    h, p = (np.load(Path(solved[name])/'fields.npz')['displacement'] for name in ('nested_source','cubic_source'))
    fields = np.load(source/'hp-fields.npz')
    coefficients, responses = fields['mode_coefficients'], fields['response_displacement']
    maps = np.load(Path(solved['prepared'])/'prolongation.npz')
    q = np.einsum('ij,ind->jnd',coefficients,h+responses)
    p = np.einsum('ij,ind->jnd',coefficients,p)
    groups = features(mesh)['groups']
    selected = np.linspace(0,len(mesh['tetrahedra'])-1,64,dtype=int)
    basis = np.vstack((-np.eye(6),np.eye(6)))
    records = []
    for order in (5,6):
        begin = time.perf_counter()
        rule = duffy_rule(order)
        energy, volume, minimum, local, volumes = integrate_batched(mesh['points'],mesh['tetrahedra'],q,
            cubic['tetrahedra'],p,maps['parents'],maps['parent_barycentric_maps'],rule)
        total = local.sum(axis=0)
        np.testing.assert_allclose(total,basis.T@energy@basis,atol=1e-12)
        reference = next(row for row in trial['rules'] if row['positive_rule_points']==order**3)
        base = np.asarray(reference['energy_matrix'])[:6,:6]
        corrected = next(row for row in audit['rules'] if row['points']==order**3)
        expected = base-np.asarray(corrected['energy_tensor'])
        identity_error = float(np.linalg.norm(total-expected)/np.linalg.norm(expected))
        assert identity_error <= .01
        scores = np.einsum('ij,tji->t',np.linalg.inv(base),local)
        assert scores.min() > -1e-12 and np.linalg.eigvalsh(local).min() > -1e-12
        grouped = integrate(mesh['points'],mesh['tetrahedra'][selected],q,cubic['tetrahedra'],p,
            maps['parents'][selected],maps['parent_barycentric_maps'][selected],rule,local_coefficients=basis)
        np.testing.assert_allclose(grouped[3],local[selected],atol=1e-12,rtol=1e-9)
        localization = {name:dict(tets=int(mask.sum()),volume_fraction=float(volumes[mask].sum()/volume),
            normalized_correction_fraction=float(scores[mask].sum()/scores.sum())) for name,mask in groups.items()}
        path = args.output.parent/f'local-{order**3}.npz'
        np.savez_compressed(path,correction_energy=local,volumes=volumes,normalized_trace=scores)
        record = dict(points=order**3,total_correction_energy=total.tolist(),galerkin_identity_relative_error=identity_error,
            localization=localization,volume=volume,minimum_sampled_jacobian=minimum,grouped_control_cells=64,
            archive_sha256=workflow.digest(path),seconds=time.perf_counter()-begin)
        records.append(record)
        workflow.save(args.output.parent/f'rule-{order**3}.json',record)
        print(f'Correction localization {order**3}: Galerkin relative error {identity_error:.9g}',flush=True)
    precision = float(np.linalg.norm(np.asarray(records[1]['total_correction_energy'])-records[0]['total_correction_energy'])/np.linalg.norm(records[1]['total_correction_energy']))
    assert precision <= .01
    workflow.save(args.output,dict(passed=True,parent25_coupled_localization=True,audit=str(args.audit),rules=records,
        correction_rule_relative_change=precision,precision_gate=.01,new_equilibrium_solve=False,reference_promoted=False,
        scope='Resolved same-domain six-load correction energy only; overlapping groups are not additive. No remaining-error, continuum, native-candidate or workflow timing claim.'))
    assert all(workflow.digest(Path(path))==expected for path,expected in inputs.items())


if __name__ == '__main__':
    main()
