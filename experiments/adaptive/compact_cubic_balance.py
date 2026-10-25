"""Normalize retained models correctly and localize the same-domain cubic correction."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import microgen_workflow as w


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    args=parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np,*_=w.imports()
    from hp_cross_energy import integrate
    from periodic import tensor_difference
    from periodic_crease_audit import features
    from quadrature_audit import duffy_rule
    from scipy.linalg import eigh
    root=w.HERE/'runs'
    parent=root/'macbook-periodic-edge-normal-lattice26-amg-v1'
    cubic=root/'macbook-periodic-edge-normal-cubic-solve-v1'
    prepared=root/'macbook-periodic-edge-normal-cubic-gates-v1'
    audit=root/'macbook-periodic-edge-normal-cubic-assessment-v1'
    audited=json.loads((audit/'result.json').read_text())
    assert audited['fresh_process_audit_passed'] and audited['retained_model_count']==24
    rows=[json.loads((p/'result.json').read_text()) for p in (parent,cubic)]
    assert all(r['qualified'] and r['completed_loads']==6 for r in rows)
    paths=[Path(__file__),audit/'result.json',audit/'manifest.json',prepared/'result.json',prepared/'mesh.npz']
    paths.extend(p/name for p in (parent,cubic) for name in ('result.json','fields.npz'))
    paths.extend(w.HERE/name for name in ('hp_cross_energy.py','periodic_crease_audit.py','cubic_reference.py','quadrature_audit.py'))
    paths.append(parent/'mesh.npz')
    hashes={str(p.resolve()):w.digest(p) for p in paths}
    for i,p in enumerate(paths):
        if p.suffix=='.py':shutil.copyfile(p,args.run/f'input-{i}-{p.name}')
    c2,c3=(np.asarray(r['stiffness']) for r in rows)
    models=[]
    for m in audited['models']:
        path=Path(m['source']);report=json.loads(path.read_text())
        matrix=np.asarray(report['models'][m['index']]['stiffness'],dtype=np.float64)
        difference=tensor_difference(c3,matrix)
        inverse=np.sort(1/(1+np.asarray(m['difference_from_cubic']['relative_energy_difference_eigenvalues']))-1)
        np.testing.assert_allclose(difference['relative_energy_difference_eigenvalues'],inverse,atol=1e-12)
        hashes[str(path.resolve())]=w.digest(path)
        models.append(dict(source=str(path),index=m['index'],difference=difference,
            matrix_sha256=hashlib.sha256(matrix.tobytes()).hexdigest()))
    worst=max(m['difference']['worst_relative_strain_energy_difference'] for m in models)
    mesh=np.load(parent/'mesh.npz');pmesh=np.load(prepared/'mesh.npz')
    offset=len(mesh['points']);count=offset+len(pmesh['points'])
    points=np.vstack((mesh['points'],pmesh['points']))
    h=np.zeros((6,count,3));p=np.zeros_like(h)
    h[:,:offset]=np.load(parent/'fields.npz')['displacement']
    p[:,offset:]=np.load(cubic/'fields.npz')['displacement']
    parents=np.arange(len(mesh['tetrahedra']))
    maps=np.repeat(np.eye(4)[None],len(parents),axis=0)
    coefficient=np.vstack((np.eye(6),-np.eye(6)))
    gain=c2-c3
    _,vectors=eigh(gain,c2);witness=vectors[:,-1]
    groups=features(mesh)['groups']
    rules=[];previous=None
    for order in (4,5):
        print(f'Cross-integrate complete same-domain fields with positive{order**3}',flush=True)
        energy,volume,minimum,local,local_volume=integrate(points,mesh['tetrahedra'],h,
            pmesh['tetrahedra']+offset,p,parents,maps,duffy_rule(order),local_coefficients=coefficient)
        correction=coefficient.T@energy@coefficient
        np.testing.assert_allclose(local.sum(axis=0),correction,atol=1e-12,rtol=1e-9)
        balance=tensor_difference(energy[:6,:6],energy[6:,6:]+correction)
        assert balance['worst_relative_strain_energy_difference']<1e-4
        relative_gain_error=float(np.linalg.norm(correction-gain)/np.linalg.norm(gain))
        assert relative_gain_error<.01 and np.linalg.eigvalsh(local).min()>-1e-12
        change=None if previous is None else float(np.linalg.norm(local-previous)/np.linalg.norm(local))
        if change is not None:assert change<.01
        scalars=dict(normalized_six_load=np.einsum('ij,tji->t',np.linalg.inv(c2),local),
            largest_gain_direction=np.einsum('i,tij,j->t',witness,local,witness))
        localization={label:{name:dict(elements=int(mask.sum()),element_fraction=float(mask.mean()),
            energy_fraction=float(scalar[mask].sum()/scalar.sum())) for name,mask in groups.items()}
            for label,scalar in scalars.items()}
        np.savez_compressed(args.run/f'local-{order**3}.npz',correction_energy=local,volume=local_volume)
        rules.append(dict(points=order**3,minimum_jacobian=minimum,volume=volume,
            galerkin_balance=balance,gain_relative_error=relative_gain_error,local_vs_previous=change,
            localization=localization))
        previous=local
    for order in (4,5):
        path=args.run/f'local-{order**3}.npz';hashes[str(path.resolve())]=w.digest(path)
    w.save(args.run/'analysis-inputs.json',hashes)
    w.save(args.run/'result.json',dict(passed=True,parent=str(parent),cubic=str(cubic),prepared=str(prepared),
        models=models,retained_model_count=24,candidate_normalized_by_each_model=True,
        worst_model_disagreement=worst,passes_one_percent_all_models=worst<=.01,rules=rules,
        correction_tensor=gain.tolist(),largest_gain_witness=witness.tolist(),inputs=hashes,
        reference_promoted=False,native_cubic_api=False,comparative_timing_claim=False,
        scope='Same-domain full-field Galerkin/localization and correct candidate-relative sensitivity scores; not remaining-error certification'))
    assert all(w.digest(Path(p))==v for p,v in hashes.items())
    print(f'All24 candidate-relative worst {worst:.9%}; no reference promoted',flush=True)


if __name__=='__main__':main()
