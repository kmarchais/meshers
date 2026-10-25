"""Replay local correction arrays on the exact audited hierarchical domain."""

import argparse
import json
import shutil
from pathlib import Path

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--audit', type=Path, required=True)
    args = parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np, *_ = w.imports()
    from scipy.linalg import eigh
    from periodic_crease_audit import features
    from hp_coupled_preflight import inherited_inputs
    audit = json.loads((args.audit/'result.json').read_text())
    assert audit['passed'] and audit['equilibrium']['fresh_operator']
    hashes = inherited_inputs(args.audit)
    for path in (Path(__file__),args.audit/'result.json',w.HERE/'periodic_crease_audit.py'):
        hashes[str(path.resolve())] = w.digest(path)
        if path.suffix == '.py':
            shutil.copyfile(path,args.run/f'input-{path.name}')
    w.save(args.run/'analysis-inputs.json',hashes)
    source = Path(audit['source'])
    mesh = np.load(source/'mesh.npz')
    groups = features(mesh)['groups']
    solved = json.loads((Path(audit['solve'])/'result.json').read_text())
    gain = np.asarray(solved['energy_gain_tensor'])
    parent = np.asarray(solved['predicted_energy_tensor'])+gain
    _, vectors = eigh(gain,parent)
    witness = vectors[:,-1]
    previous = None
    rules = []
    for count in (15,64,125):
        local = np.load(args.audit/f'local-{count}.npz')['correction_energy']
        assert local.shape == (len(mesh['tetrahedra']),6,6)
        minimum = float(np.linalg.eigvalsh(local).min())
        assert minimum >= -1e-12
        total = local.sum(axis=0)
        balance = float(np.linalg.norm(total-gain)/np.linalg.norm(gain))
        assert balance < .01
        change = None if previous is None else float(np.linalg.norm(local-previous)/np.linalg.norm(local))
        if change is not None:
            assert change < .01
        scalars = dict(normalized_six_load=np.einsum('ij,tji->t',np.linalg.inv(parent),local),
            largest_gain_direction=np.einsum('i,tij,j->t',witness,local,witness))
        localization = {}
        for label, scalar in scalars.items():
            assert scalar.min() > -1e-12 and scalar.sum() > 0
            localization[label] = {name:dict(elements=int(mask.sum()),
                element_fraction=float(mask.mean()),energy_fraction=float(scalar[mask].sum()/scalar.sum()))
                for name,mask in groups.items()}
        rules.append(dict(points=count,minimum_local_eigenvalue=minimum,
            global_gain_balance_relative=balance,local_vs_previous_relative=change,localization=localization))
        previous = local
    w.save(args.run/'result.json',dict(passed=True,audit=str(args.audit),source=str(source),rules=rules,
        largest_gain_direction=witness.tolist(),physical_nodes=len(mesh['points']),
        physical_tets=len(mesh['tetrahedra']),inputs=hashes,new_solve=False,reference_promoted=False,
        native_candidate_promoted=False,scope='Resolved correction distribution on this fixed domain, not remaining-error certification; groups overlap'))
    assert all(w.digest(Path(p)) == h for p,h in hashes.items())
    print(json.dumps(rules[-1]['localization']),flush=True)


if __name__ == '__main__':
    main()
