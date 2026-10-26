"""Necessary retained-model error floor for any tensor below a fixed parent."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from scipy.linalg import eigh, eigvalsh

import microgen_workflow as w
from cubic_refinement import MODEL_REPORTS


def screen(candidate,models):
    rows=[]
    for index,model in enumerate(models):
        values,vectors=eigh(candidate-model,model)
        load=vectors[:,0]
        observed=float(load@(candidate-model)@load/(load@model@load))
        np.testing.assert_allclose(observed,values[0],atol=1e-13,rtol=1e-12)
        rows.append(dict(model=index,minimum_relative_difference=float(values[0]),
                         maximum_relative_difference=float(values[-1]),
                         monotone_error_floor=max(0.,-float(values[0])),
                         engineering_strain=load.tolist(),rayleigh_replay=observed))
    worst=max(rows,key=lambda row:row['monotone_error_floor'])
    return dict(models=rows,necessary_error_floor=worst['monotone_error_floor'],
                worst_witness=worst,blocked_at_one_percent=worst['monotone_error_floor']>.01+1e-12,
                current_worst_disagreement=max(max(abs(r['minimum_relative_difference']),abs(r['maximum_relative_difference'])) for r in rows))


def controls():
    eye=np.eye(6)
    assert screen(1.02*eye,[eye])['necessary_error_floor']==0
    np.testing.assert_allclose(screen(.98*eye,[eye])['necessary_error_floor'],.02,atol=1e-14)
    candidate=np.diag([.98,1.02,1.,1.,1.,1.])
    assert screen(candidate,[eye])['blocked_at_one_percent']
    transform=np.random.default_rng(3).normal(size=(6,6))+3*eye
    np.testing.assert_allclose(screen(transform.T@candidate@transform,[transform.T@transform])['necessary_error_floor'],.02,atol=1e-13)
    assert not screen(.995*eye,[eye])['blocked_at_one_percent']
    return dict(too_stiff_not_blocked=True,too_soft_floor_detected=True,mixed_direction_detected=True,
                nonorthogonal_basis_invariant=True,below_target_not_blocked=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    assert not args.output.exists()
    models=[];paths=[Path(__file__),w.HERE/'cubic_refinement.py',
        w.ROOT/'handoff/macbook/monotone-screen-selection-20261026.json',
        w.ROOT/'handoff/macbook/comparison-findings-20261026.json']
    evidence=json.loads(paths[-1].read_text())
    for name,expected in evidence['lattice_stiffness']['inputs'].items():assert w.digest(Path(name))==expected
    for name in MODEL_REPORTS:
        path=w.HERE/'runs'/name;paths.append(path)
        models.extend(np.asarray(row['stiffness']) for row in json.loads(path.read_text())['models'])
    assert len(models)==24
    hashes=[hashlib.sha256(matrix.tobytes()).hexdigest() for matrix in models]
    assert hashes==evidence['lattice_stiffness']['retained_model_tensor_hashes']
    rows={}
    for row in evidence['lattice_stiffness']['rows']:
        source=Path(row['source']);paths.append(source)
        candidate=np.asarray(json.loads(source.read_text())['stiffness'])
        rows[row['route']]=screen(candidate,models)
        np.testing.assert_allclose(rows[row['route']]['current_worst_disagreement'],row['worst_model_disagreement'],atol=1e-12)
    enriched=w.HERE/'runs/macbook-lattice-longest-edge-closure-coupled-audit-v1/result.json'
    audit=json.loads(enriched.read_text());assert audit['passed'];paths.append(enriched)
    rules=[]
    for rule in audit['rules']:
        candidate=np.asarray(rule['energy_tensor'])
        rules.append(dict(points=rule['points'],**screen(candidate,models)))
    rows['complete-enriched-diagnostic']=dict(rules=rules)
    inputs={str(path):w.digest(path) for path in paths}
    shutil.copyfile(__file__,args.output.parent/'input-lattice_monotone_screen.py')
    w.save(args.output,dict(controls=controls(),rows=rows,retained_models=24,
        model_tensor_hashes=hashes,inputs=inputs,new_solve=False,reference_promoted=False,
        full_objective_complete=False,
        scope='Algebraic necessary floor for any symmetric tensor X with X<=C in Loewner order: directional too-soft errors cannot improve by lowering stiffness. Conditional fixed-parent diagnostic, not a geometry-change obstruction, attainable FE tensor or continuum bound. A floor below1% does not prove feasibility.'))


if __name__=='__main__':main()
