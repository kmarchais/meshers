"""Collect small lattice diagnostics, checking their complete recorded inputs."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh

import microgen_workflow as w
from cubic_refinement import MODEL_REPORTS


def collect():
    inputs={str(Path(__file__)):w.digest(Path(__file__))}
    def read(name):
        path=w.HERE/'runs'/name/'result.json'
        record=json.loads(path.read_text());inputs[str(path)]=w.digest(path)
        for name,expected in record['inputs'].items():
            source=Path(name)
            assert w.digest(source)==expected,name
            inputs[name]=expected
        return record
    control=read('macbook-traction-controls-v1')
    assert all(control['controls'].values())
    samples=[]
    for version in (1,2):
        sample=read(f'macbook-lattice-traction-sample-v{version}')
        dense=read(f'macbook-lattice-traction-quadrature-v{version}')
        assert dense['all_groups_quadrature_stable'] and dense['unchanged_inputs']
        assert dense['maximum_independent_stress_relative_error']<1e-7
        routes={}
        for label,groups in dense['routes'].items():
            total=sum(g['rules'][-1]['normalized_trace'] for g in groups)
            rows=[dict(group=g['group'],population=g['population'],sampled=g['sampled'],
                       score=g['rules'][-1]['normalized_trace'],share=g['rules'][-1]['normalized_trace']/total,
                       quadrature_change=g['rules'][-1]['gram_relative_change']) for g in groups]
            routes[label]=dict(total_residual_score=total,groups=rows,
                curve_incident_share=sum(g['share'] for g in rows if g['group'] in ('interior_curve','wall_curve')),
                bulk_interface_share=next(g['share'] for g in rows if g['group']=='interior_bulk'))
        samples.append(dict(seed=sample['seed'],routes=routes,
            independent_stress_relative_error=dense['maximum_independent_stress_relative_error']))
    mono=read('macbook-lattice-monotone-screen-v1')
    assert mono['retained_models']==24 and all(mono['controls'].values())
    models=[]
    for name in MODEL_REPORTS:
        path=w.HERE/'runs'/name
        models.extend(np.asarray(m['stiffness']) for m in json.loads(path.read_text())['models'])
    evidence_path=w.ROOT/'handoff/macbook/comparison-findings-20261026.json'
    evidence=json.loads(evidence_path.read_text());inputs[str(evidence_path)]=w.digest(evidence_path)
    def budget(row,candidate):
        low=min(1+m['minimum_relative_difference'] for m in row['models'])
        high=max(1+m['maximum_relative_difference'] for m in row['models'])
        alpha_low=.99/low;alpha_high=min(1.,1.01/high)
        alpha=min(1.,2/(low+high))
        worst=max(abs(alpha*low-1),abs(alpha*high-1))
        replay=max(float(abs(eigvalsh(alpha*candidate-m,m)).max()) for m in models)
        np.testing.assert_allclose(worst,replay,atol=1e-13)
        return dict(current_disagreement=row['current_worst_disagreement'],
            necessary_monotone_error_floor=row['necessary_error_floor'],
            per_model_floor=[r['monotone_error_floor'] for r in row['models']],
            witness=row['worst_witness'],uniform_softening_interval=[alpha_low,alpha_high],
            interval_feasible=alpha_low<=alpha_high,
            algebraic_balanced_scale=alpha,algebraic_balanced_disagreement=worst,
            physical_field_or_mesh_constructed=False)
    budgets={}
    for row in evidence['lattice_stiffness']['rows']:
        source=Path(row['source']);candidate=np.asarray(json.loads(source.read_text())['stiffness'])
        budgets[row['route']]=budget(mono['rows'][row['route']],candidate)
    audit_path=w.HERE/'runs/macbook-lattice-longest-edge-closure-coupled-audit-v1/result.json'
    audit=json.loads(audit_path.read_text())
    budgets['complete-enriched-diagnostic']=[dict(points=r['points'],**budget(s,np.asarray(r['energy_tensor'])))
        for s,r in zip(mono['rows']['complete-enriched-diagnostic']['rules'],audit['rules'],strict=True)]
    workers={}
    for name in ('macbook-traction-controls-v1','macbook-lattice-traction-sample-v1',
                 'macbook-lattice-traction-quadrature-v1','macbook-lattice-traction-sample-v2',
                 'macbook-lattice-traction-quadrature-v2','macbook-lattice-monotone-screen-v1'):
        path=w.HERE/'runs'/name/'worker.json';worker=json.loads(path.read_text());inputs[str(path)]=w.digest(path)
        assert worker['returncode']==0 and worker['limit_failure'] is None
        workers[name]=dict(seconds=worker['wall_s'],peak_gib=worker['peak_rss_bytes']/1024**3)
    return dict(samples=samples,monotone_budgets=budgets,resources=workers,inputs=inputs,
        model_tensor_hashes=mono['model_tensor_hashes'],all_24_models_retained=True,
        all_source_meshes_and_fields_unchanged=True,new_FE_solve=False,reference_promoted=False,
        new_native_candidate=False,full_objective_complete=False,
        conclusions=[
            'Native curve-incident share is about52.45% in both samples, leaving about47.55% elsewhere and14.65-14.91% on bulk interfaces.',
            'Actual has larger residual scores despite better retained-model agreement. Residual size is not a replacement stiffness target.',
            'Fixed-parent decreasing corrections are not obstructed by the necessary monotone floor, but no attainable correction has been constructed.',
            'Uniform algebraic softening cannot bring compact or Uniform37 under1%; the best saved and enriched matrices have a narrow feasible algebraic interval. No material, field or geometry scaling is authorized or performed.'
        ],
        next_step='Test a genuinely equilibrated stress reconstruction on a tiny control with free-wall, periodic and rigid-body compatibility before using an estimator to allocate mesh changes. Do not repeat closed bubble modes, curve-only refinement or global sizing sweeps.',
        scope='Read-only small diagnostics and algebraic sensitivity controls. Traction scores are not stiffness errors or continuum bounds; scaled matrices are not physical FE results or promoted references.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    assert not args.output.exists();w.save(args.output,collect())
