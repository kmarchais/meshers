"""Aggregate repeated mesh-only measurements and replay their actual artifacts."""

import argparse
import json
from pathlib import Path
import statistics
import time

import microgen_workflow as w


def distribution(values):
    values=[x for x in values if x is not None]
    return None if not values else dict(median=statistics.median(values),minimum=min(values),maximum=max(values))


def verify_run(run,replay_geometry=False):
    np,*_=w.imports()
    from periodic import periodic_pairs
    from facet_incidence_audit import audit
    import scorecard
    report=json.loads((run/'summary.json').read_text())
    manifest=json.loads((run/'manifest.json').read_text())
    assert manifest['no_fe_solve'] and manifest['threads']==1
    hashes=json.loads((run/'analysis-inputs.json').read_text())
    resolved=[]
    for name,expected in hashes.items():
        path=Path(name)
        if w.digest(path)!=expected and path.suffix=='.py':
            matches=[p for p in run.glob('*.py') if w.digest(p)==expected]
            assert matches,f'Missing measured source {path}'
            path=matches[0]
        assert w.digest(path)==expected,str(path)
        resolved.append(str(path))
    artifacts=0
    replays=[]
    for row in report['rows']:
        directory=Path(row['directory'])
        actual=json.loads((directory/'result.json').read_text())
        assert not actual['solve_performed'] and actual['vtk_threads']==1
        for key in ('status','qualified','config','generation_s','elements','nodes','geometry','quality_quantiles'):
            if key in row:assert actual[key]==row[key]
        worker=json.loads((directory/'worker.json').read_text())
        assert worker==row['worker']
        assert worker['numerical_thread_policy']['accelerate_blas_threading']==1
        if not (directory/'mesh.npz').exists():
            assert not row['qualified'];continue
        assert w.digest(directory/'mesh.npz')==row['mesh_sha256']
        mesh=np.load(directory/'mesh.npz')
        points,cells,surface=(mesh[k] for k in ('points','tetrahedra','surface'))
        assert len(points)==row['nodes']<=250000 and len(cells)==row['elements']<=w.CAP
        assert cells.shape[1]==(4 if row['config']['order']==1 else 10)
        q,_,_=scorecard.tet_metrics(points,cells[:,:4])
        measured=np.quantile(q,[0,.01,.05,.5,.95])
        np.testing.assert_allclose(measured,list(row['quality_quantiles'].values()),atol=1e-12)
        assert int((q<.1).sum())==row['elements_below_q01']
        assert audit(points,cells,surface)==row['geometry']['facet_consistency']
        try:
            _,_,periodic=periodic_pairs(points,surface)
        except Exception as error:
            assert 'periodic_error' in row and str(error)==row['periodic_error']
        else:assert periodic==row['periodic']
        if row['config']['arm']=='microgen' and row['config']['order']==2:
            edges=np.array(((0,1),(1,2),(2,0),(0,3),(1,3),(2,3)))
            np.testing.assert_allclose(points[cells[:,4:]],points[cells[:,edges]].mean(axis=2),atol=1e-15)
        assert row['qualified']==(w.geometry_is_admitted(row['geometry']) and 'periodic' in row)
        if replay_geometry:
            geometry=w.geometry(directory,report['case'],points,cells,surface,box_projection_fallback=True)
            for key in ('sampled_two_sided_max','relative_volume_error','q_min_corner','q_p01_corner','minimum_sampled_jacobian'):
                np.testing.assert_allclose(geometry[key],row['geometry'][key],atol=1e-12)
            assert geometry['geometry_admitted']==row['geometry']['geometry_admitted']
        artifacts+=1
        replays.append(dict(directory=str(directory),qualified=row['qualified'],mesh_sha256=row['mesh_sha256'],used_nodes=int(np.unique(cells).size),total_nodes=len(points)))
    return dict(run=str(run),input_hashes=len(hashes),artifact_replays=artifacts,full_geometry_replayed=replay_geometry,rows=replays)


def aggregate(run):
    report=json.loads((run/'summary.json').read_text())
    assert report['repetitions']==5 and report['warmups']==1
    groups={}
    for row in report['rows']:
        if row['warmup']:continue
        key=row['config']['name'];groups.setdefault(key,[]).append(row)
    rows=[]
    reference=json.loads((w.HERE/'runs/g3-microgen-defaults-v4'/(report['case']+'-specimen.json')).read_text())
    for name,values in groups.items():
        assert len(values)==5 and all(v['config']==values[0]['config'] for v in values)
        geometries=[v.get('geometry',{}) for v in values]
        row=dict(case=report['case'],name=name,config=values[0]['config'],run=str(run),measured_repetitions=5,
            qualified_repetitions=sum(v['qualified'] for v in values),statuses=[v['status'] for v in values],
            generation_s=distribution([v.get('generation_s') for v in values]),
            failed_attempt_s=distribution([v.get('failed_attempt_s') for v in values]),
            elements=distribution([v.get('elements') for v in values]),nodes=distribution([v.get('nodes') for v in values]),
            sampled_distance=distribution([g.get('sampled_two_sided_max') for g in geometries]),
            absolute_relative_volume_error=distribution([abs(g['relative_volume_error']) if 'relative_volume_error' in g else None for g in geometries]),
            quality={k:distribution([v.get('quality_quantiles',{}).get(k) for v in values]) for k in ('min','p01','p05','median','p95')},
            error_messages=sorted(set(v['error'] for v in values if 'error' in v)),
            periodic_errors=sorted(set(v['periodic_error'] for v in values if 'periodic_error' in v)),
            volume=distribution([g.get('volume') for g in geometries]),
            volume_reference=reference['volume'],volume_scramble_spread=reference['volume_scramble_spread'],
            volume_reference_kind='analytic' if report['case']=='schwarz-solid' else 'frozen four-scramble estimate',
            reference_volume_estimates=reference['volume_sobol_estimates'])
        failed=[]
        for v,g in zip(values,geometries):
            if v.get('periodic_error'):failed.append('periodicity')
            if v.get('error'):failed.append('generation')
            if g:
                if g['sampled_two_sided_max']>.002:failed.append('distance')
                if abs(g['relative_volume_error'])>.01:failed.append('volume')
                if g['q_min_corner']<.1:failed.append('quality')
                if g['projection_failed']:failed.append('projection')
                if g['minimum_sampled_jacobian']<=0:failed.append('Jacobian')
                if not g['oriented_facet_consistency_passed']:failed.append('facets')
        row['failed_gates']=sorted(set(failed))
        rows.append(row)
    ratios=[]
    for order in (1,2):
        actuals=[r for r in rows if r['config']['arm']=='microgen' and r['config']['order']==order and r['qualified_repetitions']==5]
        natives=[r for r in rows if r['config']['arm']=='meshers' and r['config']['order']==order and r['qualified_repetitions']==5]
        for a in actuals:
            for n in natives:
                quality={k:n['quality'][k]['median']>=a['quality'][k]['median'] for k in n['quality']}
                ratios.append(dict(order=order,actual=a['name'],native=n['name'],
                    generation_ratio_actual_over_native=a['generation_s']['median']/n['generation_s']['median'],
                    tet_ratio_actual_over_native=a['elements']['median']/n['elements']['median'],
                    node_ratio_actual_over_native=a['nodes']['median']/n['nodes']['median'],
                    distance_ratio_actual_over_native=a['sampled_distance']['median']/n['sampled_distance']['median'],
                    native_quality_not_worse=quality,
                    all_listed_quality_quantiles_not_worse=all(quality.values()),
                    scope='Admitted same-order comparison; not stiffness equivalence or certified surface/volume accuracy.'))
    return dict(case=report['case'],rows=rows,admitted_pair_ratios=ratios)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--runs',type=Path,nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--pilots',type=Path,nargs='*',default=[])
    args=parser.parse_args()
    assert not args.output.exists()
    started=time.perf_counter()
    reports=[aggregate(p) for p in args.runs]
    cases={r['case'] for r in reports};assert cases==set(w.CASES)
    pilots=[]
    for p in args.pilots:
        pilot=json.loads((p/'summary.json').read_text())
        pilots.append(dict(run=str(p),case=pilot['case'],summary_sha256=w.digest(p/'summary.json'),
            rows=[{k:r[k] for k in ('case','config','status','qualified','generation_s','failed_attempt_s','elements','nodes','quality_quantiles','error','periodic_error') if k in r}
                |dict(sampled_distance=r.get('geometry',{}).get('sampled_two_sided_max')) for r in pilot['rows']]))
    verification=[verify_run(p,replay_geometry=True) for p in args.runs] if args.verify else []
    inputs={str(Path(__file__)):w.digest(Path(__file__))}
    for p in args.runs:
        for name in ('summary.json','manifest.json','analysis-inputs.json'):
            inputs[str(p/name)]=w.digest(p/name)
    w.save(args.output,dict(cases=reports,pilots=pilots,verification=verification,inputs=inputs,seconds=time.perf_counter()-started,
        primary='Mesh generation, tetrahedra/nodes, same-order quality and geometry fidelity; no FE timings',
        failure_policy='No admitted speed ratio for a rejected competitor; all five statuses retained',
        volume_uncertainty='Scramble spread is empirical, not a certified confidence interval. Small volume differences are unresolved.',
        reference_promoted=False,new_equilibrium_solve=False))
    print(f'Compared {len(reports)} cases; verified {sum(r["artifact_replays"] for r in verification)} saved meshes',flush=True)


if __name__=='__main__':main()
