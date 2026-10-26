"""Verify saved comparison artifacts and generate a numeric comparison ledger."""
import argparse
import json
from pathlib import Path
import statistics
import time
import numpy as np
import microgen_workflow as w
import comparison_fixtures as fixtures
from comparison_mesh_audit import observe


def distribution(values):
    values=[v for v in values if v is not None]
    return None if not values else dict(median=statistics.median(values),minimum=min(values),maximum=max(values))


def lattice_stiffness_evidence():
    """Re-evaluate saved six-load tensors against every retained model, no solve."""
    import hashlib
    from cubic_refinement import MODEL_REPORTS
    from periodic import tensor_difference
    models=[];hashes=[];inputs={}
    for name in MODEL_REPORTS:
        path=w.HERE/'runs'/name;report=json.loads(path.read_text())
        assert report['case']=='smooth-lattice' and len(report['models'])==8 and report['stiffness_target']==.01
        inputs[str(path)]=w.digest(path)
        for model in report['models']:
            matrix=np.asarray(model['stiffness'],dtype=np.float64)
            assert matrix.shape==(6,6) and np.linalg.eigvalsh(matrix).min()>0
            models.append(matrix);hashes.append(hashlib.sha256(matrix.tobytes()).hexdigest())
    previous=json.loads((w.HERE/'runs/macbook-lattice-uniform37-assessment-v1/result.json').read_text())
    assert len(models)==len(set(hashes))==24 and previous['retained_model_tensor_hashes']==hashes
    routes={'compact-native':'macbook-cap-phase-lattice26-amg-v1','uniform37-native':'macbook-lattice-uniform37-amg-v1',
        'experimental-saved-relocation':'macbook-lattice-relocation-solve-v1','modified-actual-with-star-refinement':'macbook-lattice-actual-live-solve-v1'}
    expected={'compact-native':previous['previous_worst_model_disagreement'],'uniform37-native':previous['worst_model_disagreement']}
    for key,source in [('experimental-saved-relocation','macbook-lattice-relocation-assessment-v1'),
                       ('modified-actual-with-star-refinement','macbook-lattice-actual-live-assessment-v1')]:
        path=w.HERE/'runs'/source/'result.json';assessment=json.loads(path.read_text());inputs[str(path)]=w.digest(path)
        assert assessment['retained_model_count']==24
        if key=='experimental-saved-relocation':assert assessment['retained_model_tensor_hashes']==hashes
        expected[key]=assessment.get('worst_model_disagreement',assessment.get('worst_retained_model_disagreement'))
    rows=[]
    for key,source in routes.items():
        path=w.HERE/'runs'/source/'result.json';row=json.loads(path.read_text());inputs[str(path)]=w.digest(path)
        assert row['completed_loads']==6
        tensor=np.asarray(row['stiffness']);assert tensor.shape==(6,6) and np.linalg.eigvalsh(tensor).min()>0
        disagreements=[tensor_difference(tensor,model)['worst_relative_strain_energy_difference'] for model in models]
        error=max(disagreements)
        np.testing.assert_allclose(error,expected[key],atol=1e-12,rtol=1e-12)
        rows.append(dict(route=key,source=str(path),nodes=row['nodes'],elements=row['elements'],
            worst_model_disagreement=error,all24_model_disagreements=disagreements,passes_one_percent=error<=.01))
    return dict(retained_model_count=24,retained_model_tensor_hashes=hashes,rows=rows,inputs=inputs,
        scope='Re-evaluated saved full six-load stiffness tensors and unchanged models; no fresh FE solve or comparative workflow timing. The actual star-refined mesh differs from the mesh-only route.')


def verify(run,fixture_root,cache):
    manifest=json.loads((run/'manifest.json').read_text());summary=json.loads((run/'summary.json').read_text())
    protocol=json.loads((w.ROOT/'handoff/macbook/comparison-protocol-20261026.json').read_text())
    declared={c['name']:c for c in protocol['cases'][summary['case']]}
    if any('uniform37' in c['name'] for c in manifest['configs']):
        extra=json.loads((w.ROOT/'handoff/macbook/comparison-matched-count-selection-20261026.json').read_text())
        assert summary['case']==extra['case']
        for name,expected in extra['inputs'].items():assert w.digest(w.ROOT/name)==expected
        declared.update({c['name']:c for c in extra['configs']})
    if any('public-sampling32' in c['name'] for c in manifest['configs']):
        extra=json.loads((w.ROOT/'handoff/macbook/comparison-spin-sampling-selection-20261026.json').read_text())
        assert summary['case'] in extra['cases']
        for name,expected in extra['inputs'].items():assert w.digest(w.ROOT/name)==expected
        declared.update({c['name']:c for c in extra['configs']})
    for config in manifest['configs']:
        base=config['name'].replace('straight-tet10','tet4')
        expected=declared[base]
        if 'straight-tet10' in config['name']:expected=dict(expected,name=config['name'],order=2)
        assert config==expected,'Configuration differs from frozen matrix'
    assert manifest['native_sha256']=='2005c3e01892098a9d734567220a56ed3ea5d748b20a4f96754fbb55168060b1'
    assert manifest['node_cap']==manifest['tetrahedron_cap']==250000 and not manifest['mechanical_accuracy_claim']
    assert manifest['worker_timeout_s']==300 and manifest['memory_bytes']==8*1024**3
    for name,expected in manifest['inputs'].items():
        path=Path(name)
        if not path.exists() or w.digest(path)!=expected:
            copies=[p for p in run.glob('input-*.py') if w.digest(p)==expected]
            assert copies,f'Unresolved measured input {name}'
        else:assert w.digest(path)==expected
    replays=[]
    for row in summary['rows']:
        directory=Path(row['directory']);actual=json.loads((directory/'result.json').read_text())
        assert actual['config'] in manifest['configs']
        for source in (fixture_root/summary['case']).iterdir():
            assert w.digest(directory/source.name)==w.digest(source)
        assert row['worker']==json.loads((directory/'worker.json').read_text())
        for key,value in actual.items():assert row[key]==value,(directory,key)
        assert row['worker']['numerical_thread_policy']['accelerate_blas_threading']==1
        assert actual['vtk_backend']=='Sequential' and actual['vtk_threads']==1 and not actual['solve_performed']
        mesh=directory/'mesh.npz'
        if not mesh.exists():assert not actual['qualified'];continue
        digest=w.digest(mesh);assert digest==actual['mesh_sha256']
        arrays=np.load(mesh);points,cells,surface=(arrays[k] for k in ('points','tetrahedra','surface'))
        assert len(points)==actual['nodes']<=250000 and len(cells)==actual['elements']<=250000
        assert cells.shape[1]==(4 if actual['config']['order']==1 else 10)
        key=(summary['case'],digest)
        if key not in cache:
            cache[key]=observe(directory,summary['case'],fixtures.case(summary['case'],fixture_root),points,cells,surface,fixtures.periodic(summary['case']))
        geometry,periodic,quality=cache[key]
        assert periodic==actual['periodic']
        assert quality.keys()==actual['quality_quantiles'].keys()
        np.testing.assert_allclose(list(quality.values()),list(actual['quality_quantiles'].values()),atol=1e-13,rtol=1e-13)
        for name in ('volume','relative_volume_error','q_min_corner','q_p01_corner','minimum_sampled_jacobian',
                     'sampled_two_sided_max','sampled_reverse_max','projection_failed'):
            if name in geometry and geometry[name] is not None:
                np.testing.assert_allclose(geometry[name],actual['geometry'][name],atol=1e-11,rtol=1e-11)
        assert geometry['failed_gates']==actual['geometry']['failed_gates']
        assert geometry['quadratic_conformity']==actual['geometry']['quadratic_conformity']
        assert geometry['geometry_admitted']==actual['qualified']
        if actual['qualified']:
            assert row['worker']['returncode']==0 and row['worker']['limit_failure'] is None
            for call in actual['mmg_calls']:
                assert max(call['output_nodes'],call['output_tetrahedra'])<=w.CAP
        replay=dict(mesh=str(mesh),sha256=digest,qualified=actual['qualified'],full_distance=geometry['distance_audit_complete'])
        if actual['config']['arm']=='meshers' and actual['config']['order']==2 and summary['case'] in fixtures.CORE:
            legacy=json.loads((w.ROOT/'handoff/macbook/adaptive-interior-validation-20261025.json').read_text())
            label={'schwarz-solid':'schwarz','smooth-lattice':'lattice'}.get(summary['case'],summary['case'])
            old=w.ROOT/legacy['cases'][label]['run']/'mesh.npz'
            assert w.digest(old)==legacy['inputs'][str(old.relative_to(w.ROOT))]
            with np.load(old) as saved:
                replay['prior_field_mesh_arrays_exact']=all(np.array_equal(arrays[k],saved[k]) for k in ('points','tetrahedra','surface'))
                replay['prior_field_mesh_source']=str(old)
        replays.append(replay)
    return dict(run=str(run),inputs=len(manifest['inputs']),replays=replays)


def aggregate(runs):
    groups={};pilots={};sources=[];records={}
    for run in runs:
        summary=json.loads((run/'summary.json').read_text());sources.append(dict(run=str(run),sha256=w.digest(run/'summary.json')))
        for row in summary['rows']:
            if not row['warmup']:
                records[row['directory']]=row
                pool=groups if summary['repetitions']>=5 else pilots
                pool.setdefault((summary['case'],row['config']['name']),[]).append(row)
    for key,rows in pilots.items():
        if key not in groups:groups[key]=rows
    result=[]
    for (case,name),rows in groups.items():
        assert all(r['config']==rows[0]['config'] for r in rows)
        good=[r for r in rows if r['qualified'] and r['worker']['returncode']==0 and r['worker']['limit_failure'] is None]
        geometries=[r.get('geometry',{}) for r in rows]
        result.append(dict(case=case,name=name,config=rows[0]['config'],samples=len(rows),qualified=len(good),
            generation_s=distribution([r.get('generation_s') for r in rows]),failed_attempt_s=distribution([r.get('failed_attempt_s') for r in rows]),
            nodes=distribution([r.get('nodes') for r in rows]),elements=distribution([r.get('elements') for r in rows]),
            sampled_distance=distribution([g.get('sampled_two_sided_max') for g in geometries]),
            sampled_reverse_only=distribution([g.get('sampled_reverse_max') for g in geometries if not g.get('distance_audit_complete',False)]),
            abs_relative_volume=distribution([abs(g['relative_volume_error']) for g in geometries if 'relative_volume_error' in g]),
            quality={key:distribution([r.get('quality_quantiles',{}).get(key) for r in rows]) for key in ('min','p01','p10','median','p90','p99','max')},
            failed_gates=sorted(set(k for g in geometries for k in g.get('failed_gates',[]))),
            errors=sorted(set([r['error'] for r in rows if 'error' in r]+[
                'Worker '+r['worker']['limit_failure'] for r in rows if r['worker']['limit_failure']]+[
                f'Worker returncode {r["worker"]["returncode"]}' for r in rows if r['worker']['returncode']!=0])),
            source_rows=[r['directory'] for r in rows]))
    ratios=[]
    for case in fixtures.CASES:
        for order in (1,2):
            eligible=[r for r in result if r['case']==case and r['config']['order']==order and r['qualified']==r['samples'] and r['samples']>=5]
            for actual in (r for r in eligible if r['config']['arm']=='microgen'):
                for native in (r for r in eligible if r['config']['arm']=='meshers'):
                    a,n=actual,native
                    common=set(str(Path(p).parent) for p in a['source_rows'])&set(str(Path(p).parent) for p in n['source_rows'])
                    for run in sorted(common):
                        av=[records[p]['generation_s'] for p in a['source_rows'] if str(Path(p).parent)==run]
                        nv=[records[p]['generation_s'] for p in n['source_rows'] if str(Path(p).parent)==run]
                        assert len(av)==len(nv)==5,'A comparison requires five interleaved successful measurements'
                        ratios.append(dict(case=case,order=order,actual=a['name'],native=n['name'],run=run,paired_samples=5,
                        speed_actual_over_native=statistics.median(av)/statistics.median(nv),
                        tets_actual_over_native=a['elements']['median']/n['elements']['median'],
                        nodes_actual_over_native=a['nodes']['median']/n['nodes']['median'],
                        distance_actual_over_native=a['sampled_distance']['median']/n['sampled_distance']['median'],
                        quality_native_minus_actual={k:n['quality'][k]['median']-a['quality'][k]['median'] for k in n['quality']},
                        stiffness_matched=False))
    return dict(rows=result,admitted_same_order_ratios=ratios,sources=sources,
        pilots_retained_outside_repeated_timing_pool=[r['directory'] for key,rows in pilots.items() if key in groups and groups[key] is not rows for r in rows])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true');a=p.parse_args()
    assert not a.output.exists()
    w.imports();started=time.perf_counter();cache={}
    checks=[verify(run,a.fixtures,cache) for run in a.runs] if a.verify else []
    w.save(a.output,dict(aggregate(a.runs),verification=checks,full_mesh_replays=len(cache),seconds=time.perf_counter()-started,
        lattice_stiffness=lattice_stiffness_evidence(),
        protocol_sha256=w.digest(w.ROOT/'handoff/macbook/comparison-protocol-20261026.json'),
        scope='Finite mesh-generation comparison only; no new stiffness/local-stress/continuum certification. Failed actual routes do not support admitted speed ratios.'))
