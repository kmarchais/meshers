"""Replay post-fix meshes and retain all failures without changing old evidence."""
import argparse
import json
from pathlib import Path
import statistics
import numpy as np
import microgen_workflow as w
import comparison_fixtures as fixtures
from comparison_mesh_audit import observe
from periodic_trace import periodic_pairs


def verify(runs):
    import meshers._meshers
    native=w.digest(Path(meshers._meshers.__file__))
    assert native!='2005c3e01892098a9d734567220a56ed3ea5d748b20a4f96754fbb55168060b1'
    cache={};checks=[];groups={};pilots={};inputs={}
    for run in runs:
        manifest=json.loads((run/'manifest.json').read_text());summary=json.loads((run/'summary.json').read_text())
        assert manifest['native_sha256']==native and manifest['node_cap']==manifest['tetrahedron_cap']==250000
        assert manifest['worker_timeout_s']==300 and manifest['memory_bytes']==8*1024**3
        for name,expected in manifest['inputs'].items():
            p=Path(name)
            if not p.exists() or w.digest(p)!=expected:
                copies=[p for p in run.glob('input-*.py') if w.digest(p)==expected]
                assert copies,f'Unresolved measured source: {name}'
                inputs[str(copies[0])]=expected
            inputs[name]=expected
        command=manifest['command'];fixture_root=Path(command[command.index('--fixtures')+1])
        rows=[]
        for row in summary['rows']:
            directory=Path(row['directory']);saved=json.loads((directory/'result.json').read_text())
            for key,value in saved.items():assert row[key]==value,(directory,key)
            worker=json.loads((directory/'worker.json').read_text());assert worker==row['worker']
            assert worker['returncode']==0 and worker['limit_failure'] is None
            assert worker['peak_rss_bytes']<=worker['memory_bytes']
            assert saved['config'] in manifest['configs'] and saved['vtk_threads']==1
            if not row['warmup']:
                pool=groups if summary['repetitions']>=5 else pilots
                pool.setdefault((summary['case'],saved['config']['name']),[]).append(row)
            path=directory/'mesh.npz'
            if not path.exists():assert not saved['qualified'];continue
            assert w.digest(path)==saved['mesh_sha256']
            mesh=np.load(path);points,cells,surface=(mesh[k] for k in ('points','tetrahedra','surface'))
            assert len(points)==saved['nodes']<=250000 and len(cells)==saved['elements']<=250000
            assert cells.shape[1]==(4 if saved['config']['order']==1 else 10)
            for name in (summary['case']+'-samples.npz',summary['case']+'-specimen.json'):
                assert w.digest(directory/name)==w.digest(fixture_root/summary['case']/name)
            key=(summary['case'],saved['mesh_sha256'],bool(saved['config'].get('periodic_trace')))
            if key not in cache:
                cache[key]=observe(directory,summary['case'],fixtures.case(summary['case'],fixture_root),
                    points,cells,surface,fixtures.periodic(summary['case']),
                    periodic_auditor=periodic_pairs if saved['config'].get('periodic_trace') else None)
            geometry,periodic,quality=cache[key]
            for name in ('volume','relative_volume_error','q_min_corner','minimum_sampled_jacobian',
                         'sampled_two_sided_max','sampled_reverse_max'):
                if name in geometry:np.testing.assert_allclose(geometry[name],saved['geometry'][name],atol=1e-11,rtol=1e-11)
            assert geometry['failed_gates']==saved['geometry']['failed_gates']
            assert geometry['geometry_admitted']==saved['qualified']
            np.testing.assert_allclose(list(quality.values()),list(saved['quality_quantiles'].values()),atol=1e-13,rtol=1e-13)
            if saved['qualified'] and fixtures.periodic(summary['case']):
                assert periodic['complete_existing_node_classes'] and periodic['periodic_triangles_match']
                assert saved['native_and_independent_pairs_agree']
            rows.append(dict(directory=str(directory),mesh_sha256=saved['mesh_sha256'],admitted=saved['qualified'],
                failed_gates=geometry['failed_gates'],historical_gate=saved.get('historical_all_plane_periodic_gate')))
        checks.append(dict(run=str(run),summary_sha256=w.digest(run/'summary.json'),rows=rows))
        inputs[str(run/'summary.json')]=w.digest(run/'summary.json')
        inputs[str(run/'manifest.json')]=w.digest(run/'manifest.json')
    for key,rows in pilots.items():
        if key not in groups:groups[key]=rows
    statistics_rows=[]
    for (case,name),rows in groups.items():
        geometry=rows[0].get('geometry',{})
        times=[r['generation_s'] for r in rows if 'generation_s' in r]
        statistics_rows.append(dict(case=case,name=name,config=rows[0]['config'],samples=len(rows),
            admitted=sum(r['qualified'] for r in rows),
            generation_median_s=statistics.median(times) if times else None,
            generation_range_s=[min(times),max(times)] if times else None,
            nodes=rows[0].get('nodes'),elements=rows[0].get('elements'),geometry=geometry,
            quality_quantiles=rows[0].get('quality_quantiles'),
            errors=sorted({r['error'] for r in rows if 'error' in r}),
            directories=[r['directory'] for r in rows]))
    old=w.ROOT/'handoff/macbook/comparison-findings-20261026.json'
    historical=json.loads(old.read_text())
    for name,expected in historical['lattice_stiffness']['inputs'].items():assert w.digest(Path(name))==expected
    original=w.HERE/'runs/macbook-schwarz-d-periodic-default-regression-v1/checks/result.json'
    regression=json.loads(original.read_text())
    assert regression['box']['candidate']['binary_sha256']==native
    assert len(regression['cases'])==4 and all(r['arrays_exactly_unchanged'] for r in regression['cases'].values())
    inputs[str(old)]=w.digest(old);inputs[str(original)]=w.digest(original)
    # Replay the frozen study under its original gate, resolving the preserved
    # native binary by hash instead of rewriting any old manifest or outcome.
    from comparison_verified_report import verify as verify_historical
    historical_cache={}
    historical_checks=[verify_historical(Path(row['run']),w.HERE/'runs/macbook-comparison-fixtures-v1',historical_cache)
        for row in historical['sources']]
    return dict(native_sha256=native,baseline_native_sha256='2005c3e01892098a9d734567220a56ed3ea5d748b20a4f96754fbb55168060b1',
        checks=checks,unique_mesh_replays=len(cache),rows=statistics_rows,inputs=inputs,
        historical_runs_replayed=len(historical_checks),historical_unique_mesh_replays=len(historical_cache),
        historical_checks=historical_checks,
        original_four_mesh_arrays_unchanged=True,retained_lattice_models=24,new_FE_solve=False,
        reference_promoted=False,full_objective_complete=False,
        scope='Correct material-trace matching and linear grading coverage. Earlier all-plane rule and failures retained. No success ratio for an inadmissible competitor.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists();w.imports();w.save(a.output,verify(a.runs))
