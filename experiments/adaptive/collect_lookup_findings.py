"""Collect audited lookup and cap-fitting studies without rerunning mesh workers."""
import argparse
import json
from pathlib import Path

import numpy as np
import meshers

import microgen_workflow as w
import comparison_fixtures as fixtures
from graded_lookup_preflight import native_width


def resolve_inputs(record, audit, run):
    resolved=[]
    copies=list(audit.parent.glob('input-*.py'))+list(run.glob('input-*.py'))
    for name,expected in record['inputs'].items():
        path=Path(name)
        if not path.is_absolute():path=w.ROOT/path
        if not path.exists() or w.digest(path)!=expected:
            matches=[p for p in copies if w.digest(p)==expected]
            assert matches,f'Unresolved historical source version: {name}'
            path=matches[0]
        resolved.append(dict(measured_path=name,resolved_path=str(path),sha256=expected))
    return resolved


def study(audit_name):
    audit=w.HERE/'runs'/audit_name/'result.json';record=json.loads(audit.read_text())
    run=Path(record.get('run',record.get('checks',[{}])[0].get('run','')))
    if not run.is_absolute():run=w.ROOT/run
    summary=json.loads((run/'summary.json').read_text());manifest=json.loads((run/'manifest.json').read_text())
    assert summary['repetitions']==5 and summary['warmups']==1
    assert all(r['qualified'] and r['worker']['returncode']==0 and r['worker']['limit_failure'] is None for r in summary['rows'])
    return record,dict(audit_source=str(audit.relative_to(w.ROOT)),audit_sha256=w.digest(audit),
        run=str(run.relative_to(w.ROOT)),attempts=len(summary['rows']),configs=manifest['configs'],
        directories={c['name']:[r['directory'] for r in summary['rows'] if not r['warmup'] and r['config']==c]
                     for c in manifest['configs']},
        resolved_inputs=resolve_inputs(record,audit,run))


def collect():
    finite,finite_meta=study('macbook-density-lookup-repeated-audit-v1')
    automatic,auto_meta=study('macbook-native-interp-auto-audit-v1')
    fitting,fit_meta=study('macbook-curved-density-live-fit-repeated-audit-v1')
    assert finite['passed'] and automatic['passed']
    assert all(r['samples']==r['admitted']==5 for r in fitting['rows'])
    fit_rows={r['name']:r for r in fitting['rows']}
    def stats(row):
        return dict(samples=row['samples'],generation_median_s=row['generation_median_s'],
            generation_range_s=row['generation_range_s'],nodes=row['nodes'],elements=row['elements'],
            geometry=row['geometry'],quality=row['quality_quantiles'])
    baseline=stats(fit_rows['meshers-curved-tet10-adaptive001-auto-interp'])
    candidate=stats(fit_rows['meshers-uniform16-auto-interp-cap-fit'])
    fitting_pair=dict(order=2,baseline=baseline,prototype=candidate,
        speedup=baseline['generation_median_s']/candidate['generation_median_s'],
        tets_reduction_fraction=1-candidate['elements']/baseline['elements'],
        nodes_reduction_fraction=1-candidate['nodes']/baseline['nodes'],
        matched_stiffness=False,quality_geometry_tradeoffs_retained=True)
    fixtures_root=w.HERE/'runs/macbook-linear-density-fixtures-v1'
    target=fixtures.case('graded-gyroid-density-linear',fixtures_root)
    width=native_width(json.loads((fixtures_root/'linear-density-calibration.json').read_text()))
    def prototype(x,y,z):return fixtures.raw_gyroid(x,y,z)/width(x,y,z)
    a=meshers.compile_field(prototype,gradient=False);b=meshers.compile_field(target.function,gradient=False)
    assert a is b
    points=np.random.default_rng(20261026).uniform(-.5,.5,(8192,3))
    np.testing.assert_array_equal(a.values(points),b.values(points))
    old=w.ROOT/'handoff/macbook/comparison-findings-20261026.json'
    historical=json.loads(old.read_text())
    for name,expected in historical['lattice_stiffness']['inputs'].items():assert w.digest(Path(name))==expected
    regression=w.HERE/'runs/macbook-native-interp-originals-regression-v1/result.json'
    controls=json.loads(regression.read_text());assert controls['passed'] and len(controls['rows'])==7
    assert all(r['all_arrays_exactly_unchanged'] for r in controls['rows'])
    sources=[Path(__file__),old,regression,w.HERE/'runs/macbook-native-interp-python-full-v1/worker.log',
        w.HERE/'runs/macbook-native-interp-green-v1/worker.log',
        w.HERE/'runs/macbook-curved-density-edge-fit-v1/result.json',
        w.HERE/'runs/macbook-curved-density-native-uniform-v1/result.json']
    return dict(case='graded-gyroid-density-linear',host='Apple M4, one numerical thread',date='2026-10-26',
        native_binary_sha256=w.digest(Path(meshers._meshers.__file__)),
        expression_sha256=w.digest(Path(meshers.expression.__file__)),
        finite_difference=dict(finite_meta,pairs=finite['pairs']),
        automatic=dict(auto_meta,pairs=automatic['pairs'],compiled_field_identity=automatic['automatic_compiled_field_identity']),
        curved_fitting=dict(fit_meta,pair=fitting_pair,native_only=False,
            scope='Fresh native seed plus charged generic Python cap-edge detection/fitting. Both recipes admit but quality and geometry are not matched; no FEA comparison.'),
        compiled_graph_bridge=dict(same_cached_object=True,probes=8192,values_exactly_equal=True),
        tests=dict(full_python_and_tools_passed=208,skipped=1,focused_interp_and_Schwarz_D_passed=15,
                   seven_historical_curved_meshes_exact=True),
        retained_lattice_model_tensor_hashes=historical['lattice_stiffness']['retained_model_tensor_hashes'],
        inputs={str(p.relative_to(w.ROOT)):w.digest(p) for p in sources},
        failures_preserved=dict(audit665='Verifier syntax error; corrected replay666 passes',
            red669='Expected public API regression failures before implementation',
            uniform677='New automatic evaluator alone still rejects the fixed coarse recipe'),
        new_FE_solve=False,reference_promoted=False,actual_Microgen_speedup_claim=False,
        full_objective_complete=False,
        scope='Measured Meshers evaluator gains and an explicit hybrid geometry/cost tradeoff. Fixed specimens, gates, caps, old failures and all 24 lattice tensors remain unchanged.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    assert not args.output.exists();w.save(args.output,collect())
