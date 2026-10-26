"""Independently replay evaluator-only density comparisons and compare their meshes."""
import argparse
import io
import json
from pathlib import Path
import statistics

import numpy as np
from scipy.spatial import cKDTree

import microgen_workflow as w
import comparison_fixtures as fixtures
from comparison_mesh_audit import observe


def canonical(rows):
    rows = np.sort(rows, axis=1)
    return rows[np.lexsort(rows.T[::-1])]


def compare_meshes(a, b):
    result = dict(exact_arrays={k: np.array_equal(a[k], b[k]) for k in a.files})
    if any(a[k].shape != b[k].shape for k in a.files):
        return dict(result, same_counts=False, equivalent_up_to_numbering=False)
    distance, mapping = cKDTree(a['points']).query(b['points'])
    unique = len(np.unique(mapping)) == len(mapping)
    topology = unique and all(np.array_equal(canonical(a[k]), canonical(mapping[b[k]]))
                             for k in ('tetrahedra', 'surface'))
    return dict(result, same_counts=True, maximum_matched_point_distance=float(distance.max()),
                bijective_node_map=unique, same_cell_and_facet_node_sets=bool(topology),
                equivalent_up_to_numbering=bool(topology and distance.max() <= 1e-10),
                scope='Node-set equivalence after matching within 1e-10, with independent quadratic association checks on both meshes; not bitwise geometry or inherited FE-field certification')


def controls():
    def mesh(points, cells):
        stream=io.BytesIO()
        np.savez(stream,points=points,tetrahedra=cells,surface=np.array([[0,1,2]]))
        stream.seek(0)
        return np.load(stream)
    points=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    cells=np.array([[0,1,2,3]])
    a=mesh(points,cells)
    stream=io.BytesIO();permutation=np.array([2,0,3,1]);inverse=np.argsort(permutation)
    np.savez(stream,points=points[permutation],tetrahedra=inverse[cells[:,[1,0,3,2]]],surface=inverse[a['surface']])
    stream.seek(0)
    assert compare_meshes(a,np.load(stream))['equivalent_up_to_numbering']
    assert not compare_meshes(a,mesh(points+.01,cells))['equivalent_up_to_numbering']
    assert not compare_meshes(a,mesh(points,np.array([[0,1,2,2]])))['equivalent_up_to_numbering']


def verify(run):
    w.imports()
    manifest=json.loads((run/'manifest.json').read_text())
    summary=json.loads((run/'summary.json').read_text())
    assert summary['case']=='graded-gyroid-density-linear'
    assert manifest['node_cap']==manifest['tetrahedron_cap']==250000
    assert manifest['worker_timeout_s']==300 and manifest['memory_bytes']==8*1024**3
    inputs={str(run/'manifest.json'):w.digest(run/'manifest.json'),str(run/'summary.json'):w.digest(run/'summary.json'),
            str(Path(__file__)):w.digest(Path(__file__))}
    for name,expected in manifest['inputs'].items():
        path=Path(name)
        if not path.exists() or w.digest(path)!=expected:
            copies=[p for p in run.glob('input-*.py') if w.digest(p)==expected]
            assert copies,f'Unresolved measured input {name}'
            inputs[str(copies[0])]=expected
        else:inputs[name]=expected
    command=manifest['command'];fixture_root=Path(command[command.index('--fixtures')+1])
    compiled_identity=None
    if any(c.get('native_interp_auto') for c in manifest['configs']):
        import meshers
        expression=Path(meshers.expression.__file__)
        assert w.digest(expression)==manifest['inputs'][str(expression)]
        target=fixtures.case(summary['case'],fixture_root)
        compiled=meshers.compile_field(target.function)
        assert compiled.has_gradient
        knots=np.asarray(json.loads((fixture_root/'linear-density-calibration.json').read_text())['z'])
        boundary=np.r_[knots,np.nextafter(knots,-np.inf),np.nextafter(knots,np.inf)]
        points=np.random.default_rng(20261026).uniform(-.5,.5,(8192+len(boundary),3))
        points[-len(boundary):,2]=boundary
        expected=target.function(*points.T);actual=compiled.values(points)
        np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=1e-12)
        compiled_identity=dict(probes=len(points),maximum_field_difference=float(abs(actual-expected).max()),
            analytical_gradient_enabled=True,native_graph_nodes=compiled.node_count,
            expression_sha256=w.digest(expression),
            scope='Independent owned compiled-field check from unchanged measured inputs; the worker callable check alone does not inspect automatic compilation')
    groups={};cache={};audits=[]
    for row in summary['rows']:
        directory=Path(row['directory']);saved=json.loads((directory/'result.json').read_text())
        assert all(row[k]==v for k,v in saved.items())
        worker=json.loads((directory/'worker.json').read_text());assert worker==row['worker']
        assert worker['returncode']==0 and worker['limit_failure'] is None
        assert worker['peak_rss_bytes']<=worker['memory_bytes']
        assert worker['numerical_thread_policy']['accelerate_blas_threading']==1
        assert saved['config'] in manifest['configs'] and saved['vtk_threads']==1
        assert saved['qualified'] and not saved['solve_performed']
        native=bool(saved['config'].get('native_density_lookup') or saved['config'].get('native_interp_auto'))
        assert saved['diagnostics']['evaluator']==('compiled' if native else 'python')
        assert (saved['diagnostics']['callback_calls']==0)==native
        if saved['config'].get('native_density_lookup'):assert not saved['native_analytic_gradient_enabled']
        if saved['config'].get('native_interp_auto'):assert saved['public_auto_interp']
        for name in ('graded-gyroid-density-linear-samples.npz','graded-gyroid-density-linear-specimen.json'):
            assert w.digest(directory/name)==w.digest(fixture_root/summary['case']/name)
        mesh=directory/'mesh.npz';digest=w.digest(mesh);assert digest==saved['mesh_sha256']
        arrays=np.load(mesh);points,cells,surface=(arrays[k] for k in ('points','tetrahedra','surface'))
        assert len(points)==saved['nodes']<=250000 and len(cells)==saved['elements']<=250000
        assert cells.shape[1]==(4 if saved['config']['order']==1 else 10)
        if digest not in cache:
            cache[digest]=observe(directory,summary['case'],fixtures.case(summary['case'],fixture_root),
                                  points,cells,surface,False)
        geometry,periodic,quality=cache[digest]
        assert geometry['geometry_admitted'] and not geometry['failed_gates'] and periodic==saved['periodic']
        for key in ('volume','relative_volume_error','sampled_two_sided_max','minimum_sampled_jacobian'):
            np.testing.assert_allclose(geometry[key],saved['geometry'][key],atol=1e-12,rtol=1e-12)
        np.testing.assert_allclose(list(quality.values()),list(saved['quality_quantiles'].values()),atol=1e-13,rtol=1e-13)
        audits.append(dict(directory=str(directory),mesh_sha256=digest,admitted=True))
        if not row['warmup']:groups.setdefault((saved['config']['order'],native),[]).append(row)
    pairs=[]
    for order in (1,2):
        baseline,candidate=(groups[(order,n)] for n in (False,True))
        assert len(baseline)==len(candidate)==summary['repetitions']
        stats=[]
        for rows in (baseline,candidate):
            times=[r['generation_s'] for r in rows]
            assert len({(r['nodes'],r['elements']) for r in rows})==1
            stats.append(dict(samples=len(rows),generation_median_s=statistics.median(times),
                              generation_range_s=[min(times),max(times)],nodes=rows[0]['nodes'],elements=rows[0]['elements'],
                              callback_calls=[r['diagnostics']['callback_calls'] for r in rows],
                              geometry=rows[0]['geometry'],quality=rows[0]['quality_quantiles'],
                              adapt_phase_median_s=[statistics.median(r['diagnostics']['adapt']['seconds'][i] for r in rows) for i in range(7)]))
        comparisons=[]
        for a,b in zip(baseline,candidate,strict=True):
            with np.load(Path(a['directory'])/'mesh.npz') as ma,np.load(Path(b['directory'])/'mesh.npz') as mb:
                comparisons.append(compare_meshes(ma,mb))
        pairs.append(dict(order=order,baseline=stats[0],native_lookup=stats[1],mesh_comparisons=comparisons,
                          speedup=stats[0]['generation_median_s']/stats[1]['generation_median_s'] if len(baseline)>=5 else None,
                          pilot_ratio=stats[0]['generation_median_s']/stats[1]['generation_median_s'] if len(baseline)<5 else None,
                          accuracy_matched_FEA_claim=False))
    return dict(passed=True,run=str(run),pairs=pairs,unique_mesh_replays=len(cache),audits=audits,inputs=inputs,
                automatic_compiled_field_identity=compiled_identity,
                warmups=summary['warmups'],mechanical_accuracy_claim=False,
                actual_Microgen_comparison=False,full_objective_complete=False,
                scope='Fresh same-host evaluator comparison at fixed recipes and independent mesh admission. No actual competitor ratio or new stiffness result.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    assert not args.output.exists()
    controls()
    w.save(args.output,verify(args.run))
