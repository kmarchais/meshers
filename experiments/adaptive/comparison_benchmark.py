"""Fixed eight-case same-host mesh comparison; actual Microgen seeds only."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time
import traceback

import microgen_workflow as w
import comparison_fixtures as fixtures


def configs(name):
    from mesh_only_comparison import configs as original
    if name in fixtures.CORE:
        old=original(name,True)
        linear=next(c for c in old if c['name']==('meshers-tet4-24' if name=='gyroid' else
            'meshers-tet4-adaptive' if name=='schwarz-solid' else 'meshers-tet4-historical32'))
        curved=next(c for c in old if c['name']=='meshers-curved-known')
        tuned=next(c for c in old if c['name']=='actual-explicit-tuned')
    else:
        spin=name=='spinodoid-anisotropic'
        linear=dict(arm='meshers',cells=32 if spin else 24,order=1,adapt=spin,
            tolerance=.002 if spin else 1.,min_length=.002,adapt_tolerance=.002)
        curved=dict(arm='meshers',cells=24 if spin else 16,order=2,adapt=spin,
            tolerance=1.,min_length=.002,adapt_tolerance=.003 if spin else .002)
        tuned=dict(arm='microgen',order=1,resolution=36 if spin else 30,
            first_only=True,first_hausd=.0015 if spin else .0005,first_hmax=.06 if spin else .1,
            mmg=dict(hausd=.0005,hmax=.06,hgrad=1.3),repair=not spin)
    return [dict(name='actual-default-tet4',arm='microgen',resolution=None,order=1),
        dict(tuned,name='actual-modified-tet4'),dict(linear,name='meshers-tet4'),dict(curved,name='meshers-curved-tet10')]


def worker(run,name,config,fixture_root):
    np,_,meshers,microgen,pv,_=w.imports()
    from unittest.mock import patch
    from microgen.remesh import remesh_keeping_boundaries_for_fem
    from vtkmodules.vtkCommonCore import vtkSMPTools
    from comparison_mesh_audit import observe
    assert vtkSMPTools.GetBackend()=='Sequential' and vtkSMPTools.GetEstimatedNumberOfThreads()==1
    target=fixtures.case(name,fixture_root)
    row=dict(case=name,config=config,status='started',qualified=False,mmg_calls=[],instrumentation_s=0.,
        vtk_backend='Sequential',vtk_threads=1,periodic_required=fixtures.periodic(name),solve_performed=False)
    w.save(run/'result.json',row)
    started=time.perf_counter()
    try:
        if config['arm']=='microgen':
            shape,part=fixtures.actual_shape(name,fixture_root,config['resolution'])
            options=dict(type_part=part) if isinstance(shape,microgen.Tpms) else (
                dict(resolution=config['resolution']) if name=='smooth-lattice' and config['resolution'] is not None else {})
            grid=shape.generate_volume_mesh(**options).triangulate()
            if name.startswith('spinodoid'):grid.translate((-.5,)*3,inplace=True)
            row.update(actual_resolution=getattr(shape,'resolution',None),seed_nodes=grid.n_points,seed_elements=grid.n_cells)
            if grid.n_cells>w.CAP or grid.n_points>w.CAP:raise ValueError('Microgen seed cap exceeded')
            if config.get('material_trace_remesh'):
                import inspect
                from periodic_trace import periodic_pairs
                seed_surface=grid.extract_surface(pass_pointid=True,algorithm='dataset_surface').triangulate()
                seed_faces=np.asarray(seed_surface.point_data['vtkOriginalPointIds'])[seed_surface.faces.reshape(-1,4)[:,1:]]
                input_tolerance=inspect.signature(remesh_keeping_boundaries_for_fem).parameters['tol'].default
                # Match the public input node-check scope. MMG may replace
                # seed triangulations; the final audit always requires full caps.
                _,_,row['seed_material_trace']=periodic_pairs(np.asarray(grid.points),seed_faces,
                    tolerance=input_tolerance,require_facets=False)
                expected_bounds=np.asarray(grid.points).copy()
                row['periodic_validator_override']='Explicit material-cap/full-existing-class input/output checks; public MMG passes unchanged'
            binary=Path(importlib.metadata.distribution('mmgpy').locate_file('mmgpy/bin/mmg3d_O3'))
            os.environ['PATH']=str(binary.parent)+os.pathsep+os.environ['PATH']
            original=microgen.Mmg.mmg3d
            def record(**parameters):
                observer=time.perf_counter()
                if config.get('first_only') and parameters.get('ls'):
                    shutil.copyfile(parameters['input'],parameters['output']);row['omitted_second_pass']=True
                    return
                if config.get('first_hausd') is not None and parameters.get('nofem'):
                    parameters.update(hausd=config['first_hausd'],hmax=config['first_hmax'],hgrad=1.3)
                call={k:v for k,v in parameters.items() if k not in ('input','output')}
                row['mmg_calls'].append(call);w.save(run/'result.json',row)
                row['instrumentation_s']+=time.perf_counter()-observer
                try:result=original(**parameters)
                except Exception as error:
                    captured=getattr(error.__cause__,'output',b'')
                    if captured:(run/f'mmg-call-{len(row["mmg_calls"])}.log').write_bytes(captured if isinstance(captured,bytes) else captured.encode())
                    raise
                observer=time.perf_counter()
                stage=pv.UnstructuredGrid(str(parameters['output']))
                call.update(output_tetrahedra=int(np.count_nonzero(stage.celltypes==pv.CellType.TETRA)),output_nodes=stage.n_points)
                w.save(run/'result.json',row);row['instrumentation_s']+=time.perf_counter()-observer
                if max(call['output_tetrahedra'],call['output_nodes'])>w.CAP:raise ValueError('Intermediate MMG cap exceeded')
                return result
            with patch.object(microgen.Mmg,'mmg3d',record):
                grid=remesh_keeping_boundaries_for_fem(grid,periodic=fixtures.periodic(name) and not config.get('material_trace_remesh',False),**config.get('mmg',{}))
            if config.get('material_trace_remesh'):
                import inspect
                from microgen.remesh import _snap_to_expected_bounds
                _snap_to_expected_bounds(grid,expected_bounds,inspect.signature(remesh_keeping_boundaries_for_fem).parameters['tol'].default)
            points=np.asarray(grid.points).copy();cells=np.asarray(grid.cells_dict[pv.CellType.TETRA]).copy()
            poly=grid.extract_surface(pass_pointid=True,algorithm='dataset_surface').triangulate()
            surface=np.asarray(poly.point_data['vtkOriginalPointIds'])[poly.faces.reshape(-1,4)[:,1:]]
            cells,row['reoriented_tets']=w.reorient_linear(points,cells)
            if config.get('repair'):
                from facet_incidence_audit import audit
                from interior_relocation_trial import repair
                if audit(points,cells,surface)['oriented_facet_consistency_passed']:
                    phase=time.perf_counter();points,row['repair']=repair(points,cells,surface);row['repair_s']=time.perf_counter()-phase
                else:row['repair_skipped']='Invalid oriented source, relocation forbidden'
            if config['order']==2:
                from microgen_tet10_control import straight_elevation
                phase=time.perf_counter();lift=straight_elevation(points,cells,surface)
                if not lift['within_node_cap']:raise ValueError(f'Straight Tet10 node cap exceeded: {lift["nodes"]}')
                points,cells,surface=(lift[k] for k in ('points','tetrahedra','surface'))
                row.update(straight_elevation_s=time.perf_counter()-phase,straight_elevation_explicit=True)
        else:
            field=target.function
            if config.get('native_density_lookup'):
                assert name=='graded-gyroid-density-linear'
                from graded_lookup_preflight import native_width
                prepare=time.perf_counter()
                width=native_width(json.loads((fixture_root/'linear-density-calibration.json').read_text()))
                def native_field(x,y,z):return fixtures.raw_gyroid(x,y,z)/width(x,y,z)
                field=meshers.compile_field(native_field,gradient=False)
                assert not field.has_gradient
                row.update(field_prepare_s=time.perf_counter()-prepare,native_graph_nodes=field.node_count,
                    native_analytic_gradient_enabled=False)
            if name=='smooth-lattice':
                from gradient_tie_probe import symmetric_lattice
                lattice=target.function if config.get('lattice_field')=='original' else symmetric_lattice
                def field(x,y,z):return np.minimum(lattice(x,y,z),.02)
                row.update(positive_void_clamp=.02,fixed_ordered_union_preserved=True)
            # Include compilation for the frozen spectral field in generation time.
            if name.startswith('spinodoid'):field=meshers.compile_field(field)
            mesh=meshers.generate(field,bounds=w.BOUNDS,cells=config['cells'],band=target.band,
                periodic=(fixtures.periodic(name),)*3,threads=1,optimize_passes=config.get('optimization_passes',1),
                minimum_quality=.1,max_tetrahedra=w.CAP,order=config['order'],geometry_tolerance=config['tolerance'],
                adapt=config['adapt'],min_length=config.get('min_length',.01),max_length=.1,gradation=.3,
                adapt_iterations=4,quality_guard=.3,adapt_tolerance=config.get('adapt_tolerance',.002),
                compile=config.get('compile',True))
            points,cells,surface=np.asarray(mesh.points),np.asarray(mesh.tetrahedra),np.asarray(mesh.surface)
            row['diagnostics']=mesh.diagnostics
            if config.get('fit_cap_edges'):
                assert name=='graded-gyroid-density-linear' and config['order']==2 and not fixtures.periodic(name)
                from curved_edge_fit import fit_cap_edges
                phase=time.perf_counter();points,fit_info=fit_cap_edges(target,points,cells,surface)
                row.update(curved_cap_fit=fit_info,curved_cap_fit_s=time.perf_counter()-phase,
                    native_only=False,diagnostics_scope='Native seed before the charged experimental Python fitting pass')
            if config.get('native_interp_auto'):
                assert row['diagnostics']['evaluator']=='compiled' and row['diagnostics']['callback_calls']==0
                row['public_auto_interp']=True
        row.update(generation_wall_s=time.perf_counter()-started,nodes=len(points),elements=len(cells))
        row['generation_s']=row['generation_wall_s']-row['instrumentation_s']
        if config['arm']=='meshers':
            probes=np.random.default_rng(20261026).uniform(-.5,.5,(4096,3))
            expected=target.function(*probes.T)
            identity_field=field
            if mesh.diagnostics['evaluator']=='compiled' and not hasattr(field,'values'):
                # Automatic preparation does not replace this caller's Python
                # variable. Inspect an owned compiled evaluator, not the same
                # Python callable twice. This check is outside generation time.
                identity_field=meshers.compile_field(field)
            actual=identity_field.values(probes) if hasattr(identity_field,'values') else identity_field(*probes.T)
            np.testing.assert_allclose(actual,np.minimum(expected,.02) if name=='smooth-lattice' else expected,rtol=1e-12,atol=1e-12)
            np.testing.assert_array_equal(actual<0,expected<0)
            row['generated_field_identity_checked']=4096
        if max(len(cells),len(points))>w.CAP:raise ValueError('Final node/tetrahedron cap exceeded')
        np.savez_compressed(run/'mesh.npz',points=points,tetrahedra=cells,surface=surface)
        row['mesh_sha256']=w.digest(run/'mesh.npz')
        w.save(run/'result.json',row)
        audit=time.perf_counter()
        periodic_auditor=None
        if config.get('periodic_trace'):
            from periodic_trace import periodic_pairs as periodic_auditor
            from periodic import periodic_pairs as historical_pairs
            try:
                historical_pairs(points,surface);row['historical_all_plane_periodic_gate']=dict(passed=True)
            except ValueError as error:row['historical_all_plane_periodic_gate']=dict(passed=False,error=str(error))
        row['geometry'],row['periodic'],row['quality_quantiles']=observe(run,name,target,points,cells,surface,fixtures.periodic(name),
            periodic_auditor=periodic_auditor)
        if periodic_auditor is not None and config['arm']=='meshers':
            checked,_,_=periodic_auditor(points,surface)
            for native,independent in zip(mesh.periodic_pairs,checked,strict=True):
                assert {tuple(map(int,p)) for p in native}=={tuple(map(int,p)) for p in independent}
            row['native_and_independent_pairs_agree']=True
        row.update(audit_s=time.perf_counter()-audit,qualified=row['geometry']['geometry_admitted'])
        row['status']='admitted' if row['qualified'] else 'admission-rejected'
    except Exception as error:
        row.update(status='failed',error=str(error),traceback=traceback.format_exc(),failed_attempt_s=time.perf_counter()-started)
    row['wall_s']=time.perf_counter()-started;w.save(run/'result.json',row)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--run',type=Path,required=True);p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--case',choices=fixtures.CASES,required=True)
    p.add_argument('--worker',type=json.loads);p.add_argument('--select')
    p.add_argument('--repetitions',type=int,default=1);p.add_argument('--warmups',type=int,default=0)
    p.add_argument('--actual-tet10',action='store_true',help='Also repeat explicit straight Tet10 elevation of selected actual routes')
    p.add_argument('--lattice-uniform37',action='store_true',help='Load the separately frozen matched-count lattice selection')
    p.add_argument('--spin-default-sampling',action='store_true',help='Control public sampling32 on the frozen FFT16 spinodoid field')
    p.add_argument('--periodic-trace',action='store_true',help='Explicit corrected material-cap and complete existing node-class audit')
    p.add_argument('--material-trace-remesh',action='store_true',help='Actual control: replace coordinate-only periodic validation, retain public MMG calls')
    p.add_argument('--density-adaptive','--grading-adaptive',action='store_true',help='One declared geometry-adaptive ablation on a linear grading case')
    p.add_argument('--density-lookup',nargs='?',const='finite',choices=('finite','automatic'),
        help='Pair NumPy with native finite-difference or public automatic interpolation')
    p.add_argument('--density-auto-interp',action='store_true',help='Control the newly supported public automatic interpolation path')
    p.add_argument('--density-cap-fit',action='store_true',help='Compare admitted adaptive Tet10 with one fresh uniform16 plus generic cap-fit prototype')
    a=p.parse_args()
    if a.worker is not None:return worker(a.run,a.case,a.worker,a.fixtures)
    a.run.mkdir(parents=True,exist_ok=False)
    np,_,meshers,*_=w.imports()
    import meshers._meshers,microgen.remesh,microgen.shape.tpms,microgen.shape.spinodoid
    from job_parts import worker_part
    sys.path.insert(0,str(Path.home()/'bin'))
    from pueue_progress import Job
    assert json.loads((a.fixtures/'result.json').read_text())['passed']
    selected=configs(a.case)
    if a.lattice_uniform37:
        assert a.case=='smooth-lattice'
        extra=json.loads((w.ROOT/'handoff/macbook/comparison-matched-count-selection-20261026.json').read_text())
        selected.extend(extra['configs'])
    if a.spin_default_sampling:
        extra=json.loads((w.ROOT/'handoff/macbook/comparison-spin-sampling-selection-20261026.json').read_text())
        assert a.case in extra['cases']
        selected.extend(extra['configs'])
    if a.actual_tet10:
        assert any(c['arm']=='microgen' for c in selected)
        selected.extend(dict(c,name=c['name'].replace('tet4','straight-tet10'),order=2)
            for c in list(selected) if c['arm']=='microgen')
    if a.select:
        wanted=a.select.split(',')
        if a.actual_tet10:
            wanted.extend(name.replace('tet4','straight-tet10') for name in list(wanted)
                if name.startswith('actual-') and name.endswith('-tet4'))
        assert set(wanted)<=set(c['name'] for c in selected)
        selected=[c for c in selected if c['name'] in wanted]
    if a.periodic_trace:selected=[dict(c,periodic_trace=True) for c in selected]
    if a.material_trace_remesh:
        assert a.periodic_trace and fixtures.periodic(a.case) and all(c['arm']=='microgen' for c in selected)
        selected=[dict(c,name=c['name']+'-material-trace',material_trace_remesh=True) for c in selected]
    if a.density_adaptive:
        assert a.case in ('graded-gyroid-density-linear','graded-gyroid-linear') and all(c['arm']=='meshers' for c in selected)
        selected=[dict(c,name=c['name']+'-adaptive001',adapt=True,adapt_tolerance=.001,min_length=.002) for c in selected]
    if a.density_lookup:
        assert a.case=='graded-gyroid-density-linear' and a.density_adaptive and all(c['arm']=='meshers' for c in selected)
        native=('-native-lookup',dict(compile=False,native_density_lookup=True)) if a.density_lookup=='finite' else (
            '-auto-interp',dict(compile=True,native_interp_auto=True))
        selected=[dict(c,name=c['name']+suffix,**extra) for c in selected
            for suffix,extra in (('-numpy',dict(compile=False)),native)]
    if a.density_auto_interp:
        assert a.case=='graded-gyroid-density-linear' and a.density_adaptive and not a.density_lookup
        assert all(c['arm']=='meshers' for c in selected)
        selected=[dict(c,name=c['name']+'-auto-interp',compile=True,native_interp_auto=True) for c in selected]
    if a.density_cap_fit:
        assert a.case=='graded-gyroid-density-linear' and a.density_adaptive and not a.density_lookup and not a.density_auto_interp
        assert len(selected)==1 and selected[0]['name']=='meshers-curved-tet10-adaptive001'
        baseline=dict(selected[0],name='meshers-curved-tet10-adaptive001-auto-interp',compile=True,native_interp_auto=True)
        selected=[baseline,dict(baseline,name='meshers-uniform16-auto-interp-cap-fit',adapt=False,fit_cap_edges=True)]
    paths=[Path(__file__),Path(w.__file__),Path(fixtures.__file__),w.HERE/'comparison_mesh_audit.py',w.HERE/'mesh_only_comparison.py',
        Path(meshers._meshers.__file__),Path(meshers.__file__),Path(meshers.expression.__file__),Path(meshers.intersection.__file__),
        Path(microgen.remesh.__file__),Path(microgen.shape.tpms.__file__),Path(microgen.shape.spinodoid.__file__),
        Path(importlib.metadata.distribution('mmgpy').locate_file('mmgpy/bin/mmg3d_O3')),
        w.HERE/'scorecard.py',w.HERE/'facet_incidence_audit.py',w.HERE/'box_projection.py',w.HERE/'gradient_tie_probe.py',
        w.HERE/'interior_relocation_trial.py',w.HERE/'interior_star_probe.py',w.HERE/'microgen_tet10_control.py',
        w.HERE/'quadratic_edge_removal_probe.py',w.HERE/'nested_refinement.py',
        w.ROOT.parent/'periodic-microgen/experiments/periodic-fea/periodic.py',
        fixtures.coefficients_path('spinodoid',a.fixtures),a.fixtures/'spinodoid-anisotropic-coefficients.json',a.fixtures/'result.json']
    paths.extend((a.fixtures/a.case).iterdir())
    if a.lattice_uniform37:paths.append(w.ROOT/'handoff/macbook/comparison-matched-count-selection-20261026.json')
    if a.spin_default_sampling:paths.append(w.ROOT/'handoff/macbook/comparison-spin-sampling-selection-20261026.json')
    if a.periodic_trace:paths.extend((w.HERE/'periodic_trace.py',w.ROOT/'handoff/macbook/schwarz-d-periodic-trace-decision-20261026.json'))
    if a.case=='graded-gyroid-density-linear':paths.extend((a.fixtures/'linear-density-calibration.json',
        a.fixtures/'density-controls.json',w.HERE/'linear_density_fixture.py',w.ROOT/'handoff/macbook/linear-density-selection-20261026.json'))
    if a.density_adaptive:paths.append(w.ROOT/'handoff/macbook'/('linear-density-adaptation-selection-20261026.json'
        if a.case=='graded-gyroid-density-linear' else 'linear-thickness-adaptation-selection-20261026.json'))
    if a.density_lookup:paths.extend((w.HERE/'graded_lookup_preflight.py',
        w.ROOT/'handoff/macbook/graded-lookup-experiment-selection-20261026.json'))
    if a.density_auto_interp or a.density_lookup=='automatic':
        paths.append(w.ROOT/'handoff/macbook/native-interp-auto-selection-20261026.json')
    if a.density_cap_fit:paths.extend((w.HERE/'curved_edge_fit.py',
        w.ROOT/'handoff/macbook/curved-density-live-fit-selection-20261026.json'))
    hashes={str(path):w.digest(path) for path in paths}
    w.save(a.run/'manifest.json',dict(inputs=hashes,configs=selected,platform=platform.platform(),architecture=platform.machine(),
        native_sha256=w.digest(Path(meshers._meshers.__file__)),numerical_thread_policy=w.NUMERICAL_THREAD_POLICY,
        packages={n:importlib.metadata.version(n) for n in ('meshers','microgen','mmgpy','numpy','scipy','pyvista','vtk')},
        worker_timeout_s=300,memory_bytes=w.MEMORY,node_cap=w.CAP,tetrahedron_cap=w.CAP,
        command=sys.argv,case=a.case,repetitions=a.repetitions,warmups=a.warmups,
        timing_scope='Fresh arm-specific field/shape preparation, compilation, extraction, MMG, explicit repair/elevation and common output arrays; imports, fixture loading, archive/audit and measured MMG observer I/O excluded',
        method='Corner quality q>=.1; abs-volume<=1%; finite sampled two-sided distance<=.002; positive125 curved Jacobians; corner facets and shared quadratic-edge conformity; periodic nodes and full caps where specimen periodic',
        distance_scope='Finite observations with feasible implicit feet, not certified Hausdorff or closest distances; cheap gate failures retain sampled reverse distance only',
        curved_reverse_subdivisions=[8,16],straight_observer_reduction='Exact affine Tet10 reduction verified before geometry measurement',
        mechanical_accuracy_claim=False,linear_grading='Bounded nonperiodic case for both arms; original four remain fully periodic'))
    for index,path in enumerate(paths):
        if path.suffix=='.py':shutil.copyfile(path,a.run/f'input-{index}-{path.name}')
    tasks=[(f'repeat-{i}-{c["name"]}',c) for i in range(-a.warmups,a.repetitions)
        for c in (selected if i%2==0 else list(reversed(selected)))]
    def deadline(*_):raise TimeoutError('Whole comparison deadline exceeded')
    signal.signal(signal.SIGALRM,deadline);signal.alarm(len(tasks)*310+60)
    rows=[]
    with Job(parts=[name for name,_ in tasks]) as job:
        for name,config in tasks:
            directory=a.run/name;directory.mkdir()
            for source in (a.fixtures/a.case).iterdir():shutil.copyfile(source,directory/source.name)
            with worker_part(job,name) as part:
                record=w.bounded([sys.executable,__file__,'--run',str(directory),'--fixtures',str(a.fixtures),
                    '--case',a.case,'--worker',json.dumps(config)],directory/'worker.log',part,timeout=300)
                w.save(directory/'worker.json',record)
                row=json.loads((directory/'result.json').read_text()) if (directory/'result.json').exists() else dict(status='no-result')
                row.update(worker=record,directory=str(directory),warmup=name.startswith('repeat--'));rows.append(row)
                job.result(name,seconds=record['wall_s'],generation_s=row.get('generation_s'),
                    qualified=row.get('qualified',False),elements=row.get('elements'),nodes=row.get('nodes'))
    w.save(a.run/'summary.json',dict(case=a.case,rows=rows,repetitions=a.repetitions,warmups=a.warmups,
        primary='Tet4 versus Tet4; curved Tet10 separate; modified actual routes and failures explicit',mechanical_accuracy_claim=False))
    signal.alarm(0)
    assert all(w.digest(Path(p))==h for p,h in hashes.items()),'Inputs changed during measurement'


if __name__=='__main__':main()
