"""Fresh mesh-only measurements of actual Microgen/MMG and native Meshers."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import signal
import platform
import subprocess
import sys
import time
import traceback

import microgen_workflow as w


def configs(case, strengthened=False):
    points = [dict(name='actual-default',arm='microgen',resolution=None,order=1)]
    if case in ('gyroid','schwarz-solid'):
        points.append(dict(name='actual-resolution30',arm='microgen',resolution=30,order=1))
    for cells in (16,24):
        points.append(dict(name=f'meshers-tet4-{cells}',arm='meshers',cells=cells,order=1,adapt=False,tolerance=1.))
    cells = {'gyroid':9,'schwarz-solid':16,'spinodoid':24,'smooth-lattice':26}[case]
    points.append(dict(name='meshers-curved-known',arm='meshers',cells=cells,order=2,
        adapt=case in ('spinodoid','smooth-lattice'),tolerance=1.,min_length=.002,
        adapt_tolerance={'spinodoid':.003,'smooth-lattice':.0028}.get(case,.002)))
    if strengthened:
        strong=dict(name='actual-explicit-tuned',arm='microgen',order=1,
            resolution={'gyroid':24,'schwarz-solid':30,'spinodoid':36,'smooth-lattice':51}[case],
            mmg=dict(hausd=.001 if case in ('schwarz-solid','smooth-lattice') else .0005,hmax=.06,hgrad=1.3),
            first_hausd=.0005 if case in ('gyroid','schwarz-solid') else .0015,
            first_hmax=.1 if case in ('gyroid','schwarz-solid') else .06,
            first_only=case in ('gyroid','schwarz-solid','spinodoid'),repair=case in ('schwarz-solid','smooth-lattice'))
        cells={'gyroid':(20,24),'schwarz-solid':(16,24),'spinodoid':(24,32),'smooth-lattice':(24,32)}[case]
        linear=[dict(name=f'meshers-tet4-{n}',arm='meshers',cells=n,order=1,adapt=False,tolerance=1.) for n in cells]
        quadratic=dict(strong,name='actual-explicit-tuned-tet10',order=2)
        adaptive=dict(name='meshers-tet4-adaptive',arm='meshers',cells=24 if case=='spinodoid' else 16,
            order=1,adapt=True,tolerance=1.,min_length=.002,adapt_tolerance=.001)
        quality=dict(points[-1],name='meshers-curved-opt4',optimization_passes=4)
        default_quad=dict(points[0],name='actual-default-tet10',order=2)
        points=[points[0],default_quad,strong,quadratic,*linear,adaptive,points[-1],quality]
        if case=='smooth-lattice':
            uniform=dict(name='meshers-curved-uniform25',arm='meshers',cells=25,order=2,adapt=False,tolerance=1.,min_length=.01,adapt_tolerance=.002)
            points.extend((uniform,dict(uniform,name='meshers-curved-uniform25-opt4',optimization_passes=4)))
        if case in ('spinodoid','smooth-lattice'):
            points.append(dict(name='meshers-tet4-historical32',arm='meshers',cells=32,
                order=1,adapt=True,tolerance=.002,min_length=.002,adapt_tolerance=.002,
                lattice_field='original' if case=='smooth-lattice' else None))
    return points


def spin_field(shape,meshers,np):
    f=shape._frep
    coefficients=[(*map(float,mode/f.cell_size),float(amplitude),float(phase))
        for mode,amplitude,phase in zip(f.modes,f.amplitudes,f.phases,strict=True)]
    def field(x,y,z):
        value=float(f.threshold)
        for a,b,c,amplitude,phase in coefficients:
            value=value-amplitude*np.cos(2*np.pi*(a*(x+.5)+b*(y+.5)+c*(z+.5))+phase)
        return value
    return meshers.compile_field(field)


def worker(run,case_name,config):
    np,_,meshers,microgen,pv,scorecard=w.imports()
    from periodic import periodic_pairs
    from unittest.mock import patch
    from microgen.remesh import remesh_keeping_boundaries_for_fem
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.GetBackend()=='Sequential' and vtkSMPTools.GetEstimatedNumberOfThreads()==1
    row=dict(vtk_backend='Sequential',vtk_threads=1,case=case_name,config=config,status='started',qualified=False,mmg_calls=[],solve_performed=False,instrumentation_s=0.)
    w.save(run/'result.json',row)
    started=time.perf_counter()
    try:
        case,shape,part=w.specimen(case_name)
        if config['arm']=='microgen':
            if case_name=='spinodoid' and config['resolution'] is not None:
                shape.resolution=config['resolution']
                row['fft_realization_resolution']=16
                row['sampling_resolution']=config['resolution']
            elif case_name in ('gyroid','schwarz-solid') and config['resolution'] is not None:
                from microgen.shape import surface_functions
                shape=microgen.Tpms(surface_functions.gyroid if case_name=='gyroid' else surface_functions.schwarz_p,
                    offset=1.1 if case_name=='gyroid' else 0.,resolution=config['resolution'])
            row['actual_resolution']=getattr(shape,'resolution',None)
            grid=shape.generate_volume_mesh(**({'type_part':part} if case_name in ('gyroid','schwarz-solid') else ({'resolution':config['resolution']} if case_name=='smooth-lattice' and config['resolution'] is not None else {}))).triangulate()
            if case_name=='spinodoid':grid.translate((-.5,)*3,inplace=True)
            row.update(seed_nodes=grid.n_points,seed_elements=grid.n_cells)
            if grid.n_cells>w.CAP:raise ValueError('Seed tetrahedron cap exceeded')
            binary=Path(importlib.metadata.distribution('mmgpy').locate_file('mmgpy/bin/mmg3d_O3'))
            os.environ['PATH']=str(binary.parent)+os.pathsep+os.environ['PATH']
            original=microgen.Mmg.mmg3d
            def record(**parameters):
                observer=time.perf_counter()
                if config.get('first_only') and parameters.get('ls'):
                    shutil.copyfile(parameters['input'],parameters['output'])
                    row['omitted_second_pass']=True
                    return
                if config.get('first_hausd') is not None and parameters.get('nofem'):
                    parameters.update(hausd=config['first_hausd'],hmax=config['first_hmax'],hgrad=1.3)
                call={k:v for k,v in parameters.items() if k not in ('input','output')}
                row['mmg_calls'].append(call);w.save(run/'result.json',row)
                row['instrumentation_s']+=time.perf_counter()-observer
                try:
                    result=original(**parameters)
                except Exception as error:
                    captured=getattr(error.__cause__,'output',b'')
                    if captured:(run/f'mmg-call-{len(row["mmg_calls"])}.log').write_bytes(captured if isinstance(captured,bytes) else captured.encode())
                    raise
                observer=time.perf_counter()
                stage=pv.UnstructuredGrid(str(parameters['output']))
                count=int(np.count_nonzero(stage.celltypes==pv.CellType.TETRA))
                call['output_tetrahedra']=count
                w.save(run/'result.json',row)
                row['instrumentation_s']+=time.perf_counter()-observer
                if count>w.CAP:raise ValueError('Intermediate MMG tetrahedron cap exceeded')
                return result
            with patch.object(microgen.Mmg,'mmg3d',record):
                grid=remesh_keeping_boundaries_for_fem(grid,periodic=True,**config.get('mmg',{}))
            points=np.asarray(grid.points).copy()
            tets=np.asarray(grid.cells_dict[pv.CellType.TETRA]).copy()
            poly=grid.extract_surface(pass_pointid=True,algorithm='dataset_surface').triangulate()
            surface=np.asarray(poly.point_data['vtkOriginalPointIds'])[poly.faces.reshape(-1,4)[:,1:]]
            tets,row['reoriented_tets']=w.reorient_linear(points,tets)
            if config.get('repair'):
                from facet_incidence_audit import audit
                from interior_relocation_trial import repair
                assert audit(points,tets,surface)['oriented_facet_consistency_passed']
                phase=time.perf_counter();points,row['repair']=repair(points,tets,surface)
                row['repair_s']=time.perf_counter()-phase
            if config['order']==2:
                from microgen_tet10_control import straight_elevation
                phase=time.perf_counter();elevated=straight_elevation(points,tets,surface)
                if not elevated['within_node_cap']:raise ValueError(f'Straight Tet10 node cap exceeded: {elevated["nodes"]}')
                points,tets,surface=(elevated[k] for k in ('points','tetrahedra','surface'))
                row['straight_elevation_s']=time.perf_counter()-phase
                row['straight_elevation_explicit']=True
        else:
            field=spin_field(shape,meshers,np) if case_name=='spinodoid' else case.function
            if case_name=='smooth-lattice':
                from gradient_tie_probe import symmetric_lattice
                lattice=case.function if config.get('lattice_field')=='original' else symmetric_lattice
                def field(x,y,z):return np.minimum(lattice(x,y,z),.02)
                row['positive_void_clamp']=.02
                row['field_representation']=('Original fixed cubic smooth-union order' if config.get('lattice_field')=='original'
                    else 'Value-identical symmetric min; fixed original cubic smooth-union order')
            mesh=meshers.generate(field,bounds=w.BOUNDS,cells=config['cells'],band=case.band,
                periodic=(True,)*3,threads=1,optimize_passes=config.get('optimization_passes',1),minimum_quality=.1,max_tetrahedra=w.CAP,
                order=config['order'],geometry_tolerance=config['tolerance'],adapt=config['adapt'],
                min_length=config.get('min_length',.01),max_length=.1,gradation=.3,adapt_iterations=4,quality_guard=.3,adapt_tolerance=config.get('adapt_tolerance',.002))
            points,tets,surface=np.asarray(mesh.points),np.asarray(mesh.tetrahedra),np.asarray(mesh.surface)
            row['diagnostics']=mesh.diagnostics
        row.update(generation_wall_s=time.perf_counter()-started,nodes=len(points),elements=len(tets))
        row['generation_s']=row['generation_wall_s']-row['instrumentation_s']
        if config['arm']=='meshers' and case_name in ('spinodoid','smooth-lattice'):
            queries=np.random.default_rng(20261018).uniform(-.5,.5,(4096,3))
            expected=case.function(*queries.T)
            actual=field.values(queries) if hasattr(field,'values') else field(*queries.T)
            np.testing.assert_allclose(actual,np.minimum(expected,.02) if case_name=='smooth-lattice' else expected,rtol=1e-12,atol=1e-12)
            np.testing.assert_array_equal(actual<0,expected<0)
            row['fixed_field_identity_checked']=4096
        if len(tets)>w.CAP or len(points)>250000:raise ValueError('Final node/tetrahedron cap exceeded')
        np.savez_compressed(run/'mesh.npz',points=points,tetrahedra=tets,surface=surface)
        row['mesh_sha256']=w.digest(run/'mesh.npz')
        audit_start=time.perf_counter()
        row['geometry']=w.geometry(run,case_name,points,tets,surface,box_projection_fallback=True)
        q,angles,_=scorecard.tet_metrics(points,tets[:,:4])
        row['quality_quantiles']=dict(zip(('min','p01','p05','median','p95'),map(float,np.quantile(q,[0,.01,.05,.5,.95]))))
        row['elements_below_q01']=int((q<.1).sum())
        try:
            _,_,row['periodic']=periodic_pairs(points,surface)
        except Exception as error:row['periodic_error']=str(error)
        row['audit_s']=time.perf_counter()-audit_start
        row['qualified']=row['geometry']['geometry_admitted'] and 'periodic' in row
        row['status']='admitted' if row['qualified'] else 'geometry-or-periodicity-rejected'
    except Exception as error:
        row.update(status='failed',error=str(error),traceback=traceback.format_exc())
        row['failed_attempt_s']=time.perf_counter()-started
    row['wall_s']=time.perf_counter()-started
    w.save(run/'result.json',row)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--case',choices=w.CASES,required=True)
    parser.add_argument('--worker',type=json.loads)
    parser.add_argument('--repetitions',type=int,default=1)
    parser.add_argument('--strengthened',action='store_true')
    parser.add_argument('--warmups',type=int,default=0)
    parser.add_argument('--select',help='Comma-separated declared configuration names')
    parser.add_argument('--require-admitted',action='store_true',
        help='Fail after saving evidence if any selected mesh misses admission')
    args=parser.parse_args()
    if args.worker is not None:return worker(args.run,args.case,args.worker)
    args.run.mkdir(parents=True,exist_ok=False)
    w.imports()
    from job_parts import worker_part
    sys.path.insert(0,str(Path.home()/'bin'))
    from pueue_progress import Job
    import meshers._meshers,microgen.remesh
    import microgen.shape.spinodoid,microgen.shape.tpms
    paths_modules=[Path(microgen.shape.spinodoid.__file__),Path(microgen.shape.tpms.__file__)]
    assert not (args.run/'summary.json').exists()
    source=w.HERE/'runs/g3-microgen-defaults-v4'
    paths=[Path(__file__),Path(w.__file__),Path(meshers._meshers.__file__),Path(microgen.remesh.__file__),
        w.HERE/'scorecard.py',w.HERE/'test_mesh_only_comparison.py',w.HERE/'interior_relocation_trial.py',w.HERE/'interior_star_probe.py',w.HERE/'microgen_tet10_control.py',w.HERE/'box_projection.py',w.HERE/'gradient_tie_probe.py',
        w.HERE/'runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json']
    paths.extend(paths_modules)
    paths.extend((w.HERE/'facet_incidence_audit.py',w.ROOT.parent/'periodic-microgen/experiments/periodic-fea/periodic.py',
        Path(importlib.metadata.distribution('mmgpy').locate_file('mmgpy/bin/mmg3d_O3'))))
    paths.extend(source/(args.case+suffix) for suffix in ('-samples.npz','-specimen.json'))
    hashes={str(p):w.digest(p) for p in paths}
    w.save(args.run/'analysis-inputs.json',hashes)
    for index,path in enumerate(paths):
        if path.suffix=='.py':
            shutil.copyfile(path,args.run/f'input-{index}-{path.name}')
    w.save(args.run/'manifest.json',dict(inputs=hashes,platform=platform.platform(),architecture=platform.machine(),
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY,memory_bytes=w.MEMORY,worker_timeout_s=300,
        tetrahedron_cap=w.CAP,node_cap=250000,threads=1,packages={n:importlib.metadata.version(n) for n in ('meshers','microgen','mmgpy','numpy','scipy','pyvista','vtk')},
        no_fe_solve=True,case=args.case,repetitions=args.repetitions,warmups=args.warmups,strengthened=args.strengthened,select=args.select,command=sys.argv))
    shutil.copyfile(__file__,args.run/'analysis-driver.py')
    selected=configs(args.case,args.strengthened)
    if args.select:
        requested=args.select.split(',')
        assert set(requested)<=set(c['name'] for c in selected)
        selected=[c for c in selected if c['name'] in requested]
    tasks=[(f'repeat-{repeat}-{c["name"]}',c) for repeat in range(-args.warmups,args.repetitions) for c in (selected if repeat%2==0 else list(reversed(selected)))]
    def deadline(*_):raise TimeoutError('Whole mesh comparison deadline exceeded')
    signal.signal(signal.SIGALRM,deadline)
    signal.alarm(len(tasks)*310+60)
    rows=[]
    with Job(parts=[name for name,_ in tasks]) as job:
        for name,config in tasks:
            directory=args.run/name;directory.mkdir()
            for suffix in ('-samples.npz','-specimen.json'):
                shutil.copyfile(source/(args.case+suffix),directory/(args.case+suffix))
            with worker_part(job,name) as part:
                record=w.bounded([sys.executable,__file__,'--worker',json.dumps(config),'--run',str(directory),'--case',args.case],
                    directory/'worker.log',part,timeout=300)
                w.save(directory/'worker.json',record)
                row=json.loads((directory/'result.json').read_text()) if (directory/'result.json').exists() else dict(status='no-result')
                row.update(worker=record,directory=str(directory),warmup=name.startswith('repeat--'));rows.append(row)
                job.result(name,seconds=record['wall_s'],qualified=row.get('qualified',False),elements=row.get('elements'))
    w.save(args.run/'summary.json',dict(case=args.case,rows=rows,repetitions=args.repetitions,warmups=args.warmups,strengthened=args.strengthened,
        primary='Tet4 versus Tet4; curved Tet10 is a separate capability tier',
        timing_scope='Fresh in-process specimen/field preparation, extraction, meshing and output conversion. Imports, archive writes, validation, measured MMG instrumentation and FE excluded. Failed attempts are not successful mesh timings.',
        gates='distance.002 absolute-volume1percent q.1 positive sampled Jacobians coherent facets periodic nodes/caps;250000nodes/tets8GiB',
        geometry_auditor='box-constrained observer enabled identically for both arms; no changed geometry thresholds',
        dominance='Fewer elements with no worse generation time, quality quantiles and geometry/volume error is stronger than equal-count parity; no success ratio for rejected competitors',
        actual_default_and_explicit_variants_separate=True,new_equilibrium_solve=False))
    signal.alarm(0)
    assert all(w.digest(Path(p))==h for p,h in hashes.items())
    if args.require_admitted:
        assert all(row.get('qualified',False) for row in rows),'Mesh admission failed; read summary.json'


if __name__=='__main__':main()
