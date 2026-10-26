"""Fresh operator, six-load solve and all-model assessment of saved admitted meshes."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time

import numpy as np

import microgen_workflow as w


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--stage',choices=('resource','solve','assess'),required=True)
    parser.add_argument('--resource',type=Path);parser.add_argument('--field-audit',type=Path)
    args=parser.parse_args();assert args.run.is_dir() and not (args.run/'result.json').exists()
    w.imports();shutil.copyfile(__file__,args.run/'evaluation-driver.py')
    row=json.loads((args.source/'result.json').read_text())
    assert row['qualified'] and row['geometry']['geometry_admitted'] and row['order']==2
    assert row['nodes']<=w.CAP and row['elements']<=w.CAP
    paths=[Path(__file__),args.source/'result.json',args.source/'mesh.npz',
        w.HERE/'stress_recovery.py',w.HERE/'amg_reference.py',w.HERE/'microgen_workflow.py']
    with np.load(args.source/'mesh.npz') as mesh:
        points,cells=mesh['points'],mesh['tetrahedra']
    assert (len(points),len(cells))==(row['nodes'],row['elements'])
    if args.stage=='resource':
        from operator_resource_probe import probe
        resource=probe(points,cells,args.run/'operator.json')
        result=dict(passed=True,source=str(args.source),nodes=len(points),elements=len(cells),
            resources=resource,cg_started=False)
        paths.append(w.HERE/'operator_resource_probe.py')
    elif args.stage=='solve':
        assert args.resource is not None
        from amg_reference import AuditedAMG
        from stress_recovery import solve,vector_first_reactions
        resource=json.loads((args.resource/'operator.json').read_text())
        resource_row=json.loads((args.resource/'result.json').read_text())
        assert resource_row['passed'] and resource_row['inputs'][str(args.source/'mesh.npz')]==w.digest(args.source/'mesh.npz')
        for name in ('mesh.npz',row['case']+'-samples.npz',row['case']+'-specimen.json'):
            path=args.source/name;shutil.copyfile(path,args.run/name)
            assert w.digest(path)==w.digest(args.run/name);paths.append(path)
        result=dict(row,status='solving',completed_loads=0,geometry_source=str(args.source),
            solver_quadrature_points=15,solver='cached-elastic-block-amg-cg',equilibrium_fields_replayed=False,
            reference_promoted=False,comparative_timing_claim=False,old_approximation_space_preserved=False,
            scope='Fresh six-load fields on a saved admitted experimental mesh; no live generation timing, nesting claim or reference promotion')
        w.save(args.run/'result.json',result)
        def checkpoint(completed,partial,records):
            result.update(completed_loads=completed,partial_stiffness_columns=partial.tolist(),solver_records=records)
            w.save(args.run/'result.json',result)
        def field_checkpoint(completed,field):
            np.savez_compressed(args.run/f'load-{completed}.npz',displacement=field)
        started=time.monotonic()
        fields,tensor,records=solve(points,cells,checkpoint,n_elm_gp=15,
            solver=AuditedAMG(elastic_blocks=True),field_checkpoint=field_checkpoint,
            reaction_reader=vector_first_reactions,release_reduced_cache=True)
        assert records[0]['reduced_matrix_sha256']==resource['matrix_sha256']
        assert len(records)==6 and not records[0]['hierarchy_reused'] and sum(r['hierarchy_reused'] for r in records)==5
        assert fields.shape==(6,len(points),3) and np.isfinite(fields).all()
        np.savez_compressed(args.run/'fields.npz',displacement=fields)
        result.pop('partial_stiffness_columns',None)
        result.update(status='solved-needs-reference',stiffness=tensor.tolist(),solver_records=records,
            saved_mesh_solve_s=time.monotonic()-started,completed_loads=6,new_equilibrium_fields=True,equilibrium_solve_performed=True)
        paths.extend((args.resource/'operator.json',args.resource/'result.json'))
    else:
        assert args.field_audit is not None and row['case']=='smooth-lattice'
        from cubic_refinement import MODEL_REPORTS
        from periodic import tensor_difference
        assert row['completed_loads']==6 and len(row['solver_records'])==6
        assert max(r['relative_residual'] for r in row['solver_records'])<2e-9
        fields=np.load(args.source/'fields.npz')['displacement']
        assert fields.shape==(6,len(points),3) and np.isfinite(fields).all()
        consistency_path=args.field_audit/(args.source.name+'-consistency.json')
        quadrature_path=args.field_audit/(args.source.name+'.json')
        consistency=json.loads(consistency_path.read_text());quadrature=json.loads(quadrature_path.read_text())
        assert consistency['reaction_relative_error']<1e-8 and consistency['hill_mandel_relative_error']<1e-8
        assert consistency['periodic_jump_max']<1e-9 and consistency['periodic_mid_edge_jump_max']<1e-9
        higher=quadrature['rules'][-1];assert higher['points']==125
        assert higher['energy_vs_solved_tensor']['worst_relative_strain_energy_difference']<1e-4
        assert higher['energy_vs_previous_rule']['worst_relative_strain_energy_difference']<1e-4
        models=[];hashes=[]
        for filename in MODEL_REPORTS:
            path=w.HERE/'runs'/filename;report=json.loads(path.read_text())
            assert report['case']=='smooth-lattice' and report['stiffness_target']==.01 and len(report['models'])==8
            for index,model in enumerate(report['models']):
                matrix=np.asarray(model['stiffness'],dtype=np.float64)
                assert matrix.shape==(6,6) and np.linalg.eigvalsh(matrix).min()>0
                hashes.append(hashlib.sha256(matrix.tobytes()).hexdigest())
                models.append(dict(source=str(path),index=index,difference=tensor_difference(np.asarray(row['stiffness']),matrix)))
            paths.append(path)
        previous_path=w.HERE/'runs/macbook-lattice-quadratic-quality-assessment-v1/result.json'
        previous=json.loads(previous_path.read_text())
        assert len(models)==len(set(hashes))==24 and hashes==previous['retained_model_tensor_hashes']
        worst=max(m['difference']['worst_relative_strain_energy_difference'] for m in models)
        result=dict(passed=True,case='smooth-lattice',source=str(args.source),models=models,
            retained_model_count=24,retained_model_tensor_hashes=hashes,retained_models_unchanged=True,
            worst_model_disagreement=worst,previous_best_disagreement=previous['worst_model_disagreement'],
            passes_one_percent_all_models=worst<=.01,quadrature_sensitivity=higher['energy_vs_previous_rule']['worst_relative_strain_energy_difference'],
            field_audit_passed=True,reference_promoted=False,comparative_timing_claim=False,objective_complete=False)
        paths.extend((previous_path,consistency_path,quadrature_path,args.source/'fields.npz'))
    result['inputs']={str(p):w.digest(p) for p in paths}
    w.save(args.run/'result.json',result)


if __name__=='__main__':main()
