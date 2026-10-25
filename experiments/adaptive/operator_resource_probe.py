"""Measure an actual saved Tet10 operator and AMG setup in a fresh process."""

import argparse
import json
from pathlib import Path
import shutil
from unittest.mock import patch

import microgen_workflow as w


def probe(points,cells,output):
    import psutil
    from amg_reference import AuditedAMG,reduced_matrix_key
    from capture_reference_system import CapturedSystem
    from stress_recovery import solve
    record={}

    class Capture:
        def __init__(self):self.problem,self.records=None,[]
        def __call__(self,matrix,rhs,**kwargs):
            matrix=matrix.tocsr()
            record.update(reduced_shape=list(matrix.shape),nnz=int(matrix.nnz),
                matrix_sha256=reduced_matrix_key(matrix)[1],
                csr_bytes=matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes,
                rss_before_amg_bytes=psutil.Process().memory_info().rss,cg_started=False,
                amg_setup_passed=False)
            w.save(output,record)
            print(f'Reduced operator ready: {matrix.shape[0]} rows, {matrix.nnz} nonzeros; setup next',flush=True)
            amg=AuditedAMG(elastic_blocks=True);amg.problem=self.problem
            with patch('scipy.sparse.linalg.cg',side_effect=CapturedSystem) as intercepted:
                try:amg(matrix,rhs)
                except CapturedSystem:assert intercepted.call_count==1 and not amg.records
                else:raise AssertionError('Stop at CG boundary')
            assert amg.key[1]==record['matrix_sha256']
            record.update(levels=len(amg.hierarchy.levels),
                operator_complexity=amg.hierarchy.operator_complexity(),amg_setup_passed=True,
                rss_after_amg_bytes=psutil.Process().memory_info().rss)
            w.save(output,record)
            raise CapturedSystem()

    print('Assemble the saved quadratic mesh in a fresh process; no CG solve',flush=True)
    try:solve(points,cells,n_elm_gp=15,solver=Capture(),release_reduced_cache=True)
    except CapturedSystem:assert record['amg_setup_passed']
    else:raise AssertionError('Do not return a fake solution')
    return record


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',required=True,type=Path)
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--prepared',type=Path)
    source.add_argument('--finalize-trial',action='store_true')
    args=parser.parse_args()
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np,*_=w.imports()
    prepared=args.run if args.finalize_trial else args.prepared
    paths=[Path(__file__),prepared/'mesh.npz',prepared/'manifest.json']
    paths.extend([prepared/'trial-screen.json',prepared/'analysis-inputs.json'] if args.finalize_trial
        else [prepared/'worker.json'])
    paths.extend(w.HERE/name for name in ('stress_recovery.py','amg_reference.py','microgen_workflow.py'))
    hashes={str(p.resolve()):w.digest(p) for p in paths}
    for i,p in enumerate(paths):
        if p.suffix=='.py':shutil.copyfile(p,args.run/f'input-{i}-{p.name}')
    trial=None
    if args.finalize_trial:
        trial=json.loads((prepared/'trial-screen.json').read_text())
        assert trial['status']=='trial-verified-resource-pending' and not trial['qualified']
        assert trial['trial_signal_passed'] and trial['geometry']['geometry_admitted']
        assert all(w.digest(Path(p))==h for p,h in trial['inputs'].items())
    else:
        failed=json.loads((prepared/'worker.json').read_text())
        assert failed['limit_failure']=='memory-cap' and failed['returncode']!=0
    w.save(args.run/'resource-inputs.json',hashes)
    mesh=np.load(prepared/'mesh.npz');points,cells=mesh['points'],mesh['tetrahedra']
    assert len(points)<=250000 and len(cells)<=250000 and cells.shape[1]==10
    resource=probe(points,cells,args.run/'operator.json')
    if trial is not None:
        trial.update(status='nested-preflight-passed',qualified=True,resources=resource,
            fresh_process_resource_setup=True,inputs={**trial['inputs'],**hashes})
        w.save(args.run/'result.json',trial)
    else:
        w.save(args.run/'analysis-inputs.json',hashes)
        w.save(args.run/'result.json',dict(passed=True,prepared=str(prepared.resolve()),
        nodes=len(points),elements=len(cells),resources=resource,inputs=hashes,
        cg_started=False,admission_performed=False,reference_promoted=False,
        comparative_timing_claim=False,scope='Fresh-process resource diagnosis of the exact saved mesh; not admission or an equilibrium result'))
    assert all(w.digest(Path(p))==h for p,h in hashes.items())


if __name__=='__main__':main()
