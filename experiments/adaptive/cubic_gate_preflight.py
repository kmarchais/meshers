"""Verify existing cubic embedding and actual reduced setup before any CG solve."""

import argparse
import json
from pathlib import Path
import shutil
from unittest.mock import patch

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--prepared',type=Path,required=True)
    args = parser.parse_args()
    args.run,args.source,args.prepared = (p.resolve() for p in (args.run,args.source,args.prepared))
    assert args.run.is_dir() and not (args.run/'result.json').exists()
    np,*_ = w.imports()
    from amg_reference import AuditedAMG
    from capture_reference_system import CapturedSystem
    from cubic_reference import register_cubic, basis_checks, reference_worker
    from quadrature_audit import duffy_rule
    from stress_recovery import solve as original_solve
    assembled = json.loads((args.prepared/'result.json').read_text())
    assert assembled['assembly_only'] and assembled['no_fea_solve']
    assert assembled['nodes'] <= 250000 and assembled['elements'] <= 250000
    paths = [Path(__file__),args.prepared/'result.json',args.prepared/'mesh.npz',args.prepared/'manifest.json']
    paths.extend(args.source/name for name in ('result.json','mesh.npz','fields.npz',
        'smooth-lattice-samples.npz','smooth-lattice-specimen.json'))
    paths.extend(w.HERE/name for name in ('cubic_reference.py','stress_recovery.py','amg_reference.py',
        'tet10_integration.py','microgen_workflow.py','quadrature_audit.py'))
    hashes = {str(p.resolve()):w.digest(p) for p in paths}
    manifest = json.loads((args.prepared/'manifest.json').read_text())
    for name,expected in manifest['sources'].items():
        p = Path(name)
        if not p.is_absolute(): p = w.ROOT/p
        assert w.digest(p) == expected,str(p)
        hashes[str(p.resolve())] = expected
    for i,p in enumerate(paths):
        if p.suffix == '.py': shutil.copyfile(p,args.run/f'input-{i}-{p.name}')
    target = args.run/'mesh.npz'
    shutil.copyfile(args.prepared/'mesh.npz',target)
    assert w.digest(target) == w.digest(args.prepared/'mesh.npz')
    element = register_cubic()()
    checks = basis_checks(element)
    mesh = np.load(target)
    print('Audit actual cubic geometry with positive125 and unchanged specimen samples',flush=True)
    geometry = w.geometry(args.source,'smooth-lattice',mesh['points'],mesh['tetrahedra'],mesh['surface'],
        element_type=element.name,quadrature=duffy_rule(5))
    assert geometry['geometry_admitted']
    parent = json.loads((args.source/'result.json').read_text())
    assert abs(geometry['volume']/parent['geometry']['volume']-1) < 1e-10
    resources = {}

    class Capture:
        def __init__(self): self.problem,self.records = None,[]

        def __call__(self,matrix,rhs,**kwargs):
            matrix = matrix.tocsr()
            amg = AuditedAMG(elastic_blocks=True);amg.problem = self.problem
            with patch('scipy.sparse.linalg.cg',side_effect=CapturedSystem) as intercepted:
                try: amg(matrix,rhs)
                except CapturedSystem:
                    assert intercepted.call_count == 1 and not amg.records
                else: raise AssertionError('Stop at the CG boundary')
            resources.update(reduced_shape=list(matrix.shape),nnz=int(matrix.nnz),matrix_sha256=amg.key[1],
                csr_bytes=matrix.data.nbytes+matrix.indices.nbytes+matrix.indptr.nbytes,
                amg_levels=len(amg.hierarchy.levels),operator_complexity=amg.hierarchy.operator_complexity(),
                amg_setup_passed=True,cg_started=False)
            raise CapturedSystem()

    def capture_solve(*args,**kwargs):
        kwargs['solver'] = Capture()
        return original_solve(*args,**kwargs)

    print('Reuse full parent mapping/field embedding checks, then stop actual reduced AMG before CG',flush=True)
    try:
        with patch('stress_recovery.solve',side_effect=capture_solve):
            reference_worker(args.run,args.source,args.prepared,True,True)
    except CapturedSystem:
        assert resources and resources['amg_setup_passed']
    else: raise AssertionError('Do not return a fake cubic solution')
    embedded_path = args.run/'embedding.json'
    (args.run/'result.json').rename(embedded_path)
    embedded = json.loads(embedded_path.read_text())
    assert embedded['status'] == 'embedded-parent-verified' and 'stiffness' not in embedded
    assert embedded['exact_quadratic_mapping_embedding_error'] < 1e-12
    assert embedded['embedded_parent_tensor_error'] < 1e-7
    for p in (target,embedded_path):hashes[str(p.resolve())] = w.digest(p)
    w.save(args.run/'analysis-inputs.json',hashes)
    w.save(args.run/'result.json',dict(passed=True,qualified=True,source=str(args.source),prepared=str(args.prepared),
        nodes=assembled['nodes'],elements=assembled['elements'],order=3,element_type=element.name,
        basis_checks=checks,geometry=geometry,embedding=embedded,resources=resources,
        geometry_unchanged=True,old_quadratic_space_preserved=True,cg_started=False,
        reference_promoted=False,native_cubic_api=False,comparative_timing_claim=False,inputs=hashes,
        scope='Actual fixed-domain cubic geometry/field/reduced-setup gates only; no equilibrium solution or reference promotion'))
    assert all(w.digest(Path(p)) == h for p,h in hashes.items())
    print('Cubic geometry and reduced setup pass; no CG solve started',flush=True)


if __name__ == '__main__':
    main()
