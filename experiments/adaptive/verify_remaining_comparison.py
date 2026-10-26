"""Regenerate the three additional admitted curved fixtures after the native fix."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
import microgen_workflow as w
from comparison_benchmark import configs,worker


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True)
    p.add_argument('--all',action='store_true',help='Regenerate all seven historical admitted curved cases and run changed API controls')
    a=p.parse_args()
    cases={'spinodoid-anisotropic':'macbook-verified-spinodoid-anisotropic-repeated-v1',
        'graded-gyroid-periodic':'macbook-verified-graded-periodic-repeated-v1',
        'graded-gyroid-linear':'macbook-verified-graded-linear-repeated-v1'}
    selected_controls=False
    if a.all:
        baseline=json.loads((w.ROOT/'handoff/macbook/comparison-findings-20261026.json').read_text())
        cases={r['case']:str(Path(r['source_rows'][0]).parent) for r in baseline['rows']
               if r['name']=='meshers-curved-tet10' and r['qualified']==r['samples'] and r['samples']>=5}
        assert len(cases)==7
        subprocess.run([sys.executable,'-m','pytest','-q','--tb=short',
            str(w.ROOT/'crates/meshers-python/tests/test_interp.py'),
            str(w.ROOT/'crates/meshers-python/tests/test_periodic_zero_level.py')],check=True,timeout=60)
        selected_controls=True
    fixtures=w.HERE/'runs/macbook-comparison-fixtures-v1';rows=[]
    inputs={str(p):w.digest(p) for p in (Path(__file__),Path(worker.__code__.co_filename),Path(w.__file__),
        w.HERE/'comparison_fixtures.py',w.HERE/'comparison_mesh_audit.py')}
    for case,previous in cases.items():
        config=next(c for c in configs(case) if c['name']=='meshers-curved-tet10')
        old=(w.ROOT/previous if a.all else w.HERE/'runs'/previous)/'repeat-0-meshers-curved-tet10'
        before=json.loads((old/'result.json').read_text());assert before['qualified']
        assert w.digest(old/'mesh.npz')==before['mesh_sha256']
        inputs[str(old/'mesh.npz')]=before['mesh_sha256']
        directory=a.run/case;directory.mkdir()
        for source in (fixtures/case).iterdir():shutil.copyfile(source,directory/source.name)
        worker(directory,case,config,fixtures)
        after=json.loads((directory/'result.json').read_text());assert after['qualified'],after
        with np.load(old/'mesh.npz') as original,np.load(directory/'mesh.npz') as current:
            for key in ('points','tetrahedra','surface'):np.testing.assert_array_equal(current[key],original[key])
        rows.append(dict(case=case,source=str(old),old_mesh_sha256=before['mesh_sha256'],
            new_mesh_sha256=after['mesh_sha256'],admitted=True,all_arrays_exactly_unchanged=True))
    import meshers._meshers
    native=Path(meshers._meshers.__file__);inputs[str(native)]=w.digest(native)
    expression=Path(meshers.expression.__file__);inputs[str(expression)]=w.digest(expression)
    w.save(a.run/'result.json',dict(rows=rows,passed=True,inputs=inputs,no_new_FE_solve=True,
        changed_API_and_Schwarz_D_controls=selected_controls))
