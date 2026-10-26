"""Check cap recovery through Python and replay the four unchanged defaults."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as w

BASELINE = '19f1536993017f6f610e7b199a3cd271b9e322fb60094fa4cb7984df59f567e8'


def box_control(expect_recovery):
    import meshers
    import meshers._meshers
    import numpy as np
    rows = []
    for fine in (.2, .1):
        try:
            mesh = meshers.generate(lambda x, y, z: -1. + 0*x, cells=4,
                optimize_passes=0, minimum_quality=.1, max_tetrahedra=768,
                adapt=True, min_length=.1, max_length=.6, gradation=.5,
                adapt_iterations=1, quality_guard=.1, adapt_tolerance=.01,
                size_field=lambda x, y, z: fine + (.6-fine)*x)
        except meshers.MeshingError as error:
            assert 'exceeded max_tetrahedra' in str(error), str(error)
            assert fine == .1 or not expect_recovery, str(error)
            rows.append(dict(fine=fine, status='budget-rejected', error=str(error)))
        else:
            assert expect_recovery and fine == .2
            p, t = np.asarray(mesh.points), np.asarray(mesh.tetrahedra)
            from scorecard import tet_metrics
            q, _, volumes = tet_metrics(p, t)
            assert len(t) <= 768 and q.min() >= .1 and volumes.min() > 0
            assert abs(volumes.sum()-1.) < 1e-12
            rows.append(dict(fine=fine, status='admitted', nodes=len(p), tets=len(t),
                minimum_quality=float(q.min()), volume=float(volumes.sum())))
    print(json.dumps(dict(binary_sha256=w.digest(Path(meshers._meshers.__file__)), rows=rows)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', type=Path)
    parser.add_argument('--box-worker', choices=('baseline', 'candidate'))
    args = parser.parse_args()
    if args.box_worker:
        return box_control(args.box_worker == 'candidate')
    assert args.run is not None
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/'driver.py')
    result = dict(scope='Cap correctness and exact default-array regression; no FE or timing claim',
        inputs={str(Path(__file__)): w.digest(Path(__file__))}, box={}, cases={})
    for kind, module_path in (
        ('baseline', w.ROOT/'.local-work/native-baselines'/BASELINE),
        ('candidate', w.ROOT/'crates/meshers-python/python')):
        env = dict(os.environ, PYTHONPATH=str(module_path))
        completed = subprocess.run([sys.executable, __file__, '--box-worker', kind],
            env=env, capture_output=True, text=True, timeout=60)
        (args.run/f'{kind}-box.log').write_text(completed.stdout+completed.stderr)
        assert completed.returncode == 0, completed.stdout+completed.stderr
        result['box'][kind] = json.loads(completed.stdout)
    assert result['box']['baseline']['binary_sha256'] == BASELINE
    assert result['box']['candidate']['binary_sha256'] != BASELINE
    import mesh_only_comparison as comparison
    import numpy as np
    receipt_path = w.ROOT/'handoff/macbook/adaptive-interior-validation-20261025.json'
    receipt = json.loads(receipt_path.read_text())
    result['inputs'][str(receipt_path)] = w.digest(receipt_path)
    for label, case in zip(('gyroid', 'schwarz', 'spinodoid', 'lattice'), w.CASES, strict=True):
        old = w.ROOT/receipt['cases'][label]['run']
        previous = json.loads((old/'result.json').read_text())
        assert w.digest(old/'mesh.npz') == receipt['inputs'][str(old.relative_to(w.ROOT)/'mesh.npz')]
        run = args.run/label
        run.mkdir()
        for suffix in ('-samples.npz', '-specimen.json'):
            source = old/(case+suffix)
            shutil.copyfile(source, run/source.name)
            assert w.digest(source) == w.digest(run/source.name)
            result['inputs'][str(source)] = w.digest(source)
        comparison.worker(run, case, previous['config'])
        row = json.loads((run/'result.json').read_text())
        assert row['qualified'], row
        with np.load(old/'mesh.npz') as before, np.load(run/'mesh.npz') as after:
            for key in ('points', 'tetrahedra', 'surface'):
                np.testing.assert_array_equal(before[key], after[key])
        result['cases'][case] = dict(admitted=True, arrays_exactly_unchanged=True,
            nodes=row['nodes'], tets=row['elements'], source=str(old))
        w.save(args.run/'result.json', result)
    for source in (Path(comparison.__file__), Path(w.__file__),
            w.ROOT/'crates/meshers-core/src/tet_adapt.rs'):
        result['inputs'][str(source)] = w.digest(source)
    w.save(args.run/'result.json', result)


if __name__ == '__main__':
    main()
