"""Measure cell batching before selecting a bounded positive-rule trial."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
import microgen_workflow as workflow
from hp_parent25_probe import PREPARED, H_SOLVE, P_MESH, P_SOLVE


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--preflight', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from fedoo.lib_elements.tetrahedron import Tet10
    from hp_cross_energy import integrate
    from hp_batched_energy import integrate_batched

    preflight = json.loads((args.preflight / 'result.json').read_text())
    assert preflight['passed'] and preflight['parent25_shared_domain_preflight']
    inputs = json.loads((args.preflight / 'analysis-inputs.json').read_text())
    resolved = {}
    for name, expected in inputs.items():
        path = Path(name)
        if workflow.digest(path) != expected:
            assert path.name == 'test_hp_cross_energy.py'
            snapshot = args.output.parent / 'preflight-test_hp_cross_energy.py'
            snapshot.write_bytes(subprocess.check_output(['git', 'show', 'HEAD:'+str(path.relative_to(workflow.ROOT))]))
            path = snapshot
        assert workflow.digest(path) == expected
        resolved[str(path)] = expected
    for path in (Path(__file__), workflow.HERE/'hp_batched_energy.py', workflow.HERE/'hp_cross_energy.py',
                 workflow.HERE/'test_hp_cross_energy.py', args.preflight/'result.json'):
        resolved[str(path)] = workflow.digest(path)
        if path.suffix == '.py':
            shutil.copyfile(path, args.output.parent / path.name)
    workflow.save(args.output.parent/'analysis-inputs.json', resolved)
    mesh, cubic = (dict(np.load(directory/'mesh.npz')) for directory in (H_SOLVE, P_MESH))
    h, p = (np.load(directory/'fields.npz')['displacement'] for directory in (H_SOLVE, P_SOLVE))
    maps = dict(np.load(PREPARED/'prolongation.npz'))
    element = Tet10(4)
    arguments = (mesh['points'], mesh['tetrahedra'], h, cubic['tetrahedra'], p,
        maps['parents'], maps['parent_barycentric_maps'], (element.xi_pg, element.w_pg))
    begin = time.perf_counter()
    grouped = integrate(*arguments)
    grouped_s = time.perf_counter()-begin
    begin = time.perf_counter()
    batched = integrate_batched(*arguments)
    batched_s = time.perf_counter()-begin
    np.testing.assert_allclose(batched[0], grouped[0], atol=1e-12, rtol=1e-10)
    np.testing.assert_allclose(batched[1:3], grouped[1:3], atol=1e-12, rtol=1e-10)
    basis = np.vstack((-np.eye(6), np.eye(6)))
    np.testing.assert_allclose(batched[3].sum(axis=0), basis.T @ grouped[0] @ basis, atol=1e-12, rtol=1e-9)
    np.testing.assert_allclose(batched[4].sum(), grouped[1], atol=1e-12)
    workflow.save(args.output, dict(passed=True, preflight=str(args.preflight), batch_size=256,
        grouped_s=grouped_s, batched_s=batched_s, speedup=grouped_s/batched_s,
        full_energy_maximum_absolute_difference=float(abs(batched[0]-grouped[0]).max()),
        conservative_positive_rule_time_estimate_s=2*batched_s*(125+216)/4,
        four_point_resource_probe_only=True, new_equilibrium_solve=False, reference_promoted=False))
    assert all(workflow.digest(Path(name)) == expected for name, expected in resolved.items())


if __name__ == '__main__':
    main()
