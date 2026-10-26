"""Compare the real compact exec workflow with the unchanged six-load baseline."""

import argparse
import json
from pathlib import Path
import numpy as np
import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    w.imports()
    from periodic import tensor_difference
    row, old = [json.loads((p/'result.json').read_text()) for p in (args.source, args.baseline)]
    assert row['status'] == 'solved-needs-reference' and row['completed_loads'] == 6
    assert row['generation_process_id'] == row['solver_process_id']
    assert row['solve_process_boundary'] == 'exec-replace'
    assert row['solver_memory_after_exec_bytes'] < row['generation_memory_before_exec_bytes']
    worker = json.loads((args.source/'worker.json').read_text())
    assert worker['returncode'] == 0 and worker['limit_failure'] is None
    assert worker['memory_bytes'] == 8*1024**3
    meshes = [np.load(p/'mesh.npz') for p in (args.source, args.baseline)]
    assert all(np.array_equal(meshes[0][k], meshes[1][k]) for k in meshes[0].files)
    fields = [np.load(p/'fields.npz')['displacement'] for p in (args.source, args.baseline)]
    np.testing.assert_allclose(fields[0], fields[1], rtol=1e-7, atol=1e-8)
    difference = tensor_difference(np.asarray(row['stiffness']), np.asarray(old['stiffness']))
    assert difference['worst_relative_strain_energy_difference'] < 1e-8
    assert row['total_workflow_s'] >= row['generation_s']+row['solve_s']+row['exec_and_reload_seconds']
    assert row['solver_records'][0]['reduced_matrix_sha256'] == old['solver_records'][0]['reduced_matrix_sha256']
    paths = [Path(__file__)]+[p/n for p in (args.source, args.baseline)
        for n in ('result.json', 'mesh.npz', 'fields.npz', 'worker.json')]
    w.save(args.output, dict(passed=True, identical_mesh=True, same_pid_exec=True,
        tensor_difference=difference, maximum_field_difference=float(abs(fields[0]-fields[1]).max()),
        generation_rss_bytes=row['generation_memory_before_exec_bytes'],
        fresh_solver_rss_bytes=row['solver_memory_after_exec_bytes'],
        complete_workflow_time_includes_exec=True, memory_limit_unchanged=True,
        inputs={str(p):w.digest(p) for p in paths}))
    print(f'Exec control passes; worst tensor difference {difference["worst_relative_strain_energy_difference"]:.3g}', flush=True)


if __name__ == '__main__':
    main()
