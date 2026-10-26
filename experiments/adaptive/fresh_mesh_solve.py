"""Continue a live workflow after exec releases generation's process memory."""

import argparse
import json
import os
from pathlib import Path
import time
import traceback

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', required=True, type=Path)
    parser.add_argument('--started', required=True, type=float)
    args = parser.parse_args()
    np, *_ = w.imports()
    import psutil
    from amg_reference import AuditedAMG
    from stress_recovery import solve, vector_first_reactions

    row = json.loads((args.run/'result.json').read_text())
    assert row['solve_process_boundary'] == 'exec-replace'
    assert row['generation_process_id'] == os.getpid(), 'Only continue the same live workflow PID'
    assert 0 < args.started <= row['solve_exec_clock'] <= time.perf_counter()
    assert row['qualified'] and row['geometry']['geometry_admitted']
    assert row['order'] == 2 and row['solver_quadrature_points'] == 15
    assert row['solver'] == 'cached-elastic-block-amg-cg'
    assert w.digest(args.run/'mesh.npz') == row['solve_mesh_sha256']
    row['solver_process_id'] = os.getpid()
    row['solver_memory_after_exec_bytes'] = psutil.Process().memory_info().rss
    mesh = np.load(args.run/'mesh.npz')
    assert len(mesh['points']) == row['nodes'] <= row['solve_node_cap'] <= w.CAP
    assert len(mesh['tetrahedra']) == row['elements'] <= w.CAP
    try:
        row.update(status='solving', completed_loads=0)
        w.save(args.run/'result.json', row)
        def checkpoint(completed, partial, records):
            row.update(completed_loads=completed, partial_stiffness_columns=partial.tolist(), solver_records=records)
            w.save(args.run/'result.json', row)
        def save_field(completed, displacement):
            np.savez_compressed(args.run/f'load-{completed}.npz', displacement=displacement)
        started = time.perf_counter()
        fields, tensor, records = solve(mesh['points'], mesh['tetrahedra'], checkpoint=checkpoint,
            n_elm_gp=15, solver=AuditedAMG(elastic_blocks=True), field_checkpoint=save_field,
            reaction_reader=vector_first_reactions, release_reduced_cache=True)
        solve_seconds = time.perf_counter()-started
        np.savez_compressed(args.run/'fields.npz', displacement=fields)
        row.pop('partial_stiffness_columns', None)
        row.update(solve_s=solve_seconds, stiffness=tensor.tolist(),
            solver_records=records, status='solved-needs-reference')
    except Exception as error:
        row.update(status='failed', qualified=False, error=str(error), traceback=traceback.format_exc())
    row['total_audit_workflow_s'] = row['total_workflow_s'] = time.perf_counter()-args.started
    row['exec_and_reload_seconds'] = started-row['solve_exec_clock'] if 'started' in locals() else None
    assert w.digest(args.run/'mesh.npz') == row['solve_mesh_sha256']
    w.save(args.run/'result.json', row)


if __name__ == '__main__':
    main()
