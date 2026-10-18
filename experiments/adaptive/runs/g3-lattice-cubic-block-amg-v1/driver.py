"""Test node-block AMG on one captured reduced system, not a stiffness tensor."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import signal
import sys
import time

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    import numpy as np
    import pyamg
    from scipy.sparse import load_npz
    from scipy.sparse.linalg import cg

    source = json.loads((args.source / "result.json").read_text())
    assert source["no_solve"] and source["complete_three_component_node_blocks"]
    assert workflow.digest(args.source / "matrix.npz") == source["matrix_archive_sha256"]
    assert workflow.digest(args.source / "vectors.npz") == source["vector_archive_sha256"]
    matrix = load_npz(args.source / "matrix.npz")
    vectors = np.load(args.source / "vectors.npz")
    rhs, permutation = vectors["rhs"], vectors["permutation"]
    np.testing.assert_array_equal(vectors["components"][permutation].reshape(-1, 3),
        np.tile([0, 1, 2], (source["free_nodes"], 1)))
    row = {"status": "reordering", "converged": False, "iterations": 0,
        "source": str(args.source), "method": "SA on node-interleaved BSR3 with three translations",
        "deadline_s": 240, "cg_rtol": 1e-12, "cg_atol": 1e-14,
        "scope": "One captured first-load linear solve only, not a homogenized tensor or live workflow"}
    phase = time.perf_counter()

    def save_progress():
        row["algorithm_wall_s"] = time.perf_counter() - phase
        workflow.save(args.run / "progress.json", row)

    def expired(signum, frame):
        raise TimeoutError("Node-block solver exceeds 240-second stage")

    previous_alarm = signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, row["deadline_s"])
    try:
        save_progress()
        reordered = matrix[permutation][:, permutation]
        block_matrix = reordered.tobsr(blocksize=(3, 3))
        del reordered
        row["reorder_s"] = time.perf_counter() - phase
        row["status"] = "building-hierarchy"
        save_progress()
        np.random.seed(20261018)
        candidates = np.tile(np.eye(3), (source["free_nodes"], 1))
        started = time.perf_counter()
        hierarchy = pyamg.smoothed_aggregation_solver(block_matrix, B=candidates, symmetry="symmetric")
        row.update(hierarchy_setup_s=time.perf_counter() - started,
            levels=[level.A.shape[0] for level in hierarchy.levels],
            operator_complexity=hierarchy.operator_complexity(), status="cg")
        save_progress()

        def count_iteration(_):
            row["iterations"] += 1
            if row["iterations"] % 10 == 0:
                save_progress()

        started = time.perf_counter()
        solution, info = cg(block_matrix, rhs[permutation], M=hierarchy.aspreconditioner(cycle="V"),
            rtol=1e-12, atol=1e-14, maxiter=20000, callback=count_iteration)
        row["cg_s"] = time.perf_counter() - started
        original_solution = np.empty_like(solution)
        original_solution[permutation] = solution
        residual = float(np.linalg.norm(matrix @ original_solution - rhs) / np.linalg.norm(rhs))
        row.update(status="converged" if info == 0 and residual < 2e-9 else "residual-rejected",
            converged=bool(info == 0 and residual < 2e-9), relative_residual=residual, scipy_info=int(info))
        np.savez_compressed(args.run / "solution.npz", solution=original_solution)
    except TimeoutError:
        row["status"] = "algorithm-timeout"
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_alarm)
    row["algorithm_wall_s"] = time.perf_counter() - phase
    workflow.save(args.run / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), Path(workflow.__file__)] + [args.source / name
        for name in ("result.json", "matrix.npz", "vectors.npz")]
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": 360, "algorithm_deadline_s": 240, "threads": 1,
        "pyamg": importlib.metadata.version("pyamg"), "random_seed": 20261018})
    with Job(parts=["node-block-amg"]) as job:
        with job.part("node-block-amg") as part:
            command = [sys.executable, __file__, "--run", str(args.run), "--source", str(args.source), "--worker"]
            record = workflow.bounded(command, args.run / "worker.log", part, timeout=360)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("node-block-amg", seconds=record["wall_s"], converged=int(row["converged"]),
                iterations=row["iterations"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
