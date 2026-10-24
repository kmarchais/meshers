"""Capture one actual reduced FE system without returning a fake solution."""

import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import sys

import microgen_workflow as workflow
from job_parts import worker_part

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


class CapturedSystem(Exception):
    pass


def worker(args):
    np, *_ = workflow.imports()
    from scipy.sparse import save_npz
    from cubic_reference import register_cubic
    from stress_recovery import solve

    archive = np.load(args.prepared / "mesh.npz")
    points, cells = archive["points"], archive["tetrahedra"]
    assert len(points) <= 250000 and len(cells) <= workflow.CAP and cells.shape[1] in (10, 20)
    element_type = None
    if cells.shape[1] == 20:
        register_cubic()
        element_type = "tet20_reference"
    else:
        prepared = json.loads((args.prepared / "result.json").read_text())
        assert prepared["qualified"] and prepared["geometry"]["geometry_admitted"]

    class Capture:
        def __init__(self):
            self.problem = None
            self.records = []

        def __call__(self, matrix, rhs, **kwargs):
            matrix = matrix.tocsr()
            free = self.problem._dof_free
            assert len(free) == len(rhs) and np.all(free < self.problem.n_node_dof)
            components = free // len(points)
            nodes = free % len(points)
            permutation = np.lexsort((components, nodes))
            grouped_nodes = nodes[permutation].reshape(-1, 3)
            np.testing.assert_array_equal(grouped_nodes, np.repeat(grouped_nodes[:, :1], 3, axis=1))
            np.testing.assert_array_equal(components[permutation].reshape(-1, 3),
                np.tile(self.problem.space.get_rank_vector("Disp"), (len(grouped_nodes), 1)))
            assert np.isfinite(matrix.data).all() and (matrix.diagonal() > 0).all()
            if not args.metrics_only:
                save_npz(args.run / "matrix.npz", matrix)
                np.savez_compressed(args.run / "vectors.npz", rhs=rhs, permutation=permutation,
                    coordinates=points[nodes], components=components, free_dofs=free)
            workflow.save(args.run / "result.json", {"captured": True, "no_solve": True,
                "matrix_shape": matrix.shape, "matrix_nonzeros": matrix.nnz,
                "matrix_csr_bytes": matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes,
                "full_nodes": len(points), "free_nodes": len(grouped_nodes),
                "complete_three_component_node_blocks": True,
                "archives_written": not args.metrics_only,
                "matrix_archive_sha256": None if args.metrics_only else workflow.digest(args.run / "matrix.npz"),
                "vector_archive_sha256": None if args.metrics_only else workflow.digest(args.run / "vectors.npz"),
                "scope": "First prescribed-normal-strain reduced system only; no displacement or stiffness result"})
            raise CapturedSystem()

    try:
        solve(points, cells, n_elm_gp=15, element_type=element_type, solver=Capture())
    except CapturedSystem:
        assert (args.run / "result.json").is_file()
    else:
        raise AssertionError("Capture must stop before any solution is returned")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--metrics-only", action="store_true", help="Measure assembly without saving large matrix archives")
    args = parser.parse_args()
    args.run, args.prepared = args.run.resolve(), args.prepared.resolve()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), Path(workflow.__file__), workflow.HERE / "cubic_reference.py",
        workflow.HERE / "stress_recovery.py", args.prepared / "mesh.npz", args.prepared / "result.json"]
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "platform": platform.platform(), "architecture": platform.machine(), "python": sys.version,
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY,
        "memory_enforcement": "rlimit-as-and-process-tree-rss" if sys.platform == "linux"
            else "process-tree-rss-watchdog", "rss_poll_interval_s": .25,
        "metrics_only": args.metrics_only,
        "worker_timeout_s": 600, "threads": 1, "material": {"young": 2.5, "poisson": .25},
        "load": [1, 0, 0, 0, 0, 0], "no_fea_solution": True})
    with Job(parts=["capture-reduced-system"]) as job:
        with worker_part(job, "capture-reduced-system") as part:
            command = [sys.executable, __file__, "--run", str(args.run), "--prepared", str(args.prepared), "--worker"]
            if args.metrics_only:
                command.append("--metrics-only")
            record = workflow.bounded(command, args.run / "worker.log", part, timeout=600)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("capture-reduced-system", seconds=record["wall_s"],
                peak_gib=record["peak_rss_bytes"] / 1024**3, free_nodes=row["free_nodes"])
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
