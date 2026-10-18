"""Validate an installed AMG preconditioner without changing benchmark defaults."""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import time

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def elastic_candidates(coordinates):
    import numpy as np

    candidates = np.zeros((3 * len(coordinates), 6))
    candidates.reshape(-1, 3, 6)[:, :, :3] = np.eye(3)
    candidates[0::3, 3] = -coordinates[:, 1]
    candidates[1::3, 3] = coordinates[:, 0]
    candidates[0::3, 4] = -coordinates[:, 2]
    candidates[2::3, 4] = coordinates[:, 0]
    candidates[1::3, 5] = -coordinates[:, 2]
    candidates[2::3, 5] = coordinates[:, 1]
    return candidates


class AuditedAMG:
    def __init__(self, elastic_blocks=False):
        self.records = []
        self.problem = None
        self.hierarchy = None
        self.key = None
        self.elastic_blocks = elastic_blocks
        self.permutation = None

    def __call__(self, matrix, rhs, **kwargs):
        import numpy as np
        import pyamg
        from scipy.sparse.linalg import cg

        started = time.perf_counter()
        matrix = matrix.tocsr()
        assert (matrix.diagonal() > 0).all()
        skew = matrix - matrix.T
        assert np.linalg.norm(skew.data) / np.linalg.norm(matrix.data) < 1e-12
        del skew
        digest = hashlib.sha256()
        for array in (matrix.indptr, matrix.indices, matrix.data):
            digest.update(memoryview(array))
        key = (matrix.shape, digest.hexdigest())
        setup_seconds = 0.
        reorder_seconds = 0.
        if self.hierarchy is None:
            free = self.problem._dof_free
            count = self.problem.mesh.n_nodes
            ranks = self.problem.space.get_rank_vector("Disp")
            assert len(free) == len(rhs) and np.all(free < self.problem.n_node_dof)
            candidates = np.column_stack([free // count == rank for rank in ranks]).astype(float)
            assert np.all(candidates.sum(axis=1) == 1)
            preconditioner_matrix = matrix
            if self.elastic_blocks:
                phase = time.perf_counter()
                components, nodes = free // count, free % count
                self.permutation = np.lexsort((components, nodes))
                grouped = nodes[self.permutation].reshape(-1, 3)
                np.testing.assert_array_equal(grouped, np.repeat(grouped[:, :1], 3, axis=1))
                np.testing.assert_array_equal(components[self.permutation].reshape(-1, 3),
                    np.tile(ranks, (len(grouped), 1)))
                preconditioner_matrix = matrix[self.permutation][:, self.permutation].tobsr(blocksize=(3, 3))
                candidates = elastic_candidates(self.problem.mesh.nodes[grouped[:, 0]])
                reorder_seconds = time.perf_counter() - phase
                np.random.seed(20261018)
            phase = time.perf_counter()
            self.hierarchy = pyamg.smoothed_aggregation_solver(preconditioner_matrix,
                B=candidates, symmetry="symmetric")
            setup_seconds = time.perf_counter() - phase
            self.key = key
        else:
            assert key == self.key, "Cannot reuse hierarchy for a changed reduced system"
        iterations = 0

        def count_iteration(_):
            nonlocal iterations
            iterations += 1

        solve_matrix = self.hierarchy.levels[0].A if self.elastic_blocks else matrix
        solve_rhs = rhs[self.permutation] if self.elastic_blocks else rhs
        solution, info = cg(solve_matrix, solve_rhs, M=self.hierarchy.aspreconditioner(cycle="V"),
            rtol=1e-12, atol=1e-14, maxiter=20000, callback=count_iteration)
        if self.elastic_blocks:
            original_solution = np.empty_like(solution)
            original_solution[self.permutation] = solution
            solution = original_solution
        residual = float(np.linalg.norm(matrix @ solution - rhs) / max(np.linalg.norm(rhs), 1e-30))
        assert info == 0 and residual < 2e-9, (info, residual)
        self.records.append({"iterations": iterations, "relative_residual": residual,
            "free_dofs": len(rhs), "seconds": time.perf_counter() - started,
            "hierarchy_setup_s": setup_seconds, "hierarchy_reused": bool(len(self.records)),
            "node_reordering_s": reorder_seconds, "elastic_blocks": self.elastic_blocks,
            "levels": len(self.hierarchy.levels), "operator_complexity": self.hierarchy.operator_complexity(),
            "reduced_matrix_sha256": key[1]})
        return solution


def worker(args):
    np, _, meshers, *_ = workflow.imports()
    from stress_recovery import solve, periodic_displacement_jump_error
    from tet10_integration import integrated_tensors
    from periodic import tensor_difference
    from homogenization import MATERIAL
    from cubic_reference import register_cubic, enrich
    from quadrature_audit import duffy_rule

    if args.case == "cube":
        parent = meshers.generate(lambda coordinates_x, *_: coordinates_x * 0 - 1,
            bounds=(-.5, .5) * 3, cells=4, order=2, geometry_tolerance=1,
            minimum_quality=.1, periodic=(True, True, True), threads=1)
        element_class = register_cubic()
        points, cells, surface = enrich(parent.points, parent.tetrahedra, parent.surface, element_class())
        expected, element_type = MATERIAL, element_class.name
    else:
        archive = np.load(args.source / "mesh.npz")
        points, cells, surface = (archive[name] for name in ("points", "tetrahedra", "surface"))
        expected = np.array(json.loads((args.source / "result.json").read_text())["stiffness"])
        element_type = None
    solver = AuditedAMG(elastic_blocks=args.elastic_blocks)
    started = time.perf_counter()
    fields, tensor, records = solve(points, cells, n_elm_gp=15, element_type=element_type, solver=solver)
    seconds = time.perf_counter() - started
    stress, energy, volume = integrated_tensors(points, cells, fields, element_type=element_type)
    own_error = max(float(np.linalg.norm(value - tensor) / np.linalg.norm(tensor)) for value in (stress, energy))
    jump = periodic_displacement_jump_error(points, surface, fields)
    difference = tensor_difference(tensor, expected)
    assert max(own_error, jump, difference["worst_relative_strain_energy_difference"]) < 1e-7
    assert len(records) == 6 and sum(record["hierarchy_reused"] for record in records) == 5
    _, higher, _ = integrated_tensors(points, cells, fields, duffy_rule(5), element_type=element_type)
    workflow.save(args.run / (args.case + ".json"), {"case": args.case, "passed": True,
        "nodes": len(points), "elements": len(cells), "solve_s": seconds, "difference": difference,
        "own_rule_relative_error": own_error, "periodic_jump_max": jump, "integrated_volume": volume,
        "higher_rule_fixed_field_difference": tensor_difference(higher, tensor), "solver_records": records,
        "scope": "Solver validation only; same equations, unchanged tolerances; not a live workflow timing"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=workflow.HERE / "runs/g3-schwarz-curved8-pilot-v1")
    parser.add_argument("--case", choices=("cube", "schwarz"))
    parser.add_argument("--elastic-blocks", action="store_true")
    args = parser.parse_args()
    if args.case:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), Path(workflow.__file__), workflow.HERE / "stress_recovery.py",
        workflow.HERE / "cubic_reference.py", workflow.HERE / "tet10_integration.py",
        args.source / "mesh.npz", args.source / "result.json"]
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": 300, "threads": 1, "pyamg": importlib.metadata.version("pyamg"),
        "elastic_blocks": args.elastic_blocks,
        "method": ("SA on node-major BSR3 with six elastic candidates and V-cycle CG" if args.elastic_blocks else
            "SA on component-major CSR with three translational candidates and V-cycle CG"),
        "cg_rtol": 1e-12, "cg_atol": 1e-14})
    with Job(parts=["cube", "schwarz"]) as job:
        for case in ("cube", "schwarz"):
            with job.part(case) as part:
                command = [sys.executable, __file__, "--run", str(args.run), "--source", str(args.source), "--case", case]
                if args.elastic_blocks:
                    command.append("--elastic-blocks")
                record = workflow.bounded(command, args.run / (case + ".log"), part, timeout=300)
                workflow.save(args.run / (case + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {case}.log"
                row = json.loads((args.run / (case + ".json")).read_text())
                job.result(case, solve_s=row["solve_s"], peak_gib=record["peak_rss_bytes"] / 1024**3,
                    iterations=sum(record["iterations"] for record in row["solver_records"]))
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
