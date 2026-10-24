"""Validate an installed AMG preconditioner without changing benchmark defaults."""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time

import microgen_workflow as workflow
from job_parts import worker_part

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def reduced_matrix_key(matrix):
    digest = hashlib.sha256()
    for array in (matrix.indptr, matrix.indices, matrix.data):
        digest.update(memoryview(array))
    return matrix.shape, digest.hexdigest()


def verify_symmetric(matrix, block_rows=32768):
    import numpy as np

    assert matrix.shape[0] == matrix.shape[1] and block_rows > 0
    denominator = np.linalg.norm(matrix.data)
    assert np.isfinite(denominator) and denominator > 0
    squared_skew = 0.
    for first in range(0, matrix.shape[0], block_rows):
        last = min(first + block_rows, matrix.shape[0])
        skew = matrix[first:last] - matrix[:, first:last].T
        squared_skew += float(skew.data @ skew.data)
    assert np.sqrt(squared_skew) / denominator < 1e-12


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
        self.symmetry_key = None

    def __call__(self, matrix, rhs, **kwargs):
        import numpy as np
        import pyamg
        from scipy.sparse.linalg import cg

        started = time.perf_counter()
        matrix = matrix.tocsr()
        assert (matrix.diagonal() > 0).all()
        key = reduced_matrix_key(matrix)
        symmetry_reused = self.symmetry_key is not None
        if self.symmetry_key is None:
            verify_symmetric(matrix)
            self.symmetry_key = key
        else:
            assert self.symmetry_key == key, "Changed matrix cannot reuse its symmetry certificate"
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
            "verified_symmetry_reused": symmetry_reused,
            "levels": len(self.hierarchy.levels), "operator_complexity": self.hierarchy.operator_complexity(),
            "reduced_matrix_sha256": key[1]})
        return solution


class VerifiedPrefixAMG:
    def __init__(self, fields, elastic_blocks=True):
        self.fields = fields
        self.delegate = AuditedAMG(elastic_blocks=elastic_blocks)
        self.problem = None
        self.records = []
        self.key = None

    @property
    def problem(self):
        return self.delegate.problem

    @problem.setter
    def problem(self, problem):
        self.delegate.problem = problem

    def __call__(self, matrix, rhs, **kwargs):
        import numpy as np

        if len(self.records) >= len(self.fields):
            solution = self.delegate(matrix, rhs, **kwargs)
            assert self.key is None or self.key == self.delegate.key
            self.records.append({**self.delegate.records[-1], "replayed": False})
            return solution
        started = time.perf_counter()
        matrix = matrix.tocsr()
        key = reduced_matrix_key(matrix)
        if self.key is None:
            verify_symmetric(matrix)
            self.key = key
            self.delegate.symmetry_key = key
        else:
            assert key == self.key
        field = np.asarray(self.fields[len(self.records)])
        assert field.shape == (self.problem.mesh.n_nodes, 3) and np.isfinite(field).all()
        assert self.problem.n_node_dof == field.size and np.all(self.problem._dof_free < field.size)
        solution = field.T.ravel()[self.problem._dof_free].copy()
        residual = float(np.linalg.norm(matrix @ solution - rhs) / max(np.linalg.norm(rhs), 1e-30))
        assert residual < 2e-9, "Saved prefix does not solve the current reduced equations"
        self.records.append({"iterations": 0, "relative_residual": residual, "free_dofs": len(rhs),
            "seconds": time.perf_counter() - started, "replayed": True,
            "hierarchy_reused": False, "reduced_matrix_sha256": key[1]})
        return solution


def worker(args):
    np, _, meshers, *_ = workflow.imports()
    from stress_recovery import solve, periodic_displacement_jump_error, vector_first_reactions
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
    prefix = []
    if args.replay_prefix:
        if args.case == "cube":
            baseline_fields, _, _ = solve(points, cells, n_elm_gp=15, element_type=element_type)
            prefix = list(baseline_fields[:args.prefix_count])
        else:
            prefix = list(np.load(args.source / "fields.npz")["displacement"][:args.prefix_count])
    solver = (VerifiedPrefixAMG(prefix, elastic_blocks=args.elastic_blocks) if prefix else
        AuditedAMG(elastic_blocks=args.elastic_blocks))
    reaction_errors = []

    def checked_reactions(problem):
        native = problem.get_ext_forces("MeanStrain")
        vector_first = vector_first_reactions(problem)
        error = float(np.linalg.norm(vector_first - native) / max(np.linalg.norm(native), 1e-30))
        assert error < 1e-9
        reaction_errors.append(error)
        return vector_first

    def checked_field(load, displacement):
        if load <= len(prefix):
            np.testing.assert_allclose(displacement, prefix[load - 1], atol=1e-11, rtol=0)

    started = time.perf_counter()
    fields, tensor, records = solve(points, cells, n_elm_gp=15, element_type=element_type, solver=solver,
        reaction_reader=checked_reactions, field_checkpoint=checked_field,
        release_reduced_cache=args.release_reduced_cache)
    seconds = time.perf_counter() - started
    stress, energy, volume = integrated_tensors(points, cells, fields, element_type=element_type)
    own_error = max(float(np.linalg.norm(value - tensor) / np.linalg.norm(tensor)) for value in (stress, energy))
    jump = periodic_displacement_jump_error(points, surface, fields)
    difference = tensor_difference(tensor, expected)
    assert max(own_error, jump, difference["worst_relative_strain_energy_difference"]) < 1e-7
    assert len(records) == 6 and sum(record["hierarchy_reused"] for record in records) == 5 - len(prefix)
    assert sum(record.get("replayed", False) for record in records) == len(prefix)
    _, higher, _ = integrated_tensors(points, cells, fields, duffy_rule(5), element_type=element_type)
    workflow.save(args.run / (args.case + ".json"), {"case": args.case, "passed": True,
        "nodes": len(points), "elements": len(cells), "solve_s": seconds, "difference": difference,
        "own_rule_relative_error": own_error, "periodic_jump_max": jump, "integrated_volume": volume,
        "higher_rule_fixed_field_difference": tensor_difference(higher, tensor), "solver_records": records,
        "replayed_prefix_loads": len(prefix), "maximum_reaction_reassociation_difference": max(reaction_errors),
        "scope": "Solver validation only; same equations, unchanged tolerances; not a live workflow timing"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=workflow.HERE / "runs/g3-schwarz-curved8-pilot-v1")
    parser.add_argument("--case", choices=("cube", "schwarz"))
    parser.add_argument("--elastic-blocks", action="store_true")
    parser.add_argument("--replay-prefix", action="store_true")
    parser.add_argument("--prefix-count", type=int, choices=range(1, 6), default=2)
    parser.add_argument("--release-reduced-cache", action="store_true")
    args = parser.parse_args()
    if args.case:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), Path(workflow.__file__), workflow.HERE / "stress_recovery.py",
        workflow.HERE / "cubic_reference.py", workflow.HERE / "tet10_integration.py",
        args.source / "mesh.npz", args.source / "result.json"]
    if args.replay_prefix:
        paths.append(args.source / "fields.npz")
    hashes = {str(path): workflow.digest(path) for path in paths}
    import meshers._meshers as native
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "platform": platform.platform(), "architecture": platform.machine(), "python": sys.version,
        "native_sha256": workflow.digest(Path(native.__file__)),
        "memory_enforcement": "rlimit-as-and-process-tree-rss" if sys.platform == "linux"
            else "process-tree-rss-watchdog", "rss_poll_interval_s": .25,
        "worker_timeout_s": 300, "threads": 1, "pyamg": importlib.metadata.version("pyamg"),
        "elastic_blocks": args.elastic_blocks,
        "replayed_prefix_loads": args.prefix_count if args.replay_prefix else 0,
        "release_reduced_cache": args.release_reduced_cache,
        "method": ("SA on node-major BSR3 with six elastic candidates and V-cycle CG" if args.elastic_blocks else
            "SA on component-major CSR with three translational candidates and V-cycle CG"),
        "cg_rtol": 1e-12, "cg_atol": 1e-14})
    with Job(parts=["cube", "schwarz"]) as job:
        for case in ("cube", "schwarz"):
            with worker_part(job, case) as part:
                command = [sys.executable, __file__, "--run", str(args.run), "--source", str(args.source), "--case", case]
                if args.elastic_blocks:
                    command.append("--elastic-blocks")
                if args.replay_prefix:
                    command.append("--replay-prefix")
                    command.extend(["--prefix-count", str(args.prefix_count)])
                if args.release_reduced_cache:
                    command.append("--release-reduced-cache")
                record = workflow.bounded(command, args.run / (case + ".log"), part, timeout=300)
                workflow.save(args.run / (case + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {case}.log"
                row = json.loads((args.run / (case + ".json")).read_text())
                job.result(case, solve_s=row["solve_s"], peak_gib=record["peak_rss_bytes"] / 1024**3,
                    iterations=sum(record["iterations"] for record in row["solver_records"]))
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
