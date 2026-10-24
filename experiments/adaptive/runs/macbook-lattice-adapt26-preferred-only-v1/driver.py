"""Solve and audit one admitted fixed-domain nested Tet10 refinement level."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import microgen_workflow as workflow
from cubic_refinement import MODEL_REPORTS
from job_parts import worker_part
from nested_preflight import fixture_inputs

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    np, *_ = workflow.imports()
    from periodic import tensor_difference

    if args.worker == "solve":
        from amg_reference import AuditedAMG
        from stress_recovery import solve, vector_first_reactions

        row = json.loads((args.prepared / "result.json").read_text())
        assert row["status"] == "nested-preflight-passed" and row["qualified"]
        assert row["geometry"]["geometry_admitted"] and row["order"] == 2
        assert row["nodes"] <= 250000 and row["elements"] <= workflow.CAP
        mesh = np.load(args.run / "mesh.npz")
        row.update(status="solving", solver_quadrature_points=15,
            solver="cached-elastic-block-amg-cg", completed_loads=0, equilibrium_solve_performed=True,
            geometry_source=str(args.prepared), reference_promoted=False,
            scope="Fresh six-load solve on a saved nested mesh; no live meshing workflow timing")
        workflow.save(args.run / "result.json", row)
        def checkpoint(completed, partial, records):
            row.update(completed_loads=completed, partial_stiffness_columns=partial.tolist(), solver_records=records)
            workflow.save(args.run / "result.json", row)
        def save_field(completed, field):
            np.savez_compressed(args.run / f"load-{completed}.npz", displacement=field)
        started = time.perf_counter()
        fields, tensor, records = solve(mesh["points"], mesh["tetrahedra"], checkpoint,
            n_elm_gp=15, solver=AuditedAMG(elastic_blocks=True), field_checkpoint=save_field,
            reaction_reader=vector_first_reactions, release_reduced_cache=True)
        np.savez_compressed(args.run / "fields.npz", displacement=fields)
        row.pop("partial_stiffness_columns", None)
        row.update(status="solved-needs-reference", stiffness=tensor.tolist(), solver_records=records,
            saved_mesh_solve_s=time.perf_counter() - started)
        workflow.save(args.run / "result.json", row)
    elif args.worker == "audit":
        from quadrature_audit import worker as audit

        inputs = {str(args.run / name): workflow.digest(args.run / name)
            for name in ("result.json", "mesh.npz", "fields.npz")}
        workflow.save(args.run / "audit-inputs.json", inputs)
        audit(args.run, args.run)
        assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())
    else:
        from quadrature_audit import duffy_rule
        from tet10_integration import integrated_tensors

        inputs = json.loads((args.run / "audit-inputs.json").read_text())
        assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())
        row = json.loads((args.run / "result.json").read_text())
        assert row["completed_loads"] == 6 and row["qualified"]
        assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
        audit = json.loads((args.run / (args.run.name + ".json")).read_text())
        higher = audit["rules"][-1]
        assert higher["points"] == 125
        assert higher["energy_vs_solved_tensor"]["worst_relative_strain_energy_difference"] <= 1e-4
        assert higher["energy_vs_previous_rule"]["worst_relative_strain_energy_difference"] <= 1e-4
        mesh = np.load(args.run / "mesh.npz")
        prolonged = np.load(args.prepared / "prolongation.npz")["displacement"]
        solved = np.load(args.run / "fields.npz")["displacement"]
        _, coarse_energy, _ = integrated_tensors(mesh["points"], mesh["tetrahedra"], prolonged, duffy_rule(5))
        _, correction_energy, _ = integrated_tensors(mesh["points"], mesh["tetrahedra"], prolonged - solved, duffy_rule(5))
        fine_energy = np.asarray(higher["fixed_field_energy_tensor"])
        balance = tensor_difference(coarse_energy, fine_energy + correction_energy)
        assert balance["worst_relative_strain_energy_difference"] <= 1e-4
        assert np.linalg.eigvalsh(correction_energy).min() >= -1e-12
        parent = json.loads((Path(row["source"]) / "result.json").read_text())
        tensor = np.asarray(row["stiffness"])
        correction = tensor_difference(tensor, np.asarray(parent["stiffness"]))
        assert max(correction["relative_energy_difference_eigenvalues"]) <= 1e-4
        models = []
        for name in MODEL_REPORTS:
            path = workflow.HERE / "runs" / name
            report = json.loads(path.read_text())
            assert report["case"] == row["case"] and report["stiffness_target"] == .01
            inputs[str(path)] = workflow.digest(path)
            for index, model in enumerate(report["models"]):
                models.append({"source": str(path), "index": index,
                    "difference": tensor_difference(tensor, np.asarray(model["stiffness"]))})
        assert len(models) == 24
        worst = max(model["difference"]["worst_relative_strain_energy_difference"] for model in models)
        workflow.save(args.run / "assessment.json", {"inputs": inputs, "case": row["case"],
            "fresh_field_audit_passed": True, "positive_quadrature_sensitivity_passed": True,
            "parent_to_nested_difference": correction, "prolonged_parent_energy_tensor": coarse_energy.tolist(),
            "correction_energy_tensor": correction_energy.tolist(), "galerkin_energy_balance": balance,
            "models": models, "retained_model_count": len(models), "stiffness_target": .01,
            "worst_retained_model_disagreement": worst, "passes_target_all_models": worst <= .01,
            "reference_promoted": False, "scope": "Audited fixed-domain FE refinement correction and unchanged empirical model comparisons; no continuum convergence inferred from one correction"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--worker", choices=("solve", "audit", "assess"))
    args = parser.parse_args()
    args.run, args.prepared = args.run.resolve(), args.prepared.resolve()
    if args.worker:
        worker(args)
        return
    prepared = json.loads((args.prepared / "result.json").read_text())
    parent = Path(prepared["source"])
    assert prepared["status"] == "nested-preflight-passed" and prepared["qualified"]
    manifest = json.loads((args.prepared / "manifest.json").read_text())
    for name in ("result.json", "mesh.npz", "fields.npz"):
        assert workflow.digest(parent / name) == manifest["inputs"][str(parent / name)]
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    shutil.copyfile(args.prepared / "mesh.npz", args.run / "mesh.npz")
    fixtures = fixture_inputs(args.prepared)
    for path in fixtures:
        shutil.copyfile(path, args.run / path.name)
    paths = [args.prepared / name for name in ("result.json", "mesh.npz", "prolongation.npz", "manifest.json")]
    paths.extend(fixtures)
    paths.extend(parent / name for name in ("result.json", "mesh.npz", "fields.npz"))
    paths.extend(workflow.HERE / name for name in ("nested_solve.py", "stress_recovery.py", "amg_reference.py",
        "tet10_integration.py", "quadrature_audit.py", "microgen_workflow.py", "job_parts.py", "cubic_refinement.py",
        "nested_preflight.py"))
    paths.extend(workflow.HERE / "runs" / name for name in MODEL_REPORTS)
    inputs = {str(path): workflow.digest(path) for path in paths}
    import meshers._meshers as native
    (args.run / "requirements.txt").write_text(subprocess.check_output(
        ["uv", "pip", "freeze", "--python", sys.executable], text=True))
    workflow.save(args.run / "manifest.json", {"inputs": inputs, "memory_bytes": workflow.MEMORY,
        "platform": platform.platform(), "architecture": platform.machine(), "python": sys.version,
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY,
        "memory_enforcement": "rlimit-as-and-process-tree-rss" if sys.platform == "linux"
            else "process-tree-rss-watchdog", "rss_poll_interval_s": .25,
        "native_path": native.__file__, "native_sha256": workflow.digest(Path(native.__file__)),
        "versions": {name: importlib.metadata.version(name) for name in
            ("numpy", "scipy", "fedoo", "microgen", "mmgpy", "pyamg", "pyvista", "vtk")},
        "threads": 1, "timeouts_s": {"solve": 1200, "audit": 300, "assess": 300},
        "quadrature_sensitivity_limit": 1e-4, "galerkin_energy_balance_limit": 1e-4,
        "maximum_tetrahedra": workflow.CAP, "maximum_nodes": 250000, "stiffness_target": .01,
        "source_preflight": str(args.prepared), "no_live_workflow_timing": True})
    with Job(parts=["solve", "audit", "assess"]) as job:
        for stage in ("solve", "audit", "assess"):
            with worker_part(job, stage) as part:
                record = workflow.bounded([sys.executable, __file__, *sys.argv[1:], "--worker", stage],
                    args.run / (stage + ".log"), part, timeout=1200 if stage == "solve" else 300)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                values = {"seconds": record["wall_s"], "peak_gib": record["peak_rss_bytes"] / 1024**3}
                if stage == "assess":
                    assessed = json.loads((args.run / "assessment.json").read_text())
                    values.update(worst_disagreement=assessed["worst_retained_model_disagreement"])
                job.result(stage, **values)
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
