"""Fresh-process saved cubic field audit and retained-model comparison."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    np, *_ = workflow.imports()
    from cubic_reference import register_cubic, basis_checks
    from homogenization import MATERIAL, strain_tensor
    from periodic import tensor_difference
    from quadrature_audit import duffy_rule
    from stress_recovery import periodic_displacement_jump_error
    from tet10_integration import integrated_tensors

    row = json.loads((args.source / "result.json").read_text())
    assert row["qualified"] and row["geometry"]["geometry_admitted"]
    assert row["order"] == 3 and row["completed_loads"] == 6
    assert len(row["solver_records"]) == 6
    assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
    assert row["quadrature_sensitivity_passed"] and row["quadrature_sensitivity_limit"] <= 1e-4
    mesh = np.load(args.prepared / "mesh.npz")
    fields = np.load(args.source / "fields.npz")["displacement"]
    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    assert cells.shape == (row["elements"], 20) and len(cells) <= workflow.CAP
    assert len(points) == row["nodes"] <= 250000
    assert fields.shape == (6, len(points), 3) and np.isfinite(fields).all()
    element_class = register_cubic()
    element = element_class()
    checks = basis_checks(element)
    expected = np.asarray(row["stiffness"])
    rules = []
    previous = None
    for count in (None, 4, 5):
        quadrature = None if count is None else duffy_rule(count)
        affine = np.array([element.xi_nd @ strain_tensor(load).T for load in np.eye(6)])
        patch_stress, patch_energy, patch_volume = integrated_tensors(
            element.xi_nd, np.arange(20)[None], affine, quadrature, element_type=element_class.name)
        np.testing.assert_allclose(patch_stress, MATERIAL / 6, atol=1e-13, rtol=1e-13)
        np.testing.assert_allclose(patch_energy, MATERIAL / 6, atol=1e-13, rtol=1e-13)
        assert abs(patch_volume - 1 / 6) < 1e-14
        stress, energy, volume = integrated_tensors(points, cells, fields, quadrature,
            element_type=element_class.name)
        error = float(np.linalg.norm(stress - expected) / np.linalg.norm(expected))
        skew = float(np.linalg.norm(stress - stress.T) / np.linalg.norm(stress))
        energy_difference = tensor_difference(energy, expected)
        if count is None:
            assert max(error, skew, energy_difference["worst_relative_strain_energy_difference"]) < 1e-7
        assert abs(volume / row["geometry"]["volume"] - 1) < 1e-10
        rules.append({"points": 15 if count is None else count**3,
            "reaction_relative_error": error, "stress_skew_relative": skew,
            "fixed_field_energy_difference": energy_difference,
            "energy_vs_previous_rule": None if previous is None else tensor_difference(energy, previous),
            "volume": volume, "sampled_jacobians_positive": True})
        previous = energy
    assert rules[-1]["fixed_field_energy_difference"]["worst_relative_strain_energy_difference"] <= 1e-4
    assert rules[-1]["energy_vs_previous_rule"]["worst_relative_strain_energy_difference"] <= 1e-4
    jump = periodic_displacement_jump_error(points, surface, fields)
    assert jump < 1e-9
    models = []
    for path in args.models:
        report = json.loads(path.read_text())
        assert report["case"] == row["case"]
        for index, model in enumerate(report["models"]):
            models.append({"source": str(path), "index": index,
                "difference_from_cubic": tensor_difference(np.asarray(model["stiffness"]), expected)})
    assert len(models) == 24
    candidates = []
    for path in args.candidates:
        candidate = json.loads(path.read_text())
        assert candidate["qualified"] and candidate["geometry"]["geometry_admitted"]
        assert candidate["case"] == row["case"] and candidate["elements"] <= workflow.CAP
        tensor = np.asarray(candidate["stiffness"])
        retained = [tensor_difference(tensor, np.asarray(model["stiffness"]))[
            "worst_relative_strain_energy_difference"] for source in args.models
            for model in json.loads(source.read_text())["models"]]
        candidates.append({"source": str(path), "nodes": candidate["nodes"],
            "difference_from_cubic": tensor_difference(tensor, expected),
            "retained_model_difference_range": [min(retained), max(retained)],
            "passes_one_percent_all_retained_models": max(retained) <= .01})
    workflow.save(args.run / "result.json", {"case": row["case"], "source": str(args.source),
        "basis_checks": checks, "fresh_process_audit_passed": True, "rules": rules,
        "periodic_jump_max": jump, "retained_model_count": len(models),
        "models": models, "candidates": candidates,
        "scope": "Fixed-field cubic consistency and empirical model disagreements; no continuum certificate or live timing claim"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--models", nargs="+", type=Path, required=True)
    parser.add_argument("--candidates", nargs="+", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), args.source / "result.json", args.source / "fields.npz",
        args.prepared / "mesh.npz", *args.models, *args.candidates]
    paths.extend(workflow.HERE / name for name in ("cubic_reference.py", "microgen_workflow.py",
        "tet10_integration.py", "stress_recovery.py", "quadrature_audit.py"))
    inputs = {str(path): workflow.digest(path) for path in paths}
    source_manifest = json.loads((args.source / "manifest.json").read_text())["sources"]
    assert source_manifest[str(args.prepared / "mesh.npz")] == inputs[str(args.prepared / "mesh.npz")]
    workflow.save(args.run / "manifest.json", {"inputs": inputs, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": 300, "threads": 1, "quadrature_points": [15, 64, 125]})
    with Job(parts=["audit-and-compare"]) as job:
        with job.part("audit-and-compare") as part:
            record = workflow.bounded([sys.executable, __file__, *sys.argv[1:], "--worker"],
                args.run / "worker.log", part, timeout=300)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            job.result("audit-and-compare", seconds=record["wall_s"],
                peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == value for path, value in inputs.items())


if __name__ == "__main__":
    main()
