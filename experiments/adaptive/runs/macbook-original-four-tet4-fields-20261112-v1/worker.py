"""Audit six-load stiffness of two admitted saved Tet4 meshes, without a reference claim."""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", RAYON_NUM_THREADS="1")
import microgen_workflow as w
from job_parts import worker_part
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(source, run):
    np, *_ = w.imports()
    from amg_reference import AuditedAMG
    from stress_recovery import solve, vector_first_reactions
    from tet10_integration import audit_saved, integrated_tensors
    from quadrature_audit import duffy_rule
    from periodic import periodic_pairs, tensor_difference

    original = json.loads((source / "result.json").read_text())
    assert original["qualified"] and original["geometry"]["geometry_admitted"]
    assert original["case"] in w.CASES and original["periodic_required"]
    assert not original["geometry"]["failed_gates"]
    assert w.digest(source / "mesh.npz") == original["mesh_sha256"]
    run.mkdir(exist_ok=False)
    shutil.copyfile(source / "mesh.npz", run / "mesh.npz")
    mesh = np.load(run / "mesh.npz")
    points, cells, surface = (mesh[k] for k in ("points", "tetrahedra", "surface"))
    assert cells.shape[1] == 4 and max(len(points), len(cells)) <= 250000
    periodic_pairs(points, surface)
    row = dict(original, source=str(source), completed_loads=0, solver_quadrature_points=4,
               reference_promoted=False, equilibrium_solve_performed=True,
               scope="Saved Tet4 FE consistency and pairwise stiffness only; no continuum reference or live generation timing")
    def checkpoint(completed, partial, records):
        row.update(completed_loads=completed, partial_stiffness_columns=partial.tolist(), solver_records=records)
        w.save(run / "result.json", row)
    started = time.perf_counter()
    fields, tensor, records = solve(points, cells, checkpoint=checkpoint, n_elm_gp=4,
        solver=AuditedAMG(elastic_blocks=True), reaction_reader=vector_first_reactions,
        release_reduced_cache=True)
    row.pop("partial_stiffness_columns", None)
    row.update(stiffness=tensor.tolist(), solver_records=records,
               saved_mesh_solve_s=time.perf_counter() - started, status="solved-needs-reference")
    assert len(records) == 6 and max(r["relative_residual"] for r in records) < 2e-9
    assert np.linalg.eigvalsh(tensor).min() > 0
    np.savez_compressed(run / "fields.npz", displacement=fields)
    w.save(run / "result.json", row)
    audit_saved(run, run / "field-audit.json")
    stress, energy, volume = integrated_tensors(points, cells, fields, duffy_rule(2))
    positive = dict(points=8, positive_weights=True,
        mixed_load_energy_difference=tensor_difference(energy, tensor),
        mixed_load_stress_difference=tensor_difference((stress + stress.T) / 2, tensor),
        relative_volume_difference=float(volume / row["geometry"]["volume"] - 1))
    assert max(positive[k]["worst_relative_strain_energy_difference"] for k in
               ("mixed_load_energy_difference", "mixed_load_stress_difference")) < 1e-7
    assert abs(positive["relative_volume_difference"]) < 1e-10
    assert w.digest(run / "mesh.npz") == original["mesh_sha256"]
    w.save(run / "positive-quadrature-audit.json", positive)
    row.update(field_audit_passed=True, positive_quadrature_passed=True,
               fields_sha256=w.digest(run / "fields.npz"), assessment="FE-consistent; reference unassessed")
    w.save(run / "result.json", row)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--sources", type=Path, nargs=2)
    p.add_argument("--worker", type=Path)
    a = p.parse_args()
    if a.worker:
        worker(a.worker, a.run)
        return
    assert a.sources
    a.run.mkdir(parents=True, exist_ok=False)
    names = ("uniform", "adaptive")
    modules = ["microgen_workflow.py", "job_parts.py", "amg_reference.py", "stress_recovery.py",
               "tet10_integration.py", "quadrature_audit.py", "scorecard.py"]
    shutil.copyfile(__file__, a.run / "driver.py")
    snapshots = a.run / "audit-source"
    snapshots.mkdir()
    source_paths = [w.HERE / name for name in modules] + [w.ROOT.parent / "periodic-microgen/experiments/periodic-fea" / name for name in ("homogenization.py", "periodic.py")]
    for source in source_paths:
        shutil.copyfile(source, snapshots / source.name)
    inputs = {str(s / n): w.digest(s / n) for s in a.sources for n in ("result.json", "mesh.npz")}
    fixtures = [json.loads((s / "result.json").read_text())["fixture_hashes"] for s in a.sources]
    assert fixtures[0] == fixtures[1], "The physical specimen must be identical"
    manifest = dict(inputs=inputs, fixture_hashes=fixtures[0], source_sha256={str(p): w.digest(p) for p in [Path(__file__).resolve(), *source_paths]},
                    timeout_s=600, memory_bytes=12 * 1024**3, numerical_threads=1,
                    new_fea=True, reference_promoted=False, live_workflow_timing=False,
                    material=dict(E=2.5, nu=0.25), loads=6)
    w.save(a.run / "manifest.json", manifest)
    rows = []
    with Job(parts=list(names)) as job:
        for name, source in zip(names, a.sources, strict=True):
            with worker_part(job, name, memory_bytes=12 * 1024**3) as part:
                resource = w.bounded([sys.executable, str(Path(__file__).resolve()), "--run",
                    str((a.run / name).resolve()), "--worker", str(source.resolve())],
                    a.run / (name + ".log"), part, timeout=600, memory_bytes=12 * 1024**3)
                w.save(a.run / (name + "-resource.json"), resource)
                assert resource["returncode"] == 0, (name, "Read the saved worker log")
                row = json.loads((a.run / name / "result.json").read_text())
                assert row["field_audit_passed"] and row["positive_quadrature_passed"]
                rows.append(row)
                job.result(name, loads=row["completed_loads"], seconds=resource["wall_s"],
                           peak_gib=resource["peak_rss_bytes"] / 1024**3)
    assert all(w.digest(Path(path)) == digest for path, digest in inputs.items())
    assert all(w.digest(Path(p)) == digest for p, digest in manifest["source_sha256"].items())
    w.imports()
    from periodic import tensor_difference
    import numpy as np
    comparison = tensor_difference(np.asarray(rows[1]["stiffness"]), np.asarray(rows[0]["stiffness"]))
    w.save(a.run / "assessment.json", dict(complete=True, cases=names, pairwise_mixed_load_difference=comparison,
        source_inputs=inputs, new_fea=True, reference_promoted=False, reference_accuracy_assessed=False,
        scope="Different admitted polyhedral domains and non-nested spaces; pairwise agreement is not convergence, accuracy equivalence or a continuum bound"))


if __name__ == "__main__":
    main()
