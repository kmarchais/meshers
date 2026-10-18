"""Measure recovery localization on an immutable, audited Tet4 solution."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(source, run):
    np, *_ = workflow.imports()
    from sizing_pilot import size_from_recovery
    from tet10_integration import audit_saved

    audit_saved(source, run / "field-consistency.json")
    row = json.loads((source / "result.json").read_text())
    archive = np.load(source / "mesh.npz")
    mesh = SimpleNamespace(**{name: archive[name] for name in ("points", "tetrahedra", "surface")})
    assert mesh.tetrahedra.shape[1] == 4
    fields = np.load(source / "fields.npz")["displacement"]
    _, targets, mass = size_from_recovery(mesh, fields, np.asarray(row["stiffness"]), 0.065)
    volumes = np.linalg.det(mesh.points[mesh.tetrahedra[:, 1:]] - mesh.points[mesh.tetrahedra[:, :1]]) / 6
    assert (volumes > 0).all() and np.isfinite(mass).all() and (mass >= 0).all() and mass.sum() > 0
    fractions = []
    for fraction in (0.05, 0.1, 0.2, 0.5):
        selected = np.argsort(mass / volumes)[-int(np.ceil(fraction * len(mass))):]
        fractions.append({"element_fraction": fraction, "volume_fraction": float(volumes[selected].sum() / volumes.sum()),
            "recovery_mass_fraction": float(mass[selected].sum() / mass.sum())})
    workflow.save(run / "result.json", {"source": str(source), "elements": len(mass),
        "recovery_trace_indicator": float(mass.sum()), "density_ranked_capture": fractions,
        "target_length_quantiles": np.quantile(targets, [0, 0.1, 0.5, 0.9, 1]).tolist(),
        "scope": "Recovery localization only, not observed error, stiffness error bounds or local stress certification",
        "timing_scope": "Saved-field diagnostic, not a live adaptive workflow"})
    np.savez_compressed(run / "indicators.npz", mass=mass, volumes=volumes, target_length=targets)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.source, args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    files = [args.source / name for name in ("result.json", "mesh.npz", "fields.npz")]
    files += [Path(__file__), workflow.HERE / "sizing_pilot.py", workflow.HERE / "stress_recovery.py",
        workflow.HERE / "tet10_integration.py", Path(workflow.__file__)]
    hashes = {str(path): workflow.digest(path) for path in files}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 300})
    with Job(parts=["recovery"]) as job:
        with job.part("recovery") as part:
            record = workflow.bounded([sys.executable, __file__, "--worker", "--source", str(args.source),
                "--run", str(args.run)], args.run / "worker.log", part, timeout=300)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("recovery", seconds=record["wall_s"], elements=row["elements"],
                top_20_capture=row["density_ranked_capture"][2]["recovery_mass_fraction"],
                peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
