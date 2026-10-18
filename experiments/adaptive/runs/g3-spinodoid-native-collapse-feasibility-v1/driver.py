"""Exercise native collapse predicates on the saved rejected spinodoid."""

import json
import os
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def main():
    import numpy as np
    run = Path(__file__).resolve().parent
    source = workflow.HERE / "runs/g3-spinodoid-coarse16-opt4-v1"
    stars = workflow.HERE / "runs/g3-spinodoid-coarse16-projection-stars-v1/result.json"
    probe = ROOT / "crates/meshers-core/src/tet_adapt/collapse_probe.rs"
    native = ROOT / "crates/meshers-core/src/tet_adapt.rs"
    assert not (run / "manifest.json").exists()
    paths = [source / "mesh.npz", source / "fourier-coefficients.json", stars, probe, native, Path(__file__)]
    hashes = {str(path): workflow.digest(path) for path in paths}
    shutil.copyfile(probe, run / "probe.rs")
    shutil.copyfile(native, run / "tet_adapt-snapshot.rs")
    mesh = np.load(source / "mesh.npz")
    coefficients = json.loads((source / "fourier-coefficients.json").read_text())
    payload = {name: mesh[name].tolist() for name in ("points", "tetrahedra", "surface")}
    payload.update(coefficients)
    payload["origins"] = [row["origin"] for row in json.loads(stars.read_text())["vertices"]]
    workflow.save(run / "input.json", payload)
    os.environ.update(RUSTUP_TOOLCHAIN="1.95.0", CARGO_BUILD_JOBS="2", OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1", RAYON_NUM_THREADS="1", MESHERS_COLLAPSE_INPUT=str(run / "input.json"),
        MESHERS_COLLAPSE_OUTPUT=str(run / "result.json"),
        PATH=str(Path.home() / ".cargo/bin") + os.pathsep + os.environ["PATH"])
    workflow.save(run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "timeout_s": 600, "build_threads": 2, "numeric_threads": 1,
        "temporary_test_only_probe": True, "production_extension_rebuilt": False})
    with Job(parts=["native-collapse-predicate"]) as job:
        with job.part("native-collapse-predicate") as part:
            command = ["cargo", "test", "--release", "--locked", "-p", "meshers-core", "--lib",
                "inspect_saved_pending_collapse", "--", "--ignored", "--test-threads=1"]
            record = workflow.bounded(command, run / "worker.log", part, timeout=600)
            workflow.save(run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            rows = json.loads((run / "result.json").read_text())["rows"]
            job.result("native-collapse-predicate", seconds=record["wall_s"],
                accepted=sum(len(row["native_accepted_collapses"]) for row in rows))
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
