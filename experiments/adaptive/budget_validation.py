"""Rebuild and validate the adaptive tetrahedron-budget failure on G3."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
        RAYON_NUM_THREADS="1", CARGO_BUILD_JOBS="2", RUSTUP_TOOLCHAIN="1.95.0",
        VIRTUAL_ENV=str(workflow.ROOT / ".venv"),
        PATH=str(Path.home() / ".cargo/bin") + os.pathsep + os.environ["PATH"])
    paths = [workflow.ROOT / name for name in (
        "crates/meshers-core/src/tet_adapt.rs", "crates/meshers-core/src/implicit.rs",
        "crates/meshers-python/tests/test_size_field.py", "experiments/adaptive/sizing_pilot.py",
        "experiments/adaptive/stress_recovery.py", "experiments/adaptive/budget_validation.py")]
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"source_sha256": hashes,
        "memory_bytes": workflow.MEMORY, "build_threads": 2, "numeric_threads": 1,
        "tetrahedra_cap": workflow.CAP, "rust": "1.95.0"})
    replay = args.run / "replay"
    replay.mkdir()
    tasks = [
        ("build", [str(workflow.ROOT / ".venv/bin/maturin"), "develop", "--release",
            "--locked", "--manifest-path", "crates/meshers-python/Cargo.toml",
            "--features", "experimental-surfaces"], 1200),
        ("core-tests", ["cargo", "test", "--release", "--locked", "-p", "meshers-core",
            "--lib", "--", "--test-threads=1"], 900),
        ("python-tests", [sys.executable, "-m", "pytest", "-q",
            "crates/meshers-python/tests/test_size_field.py",
            "crates/meshers-python/tests/test_api.py",
            "crates/meshers-python/tests/test_adaptive.py",
            "crates/meshers-python/tests/test_fedoo.py",
            "experiments/adaptive/test_reference_admission.py"], 600),
        ("replay", [sys.executable, "experiments/adaptive/sizing_pilot.py", "--worker", "uniform",
            "--run", str(replay), "--case", "smooth-lattice", "--cells", "24",
            "--coarse-cells", "24", "--min-length", ".002", "--nominal",
            "0.02247651505130419"], 180)]
    workers = {}
    with Job(parts=[name for name, _, _ in tasks]) as job:
        for name, command, timeout in tasks:
            with job.part(name) as part:
                workers[name] = workflow.bounded(command, args.run / (name + ".log"), part,
                    timeout=timeout)
                workflow.save(args.run / "workers.json", workers)
                assert workers[name]["returncode"] == 0, f"Read {name}.log"
                job.result(name, seconds=workers[name]["wall_s"],
                    peak_gib=workers[name]["peak_rss_bytes"] / 1024**3)
    row = json.loads((replay / "uniform.json").read_text())
    assert row["status"] == "failed" and not row["qualified"]
    assert "adaptation exceeded max_tetrahedra" in row["error"]
    assert "stiffness" not in row and "solve_s" not in row
    assert not (replay / "uniform-fields.npz").exists()
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())
    import meshers._meshers
    workflow.save(args.run / "result.json", {"budget_failure_reproduced_and_stopped": True,
        "no_final_solve": True, "sources_unchanged": True,
        "native_sha256": workflow.digest(Path(meshers._meshers.__file__)), "workers": workers,
        "replay_error": row["error"]})


if __name__ == "__main__":
    main()
