"""Validate corrected tie gradients with the full Python suite and actual specimens."""

import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job
import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", RAYON_NUM_THREADS="1")
run = Path(__file__).parent
assert not (run / "python-tests.json").exists()
import meshers._meshers

sources = [Path(__file__), ROOT / "crates/meshers-python/python/meshers/expression.py",
    ROOT / "crates/meshers-python/tests/test_compiled.py", Path(meshers._meshers.__file__),
    ROOT / "experiments/adaptive/gradient_tie_probe.py", ROOT / "experiments/adaptive/settle_controls.py",
    ROOT / "experiments/adaptive/admission_pilot.py", Path(workflow.__file__)]
hashes = {str(path): workflow.digest(path) for path in sources}
workflow.save(run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
    "worker_timeouts_s": {"python-tests": 1200, "lattice-gradients": 120, "geometry-controls": 600},
    "numerical_thread_environment": 1, "test_threads": "Existing tests retain their explicit concurrency checks",
    "native_rebuild_required": False, "change_scope": "Python derivative graph; same embedded native compiler"})
commands = [
    ("python-tests", [sys.executable, "-m", "pytest", "-q", str(ROOT / "crates/meshers-python/tests")], 1200),
    ("lattice-gradients", [sys.executable, str(ROOT / "experiments/adaptive/gradient_tie_probe.py"),
        "--output", str(run / "lattice-gradients-result.json")], 120),
    ("geometry-controls", [sys.executable, str(ROOT / "experiments/adaptive/settle_controls.py"),
        "--run", str(run / "controls")], 600),
]
with Job(parts=[name for name, _, _ in commands]) as job:
    for name, command, timeout in commands:
        with job.part(name) as part:
            progress_path = str(run / (name + "-progress.json"))
            isolated = ["env", "PUEUE_PROGRESS_FILE=" + progress_path, "GPU_PROGRESS_FILE=" + progress_path, *command]
            record = workflow.bounded(isolated, run / (name + ".log"), part, timeout=timeout)
            workflow.save(run / (name + ".json"), record)
            assert record["returncode"] == 0, f"Read {name}.log"
            job.result(name, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())
