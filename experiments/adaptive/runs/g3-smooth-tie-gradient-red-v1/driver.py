"""Require the smooth-composition gradient regression to fail before correction."""

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
assert not (run / "result.json").exists()
import meshers._meshers

sources = [Path(__file__), ROOT / "crates/meshers-python/python/meshers/expression.py",
           ROOT / "crates/meshers-python/tests/test_compiled.py", Path(meshers._meshers.__file__)]
hashes = {str(path): workflow.digest(path) for path in sources}
workflow.save(run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
    "worker_timeout_s": 300, "threads": 1, "expected_failure": True})
with Job(parts=["gradient-red"]) as job:
    with job.part("gradient-red") as part:
        record = workflow.bounded([sys.executable, "-m", "pytest", "-q",
            str(ROOT / "crates/meshers-python/tests/test_compiled.py"), "-k", "tied_smooth_composition"],
            run / "worker.log", part, timeout=300)
        log = (run / "worker.log").read_text()
        expected = (record["returncode"] == 1 and record["limit_failure"] is None
            and "smooth-composition tie gradient mismatch" in log and "2 failed" in log)
        workflow.save(run / "result.json", {"worker": record, "regression_failed_as_expected": expected})
        assert expected, "Expected the two gradient assertions, not a harness or resource failure"
        job.result("gradient-red", expected_failure=int(expected), seconds=record["wall_s"])
assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())
