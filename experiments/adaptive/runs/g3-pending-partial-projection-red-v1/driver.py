"""Prove the focused pending-wall regression fails before changing the algorithm."""

import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job
import microgen_workflow as workflow

os.environ.update(CARGO_BUILD_JOBS="2", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
run = Path(__file__).parent
assert not (run / "result.json").exists()
sources = [Path(__file__), ROOT / "crates/meshers-core/src/tet_adapt.rs", ROOT / "Cargo.lock"]
hashes = {str(path): workflow.digest(path) for path in sources}
workflow.save(run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
    "worker_timeout_s": 1200, "build_threads": 2, "test_threads": 1, "expected_failure": True})
with Job(parts=["regression-red"]) as job:
    with job.part("regression-red") as part:
        record = workflow.bounded(["cargo", "+1.95.0", "test", "--release", "--locked",
            "-p", "meshers-core", "pending_wall_projection_makes_safe_partial_progress",
            "--", "--test-threads=1"], run / "worker.log", part, timeout=1200)
        log = (run / "worker.log").read_text()
        expected = (record["returncode"] != 0 and record["limit_failure"] is None
            and "pending wall vertex did not make safe partial progress" in log
            and "test result: FAILED" in log)
        workflow.save(run / "result.json", {"worker": record, "regression_failed_as_expected": expected})
        assert expected, "Expected the focused assertion failure, not a build or resource failure"
        job.result("regression-red", expected_failure=int(expected), seconds=record["wall_s"])
assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())
