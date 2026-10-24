"""Exercise real worker exits and resource stops, without numerical dependencies."""

import sys
import os

import pytest

import microgen_workflow as workflow


class Part:
    def update(self, *args, **kwargs):
        pass


@pytest.mark.skipif(sys.platform != "darwin", reason="Apple Accelerate control")
def test_accelerate_runs_with_verified_single_thread_policy():
    import ctypes
    import numpy as np

    accelerate = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
    assert accelerate.BLASGetThreading() == 1
    assert workflow.NUMERICAL_THREAD_POLICY["accelerate_blas_threading"] == 1
    assert os.environ["VECLIB_MAXIMUM_THREADS"] == "1"
    matrix = np.eye(8) * 2 + np.ones((8, 8))
    np.testing.assert_allclose(matrix @ np.linalg.solve(matrix, np.eye(8)), np.eye(8), atol=1e-14)
    assert accelerate.BLASGetThreading() == 1


@pytest.mark.parametrize("mode", ["success", "timeout", "memory"])
def test_real_bounded_worker(tmp_path, monkeypatch, mode):
    # A low RSS ceiling is a controlled watchdog test, not a scientific cap change.
    monkeypatch.setattr(workflow, "MEMORY", 64 * 1024**2 if mode == "memory" else 8 * 1024**3)
    source = {
        "success": "print('worker complete')",
        "timeout": "import time; time.sleep(10)",
        "memory": "import time; data=bytearray(96*1024**2); time.sleep(10)",
    }[mode]
    # Linux RLIMIT_AS can reject the allocation before the RSS watchdog sees it.
    if mode == "memory":
        monkeypatch.setattr(sys, "platform", "darwin")
    record = workflow.bounded([sys.executable, "-c", source], tmp_path / "worker.log",
        Part(), timeout=.5 if mode == "timeout" else 5)
    assert record["limit_failure"] == {"success": None, "timeout": "timeout", "memory": "memory-cap"}[mode]
    assert (record["returncode"] == 0) == (mode == "success")
    assert record["wall_s"] < 6
    assert record["memory_bytes"] == workflow.MEMORY
    if mode == "success":
        assert (tmp_path / "worker.log").read_text().strip() == "worker complete"
    if mode == "memory":
        assert record["peak_rss_bytes"] > record["memory_bytes"]
