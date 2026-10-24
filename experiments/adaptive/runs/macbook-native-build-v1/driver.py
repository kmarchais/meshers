"""Record one bounded native setup stage through the shared Job helper."""

import argparse
import importlib.metadata
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

import microgen_workflow as workflow
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--timeout", type=int, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    assert command and args.timeout > 0
    assert sys.platform == "darwin" and platform.machine() == "arm64"
    os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
        RAYON_NUM_THREADS="1", CARGO_BUILD_JOBS="2", RUSTUP_TOOLCHAIN="1.95.0",
        VIRTUAL_ENV=str(Path(sys.executable).parent.parent))
    os.environ.pop("LD_LIBRARY_PATH", None)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    requirements = subprocess.check_output(["uv", "pip", "freeze", "--python", sys.executable], text=True)
    (args.run / "requirements.txt").write_text(requirements)
    manifest = {"command": command, "platform": platform.platform(), "architecture": platform.machine(),
        "python": sys.version, "executable": sys.executable, "memory_bytes": workflow.MEMORY,
        "memory_enforcement": "process-tree-rss-watchdog", "rss_poll_interval_s": .25,
        "timeout_s": args.timeout, "numerical_threads": 1, "build_jobs": 2,
        "rust": subprocess.check_output(["rustc", "+1.95.0", "--version"], text=True).strip(),
        "inputs": {str(path): workflow.digest(path) for path in
            (Path(__file__), Path(workflow.__file__), workflow.ROOT / "Cargo.lock")},
        "versions": {name: importlib.metadata.version(name) for name in
            ("numpy", "scipy", "fedoo", "microgen", "mmgpy", "pyamg", "pyvista", "vtk", "maturin", "psutil")}}
    workflow.save(args.run / "manifest.json", manifest)
    with Job(parts=["worker"]) as job:
        with worker_part(job, "worker") as part:
            record = workflow.bounded(command, args.run / "worker.log", part, timeout=args.timeout)
            workflow.save(args.run / "worker.json", record)
            job.result("worker", seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3,
                returncode=record["returncode"])
            assert record["returncode"] == 0, "Read worker.log and worker.json"


if __name__ == "__main__":
    main()
