"""Archive exact native inputs, run existing bounded-build checks, then seal binary.

Invoke inside macbook_stage.py, which owns the Job progress and process-tree
timeout/RSS watchdog. This wrapper does not submit work or change build policy.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--build-run", type=Path, required=True)
    args = parser.parse_args()
    archive = args.archive.resolve()
    archive.mkdir(parents=True, exist_ok=False)
    core = w.ROOT/"crates/meshers-core"
    sources = [w.ROOT/"Cargo.toml", w.ROOT/"Cargo.lock", core/"Cargo.toml",
               *sorted(core.rglob("*.rs")), Path(__file__).resolve(), w.HERE/"pmy_background_build.py"]
    inputs = {str(path.relative_to(w.ROOT)): w.digest(path) for path in sources}
    for path in sources:
        destination = archive/"sources"/path.relative_to(w.ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
    w.save(archive/"source-manifest.json", dict(inputs=inputs, status="snapshotted-before-build"))
    with (archive/"build.log").open("w") as stream:
        finished = subprocess.run([sys.executable, str(w.HERE/"pmy_background_build.py")],
            cwd=w.ROOT, stdout=stream, stderr=subprocess.STDOUT, check=False)
    assert finished.returncode == 0, f"Read {archive/'build.log'}; no binary receipt created"
    assert all(w.digest(w.ROOT/name) == expected == w.digest(archive/"sources"/name)
               for name, expected in inputs.items()), "Native source changed during build"
    executable = archive/"pmy_background"
    shutil.copy2(w.ROOT/"target/release/examples/pmy_background", executable)
    receipt = dict(binary_sha256=w.digest(executable), inputs=inputs,
        build_job=os.environ.get("PUEUE_TASK_ID"), build_run=str(args.build_run.resolve()),
        build_log_sha256=w.digest(archive/"build.log"),
        source_snapshot_manifest_sha256=w.digest(archive/"source-manifest.json"),
        rust=subprocess.check_output(["rustc", "+1.95.0", "--version"], text=True).strip(),
        scope="Experimental supplied-background observer snapshots; identical meshing guards and defaults; installed Python package unchanged")
    w.save(archive/"receipt.json", receipt)
    print(f"Native archive verified: {len(inputs)} inputs; binary {receipt['binary_sha256']}", flush=True)


if __name__ == "__main__":
    main()
