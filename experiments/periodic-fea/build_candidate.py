"""Build and test an isolated optimizer candidate without replacing installed packages."""

import argparse
import json
import os
import subprocess
import sys
import threading
import time
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report
from run import digest, write_json


def build(worktree, output, target_directory=None):
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    live = {"done": 0, "message": "Verify committed candidate source"}
    stop = threading.Event()

    def heartbeat():
        while not stop.is_set():
            report(live["done"] / 3, live["message"], eta_s=max(0, 240 - (time.monotonic() - started)),
                   seconds=round(time.monotonic() - started))
            stop.wait(2)

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        def git(*args):
            return subprocess.check_output(["git", *args], cwd=worktree, text=True).strip()

        if git("status", "--porcelain", "--untracked-files=no"):
            raise ValueError("Commit candidate changes before building")
        commit = git("rev-parse", "HEAD")
        files = git("ls-files", "crates/meshers-core", "crates/meshers-python", "Cargo.toml", "Cargo.lock").splitlines()
        hashes = {name: digest(worktree / name) for name in files}
        env = os.environ.copy()
        env["CARGO_BUILD_JOBS"] = "2"
        env["CARGO_TARGET_DIR"] = str(target_directory or output / "target")
        commands = [
            ["cargo", "test", "-p", "meshers-core", "--test", "implicit", "--release"],
            [str(Path.home() / "dev/meshers/.venv/Scripts/maturin.exe"), "build", "--release",
             "--manifest-path", "crates/meshers-python/Cargo.toml", "--interpreter", sys.executable,
             "--out", str(output / "wheels")],
        ]
        records = []
        for index, command in enumerate(commands):
            live["message"] = ["Compile and run implicit core checks", "Build isolated experimental Python wheel"][index]
            with (output / f"build-{index}.log").open("w", encoding="utf-8") as log:
                completed = subprocess.run(command, cwd=worktree, env=env, stdout=log,
                                           stderr=subprocess.STDOUT, timeout=1200)
            records.append({"command": command, "exit_code": completed.returncode})
            if completed.returncode:
                print((output / f"build-{index}.log").read_text()[-4000:], flush=True)
                raise RuntimeError(f"Candidate build stage {index} failed")
            live["done"] += 1
        for name, expected in hashes.items():
            if digest(worktree / name) != expected:
                raise ValueError(f"Source changed while building: {name}")
        wheel = next((output / "wheels").glob("*.whl"))
        with zipfile.ZipFile(wheel) as archive:
            for name in archive.namelist():
                if name.startswith("meshers/"):
                    archive.extract(name, output / "package")
        native = next((output / "package/meshers").glob("*.pyd"))
        result = {"state": "complete", "source_commit": commit, "source_sha256": hashes,
                  "native_sha256": digest(native), "wheel_sha256": digest(wheel), "commands": records,
                  "worktree": str(worktree), "seconds": time.monotonic() - started}
        write_json(output / "build.json", result)
        live["done"] = 3
        report(1, "Experimental optimizer wheel and core checks complete", eta_s=0)
        print(json.dumps({k: result[k] for k in ["source_commit", "native_sha256", "seconds"]}), flush=True)
    finally:
        stop.set()
        thread.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worktree", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target-directory", type=Path, help="Reuse this experiment's serial build cache")
    args = parser.parse_args()
    build(args.worktree.resolve(), args.output.resolve(),
          args.target_directory.resolve() if args.target_directory else None)
