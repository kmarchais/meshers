"""Archive and build an isolated native wheel inside macbook_stage.py.

Queue macbook_stage.py with a 1200 s timeout and 8 GiB RSS ceiling. It owns
standard Job progress and the process-tree watchdog. The installed wheel is
left unchanged for the old-binary BCC control.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--scope", default="Experimental linear snapshots only; ordinary scalar generation and admission unchanged")
    args = parser.parse_args()
    archive = args.archive.resolve()
    archive.mkdir(parents=True, exist_ok=False)
    core, binding = [w.ROOT/"crates"/name for name in ("meshers-core", "meshers-python")]
    sources = {w.ROOT/"Cargo.toml", w.ROOT/"Cargo.lock", core/"Cargo.toml",
        binding/"Cargo.toml", binding/"Cargo.lock", binding/"pyproject.toml",
        binding/"README.md", binding/"LICENSE", binding/"THIRDPARTY.json",
        Path(__file__).resolve(), w.HERE/"cad_bcc_linear_capture.py"}
    sources.update(core.rglob("*.rs"))
    sources.update((binding/"src").rglob("*.rs"))
    sources.update((binding/"python").rglob("*.py"))
    sources.update((binding/"tests").rglob("*.py"))
    inputs = {str(path.relative_to(w.ROOT)): w.digest(path) for path in sorted(sources)}
    for path in sorted(sources):
        copied = archive/"sources"/path.relative_to(w.ROOT)
        copied.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, copied)
    w.save(archive/"source-manifest.json", dict(inputs=inputs, status="snapshotted-before-build"))
    manifest = str(binding/"Cargo.toml")
    commands = [
        ["cargo", "+1.95.0", "test", "--locked", "--release", "-p", "meshers-core",
         "--features", "experimental-quads"],
        ["cargo", "+1.95.0", "clippy", "--locked", "--release", "-p", "meshers-core", "--all-targets",
         "--features", "experimental-quads", "--", "-D", "warnings"],
        ["cargo", "+1.95.0", "check", "--locked", "--release", "--workspace", "--all-targets"],
        ["cargo", "+1.95.0", "check", "--locked", "--release", "--manifest-path", manifest],
        ["cargo", "+1.95.0", "clippy", "--locked", "--release", "--manifest-path", manifest,
         "--features", "experimental-surfaces", "--", "-D", "warnings"],
        [sys.executable, "-m", "maturin", "build", "--locked", "--release", "--manifest-path", manifest,
         "--features", "experimental-surfaces", "--interpreter", sys.executable, "--out", str(archive/"wheels")],
    ]
    for index, command in enumerate(commands):
        print(f"Native validation {index+1}/{len(commands)}", flush=True)
        with (archive/f"build-{index}.log").open("w") as stream:
            subprocess.run(command, cwd=w.ROOT, stdout=stream, stderr=subprocess.STDOUT, check=True)
    wheels = list((archive/"wheels").glob("*.whl"))
    assert len(wheels) == 1
    with zipfile.ZipFile(wheels[0]) as wheel:
        wheel.extractall(archive/"python")
    environment = dict(os.environ, PYTHONPATH=str(archive/"python"))
    test_code = (
        "from pathlib import Path; import meshers, pytest; "
        f"assert Path(meshers.__file__).resolve().is_relative_to(Path({str(archive/'python')!r})); "
        "assert meshers._meshers._LINEAR_OBSERVATIONS; "
        f"raise SystemExit(pytest.main(['-q', {str(binding/'tests/test_observations.py')!r}]))"
    )
    with (archive/"python-tests.log").open("w") as stream:
        subprocess.run([sys.executable, "-c", test_code],
            cwd=w.ROOT, env=environment, stdout=stream, stderr=subprocess.STDOUT, check=True)
    assert all(w.digest(w.ROOT/name) == digest == w.digest(archive/"sources"/name)
               for name, digest in inputs.items()), "Source changed during build"
    binaries = list((archive/"python/meshers").glob("_meshers*.so"))
    assert len(binaries) == 1
    w.save(archive/"receipt.json", dict(inputs=inputs, wheel_sha256=w.digest(wheels[0]),
        native_binary_sha256=w.digest(binaries[0]), package=str(archive/"python"),
        build_job=os.environ.get("PUEUE_TASK_ID"), installed_package_changed=False,
        scope=args.scope))
    print("Isolated observer wheel and off/on controls completed", flush=True)


if __name__ == "__main__":
    main()
