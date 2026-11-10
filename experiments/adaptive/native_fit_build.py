"""Fail-fast native fitting validation and editable arm64 rebuild, inside a queued stage."""

import subprocess
import sys

commands = [
    [
        "cargo",
        "+1.95.0",
        "clippy",
        "--locked",
        "--release",
        "-p",
        "meshers-core",
        "--all-targets",
        "--features",
        "experimental-quads",
        "--",
        "-D",
        "warnings",
    ],
    [
        "cargo",
        "+1.95.0",
        "test",
        "--locked",
        "--release",
        "-p",
        "meshers-core",
        "--features",
        "experimental-quads",
        "curved_fit",
        "--",
        "--nocapture",
    ],
    [
        sys.executable,
        "-m",
        "maturin",
        "develop",
        "--release",
        "--manifest-path",
        "crates/meshers-python/Cargo.toml",
        "--features",
        "experimental-surfaces",
    ],
    [sys.executable, "-m", "pytest", "-q", "crates/meshers-python/tests/test_curved_fit.py"],
]
for command in commands:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, check=True)
