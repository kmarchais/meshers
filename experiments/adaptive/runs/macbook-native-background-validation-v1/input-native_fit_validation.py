"""Full regression checks for the native fitting change, inside one bounded queue stage."""

import subprocess
import sys

commands = [
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
    ],
    [
        "cargo",
        "+1.95.0",
        "clippy",
        "--locked",
        "--release",
        "--manifest-path",
        "crates/meshers-python/Cargo.toml",
        "--all-targets",
        "--features",
        "experimental-surfaces",
        "--",
        "-D",
        "warnings",
    ],
    [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "--tb=short",
        "crates/meshers-python/tests",
        "tools/tests",
    ],
]
for command in commands:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, check=True)
