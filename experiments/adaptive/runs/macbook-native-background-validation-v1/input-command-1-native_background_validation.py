"""Check the experimental background seam and rebuild the normal Python binding."""

import subprocess
import sys


commands = [
    ["cargo", "+1.95.0", "check", "--locked", "-p", "meshers-core", "--all-targets"],
    [
        sys.executable, "-m", "maturin", "develop", "--release", "--manifest-path",
        "crates/meshers-python/Cargo.toml", "--features", "experimental-surfaces",
    ],
    [sys.executable, "experiments/adaptive/native_fit_validation.py"],
]
for command in commands:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, check=True)
