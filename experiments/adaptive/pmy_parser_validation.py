"""Check native rollback and build the independent input parser inside qrun."""

import subprocess


commands = [
    ["cargo", "+1.95.0", "test", "--locked", "--release", "-p", "meshers-core",
     "--features", "experimental-quads"],
    ["cargo", "+1.95.0", "build", "--locked", "--release", "-p", "meshers-core",
     "--features", "experimental-quads", "--example", "pmy_input_roundtrip"],
    ["cargo", "+1.95.0", "clippy", "--locked", "--release", "--manifest-path",
     "crates/meshers-python/Cargo.toml", "--all-targets", "--features",
     "experimental-surfaces", "--", "-D", "warnings"],
]
for command in commands:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, check=True)
