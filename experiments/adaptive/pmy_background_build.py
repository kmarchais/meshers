"""Build and check the experimental background path inside a bounded queue stage."""

import subprocess


commands = [
    [
        "cargo", "+1.95.0", "test", "--locked", "--release", "--lib", "-p", "meshers-core",
        "--features", "experimental-quads",
        "implicit::topology_tests::",
        "--", "--nocapture",
    ],
    [
        "cargo", "+1.95.0", "build", "--locked", "--release", "-p", "meshers-core",
        "--features", "experimental-quads", "--example", "pmy_background",
    ],
    [
        "cargo", "+1.95.0", "clippy", "--locked", "--release", "-p", "meshers-core",
        "--all-targets", "--features", "experimental-quads", "--", "-D", "warnings",
    ],
]
for command in commands:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, check=True)
