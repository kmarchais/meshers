# Copyright (c) 2026 Kevin Marchais

"""Run a binding test in an isolated interpreter and show both streams on failure."""

import subprocess  # ruff: ignore[suspicious-subprocess-import] - tests execute the current interpreter
import sys
from pathlib import Path

import pytest


def run_python(
    *args: str,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    timeout: int = 90,
) -> str:
    """Run a test script in the active Python environment.

    Returns:
        Standard output from the successful child process.
    """
    result = subprocess.run(  # ruff: ignore[subprocess-without-shell-equals-true] - fixed executable, no shell
        [sys.executable, *args],
        cwd=cwd,
        env=env,
        timeout=timeout,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode:
        pytest.fail(f"Child Python exited with {result.returncode}\n{result.stdout}{result.stderr}")
    return result.stdout
