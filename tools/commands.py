# Copyright (c) 2026 Kevin Marchais

"""Run synchronous maintainer commands with an argument list and an explicit executable."""

from __future__ import annotations

import shutil
import subprocess  # ruff: ignore[suspicious-subprocess-import] - maintainer tools launch child processes.
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping


CalledProcessError = subprocess.CalledProcessError
TimeoutExpired = subprocess.TimeoutExpired


def run_command(
    executable: str | Path,
    *args: str,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    timeout: int | None = None,
    stdout_encoding: str | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run a tool, inheriting stderr and collecting stdout when an encoding is supplied.

    Returns:
        The completed process with decoded output when requested.

    Raises:
        FileNotFoundError: The requested executable could not be located.
    """
    path = shutil.which(str(executable))
    if path is None:
        msg = f"Executable not found: {executable}"
        raise FileNotFoundError(msg)
    # Resolving a venv interpreter symlink would switch to the base installation.
    return subprocess.run(  # ruff: ignore[subprocess-without-shell-equals-true] - resolved executable and separate arguments.
        [str(Path(path).absolute()), *args],
        cwd=cwd,
        env=env,
        timeout=timeout,
        stdout=subprocess.PIPE if stdout_encoding else None,
        text=True,
        encoding=stdout_encoding,
        check=True,
    )
