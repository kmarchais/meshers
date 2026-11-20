# Copyright (c) 2026 Kevin Marchais

"""Describe a benchmark configuration using public hardware and OS specifications."""

import json
import os
import platform
from pathlib import Path

if __package__:
    from .commands import run_command
else:
    from commands import run_command


def hardware_snapshot() -> dict[str, str | int | None]:
    """Return the CPU, OS-visible physical memory and operating system at collection time."""
    system = platform.system()
    cpu = platform.processor()
    ram_bytes = None
    if system == "Windows":
        query = (
            "[Console]::OutputEncoding = [Text.Encoding]::UTF8; "
            "@{cpu=(Get-CimInstance Win32_Processor | Select-Object -First 1).Name;"
            "ram=(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory} "
            "| ConvertTo-Json -Compress"
        )
        powershell = (
            Path(os.environ.get("SYSTEMROOT") or os.environ.get("WINDIR") or r"C:\Windows")
            / "System32/WindowsPowerShell/v1.0/powershell.exe"
        )
        specs = json.loads(
            run_command(
                powershell,
                "-NoProfile",
                "-NonInteractive",
                "-Command",
                query,
                stdout_encoding="utf-8",
            ).stdout,
        )
        cpu, ram_bytes = specs["cpu"], int(specs["ram"])
    elif system == "Darwin":
        cpu = run_command(
            "/usr/sbin/sysctl",
            "-n",
            "machdep.cpu.brand_string",
            stdout_encoding="utf-8",
        ).stdout.strip()
        ram_bytes = int(
            run_command("/usr/sbin/sysctl", "-n", "hw.memsize", stdout_encoding="utf-8").stdout,
        )
    elif system == "Linux":
        lines = Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines()
        cpu = next(
            (line.partition(":")[2].strip() for line in lines if line.startswith("model name")),
            cpu,
        )
        ram_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    return {
        "cpu_model": cpu or "unknown",
        "cpu_count": os.cpu_count(),
        "ram_bytes": ram_bytes,
        "machine": platform.machine(),
        "platform": platform.platform(),
    }
