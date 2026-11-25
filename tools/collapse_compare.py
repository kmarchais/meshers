# Copyright (c) 2026 Kevin Marchais

"""Replay the collapse guard comparison with two unpacked release wheels."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from tools.benchmark_hardware import hardware_snapshot
from tools.collapse_sample import CASES
from tools.commands import run_command

if TYPE_CHECKING:
    from collections.abc import Callable


def compare(root: Path, progress: Callable[[float, str], None]) -> None:
    """Alternate wheel order per pair and require exact agreement across all samples.

    Raises:
        ValueError: The native binaries coincide or a fixture's validation or arrays differ.
    """
    packages = {side: root / f"{side}-packages" for side in ("base", "candidate")}
    hashes = {
        side: hashlib.sha256(
            next(
                p for p in (path / "meshers").glob("_meshers.*") if p.suffix in {".so", ".pyd"}
            ).read_bytes(),
        ).hexdigest()
        for side, path in packages.items()
    }
    if hashes["base"] == hashes["candidate"]:
        msg = "Base and candidate load the same native binary"
        raise ValueError(msg)
    steps = [
        (other, block, case, side)
        for other, pairs in (("base", 7), ("candidate", 13))
        for block in range(-1, pairs)
        for case in CASES
        for side in ("ab" if block < 0 or block % 2 == 0 else "ba")
    ]
    references = {}
    sources = {}
    runs = []
    for i, (other, block, case, side) in enumerate(steps):
        stage = "base" if side == "a" else other
        progress(i / len(steps), f"{case} {stage}, pair {block}")
        result = run_command(
            sys.executable,
            "-m",
            "tools.collapse_sample",
            case,
            hashes[stage],
            env=dict(os.environ, PYTHONPATH=str(packages[stage])),
            stdout_encoding="utf-8",
        )
        row = json.loads(result.stdout)
        reference = {key: row[key] for key in ("validation", "hashes")}
        if references.setdefault(case, reference) != reference:
            msg = f"Mesh or validation changed for {case} in {stage} pair {block}"
            raise ValueError(msg)
        sources[stage] = row["identity"]
        runs.append({
            "case": case,
            "side": side,
            "block": block,
            "other": other,
            "warmup": block < 0,
            "metrics": row["metrics"],
        })
    record = {
        "hardware": hardware_snapshot(),
        "sources": sources,
        "references": references,
        "runs": runs,
    }
    (root / "compare.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    """Run from the candidate checkout using the active numerical environment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="Contains base-packages and candidate-packages")
    args = parser.parse_args()
    compare(args.root.resolve(), lambda _, label: sys.stderr.write(label + "\n"))


if __name__ == "__main__":
    main()
