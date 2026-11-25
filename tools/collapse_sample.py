# Copyright (c) 2026 Kevin Marchais

"""Replay a release fixture and fingerprint its returned arrays outside generation timing."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from importlib import import_module
from typing import TYPE_CHECKING

from tools.performance_gate import configure_worker

if TYPE_CHECKING:
    import numpy as np
    from meshers import Mesh

CASES = ("adaptive_spinodoid", "spinodoid", "adaptive_lattice")


def array_hash(array: np.ndarray) -> str:
    """Return the SHA-256 of type, shape and C-order bytes."""
    header = json.dumps([array.dtype.str, array.shape]).encode("ascii")
    return hashlib.sha256(header + b"\0" + array.tobytes()).hexdigest()


def mesh_hashes(mesh: Mesh) -> dict[str, str]:
    """Return separate digests for coordinates, connectivity, tags and periodic pairs."""
    return {
        **{
            name: array_hash(getattr(mesh, name))
            for name in ("points", "tetrahedra", "surface", "boundary_tags")
        },
        **{f"pairs_{axis}": array_hash(pairs) for axis, pairs in enumerate(mesh.periodic_pairs)},
    }


def measure(case: str, native_sha256: str) -> dict:
    """Generate once, then run the release validator on precisely those arrays.

    Returns:
        Costs, mesh digests, release checks.

    Raises:
        ValueError: The loaded extension differs from the requested artifact.
    """
    configure_worker()
    # Numerical imports must follow thread configuration in each fresh process.
    worker = import_module("tools.performance_worker")

    process_start = time.process_time()
    mesh, elapsed, peak = worker._measure_generation(case)  # ruff: ignore[private-member-access] - pinned release fixture.
    cpu_seconds = time.process_time() - process_start
    identity = worker.runtime_identity()
    if identity["native_sha256"] != native_sha256:
        msg = "Loaded native extension differs from the pinned artifact"
        raise ValueError(msg)
    # sample() owns the release checks. Hand it the same mesh and stop its second generation.
    patch = import_module("unittest.mock").patch
    with patch.object(worker, "_measure_generation", return_value=(mesh, elapsed, peak)):
        validation = worker.sample(case)
    del validation["seconds"], validation["peak_rss_bytes"]
    return {
        "identity": identity,
        "metrics": {"wall_s": elapsed, "cpu_s": cpu_seconds, "peak_bytes": peak},
        "hashes": mesh_hashes(mesh),
        "validation": validation,
    }


def main() -> None:
    """Emit one JSON row on stdout for a fresh worker."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("native_sha256")
    args = parser.parse_args()
    json.dump(
        measure(args.case, args.native_sha256),
        sys.stdout,
        allow_nan=False,
    )


if __name__ == "__main__":
    main()
