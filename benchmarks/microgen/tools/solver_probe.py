"""Replay a saved mesh with another linear solver; keep the FEA problem unchanged.

Queue this script. Meshing and the manufactured elasticity problem stay unchanged.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from microgen_benchmark import elasticity
from microgen_benchmark.geometry import boundary
from microgen_benchmark.runner import progress, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--method", choices=["direct", "amg"], default="amg")
    parser.add_argument("--scale", type=float, default=1)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a new output path")
    mesh_path = args.artifact / "mesh.npz"
    data = np.load(mesh_path)
    points, tets = data["points"], data["tetrahedra"]
    reference = json.loads((args.artifact / "result.json").read_text(encoding="utf-8"))
    start = time.monotonic()
    progress(0, 24, "Preparing solver replay", start, 0)

    def checkpoint(**values):
        load = values.get("fea_rhs", 0)
        progress(
            load,
            24,
            f"Solver load {load}/24, iteration {values.get('fea_cg_iterations', 0)}",
            start,
            0,
        )

    metrics, fields = elasticity.solve(
        points,
        tets,
        boundary(points, tets),
        args.scale,
        linear_solver=args.method,
        checkpoint=checkpoint,
    )
    reference_fields = np.load(args.artifact / "fields.npz")
    difference = fields["displacements"] - reference_fields["displacements"]
    relative_difference = float(
        np.linalg.norm(difference) / np.linalg.norm(reference_fields["displacements"])
    )
    error_difference = abs(metrics["worst_energy_error"] - reference["worst_energy_error"])
    if relative_difference > 1e-6 or error_difference > 1e-7:
        raise ValueError(
            f"Replay disagrees with direct solution: field {relative_difference}, energy {error_difference}"
        )
    report = {
        "method": args.method,
        "metrics": metrics,
        "elapsed_s": time.monotonic() - start,
        "field_relative_difference": relative_difference,
        "energy_error_absolute_difference": error_difference,
        "mesh_sha256": hashlib.sha256(mesh_path.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(elasticity.__file__).read_bytes()).hexdigest(),
    }
    write_json(args.output, report)
    progress(24, 24, "Solver replay verified against saved direct solution", start, 0)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
