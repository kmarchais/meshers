"""Compare the same candidates against retained volume reference topologies."""

import argparse
import json
import sys
from pathlib import Path

import accuracy_metrics
from run import digest, write_json

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report


def compare(references, candidates, output, geometry_gate, power=17):
    # Read-only inputs and atomic fixed-name outputs make an interrupted run restartable.
    output.mkdir(parents=True, exist_ok=True)
    boundary_checks = []
    for path in candidates:
        row = json.loads(path.read_text())
        if row["workflow"] not in {"meshers-body-refine", "meshers-native-body"}:
            continue
        import numpy as np
        from body_refine import verify_preserved_boundary
        seed_dir = path.parent.parent / "source/body-seed"
        seed_row = json.loads((seed_dir / "result.json").read_text())
        if digest(seed_dir / "mesh.npz") != seed_row["mesh_sha256"] or digest(path.parent / "mesh.npz") != row["mesh_sha256"]:
            raise ValueError("Body-refinement boundary-check input changed")
        with np.load(seed_dir / "mesh.npz") as seed, np.load(path.parent / "mesh.npz") as mesh:
            boundary_checks.append({"artifact": str(path), **verify_preserved_boundary(
                seed["points"], seed["surface"], mesh["points"], mesh["surface"])})
    rows = []
    for index, reference in enumerate(references):
        accuracy_metrics.report = lambda fraction, message, **kwargs: report(
            (index + fraction) / len(references), f"Reference {index + 1}/{len(references)}: {message}", **kwargs)
        rows.append(accuracy_metrics.compare(reference, candidates, output / f"reference-{index}.json",
                                             power=power, geometry_gate=geometry_gate))
    write_json(output / "comparisons.json", {"state": "complete", "source_sha256": digest(Path(__file__)),
                                             "common_geometry_gate": geometry_gate, "comparisons": rows,
                                             "body_boundary_verification": boundary_checks})
    report(1, "Accuracy comparisons against all reference topologies complete", eta_s=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, nargs="+", required=True)
    parser.add_argument("--candidates", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--common-geometry-gate", type=float, default=.002)
    parser.add_argument("--power", type=int, default=17)
    args = parser.parse_args()
    if not 8 <= args.power <= 21:
        parser.error("--power must be between 8 and 21")
    compare(args.references, args.candidates, args.output, args.common_geometry_gate, args.power)
