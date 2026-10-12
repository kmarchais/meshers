"""Compare retained periodic tensors with an observed fine-grid reference."""

import argparse
import hashlib
import json
from pathlib import Path

from periodic import tensor_difference


def compare(reference_path, run_paths):
    reference = json.loads(reference_path.read_text())
    if not all(reference.get(key) for key in (
        "homogenization_audited", "geometry_eligible", "quality_eligible", "periodic_eligible",
    )):
        raise ValueError("Reference has not passed all admission and physical audits")
    tensor = reference["homogenization"]["stiffness_engineering_voigt"]
    rows = []
    for run_path in run_paths:
        document = json.loads((run_path / "results.json").read_text())
        if document["state"] != "complete":
            raise ValueError(f"Run is not independently verified: {run_path}")
        for row in document["rows"]:
            if row["case"] != reference["case"] or not row.get("homogenization_audited"):
                continue
            artifact = run_path / row["artifact"] / "result.json"
            row = json.loads(artifact.read_text())
            rows.append({
                "artifact": str(artifact),
                "artifact_sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
                **{key: row[key] for key in (
                    "workflow", "cells", "tetrahedra", "mesh_s", "q_min",
                    "geometry_eligible", "quality_eligible", "periodic_eligible",
                )},
                "audited_homogenization_s": row["homogenization"]["seconds"],
                **tensor_difference(row["homogenization"]["stiffness_engineering_voigt"], tensor),
            })
    return {
        "scope": "Observed fine-grid tensor disagreement; no continuum error certificate",
        "timing_scope": "Single samples without warmup or repeats",
        "reference": str(reference_path),
        "reference_sha256": hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        "rows": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("runs", type=Path, nargs="+")
    args = parser.parse_args()
    print(json.dumps(compare(args.reference, args.runs), indent=2))
