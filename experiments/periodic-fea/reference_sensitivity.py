"""Compare two quadratic stiffness extrapolations; no continuum error certificate."""

import argparse
import json
from pathlib import Path

import numpy as np

from periodic import tensor_difference
from run import digest


def read(path, geometry_gate=None):
    result = json.loads(path.read_text())
    geometry_ok = result.get("geometry_eligible")
    if geometry_gate is not None:
        geometry_ok = all(result.get(k) is not None and result[k] <= geometry_gate for k in
                          ["sampled_geometry_distance", "sampled_reverse_distance", "sampled_cap_overshoot"])
        geometry_ok &= result.get("sampled_outward_orientation", False)
        geometry_ok &= result.get("relative_volume_error_with_uncertainty", float("inf")) <= .01
    if not geometry_ok or not all(result.get(k) for k in ["homogenization_audited",
                                                         "quality_eligible", "periodic_eligible"]):
        raise ValueError(f"Unadmitted artifact: {path}")
    parent = json.loads((path.parent.parent / "results.json").read_text())
    if parent["state"] != "complete" or not any(
        r["artifact"] == path.parent.name for r in parent["artifact_verification"]["rows"]
    ):
        raise ValueError(f"Saved fields not independently verified: {path}")
    parent_row = next(row for row in parent["rows"] if row["artifact"] == path.parent.name)
    if parent_row["homogenization"] != result["homogenization"] or parent_row["mesh_sha256"] != result["mesh_sha256"]:
        raise ValueError(f"Artifact disagrees with independently audited parent run: {path}")
    if digest(path.parent / "mesh.npz") != result["mesh_sha256"]:
        raise ValueError(f"Saved mesh changed: {path}")
    tensor = np.asarray(result["homogenization"]["stiffness_engineering_voigt"])
    return result, tensor, {"path": str(path), "result_sha256": digest(path),
                            "mesh_sha256": result["mesh_sha256"], "native_sha256": result["native_sha256"],
                            "common_geometry_gate_override": geometry_gate}


def sensitivity(reference_paths, candidate_paths):
    references = [read(p) for p in reference_paths]
    case = references[0][0]["case"]
    binary = references[0][0]["native_sha256"]
    cells = [r[0]["cells"] for r in references]
    if not cells[0] < cells[1] < cells[2] or any(
        r[0]["case"] != case or r[0]["workflow"] != "meshers-pass1"
        or r[0]["native_sha256"] != binary for r in references
    ):
        raise ValueError("Need three increasing grids of the same admitted native case/build")
    fits = []
    for coarse, fine in zip(references[:-1], references[1:]):
        ratio = (coarse[0]["cells"] / fine[0]["cells"]) ** 2
        tensor = (fine[1] - ratio * coarse[1]) / (1 - ratio)
        fits.append({"grids": [coarse[0]["cells"], fine[0]["cells"]],
                     "extrapolated_tensor": tensor.tolist()})
    rows = []
    for path in candidate_paths:
        result, tensor, provenance = read(path)
        if result["case"] != case:
            raise ValueError("Candidate geometry differs from the reference sequence")
        rows.append({**provenance, "native_sha256": result["native_sha256"],
                     "workflow": result["workflow"], "cells": result["cells"],
                     "difference_to_each_extrapolation": [
                         tensor_difference(tensor, fit["extrapolated_tensor"]) for fit in fits]})
    return {
        "scope": "Quadratic convergence sensitivity only; no continuum error bound or certification",
        "assumption": "Leading full-tensor error proportional to 1/cells^2, including geometry/discretization effects",
        "case": case, "references": [r[2] for r in references], "fits": fits,
        "finest_observed_to_last_extrapolation": tensor_difference(
            references[-1][1], fits[-1]["extrapolated_tensor"]),
        "difference_between_extrapolations": tensor_difference(
            fits[0]["extrapolated_tensor"], fits[1]["extrapolated_tensor"]),
        "candidates": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, nargs=3, required=True)
    parser.add_argument("--candidates", type=Path, nargs="+", required=True)
    args = parser.parse_args()
    print(json.dumps(sensitivity(args.references, args.candidates), indent=2))
