"""Report warmed generation at an observed/model target, never certified continuum error."""

import argparse
import json
import statistics
from pathlib import Path

from periodic import tensor_difference
from reference_sensitivity import sensitivity
from run import digest


def qualify(reference_paths, native_path, mmg_path, repeat_path):
    model = sensitivity(reference_paths, [native_path, mmg_path])
    artifacts = [json.loads(p.read_text()) for p in [native_path, mmg_path]]
    if [r["workflow"] for r in artifacts] != ["meshers-pass1", "mmg-native-seed"]:
        raise ValueError("Expected direct native and standard native-seed MMG artifacts")
    fine = json.loads(reference_paths[-1].read_text())["homogenization"]["stiffness_engineering_voigt"]
    repeated = json.loads((repeat_path / "generation-repeats.json").read_text())
    manifest = json.loads((repeat_path / "manifest.json").read_text())
    if repeated["state"] != "complete" or repeated["manifest_sha256"] != digest(repeat_path / "manifest.json"):
        raise ValueError("Repeated generation is incomplete or its manifest changed")
    if len(repeated["rows"]) != 8 or manifest["native_sha256"] != artifacts[0]["native_sha256"]:
        raise ValueError("Need eight repetitions of the same native build")
    rows = []
    for path, artifact, modeled in zip([native_path, mmg_path], artifacts, model["candidates"]):
        matches = [r for r in repeated["rows"] if r["workflow"] == artifact["workflow"]]
        warmups = [r for r in matches if r["warmup"]]
        samples = [r for r in matches if not r["warmup"]]
        if len(warmups) != 1 or len(samples) != 3 or any(
            not r["exact_mesh_identity"] or r["cells"] != artifact["cells"]
            or r["tetrahedra"] != artifact["tetrahedra"]
            or r["audited_mesh_sha256"] != artifact["mesh_sha256"]
            or Path(r["audited_artifact"]).resolve() != path.parent.resolve() for r in matches
        ):
            raise ValueError("Repeated outputs do not identify the selected audited artifact")
        median = statistics.median(r["mesh_s"] for r in samples)
        if median <= 0 or median != repeated["medians_s"][artifact["workflow"]]:
            raise ValueError("Stored median does not match the actual samples")
        observed = tensor_difference(artifact["homogenization"]["stiffness_engineering_voigt"], fine)
        modeled_scores = [r["worst_relative_strain_energy_difference"]
                          for r in modeled["difference_to_each_extrapolation"]]
        rows.append({
            "artifact": str(path), "result_sha256": digest(path),
            "mesh_sha256": artifact["mesh_sha256"],
            **{k: artifact[k] for k in ["case", "workflow", "cells", "tetrahedra", "nodes", "q_min"]},
            "observed_worst_relative_energy_difference": observed["worst_relative_strain_energy_difference"],
            "modeled_worst_relative_energy_differences": modeled_scores,
            "observed_and_modeled_target_met": observed["worst_relative_strain_energy_difference"] <= .01
                and max(modeled_scores) <= .01,
            "generation_samples_s": [r["mesh_s"] for r in samples],
            "generation_median_s": median,
            "audited_fea_single_sample_s": artifact["homogenization"]["seconds"],
        })
    selected = all(r["observed_and_modeled_target_met"] for r in rows)
    return {
        "scope": "Admitted meshes at observed and quadratic-model 1% targets; no certified continuum error",
        "timing_scope": "One warmup plus three alternating generation samples per workflow; FEA not repeated",
        "target": .01, "selected_pair": selected,
        "generation_ratio_at_target": rows[1]["generation_median_s"] / rows[0]["generation_median_s"]
            if selected else None,
        "repeat_result": str(repeat_path / "generation-repeats.json"),
        "repeat_result_sha256": digest(repeat_path / "generation-repeats.json"),
        "reference_sensitivity": model, "rows": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, nargs=3, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--mmg", type=Path, required=True)
    parser.add_argument("--repeats", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(qualify(args.references, args.native, args.mmg, args.repeats), indent=2))
