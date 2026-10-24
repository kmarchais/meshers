"""Recheck transferred evidence, measured source versions, and full model scores."""

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow
from verify_schwarz_tet10_repeated import verify_hash

import numpy as np
from scipy.linalg import eigvalsh


def main():
    output = Path(sys.argv[1])
    assert not output.exists(), "Keep verification receipts immutable"
    transferred = json.loads((ROOT.parent / "meshers-handoff-20261024/data-manifest.json").read_text())
    files, total_bytes = 0, 0
    for repository in transferred["repositories"]:
        for item in repository["files"]:
            path = ROOT.parent / repository["name"] / item["path"]
            with path.open("rb") as stream:
                assert hashlib.file_digest(stream, "sha256").hexdigest() == item["sha256"], str(path)
            assert path.stat().st_size == item["bytes"]
            files += 1
            total_bytes += item["bytes"]
    checked, archived, comparisons, scores = 0, 0, 0, {}
    for run in sorted((workflow.HERE / "runs").glob("macbook-*")):
        if not run.is_dir() or not (run / "manifest.json").exists():
            continue
        manifest = json.loads((run / "manifest.json").read_text())
        source_ref = "2779c4a" if "numerical_thread_policy" in manifest else "2ad9895"
        for name, expected in manifest["inputs"].items():
            # Frozen Linux prefixes are resolved for reading, never rewritten on disk.
            path = Path(name.replace("/home/kmarchais/dev/", str(ROOT.parent) + "/", 1))
            path = path if path.is_absolute() else ROOT / path
            archived += verify_hash(path, expected, run, source_ref)
            checked += 1
        if not (run / "assessment.json").exists():
            continue
        assessment = json.loads((run / "assessment.json").read_text())
        row = json.loads((run / "result.json").read_text())
        tensor = np.asarray(row["stiffness"])
        assert row["completed_loads"] == 6 and assessment["retained_model_count"] == 24
        assert not assessment["reference_promoted"] and len(assessment["models"]) == 24
        errors = []
        for model in assessment["models"]:
            source = Path(model["source"].replace("/home/kmarchais/dev/", str(ROOT.parent) + "/", 1))
            source = source if source.is_absolute() else ROOT / source
            matrix = np.asarray(json.loads(source.read_text())["models"][model["index"]]["stiffness"])
            error = float(abs(eigvalsh(tensor - matrix, matrix)).max())
            assert abs(error - model["difference"]["worst_relative_strain_energy_difference"]) < 1e-12
            errors.append(error)
        assert abs(max(errors) - assessment["worst_retained_model_disagreement"]) < 1e-12
        assert assessment["passes_target_all_models"] == (max(errors) <= .01)
        comparisons += len(errors)
        scores[run.name] = max(errors)
    workflow.save(output, {"passed": True, "transferred_files_checked": files,
        "transferred_bytes_checked": total_bytes, "mac_input_hashes_checked": checked,
        "sources_verified_from_frozen_driver_or_git": archived, "mixed_load_comparisons": comparisons,
        "worst_retained_model_scores": scores, "reference_promoted": False})
    print(f"Verified {files} transferred files, {checked} input hashes and {comparisons} full model comparisons")


if __name__ == "__main__":
    main()
