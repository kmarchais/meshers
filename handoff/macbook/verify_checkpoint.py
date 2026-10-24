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
    seed_replays, failed_runs, coupled_replays = [], [], []
    for run in sorted((workflow.HERE / "runs").glob("macbook-*")):
        if not run.is_dir() or not (run / "manifest.json").exists():
            continue
        manifest = json.loads((run / "manifest.json").read_text())
        source_ref = "2779c4a" if "numerical_thread_policy" in manifest else "2ad9895"
        inputs = dict(manifest["inputs"])
        if (run / "analysis-inputs.json").exists():
            inputs.update(json.loads((run / "analysis-inputs.json").read_text()))
        for name, expected in inputs.items():
            # Frozen Linux prefixes are resolved for reading, never rewritten on disk.
            path = Path(name.replace("/home/kmarchais/dev/", str(ROOT.parent) + "/", 1))
            path = path if path.is_absolute() else ROOT / path
            snapshots = [candidate for candidate in run.glob("*.py")
                if path.suffix == ".py" and workflow.digest(candidate) == expected]
            if workflow.digest(path) != expected and snapshots:
                archived += 1
            else:
                archived += verify_hash(path, expected, run, source_ref)
            checked += 1
        if (run / "worker.json").exists() and json.loads((run / "worker.json").read_text())["returncode"] != 0:
            failed_runs.append(run.name)
        if (run / "result.json").exists():
            result = json.loads((run / "result.json").read_text())
            if result.get("unknowns") == 4 and "attempts" in result:
                from cap_quality_bound import bounds, caps
                from mmg_quality_localization import quality
                mesh = np.load(Path(result["source"]) / "mesh.npz")
                original, cells = mesh["points"], mesh["tetrahedra"]
                nodes = np.concatenate(result["orbits"])
                affected = cells[np.isin(cells, nodes).any(axis=1)]
                faces, _ = caps(mesh)
                faces = faces[np.isin(faces, nodes).any(axis=1)]
                for attempt in result["attempts"]:
                    moved = original.copy()
                    moved[nodes] = attempt["candidate_coordinates"]
                    cap_quality = bounds(moved[faces])[0]
                    volume_quality = quality(moved, affected)
                    determinant = np.linalg.det(moved[affected[:, 1:]] - moved[affected[:, :1]])
                    np.testing.assert_allclose(cap_quality.min(), attempt["minimum_cap_quality"], rtol=1e-12)
                    np.testing.assert_allclose(volume_quality.min(), attempt["minimum_volume_quality"], rtol=1e-12)
                    np.testing.assert_allclose(determinant.min(), attempt["minimum_volume_determinant"], rtol=1e-12)
                    for orbit in result["orbits"]:
                        np.testing.assert_allclose(moved[orbit[1]] - moved[orbit[0]],
                            original[orbit[1]] - original[orbit[0]], atol=1e-15)
                    assert not attempt["feasible"] and cap_quality.min() < result["target_cap_quality"]
                assert not result["candidate_found"] and not result["infeasibility_proved"]
                coupled_replays.append({"run": run.name, "attempts": len(result["attempts"]),
                    "incident_cap_triangles": len(faces), "incident_volume_tets": len(affected)})
        if (run / "summary.json").exists():
            summary = json.loads((run / "summary.json").read_text())
            if summary.get("sampled_field_identical"):
                workflow.imports()
                from cap_quality_bound import bounds, caps
                from facet_incidence_audit import audit
                from periodic import periodic_pairs
                rows = summary["rows"]
                assert rows["unmodified"]["sampled_inputs"] == rows["pretriangulated"]["sampled_inputs"]
                for mode, row in rows.items():
                    mesh = np.load(run / mode / "mesh.npz")
                    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
                    assert len(points) == row["nodes"] and len(cells) == row["elements"]
                    assert audit(points, cells, surface) == row["facet_consistency"]
                    _, triangles = caps(mesh)
                    assert int((bounds(triangles)[0] < .1).sum()) == row["cap_faces_below_quality_gate"]
                    _, _, periodic = periodic_pairs(points, surface)
                    assert periodic == row["periodic"]
                    assert row["mmg_calls"] == [] and not row["stiffness_solve_performed"]
                    seed_replays.append({"run": run.name, "mode": mode, "nodes": len(points),
                        "elements": len(cells), "mesh_sha256": workflow.digest(run / mode / "mesh.npz")})
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
        "worst_retained_model_scores": scores, "seed_artifact_replays": seed_replays,
        "failed_runs_preserved": failed_runs, "coupled_candidate_replays": coupled_replays, "reference_promoted": False})
    print(f"Verified {files} transferred files, {checked} input hashes and {comparisons} full model comparisons")


if __name__ == "__main__":
    main()
