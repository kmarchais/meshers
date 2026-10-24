"""Recheck transferred evidence, measured source versions, and full model scores."""

import hashlib
import json
from pathlib import Path
import sys
import subprocess
from functools import lru_cache

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow
from verify_schwarz_tet10_repeated import verify_hash

import numpy as np
from scipy.linalg import eigvalsh


@lru_cache(maxsize=None)
def committed_python_ref(relative, expected, legacy_ref):
    """Find exact measured Python bytes on this branch; never substitute data."""
    assert relative.endswith(".py")
    tried = set()
    def matches(ref):
        if ref in tried:
            return False
        tried.add(ref)
        result = subprocess.run(["git", "show", f"{ref}:{relative}"], cwd=ROOT, capture_output=True)
        return result.returncode == 0 and hashlib.sha256(result.stdout).hexdigest() == expected
    for ref in ("HEAD", legacy_ref):
        if matches(ref):
            return ref
    history = subprocess.check_output(["git", "log", "--format=%H", "--", relative], cwd=ROOT).decode().splitlines()
    return next((ref for ref in history if matches(ref)), None)


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
    seed_replays, failed_runs, coupled_replays, geometry_replays, workflow_replays, hp_replays = [], [], [], [], [], []
    for run in sorted((workflow.HERE / "runs").glob("macbook-*")):
        if not run.is_dir() or not (run / "manifest.json").exists():
            continue
        manifest = json.loads((run / "manifest.json").read_text())
        source_ref = "2779c4a" if "numerical_thread_policy" in manifest else "2ad9895"
        inputs = dict(manifest["inputs"])
        inputs.update(manifest.get("source_sha256", {}))
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
            elif workflow.digest(path) != expected and path.suffix == ".py" and path.resolve().is_relative_to(ROOT):
                relative = path.resolve().relative_to(ROOT).as_posix()
                if committed_python_ref(relative, expected, source_ref) is not None:
                    archived += 1
                else:
                    archived += verify_hash(path, expected, run, source_ref)
            else:
                archived += verify_hash(path, expected, run, source_ref)
            checked += 1
        if (run / "worker.json").exists() and json.loads((run / "worker.json").read_text())["returncode"] != 0:
            failed_runs.append(run.name)
        if (run / "result.json").exists():
            result = json.loads((run / "result.json").read_text())
            if "complete_field_reduced_residuals" in result:
                assert result["passed"] and result["retained_model_count"] == 24
                assert max(result["auxiliary_reduced_residuals"]) < 2e-9
                assert max(result["complete_field_reduced_residuals"]) < 2e-9
                assert result["mode_equilibrium_relative_residual"] < 2e-9
                assert result["positive_rule_energy_difference"] < 1e-4
                tensor = np.asarray(result["rules"][-1]["energy_tensor"])
                assert result["rules"][-1]["points"] == 216
                errors = []
                for model in result["models"]:
                    path = Path(model["source"])
                    reference = np.asarray(json.loads(path.read_text())["models"][model["index"]]["stiffness"])
                    values = eigvalsh(tensor - reference, reference)
                    np.testing.assert_allclose(values, model["signed_relative_energy_eigenvalues"], atol=1e-12)
                    errors.append(float(abs(values).max()))
                assert len(errors) == 24
                np.testing.assert_allclose(max(errors), result["worst_retained_model_disagreement"], atol=1e-12)
                assert result["passes_target_all_models"] == (max(errors) <= .01)
                assert not result["reference_promoted"] and not result["native_candidate_promoted"]
                comparisons += len(errors)
                scores[run.name] = max(errors)
                hp_replays.append({"run": run.name, "model_count": len(errors), "worst_disagreement": max(errors)})
            if "worst_reverse_samples" in result:
                import scorecard
                source_mesh = np.load(Path(result["source"]) / "mesh.npz")
                queries = np.asarray([row["query"] for row in result["worst_reverse_samples"]])
                distances = scorecard.distance_to_mesh(source_mesh["points"], source_mesh["surface"], queries)
                np.testing.assert_allclose(distances, [row["distance"] for row in result["worst_reverse_samples"]], atol=1e-14)
                np.testing.assert_allclose(distances.max(), result["reverse_maximum"], atol=1e-14)
                geometry_replays.append({"run": run.name, "reverse_witnesses": len(queries),
                    "maximum_distance": float(distances.max())})
            if "repair" in result and (run / "mesh.npz").exists():
                from mmg_quality_localization import quality
                original = np.load(Path(result["source"]) / "mesh.npz")
                repaired = np.load(run / "mesh.npz")
                for key in ("tetrahedra", "surface"):
                    np.testing.assert_array_equal(original[key], repaired[key])
                boundary = np.unique(original["surface"])
                np.testing.assert_array_equal(original["points"][boundary], repaired["points"][boundary])
                cells, points = repaired["tetrahedra"], repaired["points"]
                measured = quality(points, cells)
                signed = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])
                assert measured.min() >= .1 and signed.min() > 0
                np.testing.assert_allclose(measured.min(), result["repair"]["final_minimum_quality"], atol=1e-14)
                geometry_replays.append({"run": run.name, "minimum_quality": float(measured.min()),
                    "unchanged_boundary_and_connectivity": True})
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
            if summary.get("case") == "schwarz-solid" and summary.get("live_workflow_timing"):
                from verify_schwarz_tet10_repeated import model_errors
                for path, expected in summary["artifact_sha256"].items():
                    assert workflow.digest(Path(path)) == expected
                baseline = json.loads((workflow.HERE / "runs/g3-microgen-defaults-v4/summary.json").read_text())
                levels = {row["cells"]: row for row in baseline["rows"]
                    if row["case"] == "schwarz-solid" and row["route"] == "meshers-curved"}
                coarse, fine = (np.asarray(levels[count]["stiffness"]) for count in (16, 20))
                ratios = [16 / 20, (levels[16]["elements"] / levels[20]["elements"])**(1 / 3)]
                parent = json.loads((run / "microgen-parent/result.json").read_text())
                rows, ranges = {}, {}
                for arm in ("microgen", "meshers"):
                    row = json.loads((run / f"{arm}-tet10/result.json").read_text())
                    rows[arm] = row
                    assert row["qualified"] and row["geometry"]["geometry_admitted"]
                    assert row["nodes"] <= 250000 and row["elements"] <= 250000
                    assert len(row["solver_records"]) == 6
                    assert max(load["relative_residual"] for load in row["solver_records"]) < 2e-9
                    tensor = np.asarray(row["stiffness"])
                    errors = model_errors((tensor + tensor.T) / 2, coarse, fine, ratios)
                    reference = json.loads((run / f"{arm}-reference.json").read_text())
                    np.testing.assert_allclose([min(errors), max(errors)], reference["declared_order_difference_range"], atol=1e-12)
                    assert max(errors) <= .001
                    ranges[arm] = [min(errors), max(errors)]
                    field_audit = json.loads((run / f"audits/{arm}-tet10.json").read_text())
                    assert max(field_audit["reaction_relative_error"], field_audit["hill_mandel_relative_error"]) < 1e-7
                    assert field_audit["periodic_jump_max"] < 1e-9
                repair = json.loads((run / "microgen-repair/result.json").read_text()) if summary["explicit_interior_repair"] else None
                totals = {"microgen": parent["total_workflow_s"] + rows["microgen"]["total_saved_mesh_study_s"],
                    "meshers": rows["meshers"]["total_workflow_s"]}
                if repair:
                    totals["microgen"] += repair["total_saved_mesh_study_s"]
                assert totals == summary["science_audited_workflow_s"]
                assert summary["both_pass_declared_reference_models"] and summary["timing_clock_consistent"]
                np.testing.assert_allclose(totals["microgen"] / totals["meshers"],
                    summary["single_pair_workflow_ratio_microgen_over_meshers"], atol=1e-12)
                workflow_replays.append({"run": run.name, "model_difference_ranges": ranges,
                    "science_audited_workflow_s": totals, "repetitions": 1,
                    "explicit_interior_repair": repair is not None})
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
    workflow.save(output, {"passed": True, "verifier_sha256": workflow.digest(Path(__file__)),
        "repository_head_at_verification": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "transferred_files_checked": files,
        "transferred_bytes_checked": total_bytes, "mac_input_hashes_checked": checked,
        "sources_verified_from_frozen_driver_or_git": archived, "mixed_load_comparisons": comparisons,
        "worst_retained_model_scores": scores, "seed_artifact_replays": seed_replays,
        "failed_runs_preserved": failed_runs, "coupled_candidate_replays": coupled_replays,
        "geometry_witness_replays": geometry_replays, "same_host_workflow_replays": workflow_replays,
        "coupled_modal_field_replays": hp_replays,
        "reference_promoted": False})
    print(f"Verified {files} transferred files, {checked} input hashes and {comparisons} full model comparisons")


if __name__ == "__main__":
    main()
