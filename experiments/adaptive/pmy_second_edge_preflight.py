"""One second field-driven PMY edge pass, background preflight only.

Recompute the unchanged field-only criterion on job769's exact background.
Freeze the same deterministic quality-guarded periodic edge-star decisions before
observing held-out queries. Empirical Tet10 count predictions are not bounds.
"""

import argparse
from itertools import combinations
import json
import math
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
TIMEOUT = 300


def count_features(case, points, cells):
    import numpy as np

    edges = np.unique(np.sort(cells[:, np.asarray(list(combinations(range(4), 2)))].reshape(-1, 2), axis=1), axis=0)
    values = np.asarray(case.function(*points.T))[cells]
    below, above = (values < case.levels[0]).all(1), (values > case.levels[1]).all(1)
    inside = ((values >= case.levels[0]) & (values <= case.levels[1])).all(1)
    return dict(nodes=len(points), tets=len(cells), unique_edges=len(edges),
        full_background_tet10_nodes=len(points)+len(edges), active_tets=int((~(below | above)).sum()),
        wholly_in_band_tets=int(inside.sum()), crossing_tets=int((~(below | above | inside)).sum()))


def node_risk(case, before, after, native):
    """Four explicitly empirical count projections, never a resource bound."""
    import numpy as np

    current = json.loads((native/"edge/result.json").read_text())
    previous = json.loads((native/"baseline/result.json").read_text())
    old_input = json.loads((native/"baseline/input.json").read_text())
    older = count_features(case, np.asarray(old_input["points"]), np.asarray(old_input["tetrahedra"]))
    models = []
    for feature in ("active_tets", "tets", "full_background_tet10_nodes"):
        factor = after[feature]/before[feature]
        models.append(dict(model="scale_observed_native_by_"+feature,
            predicted_nodes=math.ceil(current["nodes"]*factor), predicted_tets=math.ceil(current["elements"]*factor), factor=factor))
    growth = before["tets"]-older["tets"]
    assert growth > 0
    node_slope = (current["nodes"]-previous["nodes"])/growth
    tet_slope = (current["elements"]-previous["elements"])/growth
    models.append(dict(model="extrapolate_previous_increment_per_background_tet",
        predicted_nodes=math.ceil(current["nodes"]+node_slope*(after["tets"]-before["tets"])),
        predicted_tets=math.ceil(current["elements"]+tet_slope*(after["tets"]-before["tets"])),
        observed_node_increment_per_background_tet=node_slope, observed_tet_increment_per_background_tet=tet_slope))
    for row in models:
        row["node_headroom_to_cap"] = w.CAP-row["predicted_nodes"]
        row["tet_headroom_to_cap"] = w.CAP-row["predicted_tets"]
    return dict(observed_prior_background=older, observed_current_background=before, proposed_background=after,
        observed_prior_native=dict(nodes=previous["nodes"], tets=previous["elements"]),
        observed_current_native=dict(nodes=current["nodes"], tets=current["elements"], status=current["status"],
            forward=current["geometry"]["sampled_forward_max"], reverse=current["geometry"]["sampled_reverse_max"]),
        models=models, minimum_model_node_prediction=min(row["predicted_nodes"] for row in models),
        maximum_model_node_prediction=max(row["predicted_nodes"] for row in models),
        any_prediction_over_cap=any(max(row["predicted_nodes"], row["predicted_tets"]) > w.CAP for row in models),
        is_bound=False, clipped_tet10_admission_unassessed=True,
        limitations="Only two previous backgrounds; clipping, snap restoration and quality repair may change count growth. The model spread is not an uncertainty interval. Every future native result still requires actual 250000 node/tet caps.")


def provenance(a):
    """Resolve the native measurement's source from its immutable archive."""
    archive = w.ROOT/".local-work/macbook-pmy-background-native-v3"
    receipt = json.loads((archive/"receipt.json").read_text())
    assert w.digest(archive/"pmy_background") == receipt["binary_sha256"]
    records, files = [], [archive/"receipt.json", archive/"pmy_background"]
    for relative, expected in receipt["inputs"].items():
        saved = archive/"sources"/relative
        assert w.digest(saved) == expected, f"Frozen native source archive mismatch: {relative}"
        working = w.ROOT/relative
        current_hash = w.digest(working) if working.exists() else None
        records.append(dict(original_path=relative, resolved_source=str(saved), sha256=expected,
            working_copy_sha256_at_manifest=current_hash, working_copy_differs=current_hash != expected))
        files.append(saved)
    return dict(binary_sha256=receipt["binary_sha256"], receipt_sha256=w.digest(archive/"receipt.json"),
        source_resolution="Frozen archive for every native-v3 source, never current modified source",
        sources=records), files


def worker(a):
    import numpy as np
    from vtkmodules.vtkCommonCore import vtkSMPTools
    from microgen.shape import surface_functions
    import pmy_edge_sampling_control as edge
    import pmy_global_preclip_preflight as preflight
    import pmy_interpolation_control as interpolation
    import pmy_neck_control as neck
    import pmy_residual_probe as residual
    import tpms_family_fixtures as fixtures
    from periodic_trace import periodic_pairs

    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    prior = json.loads((a.source/"result.json").read_text())
    assert prior["complete"] and prior["background"]["passed"] and not prior["native_clipping_performed"]
    for name, expected in prior["frozen_decision_hashes"].items():
        assert w.digest(a.source/name) == expected
    # Query content in the saved source report is never passed to either rule.
    previous_criterion = prior["criterion"]["after"]
    prior_metrics, prior_hashes = prior["background"]["after"], prior["inputs"]
    del prior
    trace_proof = json.loads((a.trace_proof/"result.json").read_text())
    assert trace_proof["complete"] and trace_proof["trace_contract_verified"]
    case = fixtures.case(NAME, a.fixtures)
    specimen = list((a.fixtures/NAME).iterdir())
    for path in [Path(surface_functions.__file__), Path(neck.__file__), *specimen]:
        assert w.digest(path) == prior_hashes[str(path)]
    native_source, archive_files = provenance(a)
    native_rows = []
    for arm in ("baseline", "edge"):
        row = json.loads((a.native/arm/"result.json").read_text())
        assert row["complete"] and row["status"] == "admission-rejected"
        assert w.digest(a.native/arm/"mesh.npz") == row["mesh_sha256"]
        native_rows += [a.native/arm/name for name in ("result.json", "input.json", "mesh.npz")]
    files = [Path(__file__), Path(edge.__file__), Path(preflight.__file__), Path(interpolation.__file__),
        Path(neck.__file__), Path(residual.__file__), Path(fixtures.__file__), Path(surface_functions.__file__),
        a.source/"result.json", a.source/"manifest.json", a.source/"background-edge-round.npz",
        a.source/"criterion-observation.npz", a.native/"manifest.json", a.native/"baseline-equivalence.json",
        a.trace_proof/"result.json", a.witnesses/"result.json", a.fixtures/"result.json", *specimen, *native_rows, *archive_files]
    files += [w.HERE/name for name in ("microgen_workflow.py", "scorecard.py", "periodic_trace.py",
        "mmg_quality_localization.py", "pmy_barycenter_control.py", "pmy_failure_probe.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs", ".toml", ".lock"):
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, vtk_backend="Sequential", node_cap=w.CAP, tetrahedron_cap=w.CAP,
        native_measurement_provenance=native_source, source_background=str(a.source),
        selection="Recompute unchanged class mismatch plus .002 normalized affine-wall residual; constrained cap gradient",
        subdivision="Same deterministic longest-edge periodic-orbit disjoint-star pass, every child positive and q>=.1",
        witness_use="Observer only after selections, edge decisions and candidate background are saved; never passed to either selection function",
        count_prediction="Empirical Tet10 node estimates only, not a bound or admission",
        native_clipping_performed=False, candidate_generated=False, solve_performed=False))
    with np.load(a.source/"background-edge-round.npz") as saved:
        points, cells, root_owner = saved["points"], saved["tetrahedra"], saved["original_tet"]
    result = dict(inputs=hashes, admitted=False, candidate_generated=False, native_clipping_performed=False,
        admission_gates_unchanged=True, source=str(a.source), native_measurement_provenance=native_source)
    with Job(parts={"controls": 1, "freeze-criterion": 2, "edge-stars": 5, "background-audit": 2,
                    "count-risk": 1, "residual-observers": 2}) as job:
        with job.part("controls") as part:
            result["controls"] = edge.controls()
            before, original_det, boundary = preflight.metrics(points, cells)
            assert before == prior_metrics
            _, classes, _ = periodic_pairs(points, boundary, declared_bounds=([-.5]*3, [.5]*3))
            part.update(1, "Verified exact source background and unchanged edge controls", nodes=len(points), tets=len(cells), minimum_q=before["minimum_quality"])
        with job.part("freeze-criterion") as part:
            criterion_before, selection = edge.criterion(case, points, cells, neck.raw_gradient)
            assert criterion_before == previous_criterion
            with np.load(a.source/"criterion-observation.npz") as previous:
                for key, value in selection.items():
                    np.testing.assert_array_equal(value, previous[key])
            np.savez_compressed(a.run/"selection.npz", **selection)
            result["selection"] = criterion_before
            result["selection_sha256"] = w.digest(a.run/"selection.npz")
            part.update(1, "Frozen recomputed same-rule marks; exact prior observer replay", faces=int(selection["selected_faces"].sum()), tets=int(selection["selected_tets"].sum()))
        with job.part("edge-stars") as part:
            output = edge.edge_pass(points, cells, selection["faces"], selection["selected_faces"], selection["selected_tets"], classes,
                lambda fraction, completed, total, counts, nodes, tets: part.update(fraction,
                    "Testing second-pass complete periodic edge stars", completed=completed, orbits=total,
                    accepted=counts["accepted"], overlap=counts["overlap"], quality_rejected=counts["quality"], nodes=nodes, tets=tets))
            refined, tets, owners, records, failed, decisions, coverage, arrays = output
            w.save(a.run/"edge-decisions.json", records)
            w.save(a.run/"quality-failure-proof.json", failed)
            np.savez_compressed(a.run/"edge-plan.npz", **arrays)
            np.savez_compressed(a.run/"background-edge-round.npz", points=refined, tetrahedra=tets,
                                original_tet=owners, root_tet=root_owner[owners])
            result.update(decisions=decisions, coverage=coverage,
                all_quality_infeasible_orbits=sum(not row["quality_ok"] for row in records),
                overlapping_quality_failures=sum(not row["quality_ok"] and row["overlap"] for row in records))
            result["frozen_decision_hashes"] = {name: w.digest(a.run/name) for name in (
                "selection.npz", "edge-decisions.json", "quality-failure-proof.json", "edge-plan.npz", "background-edge-round.npz")}
            part.update(1, "Saved complete second-pass decisions and failure proofs", accepted=decisions["accepted"], nodes=len(refined), tets=len(tets))
        with job.part("background-audit") as part:
            after, determinant, _ = preflight.metrics(refined, tets)
            error = float(np.max(abs(np.bincount(owners, weights=determinant)-original_det)/original_det))
            assert error < 1e-10
            np.testing.assert_array_equal(refined[:len(points)], points)
            result["background"] = dict(before=before, after=after, maximum_parent_volume_error=error,
                original_points_bitwise_preserved=True, passed=bool(after["quality_floor_met"] and after["within_caps"]))
            part.update(1, "Verified quality, original points, parent volumes and exact periodic facets", minimum_q=after["minimum_quality"], nodes=len(refined), tets=len(tets))
        with job.part("count-risk") as part:
            before_counts = count_features(case, points, cells)
            after_counts = count_features(case, refined, tets)
            result["tet10_count_risk"] = node_risk(case, before_counts, after_counts, a.native)
            part.update(1, "Computed explicit empirical clipped-Tet10 count projections", prediction_min=result["tet10_count_risk"]["minimum_model_node_prediction"], prediction_max=result["tet10_count_risk"]["maximum_model_node_prediction"], cap=w.CAP)
        with job.part("residual-observers") as part:
            criterion_after, after_arrays = edge.criterion(case, refined, tets, neck.raw_gradient)
            np.savez_compressed(a.run/"criterion-observation.npz", **after_arrays)
            result["criterion"] = dict(before=criterion_before, after=criterion_after)
            witnesses = json.loads((a.witnesses/"result.json").read_text())["witnesses"]
            source_witnesses = {row["query_index"]: row for row in json.loads((a.source/"result.json").read_text())["witnesses"]}
            original = residual.AffineBackground(points, cells, np.asarray(case.function(*points.T)))
            updated = residual.AffineBackground(refined, tets, np.asarray(case.function(*refined.T)))
            rows = []
            for index, witness in enumerate(witnesses):
                query = np.asarray(witness["probe"])
                old = residual.affine_boundary_distance(original, query, case.levels)
                saved = source_witnesses[witness["query_index"]]["after"]
                np.testing.assert_allclose(old["nearest_distance"], saved["nearest_distance"], atol=1e-12, rtol=0)
                new = residual.affine_boundary_distance(updated, query, case.levels)
                containing = original.interpolation(query, witness["wall_level"])["containing_tets"]
                rows.append(dict(query_index=witness["query_index"], probe=witness["probe"], before=old, after=new,
                    containing_tet_split=bool((arrays["original_tet_edge"][containing] >= 0).any())))
                part.update((index+1)/len(witnesses), "Observed unchanged held-out affine queries after all choices", completed=index+1, witnesses=len(witnesses), distance=new["nearest_distance"])
            result["witnesses"] = rows
            result["witness_summary"] = dict(count=len(rows),
                within_distance_gate=sum(row["after"]["nearest_is_global_for_affine_representation"] and row["after"]["nearest_distance"] <= .002 for row in rows),
                improved=sum(row["after"]["nearest_distance"] is not None and row["after"]["nearest_distance"] < row["before"]["nearest_distance"]-1e-12 for row in rows),
                worsened=sum(row["after"]["nearest_distance"] is not None and row["after"]["nearest_distance"] > row["before"]["nearest_distance"]+1e-12 for row in rows),
                worst_after=max(row["after"]["distance_lower_bound"] for row in rows),
                scope="Reverse affine observers only; improvements cannot establish native forward geometry or final admission")
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    assert all(w.digest(a.run/name) == expected for name, expected in result["frozen_decision_hashes"].items())
    result.update(complete=True, input_hashes_verified=True,
        interpretation="One second background edge pass, no native generation. Existing native forward error remains .01449. Parent must inspect background gates/count risk before deciding on a new native control.")
    w.save(a.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-edge-sampling-v1")
    parser.add_argument("--native", type=Path, default=w.HERE/"runs/macbook-pmy-rollback-native-v1")
    parser.add_argument("--trace-proof", type=Path, default=w.HERE/"runs/macbook-pmy-saved-trace-replay-v1")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--witnesses", type=Path, default=w.HERE/"runs/macbook-pmy-residual-probe-v1")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "source", "native", "trace_proof", "fixtures", "witnesses"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-pmy_second_edge_preflight.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **kw: None), timeout=TIMEOUT)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    findings = json.loads((args.run/"result.json").read_text())
    result("PMY second edge decisions", **findings["decisions"])
    after = findings["background"]["after"]
    result("PMY second background preflight", nodes=after["nodes"], tets=after["tets"], minimum_q=after["minimum_quality"], passed=findings["background"]["passed"])
    result("PMY empirical Tet10 node risk", predicted_min=findings["tet10_count_risk"]["minimum_model_node_prediction"], predicted_max=findings["tet10_count_risk"]["maximum_model_node_prediction"], is_bound=False)
    result("PMY held-out affine witnesses", **{key: value for key, value in findings["witness_summary"].items() if key != "scope"})


if __name__ == "__main__":
    main()
