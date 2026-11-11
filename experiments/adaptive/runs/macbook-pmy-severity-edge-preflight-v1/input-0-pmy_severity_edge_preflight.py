"""Materialize exactly the frozen job784 severity-only edge schedule.

Start again from job769, change ordering only, and retain every original edge
choice, field criterion, periodic orbit and admission guard. No native meshing.
"""

import argparse
from itertools import combinations
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
TIMEOUT = 300


def materialize(points, cells, edges, records, outcomes, classes, progress):
    import numpy as np
    import pmy_edge_sampling_control as edge

    lookup = {tuple(nodes): i for i, nodes in enumerate(edges)}
    incidence = [[] for _ in edges]
    groups = {}
    for e, nodes in enumerate(edges):
        groups.setdefault(tuple(sorted(map(int, classes[nodes]))), []).append(e)
    for parent, tet in enumerate(cells):
        for offset, (i, j) in enumerate(edge.EDGES):
            edge_id = lookup[tuple(sorted((tet[i], tet[j])))]
            incidence[edge_id].append((parent, offset))
    nodes = points.tolist()
    midpoint = np.full(len(edges), -1, dtype=np.int64)
    selected = np.full(len(cells), -1, dtype=np.int64)
    accepted = [row["orbit"] for row in outcomes if row["reason"] == "accepted"]
    live, minimum_checked = len(cells), 1.
    for completed, orbit in enumerate(accepted):
        row = records[orbit]
        group = row["edges"]
        assert groups[tuple(row["key"])] == group
        incident = [(e, parent, offset) for e in group for parent, offset in incidence[e]]
        parents = [parent for _, parent, _ in incident]
        assert parents == row["incident_tets"] and len(set(parents)) == len(parents)
        assert (selected[parents] < 0).all()
        mids = edge.orbit_midpoints(points, edges, group, classes)
        by_edge = dict(zip(group, mids, strict=True))
        minimum_q, minimum_det = 1., float("inf")
        for e, parent, offset in incident:
            _, det, q = edge.child_metrics(points[cells[parent]], edge.EDGES[offset], by_edge[e])
            minimum_q, minimum_det = min(minimum_q, float(q.min())), min(minimum_det, float(det.min()))
        assert minimum_q >= .1 and minimum_det > 0
        assert minimum_q == row["minimum_child_q"] and minimum_det == row["minimum_child_determinant"]
        assert len(nodes)+len(group) <= w.CAP and live+len(parents) <= w.CAP
        for e, value in zip(group, mids, strict=True):
            assert midpoint[e] < 0
            midpoint[e] = len(nodes)
            nodes.append(value.tolist())
        for e, parent, _ in incident:
            selected[parent] = e
        live += len(parents)
        minimum_checked = min(minimum_checked, minimum_q)
        if completed % 128 == 0 or completed+1 == len(accepted):
            progress((completed+1)/len(accepted), completed+1, len(accepted), len(nodes), live)
    refined = np.asarray(nodes)
    assert len(np.unique(refined, axis=0)) == len(refined)
    np.testing.assert_array_equal(refined[:len(points)], points)
    children, owners = [], []
    for parent, tet in enumerate(cells):
        e = selected[parent]
        if e < 0:
            children.append(tet)
            owners.append(parent)
        else:
            for endpoint in edges[e]:
                child = tet.copy()
                child[child == endpoint] = midpoint[e]
                children.append(child)
                owners.append(parent)
    assert len(children) == live
    return refined, np.asarray(children), np.asarray(owners), midpoint, selected, minimum_checked


def worker(a):
    import numpy as np
    from vtkmodules.vtkCommonCore import vtkSMPTools
    from microgen.shape import surface_functions
    import pmy_edge_sampling_control as edge
    import pmy_edge_causal_audit as causal
    import pmy_second_edge_preflight as second
    import pmy_global_preclip_preflight as preflight
    import pmy_neck_control as neck
    import pmy_residual_probe as residual
    import tpms_family_fixtures as fixtures
    from periodic_trace import periodic_pairs

    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    audit = json.loads((a.audit/"result.json").read_text())
    assert audit["complete"] and not audit["global_candidate_generated"]
    frozen = json.loads((a.audit/"severity-schedule-dry-run.json").read_text())
    original = json.loads((a.original/"result.json").read_text())
    parent = json.loads((a.parent/"result.json").read_text())
    for path in (a.original/"selection.npz", a.original/"edge-plan.npz", a.original/"edge-decisions.json",
                 a.parent/"background-edge-round.npz"):
        assert w.digest(path) == audit["inputs"][str(path)]
    assert frozen["nodes"] == audit["severity_dry_run"]["nodes"] == 16905
    assert frozen["tets"] == audit["severity_dry_run"]["tets"] == 85127
    case = fixtures.case(NAME, a.fixtures)
    specimen = list((a.fixtures/NAME).iterdir())
    for path in [Path(surface_functions.__file__), Path(neck.__file__), *specimen]:
        assert w.digest(path) == original["inputs"][str(path)]
    native_source, archive_files = second.provenance(a)
    files = [Path(__file__), Path(edge.__file__), Path(causal.__file__), Path(second.__file__),
        Path(preflight.__file__), Path(neck.__file__), Path(residual.__file__), Path(fixtures.__file__),
        Path(surface_functions.__file__), a.audit/"result.json", a.audit/"severity-schedule-dry-run.json",
        a.original/"result.json", a.original/"selection.npz", a.original/"edge-plan.npz", a.original/"edge-decisions.json",
        a.parent/"result.json", a.parent/"background-edge-round.npz", a.witnesses/"result.json", *specimen, *archive_files]
    files += [a.native/arm/name for arm in ("baseline", "edge") for name in ("result.json", "input.json")]
    files += [w.HERE/name for name in ("microgen_workflow.py", "scorecard.py", "periodic_trace.py",
        "mmg_quality_localization.py", "pmy_barycenter_control.py", "pmy_interpolation_control.py", "pmy_failure_probe.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs", ".toml", ".lock"):
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, vtk_backend="Sequential", node_cap=w.CAP, tetrahedron_cap=w.CAP,
        native_measurement_provenance=native_source, only_change="Frozen job784 severity-only ordering on job769 background",
        selection="Exact job783 marks and edge choices; no alternative-edge substitution",
        schedule_sha256=w.digest(a.audit/"severity-schedule-dry-run.json"),
        guards="All periodic star children positive, q>=.1, count caps, exact old points and parent volumes",
        heldout_use="Observers only after the schedule is materialized and saved",
        known_estimator_limit="Some held-out affine gaps exceed .004 while containing-tet residual is .001832; ordering alone does not address omitted marks",
        native_generation_performed=False, candidate_generated=False, solve_performed=False))
    with np.load(a.parent/"background-edge-round.npz") as saved:
        points, cells, root_owner = saved["points"], saved["tetrahedra"], saved["original_tet"]
    marks, plan = np.load(a.original/"selection.npz"), np.load(a.original/"edge-plan.npz")
    records = json.loads((a.original/"edge-decisions.json").read_text())
    result = dict(inputs=hashes, admitted=False, candidate_generated=False, native_generation_performed=False,
        native_clipping_performed=False, admission_gates_unchanged=True, native_measurement_provenance=native_source)
    with Job(parts={"schedule-replay": 1, "materialize": 4, "background-audit": 2,
                    "field-and-count-audit": 2, "heldout-observers": 1}) as job:
        with job.part("schedule-replay") as part:
            faces, _, _, boundary = preflight.inventory(cells)
            np.testing.assert_array_equal(faces, marks["faces"])
            _, classes, _ = periodic_pairs(points, boundary, declared_bounds=([-.5]*3, [.5]*3))
            accepted, occupied, replay = causal.schedule(records, [row["orbit"] for row in frozen["outcomes"]], len(points), len(cells))
            assert replay["outcomes"] == frozen["outcomes"] and replay["counts"] == frozen["counts"]
            assert replay["nodes"] == frozen["nodes"] and replay["tets"] == frozen["tets"]
            shutil.copyfile(a.original/"selection.npz", a.run/"selection.npz")
            shutil.copyfile(a.audit/"severity-schedule-dry-run.json", a.run/"frozen-schedule.json")
            result["decisions"] = replay["counts"]
            part.update(1, "Exact frozen schedule and guards reproduced", accepted=int(accepted.sum()), nodes=replay["nodes"], tets=replay["tets"])
        with job.part("materialize") as part:
            refined, tets, owners, midpoint, selected, minimum_child = materialize(
                points, cells, plan["edges"], records, frozen["outcomes"], classes,
                lambda fraction, done, total, nodes, tets: part.update(fraction,
                    "Materializing only frozen accepted periodic stars", accepted=done, total=total, nodes=nodes, tets=tets))
            assert len(refined) == 16905 and len(tets) == 85127
            np.testing.assert_array_equal(selected >= 0, occupied >= 0)
            lookup = {tuple(nodes): i for i, nodes in enumerate(plan["edges"])}
            split_faces = np.array([any(midpoint[lookup[tuple(sorted(e))]] >= 0 for e in combinations(face, 2)) for face in faces])
            np.savez_compressed(a.run/"background-edge-round.npz", points=refined, tetrahedra=tets,
                                original_tet=owners, root_tet=root_owner[owners])
            np.savez_compressed(a.run/"edge-plan.npz", edges=plan["edges"], edge_midpoint=midpoint,
                original_tet_edge=selected, preferred_face_edge=plan["preferred_face_edge"], preferred_tet_edge=plan["preferred_tet_edge"],
                marked_faces=marks["selected_faces"], marked_tets=marks["selected_tets"],
                remaining_faces=marks["selected_faces"] & ~split_faces,
                remaining_tets=marks["selected_tets"] & (selected < 0))
            result["coverage"] = dict(marked_faces_split=int((marks["selected_faces"] & split_faces).sum()),
                remaining_unsplit_faces=int((marks["selected_faces"] & ~split_faces).sum()),
                marked_tets_split=int((marks["selected_tets"] & (selected >= 0)).sum()),
                remaining_unsplit_tets=int((marks["selected_tets"] & (selected < 0)).sum()))
            result["frozen_decision_hashes"] = {name: w.digest(a.run/name) for name in (
                "selection.npz", "frozen-schedule.json", "background-edge-round.npz", "edge-plan.npz")}
            part.update(1, "Materialized counts match frozen dry run exactly", nodes=len(refined), tets=len(tets), minimum_child_q=minimum_child)
        with job.part("background-audit") as part:
            before, original_det, _ = preflight.metrics(points, cells)
            after, determinant, _ = preflight.metrics(refined, tets)
            assert before == parent["background"]["after"]
            error = float(np.max(abs(np.bincount(owners, weights=determinant)-original_det)/original_det))
            assert error < 1e-10
            np.testing.assert_array_equal(refined[:len(points)], points)
            result["background"] = dict(before=before, after=after, maximum_parent_volume_error=error,
                original_points_bitwise_preserved=True, passed=bool(after["quality_floor_met"] and after["within_caps"]))
            top = []
            outcome_by_orbit = {row["orbit"]: row for row in replay["outcomes"]}
            for entity in audit["top_field_entities"]:
                index, kind, orbit = entity["index"], entity["kind"], entity["preferred_orbit"]
                top.append(dict(kind=kind, index=index, severity_rank=entity["severity_rank"], residual=entity["residual"],
                    original_preferred_status=entity["preferred_status"], new_preferred_status=outcome_by_orbit[orbit]["reason"],
                    preferred_quality_feasible=records[orbit]["quality_ok"],
                    entity_split=bool(split_faces[index] if kind == "face" else selected[index] >= 0)))
            result["top_field_entity_coverage"] = top
            result["remaining_quality_blocks"] = dict(
                formerly_quality_rejected_top_entities=sum(row["original_preferred_status"] == "quality" for row in top),
                still_infeasible_preferred_top_entities=sum(not row["preferred_quality_feasible"] for row in top),
                quality_infeasible_top_entities_still_unsplit=sum(not row["preferred_quality_feasible"] and not row["entity_split"] for row in top),
                alternative_edges_not_enabled=True)
            part.update(1, "Verified unchanged quality, volume and exact periodic gates", minimum_q=after["minimum_quality"], unsplit_top=sum(not row["entity_split"] for row in top))
        with job.part("field-and-count-audit") as part:
            before_criterion, before_arrays = edge.criterion(case, points, cells, neck.raw_gradient)
            for key, value in before_arrays.items():
                np.testing.assert_array_equal(value, marks[key])
            after_criterion, after_arrays = edge.criterion(case, refined, tets, neck.raw_gradient)
            np.savez_compressed(a.run/"criterion-observation.npz", **after_arrays)
            result["criterion"] = dict(before=before_criterion, after=after_criterion,
                longest_first_control=original["criterion"]["after"])
            result["tet10_count_risk"] = second.node_risk(case,
                second.count_features(case, points, cells), second.count_features(case, refined, tets), a.native)
            part.update(1, "Measured actual field residuals and fresh empirical count risk", max_face=after_criterion["faces"]["maximum_normalized_residual"], max_tet=after_criterion["tets"]["maximum_normalized_residual"], predicted_nodes=result["tet10_count_risk"]["maximum_model_node_prediction"])
        with job.part("heldout-observers") as part:
            witnesses = json.loads((a.witnesses/"result.json").read_text())["witnesses"]
            source_witnesses = {row["query_index"]: row for row in parent["witnesses"]}
            original_field = residual.AffineBackground(points, cells, np.asarray(case.function(*points.T)))
            refined_field = residual.AffineBackground(refined, tets, np.asarray(case.function(*refined.T)))
            rows = []
            for index, witness in enumerate(witnesses):
                query = np.asarray(witness["probe"])
                old = residual.affine_boundary_distance(original_field, query, case.levels)
                np.testing.assert_allclose(old["nearest_distance"], source_witnesses[witness["query_index"]]["after"]["nearest_distance"], atol=1e-12, rtol=0)
                new = residual.affine_boundary_distance(refined_field, query, case.levels)
                containing = original_field.interpolation(query, witness["wall_level"])["containing_tets"]
                rows.append(dict(query_index=witness["query_index"], probe=witness["probe"], before=old, after=new,
                    containing_tet_split=bool((selected[containing] >= 0).any()),
                    containing_tet_marked=bool(marks["selected_tets"][containing].any()),
                    containing_tet_field_estimate=marks["tet_estimate"][containing].tolist()))
                part.update((index+1)/len(witnesses), "Observed held-out distances without altering marks", completed=index+1, witnesses=len(witnesses), distance=new["nearest_distance"])
            result["witnesses"] = rows
            result["witness_summary"] = dict(count=len(rows),
                within_distance_gate=sum(row["after"]["nearest_is_global_for_affine_representation"] and row["after"]["nearest_distance"] <= .002 for row in rows),
                improved=sum(row["after"]["nearest_distance"] is not None and row["after"]["nearest_distance"] < row["before"]["nearest_distance"]-1e-12 for row in rows),
                worsened=sum(row["after"]["nearest_distance"] is not None and row["after"]["nearest_distance"] > row["before"]["nearest_distance"]+1e-12 for row in rows),
                worst_after=max(row["after"]["distance_lower_bound"] for row in rows))
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    assert all(w.digest(a.run/name) == expected for name, expected in result["frozen_decision_hashes"].items())
    result.update(complete=True, input_hashes_verified=True,
        interpretation="Severity-only background control, not a native candidate. Quality-blocked preferred edges and omitted geometric marks remain visible. No claim that ordering resolves geometry globally.")
    w.save(a.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--audit", type=Path, default=w.HERE/"runs/macbook-pmy-edge-causal-audit-v1")
    parser.add_argument("--original", type=Path, default=w.HERE/"runs/macbook-pmy-second-edge-preflight-v1")
    parser.add_argument("--parent", type=Path, default=w.HERE/"runs/macbook-pmy-edge-sampling-v1")
    parser.add_argument("--native", type=Path, default=w.HERE/"runs/macbook-pmy-rollback-native-v1")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--witnesses", type=Path, default=w.HERE/"runs/macbook-pmy-residual-probe-v1")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "audit", "original", "parent", "native", "fixtures", "witnesses"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-pmy_severity_edge_preflight.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **kw: None), timeout=TIMEOUT)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    findings = json.loads((args.run/"result.json").read_text())
    after = findings["background"]["after"]
    result("PMY severity background", nodes=after["nodes"], tets=after["tets"], minimum_q=after["minimum_quality"], passed=findings["background"]["passed"])
    result("PMY severity field residual", max_face=findings["criterion"]["after"]["faces"]["maximum_normalized_residual"], max_tet=findings["criterion"]["after"]["tets"]["maximum_normalized_residual"])
    result("PMY remaining quality blocks", **findings["remaining_quality_blocks"])
    result("PMY held-out affine witnesses", **findings["witness_summary"])


if __name__ == "__main__":
    main()
