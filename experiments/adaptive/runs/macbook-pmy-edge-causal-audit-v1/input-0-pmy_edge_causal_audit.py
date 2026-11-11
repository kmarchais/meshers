"""Read-only audit of stalled PMY field-error refinement after job783.

Inspect the fixed top 32 marked faces and 32 marked tetrahedra by field residual.
Replay existing decisions, dry-run severity ordering of the same candidate orbits,
and measure hypothetical local alternative-edge children. No global background
candidate or native mesh is generated. Held-out queries are final observers only.
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
TOP_COUNT = 32
TIMEOUT = 300


def schedule(records, order, nodes, tets):
    import numpy as np

    occupied = np.full(tets, -1, dtype=np.int64)
    selected = np.zeros(len(records), dtype=bool)
    counts = dict(accepted=0, overlap=0, quality=0, budget=0)
    outcomes = []
    for orbit in order:
        row = records[orbit]
        parents = np.asarray(row["incident_tets"], dtype=np.int64)
        overlap = len(np.unique(parents)) != len(parents) or bool((occupied[parents] >= 0).any())
        within = nodes+len(row["edges"]) <= w.CAP and tets+len(parents) <= w.CAP
        reason = "overlap" if overlap else "quality" if not row["quality_ok"] else "budget" if not within else "accepted"
        counts[reason] += 1
        blockers = sorted(set(map(int, occupied[parents][occupied[parents] >= 0])))
        outcomes.append(dict(orbit=int(orbit), reason=reason, earlier_blocking_orbits=blockers))
        if reason == "accepted":
            selected[orbit] = True
            occupied[parents] = orbit
            nodes += len(row["edges"])
            tets += len(parents)
    return selected, occupied, dict(counts=counts, nodes=nodes, tets=tets, outcomes=outcomes)


def measure(case, points, entities, *, faces):
    import numpy as np
    from pmy_interpolation_control import wall_estimates
    from pmy_neck_control import raw_gradient
    from pmy_barycenter_control import classify

    estimates = wall_estimates(case, points, entities, raw_gradient, faces=faces)
    actual = case.function(*points[entities].mean(1).T)
    affine = case.function(*points.T)[entities].mean(1)
    return dict(maximum=float(estimates["maximum"].max(initial=0)),
        geometry_selected=int(estimates["selected"].sum()), unresolved=int(estimates["unresolved"].sum()),
        active=int(estimates["active"].sum()),
        class_mismatches=int((classify(actual, case.levels) != classify(affine, case.levels)).sum()))


def local_alternative(case, points, cells, edges, group, incidence, classes):
    import numpy as np
    import pmy_edge_sampling_control as edge

    incident = [(e, parent, offset) for e in group for parent, offset in incidence[e]]
    parents = [parent for _, parent, _ in incident]
    mids = edge.orbit_midpoints(points, edges, group, classes)
    mid_by_edge = dict(zip(group, mids, strict=True))
    minimum_q, minimum_det = 1., float("inf")
    for e, parent, offset in incident:
        _, determinant, q = edge.child_metrics(points[cells[parent]], edge.EDGES[offset], mid_by_edge[e])
        minimum_q = min(minimum_q, float(q.min()))
        minimum_det = min(minimum_det, float(determinant.min()))
    summary = dict(physical_edges=len(group), complete_star_tets=len(parents),
        minimum_child_q=minimum_q, minimum_child_determinant=minimum_det,
        quality_feasible=minimum_q >= .1 and minimum_det > 0,
        intrinsic_star_overlap=len(set(parents)) != len(parents), added_nodes=len(group), added_tets=len(parents))
    if summary["intrinsic_star_overlap"]:
        return summary, None
    new_ids = {e: len(points)+i for i, e in enumerate(group)}
    children = []
    for e, parent, _ in incident:
        for endpoint in edges[e]:
            child = cells[parent].copy()
            child[child == endpoint] = new_ids[e]
            children.append(child)
    children = np.asarray(children)
    node_ids = np.unique(np.r_[cells[parents].ravel(), children.ravel()])
    inverse = {int(node): i for i, node in enumerate(node_ids)}
    local_points = np.array([points[node] if node < len(points) else mids[node-len(points)] for node in node_ids])
    before = np.array([[inverse[int(node)] for node in tet] for tet in cells[parents]])
    after = np.array([[inverse[int(node)] for node in tet] for tet in children])
    face_indices = np.asarray(list(combinations(range(4), 3)))
    before_faces = np.unique(np.sort(before[:, face_indices].reshape(-1, 3), axis=1), axis=0)
    after_faces = np.unique(np.sort(after[:, face_indices].reshape(-1, 3), axis=1), axis=0)
    summary["whole_star_before"] = dict(faces=measure(case, local_points, before_faces, faces=True), tets=measure(case, local_points, before, faces=False))
    summary["whole_star_after"] = dict(faces=measure(case, local_points, after_faces, faces=True), tets=measure(case, local_points, after, faces=False))
    return summary, (local_points, inverse, new_ids)


def worker(a):
    import numpy as np
    from vtkmodules.vtkCommonCore import vtkSMPTools
    import pmy_edge_sampling_control as edge
    import pmy_global_preclip_preflight as preflight
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
    parent = json.loads((a.parent/"result.json").read_text())
    assert w.digest(a.parent/"background-edge-round.npz") == parent["frozen_decision_hashes"]["background-edge-round.npz"]
    case = fixtures.case(NAME, a.fixtures)
    files = [Path(__file__), a.source/"result.json", a.source/"manifest.json", a.source/"selection.npz",
        a.source/"edge-plan.npz", a.source/"edge-decisions.json", a.source/"background-edge-round.npz",
        a.source/"criterion-observation.npz", a.parent/"result.json", a.parent/"background-edge-round.npz",
        a.witnesses/"result.json", *list((a.fixtures/NAME).iterdir())]
    files += [w.HERE/name for name in ("pmy_edge_sampling_control.py", "pmy_interpolation_control.py",
        "pmy_global_preclip_preflight.py", "pmy_barycenter_control.py", "pmy_neck_control.py",
        "pmy_residual_probe.py", "tpms_family_fixtures.py", "microgen_workflow.py", "mmg_quality_localization.py",
        "periodic_trace.py", "scorecard.py")]
    from microgen.shape import surface_functions
    files.append(Path(surface_functions.__file__))
    for path in [Path(surface_functions.__file__), w.HERE/"pmy_neck_control.py", *list((a.fixtures/NAME).iterdir())]:
        assert w.digest(path) == prior["inputs"][str(path)]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix == ".py":
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, vtk_backend="Sequential", top_faces=TOP_COUNT, top_tets=TOP_COUNT,
        diagnostic_selection="Highest pre-pass normalized residual among frozen marked entities; descending value then entity index",
        severity_priority="Maximum field residual pooled over all marked faces/tets requesting a periodic edge orbit",
        schedule_control="Same candidate orbits, exact stored child-quality feasibility, disjoint-star and 250000 caps; change ordering only",
        alternatives="All existing edges of each top entity; whole periodic star q and same physical-field/cap residuals on temporary children",
        local_children_only=True, global_candidate_generated=False, native_generation_performed=False,
        heldout_use="Final schedule-coverage observer only; no query coordinates enter ranking, alternative choices or severity scores"))
    with np.load(a.parent/"background-edge-round.npz") as source:
        points, cells = source["points"], source["tetrahedra"]
    marks, plan = np.load(a.source/"selection.npz"), np.load(a.source/"edge-plan.npz")
    records = json.loads((a.source/"edge-decisions.json").read_text())
    edges = plan["edges"]
    result = dict(inputs=hashes, global_candidate_generated=False, native_generation_performed=False, admitted=False)
    with Job(parts={"frozen-schedule-replay": 2, "field-severity-dry-run": 1,
                    "local-alternative-audit": 6, "heldout-observers": 1}) as job:
        with job.part("frozen-schedule-replay") as part:
            faces, counts, _, boundary = preflight.inventory(cells)
            np.testing.assert_array_equal(faces, marks["faces"])
            _, classes, _ = periodic_pairs(points, boundary, declared_bounds=([-.5]*3, [.5]*3))
            edge_lookup = {tuple(e): i for i, e in enumerate(edges)}
            incidence = [[] for _ in edges]
            for parent_id, tet in enumerate(cells):
                for offset, (i, j) in enumerate(edge.EDGES):
                    edge_id = edge_lookup[tuple(sorted((tet[i], tet[j])))]
                    incidence[edge_id].append((parent_id, offset))
            edge_record = {e: i for i, row in enumerate(records) for e in row["edges"]}
            groups = {}
            for e, nodes in enumerate(edges):
                groups.setdefault(tuple(sorted(map(int, classes[nodes]))), []).append(e)
            actual, occupied, replay = schedule(records, range(len(records)), len(points), len(cells))
            assert replay["counts"] == prior["decisions"]
            assert replay["nodes"] == prior["background"]["after"]["nodes"] and replay["tets"] == prior["background"]["after"]["tets"]
            assert all(row["reason"] == records[row["orbit"]]["reason"] for row in replay["outcomes"])
            result["original_schedule"] = {key: value for key, value in replay.items() if key != "outcomes"}
            w.save(a.run/"original-schedule-replay.json", replay)
            part.update(1, "Reproduced all original accept/reject decisions and exact counts", orbits=len(records), accepted=replay["counts"]["accepted"])
        with job.part("field-severity-dry-run") as part:
            scores = np.zeros(len(records))
            targets = []
            for kind, entities, selected, estimates, preferred in (
                ("face", faces, marks["selected_faces"], marks["face_estimate"], plan["preferred_face_edge"]),
                ("tet", cells, marks["selected_tets"], marks["tet_estimate"], plan["preferred_tet_edge"])):
                ranked = sorted(map(int, np.flatnonzero(selected)), key=lambda i: (-float(estimates[i]), i))
                for index in ranked:
                    orbit = edge_record[int(preferred[index])]
                    scores[orbit] = max(scores[orbit], float(estimates[index]))
                targets.extend(dict(kind=kind, index=index, severity_rank=rank+1, residual=float(estimates[index]))
                               for rank, index in enumerate(ranked[:TOP_COUNT]))
            lengths = ((points[edges[:, 1]]-points[edges[:, 0]])**2).sum(1)
            order = sorted(range(len(records)), key=lambda i: (-scores[i], -float(lengths[records[i]["edges"]].max()),
                                                                tuple(records[i]["key"]), records[i]["edges"][0]))
            severity, severity_occupied, alternative = schedule(records, order, len(points), len(cells))
            alternative["orbit_severity"] = scores.tolist()
            w.save(a.run/"severity-schedule-dry-run.json", alternative)

            def coverage(accepted, parent_occupied):
                accepted_edges = np.zeros(len(edges), dtype=bool)
                for orbit in np.flatnonzero(accepted):
                    accepted_edges[records[orbit]["edges"]] = True
                split_faces = np.array([any(accepted_edges[edge_lookup[tuple(sorted(e))]] for e in combinations(face, 2)) for face in faces])
                rows = {}
                for kind, marked, estimates, split in (("faces", marks["selected_faces"], marks["face_estimate"], split_faces),
                                                      ("tets", marks["selected_tets"], marks["tet_estimate"], parent_occupied >= 0)):
                    remaining = marked & ~split
                    rows[kind] = dict(marked=int(marked.sum()), split=int((marked & split).sum()), unsplit=int(remaining.sum()),
                        maximum_unsplit_residual=float(estimates[remaining].max(initial=0)),
                        sum_split_mark_residual=float(estimates[marked & split].sum()))
                return rows

            result["severity_dry_run"] = dict(counts=alternative["counts"], nodes=alternative["nodes"], tets=alternative["tets"],
                original_coverage=coverage(actual, occupied), severity_coverage=coverage(severity, severity_occupied),
                scope="Coverage prediction only; local/global field error after this schedule is not measured")
            w.save(a.run/"diagnostic-targets.json", targets)
            part.update(1, "Frozen global field-ranked diagnostics and severity-only dry run", faces=TOP_COUNT, tets=TOP_COUNT, accepted=alternative["counts"]["accepted"])
        with job.part("local-alternative-audit") as part:
            cache, rows = {}, []
            for done, target in enumerate(targets):
                kind, index = target["kind"], target["index"]
                entity = faces[index] if kind == "face" else cells[index]
                preferred = int(plan["preferred_face_edge" if kind == "face" else "preferred_tet_edge"][index])
                preferred_orbit = edge_record[preferred]
                prior_outcome = replay["outcomes"][preferred_orbit]
                row = dict(target, vertices=entity.tolist(), periodic_entity_key=sorted(map(int, classes[entity])),
                    preferred_edge=preferred, preferred_orbit=preferred_orbit,
                    preferred_status=records[preferred_orbit]["reason"], orbit_severity=float(scores[preferred_orbit]),
                    severity_schedule_preferred_accepted=bool(severity[preferred_orbit]),
                    earlier_blockers=[dict(orbit=i, pooled_severity=float(scores[i]), minimum_child_q=records[i]["minimum_child_q"])
                                      for i in prior_outcome["earlier_blocking_orbits"]], alternatives=[])
                for vertices in combinations(entity, 2):
                    e = edge_lookup[tuple(sorted(vertices))]
                    key = tuple(sorted(map(int, classes[edges[e]])))
                    group = groups[key]
                    if key not in cache:
                        cache[key] = local_alternative(case, points, cells, edges, group, incidence, classes)
                    summary, local = cache[key]
                    parent_ids = sorted({parent for member in group for parent, _ in incidence[member]})
                    existing = edge_record.get(e)
                    observation = dict(summary, edge=e, vertices=edges[e].tolist(), is_preferred=e == preferred,
                        existing_candidate_orbit=existing,
                        current_status=records[existing]["reason"] if existing is not None else "not-requested-by-any-mark",
                        conflicts_with_final_accepted_orbits=sorted(set(map(int, occupied[parent_ids][occupied[parent_ids] >= 0]))))
                    if local is not None:
                        local_points, inverse, new_ids = local
                        original_entity = np.array([[inverse[int(node)] for node in entity]])
                        children = []
                        for endpoint in edges[e]:
                            children.append([inverse[new_ids[e]] if node == endpoint else inverse[int(node)] for node in entity])
                        observation["entity_before"] = measure(case, local_points, original_entity, faces=kind == "face")
                        observation["entity_after"] = measure(case, local_points, np.asarray(children), faces=kind == "face")
                        observation["entity_residual_improves"] = bool(observation["entity_after"]["unresolved"] == 0 and
                            observation["entity_after"]["maximum"] < observation["entity_before"]["maximum"]-1e-12)
                    row["alternatives"].append(observation)
                rows.append(row)
                part.update((done+1)/len(targets), "Measured all existing edge alternatives for top field entities", completed=done+1, entities=len(targets), complete_orbits=len(cache))
            result["top_field_entities"] = rows
            with np.load(a.source/"criterion-observation.npz") as after, np.load(a.source/"background-edge-round.npz") as output:
                persistence = {}
                for kind, before_entities, after_entities, estimates in (
                    ("faces", faces, after["faces"], after["face_estimate"]),
                    ("tets", cells, output["tetrahedra"], after["tet_estimate"])):
                    lookup = {tuple(sorted(map(int, entity))): index for index, entity in enumerate(before_entities)}
                    ranked = sorted(range(len(estimates)), key=lambda i: (-float(estimates[i]), i))[:TOP_COUNT]
                    persistence[kind] = [dict(after_index=i, residual=float(estimates[i]),
                        unchanged_original_entity=lookup.get(tuple(sorted(map(int, after_entities[i]))))) for i in ranked]
                result["highest_post_pass_field_entities"] = persistence
            result["summary"] = dict(
                top_entities=len(rows), top_periodic_entity_classes=len({(row["kind"], tuple(row["periodic_entity_key"])) for row in rows}),
                preferred_overlap=sum(row["preferred_status"] == "overlap" for row in rows),
                preferred_quality_reject=sum(row["preferred_status"] == "quality" for row in rows),
                preferred_accepted=sum(row["preferred_status"] == "accepted" for row in rows),
                preferred_severity_accepted=sum(row["severity_schedule_preferred_accepted"] for row in rows),
                entities_with_quality_feasible_improving_edge=sum(any(alt["quality_feasible"] and not alt["intrinsic_star_overlap"] and alt.get("entity_residual_improves", False) for alt in row["alternatives"]) for row in rows),
                entities_with_nonpreferred_quality_feasible_improvement=sum(any(not alt["is_preferred"] and alt["quality_feasible"] and not alt["intrinsic_star_overlap"] and alt.get("entity_residual_improves", False) for alt in row["alternatives"]) for row in rows))
        with job.part("heldout-observers") as part:
            queries = json.loads((a.witnesses/"result.json").read_text())["witnesses"]
            background = residual.AffineBackground(points, cells, np.asarray(case.function(*points.T)))
            observers = []
            for index, witness in enumerate(queries):
                containing = background.interpolation(np.asarray(witness["probe"]), witness["wall_level"])["containing_tets"]
                observers.append(dict(query_index=witness["query_index"], containing_tets=containing,
                    containing_marked=bool(marks["selected_tets"][containing].any()),
                    original_schedule_split=bool((occupied[containing] >= 0).any()),
                    severity_schedule_would_split=bool((severity_occupied[containing] >= 0).any()),
                    field_residuals=marks["tet_estimate"][containing].tolist()))
                part.update((index+1)/len(queries), "Observed held-out coverage after all field-only diagnostics", completed=index+1, queries=len(queries))
            result["heldout_coverage_observers"] = observers
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    result.update(complete=True, input_hashes_verified=True,
        interpretation="Read-only cause audit. Temporary local child estimates and schedule coverage are not new meshes, global error reductions or native admission. Selection changes require a separate parent decision.")
    w.save(a.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-second-edge-preflight-v1")
    parser.add_argument("--parent", type=Path, default=w.HERE/"runs/macbook-pmy-edge-sampling-v1")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--witnesses", type=Path, default=w.HERE/"runs/macbook-pmy-residual-probe-v1")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "source", "parent", "fixtures", "witnesses"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-pmy_edge_causal_audit.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **kw: None), timeout=TIMEOUT)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    findings = json.loads((args.run/"result.json").read_text())
    result("PMY frozen field-entity causes", **findings["summary"])
    dry = findings["severity_dry_run"]
    result("PMY severity schedule dry run", nodes=dry["nodes"], tets=dry["tets"], **dry["counts"])


if __name__ == "__main__":
    main()
