"""One quality-guarded periodic edge-star pass using the frozen job764 marks.

This changes only the subdivision operation. Original coordinates, physical
field, specimen, marks, quality floor and count caps remain fixed. No clipping,
curving or FEA is performed. Witness coordinates are opened only after the
complete edge-selection and admission decisions have been saved.
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
QUALITY = .1
EDGES = tuple(combinations(range(4), 2))


def child_metrics(coordinates, endpoints, midpoint):
    import numpy as np
    from mmg_quality_localization import quality

    children = np.repeat(coordinates[None], 2, axis=0)
    children[0, endpoints[0]] = midpoint
    children[1, endpoints[1]] = midpoint
    determinant = np.linalg.det(children[:, 1:]-children[:, :1])
    q = quality(children.reshape(-1, 3), np.arange(8).reshape(2, 4))
    assert np.isfinite(determinant).all() and np.isfinite(q).all()
    return children, determinant, q


def orbit_midpoints(points, edges, group, classes):
    """Translate one midpoint exactly to each periodic copy of a short edge."""
    import numpy as np

    ordered = [edge[np.argsort(classes[edge])] for edge in edges[group]]
    reference = points[ordered[0]]
    midpoint = reference.mean(0)
    result = []
    for edge in ordered:
        delta = points[edge]-reference
        shift = np.rint(delta[0])
        np.testing.assert_allclose(delta, np.tile(shift, (2, 1)), atol=1e-14, rtol=0)
        assert (abs(shift) <= 1).all()
        value = midpoint+shift
        for axis in range(3):
            for side in (-.5, .5):
                if (points[edge, axis] == side).all():
                    value[axis] = side
        np.testing.assert_allclose(value, points[edge].mean(0), atol=1e-14, rtol=0)
        result.append(value)
    return np.asarray(result)


def controls():
    import numpy as np

    regular = np.array([[0., 0., 0.], [1., 0., 0.], [.5, np.sqrt(3)/2, 0.],
                        [.5, np.sqrt(3)/6, np.sqrt(2/3)]])
    children, determinant, q = child_metrics(regular, (0, 1), regular[:2].mean(0))
    assert (determinant > 0).all() and (q >= QUALITY).all()
    original = np.linalg.det(regular[1:]-regular[0])
    np.testing.assert_allclose(determinant.sum(), original, atol=1e-14, rtol=0)
    inverted = regular[[1, 0, 2, 3]]
    _, bad_det, _ = child_metrics(inverted, (0, 1), inverted[:2].mean(0))
    assert (bad_det < 0).all()
    flat = regular.copy()
    flat[3, 2] *= .001
    _, positive, low = child_metrics(flat, (0, 1), flat[:2].mean(0))
    assert (positive > 0).all() and (low < QUALITY).all()
    points = np.array([[-.5, .125, .25], [-.5, .375, .25], [.5, .125, .25], [.5, .375, .25]])
    mids = orbit_midpoints(points, np.array([[0, 1], [2, 3]]), [0, 1], np.array([0, 1, 0, 1]))
    np.testing.assert_array_equal(mids[1]-mids[0], [1., 0., 0.])
    assert len(np.unique(children.reshape(-1, 3), axis=0)) == 5
    faces = np.array(list(combinations(range(4), 3)))
    split = edge_pass(regular, np.arange(4)[None], faces, np.ones(4, dtype=bool),
                     np.ones(1, dtype=bool), np.arange(4), lambda *args: None)
    refined, cells, owners, _, _, counts, _, _ = split
    assert len(refined) == 5 and len(cells) == 2 and (owners == 0).all()
    assert counts["accepted"] == 1 and counts["overlap"] >= 1
    np.testing.assert_array_equal(refined[:4], regular)
    return dict(regular_children_accepted=True, positive_parent_volume_preserved=True,
                inverted_children_rejected=True, positive_low_quality_children_rejected=True,
                periodic_midpoints_exact=True, overlapping_stars_rejected=True,
                accepted_star_preserves_original_vertices=True)


def edge_pass(points, cells, faces, marked_faces, marked_tets, classes, progress):
    import numpy as np

    raw = np.sort(cells[:, np.asarray(EDGES)].reshape(-1, 2), axis=1)
    edges, inverse = np.unique(raw, axis=0, return_inverse=True)
    lookup = {tuple(edge): i for i, edge in enumerate(edges)}
    incidence = [[] for _ in edges]
    for row, edge in enumerate(inverse):
        incidence[edge].append(divmod(row, len(EDGES)))
    keys = [tuple(sorted(map(int, classes[edge]))) for edge in edges]
    assert all(a != b for a, b in keys), "Unexpected edge spanning a whole periodic cell"
    groups = {}
    for index, key in enumerate(keys):
        groups.setdefault(key, []).append(index)
    squared = ((points[edges[:, 1]]-points[edges[:, 0]])**2).sum(1)

    def preferred(entities, marked):
        selected = np.full(len(entities), -1, dtype=np.int64)
        for index in np.flatnonzero(marked):
            local = [lookup[tuple(sorted(edge))] for edge in combinations(entities[index], 2)]
            # Exact longest length, then periodic node-class key and edge index.
            selected[index] = min(local, key=lambda edge: (-squared[edge], keys[edge], edge))
        return selected

    face_edge, tet_edge = preferred(faces, marked_faces), preferred(cells, marked_tets)
    requested = set(keys[i] for i in np.concatenate((face_edge[marked_faces], tet_edge[marked_tets])))
    ordered = sorted((groups[key] for key in requested),
                     key=lambda group: (-float(squared[group].max()), keys[group[0]], group[0]))
    selected = np.full(len(cells), -1, dtype=np.int64)
    edge_midpoint = np.full(len(edges), -1, dtype=np.int64)
    nodes = points.tolist()
    live = len(cells)
    records, failed = [], []
    accepted_keys = set()
    counts = dict(accepted=0, overlap=0, quality=0, budget=0)
    for index, group in enumerate(ordered):
        incident = [(int(edge), int(parent), int(offset)) for edge in group for parent, offset in incidence[edge]]
        parents = np.array([parent for _, parent, _ in incident], dtype=np.int64)
        mids = orbit_midpoints(points, edges, group, classes)
        midpoint_by_edge = dict(zip(group, mids, strict=True))
        # Measure every candidate even if overlap wins the rejection order.
        minimum_q, minimum_det = 1., float("inf")
        worst = None
        for edge, parent, offset in incident:
            children, determinant, q = child_metrics(points[cells[parent]], EDGES[offset], midpoint_by_edge[edge])
            if float(q.min()) < minimum_q:
                minimum_q = float(q.min())
                worst = dict(edge=edge, parent=parent, local_edge=list(EDGES[offset]),
                             children=children.tolist(), determinant=determinant.tolist(), quality=q.tolist())
            minimum_det = min(minimum_det, float(determinant.min()))
        overlap = len(np.unique(parents)) != len(parents) or bool((selected[parents] >= 0).any())
        quality_ok = minimum_det > 0 and minimum_q >= QUALITY
        budget_ok = len(nodes)+len(group) <= w.CAP and live+len(incident) <= w.CAP
        reason = "overlap" if overlap else "quality" if not quality_ok else "budget" if not budget_ok else "accepted"
        counts[reason] += 1
        record = dict(orbit=index, key=list(keys[group[0]]), edges=list(map(int, group)),
                      incident_tets=parents.tolist(), reason=reason, minimum_child_q=minimum_q,
                      minimum_child_determinant=minimum_det, quality_ok=quality_ok,
                      overlap=overlap, budget_ok=budget_ok)
        records.append(record)
        if not quality_ok:
            failed.append(dict(orbit=index, rejection_reason=reason, worst_child=worst))
        if reason == "accepted":
            accepted_keys.add(keys[group[0]])
            for edge, midpoint in zip(group, mids, strict=True):
                edge_midpoint[edge] = len(nodes)
                nodes.append(midpoint.tolist())
            for edge, parent, offset in incident:
                assert selected[parent] < 0
                selected[parent] = edge
            live += len(incident)
        if index % 128 == 0 or index+1 == len(ordered):
            progress((index+1)/max(len(ordered), 1), index+1, len(ordered), counts, len(nodes), live)
    refined = np.asarray(nodes)
    assert len(np.unique(refined, axis=0)) == len(refined), "Midpoint duplicates an existing vertex"
    np.testing.assert_array_equal(refined[:len(points)], points)
    tets, owners = [], []
    for parent, tet in enumerate(cells):
        edge = selected[parent]
        if edge < 0:
            tets.append(tet.tolist())
            owners.append(parent)
            continue
        endpoints = edges[edge]
        for endpoint in endpoints:
            child = tet.copy()
            child[child == endpoint] = edge_midpoint[edge]
            tets.append(child.tolist())
            owners.append(parent)
    tets, owners = np.asarray(tets), np.asarray(owners)
    assert len(tets) == live
    face_preferred = np.array([edge >= 0 and keys[edge] in accepted_keys for edge in face_edge])
    tet_preferred = np.array([edge >= 0 and keys[edge] in accepted_keys for edge in tet_edge])
    face_split = np.array([any(edge_midpoint[lookup[tuple(sorted(edge))]] >= 0
                             for edge in combinations(face, 2)) for face in faces])
    coverage = dict(selected_faces=int(marked_faces.sum()), selected_tets=int(marked_tets.sum()),
                    preferred_faces_accepted=int((marked_faces & face_preferred).sum()),
                    preferred_tets_accepted=int((marked_tets & tet_preferred).sum()),
                    marked_faces_split=int((marked_faces & face_split).sum()),
                    marked_tets_split=int((marked_tets & (selected >= 0)).sum()),
                    remaining_unsplit_faces=int((marked_faces & ~face_split).sum()),
                    remaining_unsplit_tets=int((marked_tets & (selected < 0)).sum()))
    arrays = dict(edges=edges, edge_midpoint=edge_midpoint, original_tet_edge=selected,
                  preferred_face_edge=face_edge, preferred_tet_edge=tet_edge,
                  marked_faces=marked_faces, marked_tets=marked_tets,
                  remaining_faces=marked_faces & ~face_split, remaining_tets=marked_tets & (selected < 0))
    return refined, tets, owners, records, failed, counts, coverage, arrays


def criterion(case, points, cells, gradient):
    import numpy as np
    import pmy_global_preclip_preflight as preflight
    from pmy_interpolation_control import wall_estimates

    faces, counts, _, _ = preflight.inventory(cells)
    prior = preflight.selection(case, points, cells, faces, counts)
    face = wall_estimates(case, points, faces, gradient, faces=True)
    tet = wall_estimates(case, points, cells, gradient, faces=False)
    selected_faces, _, _, _, _ = preflight.paired_selection(
        points, faces, counts, prior["selected_faces"] | face["selected"], prior["face_centers"])

    def summary(old, current, selected):
        return dict(class_selected=int(old.sum()), geometry_selected=int(current["selected"].sum()),
                    combined=int(selected.sum()), wall_active=int(current["active"].sum()),
                    unresolved_gradient=int(current["unresolved"].sum()),
                    maximum_normalized_residual=float(current["maximum"].max()),
                    active_residual_quantiles=np.quantile(current["maximum"][current["active"]], [.5, .9, .99, 1]).tolist())
    return dict(faces=summary(prior["selected_faces"], face, selected_faces),
                tets=summary(prior["selected_tets"], tet, prior["selected_tets"] | tet["selected"])), dict(
        faces=faces, selected_faces=selected_faces, selected_tets=prior["selected_tets"] | tet["selected"],
        face_estimate=face["maximum"], tet_estimate=tet["maximum"])


def worker(a):
    import numpy as np
    from vtkmodules.vtkCommonCore import vtkSMPTools
    from microgen.shape import surface_functions
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

    source = a.source/"input.json"
    prepared = json.loads((a.source.parent/"prepare-result.json").read_text())
    prior = json.loads((a.selection/"result.json").read_text())
    prior_manifest = json.loads((a.selection/"manifest.json").read_text())
    native_manifest = json.loads((a.source.parent/"manifest.json").read_text())
    assert w.digest(source) == prepared["generated_input_hashes"]["one-round"]
    assert w.digest(source) == prior_manifest["inputs"][str(source)]
    assert w.digest(a.selection/"selection.npz") == prior["selection_sha256"]
    for path in (Path(surface_functions.__file__), Path(neck.__file__)):
        assert w.digest(path) == native_manifest["inputs"][str(path.resolve())]
    case = fixtures.case(NAME, a.fixtures)
    frozen = json.loads((a.fixtures/"result.json").read_text())
    known = {Path(path).resolve(): digest for path, digest in frozen["inputs"].items()}
    specimen = list((a.fixtures/NAME).iterdir())
    assert all(w.digest(path) == known[path.resolve()] for path in specimen)
    files = [Path(__file__), Path(preflight.__file__), Path(interpolation.__file__), Path(neck.__file__),
             Path(residual.__file__), Path(fixtures.__file__), Path(surface_functions.__file__), source,
             a.source.parent/"prepare-result.json", a.source.parent/"manifest.json",
             a.selection/"result.json", a.selection/"manifest.json", a.selection/"selection.npz",
             a.witnesses/"result.json", a.fixtures/"result.json", *specimen]
    files += [w.HERE/name for name in ("microgen_workflow.py", "scorecard.py", "periodic_trace.py",
                                      "mmg_quality_localization.py", "pmy_barycenter_control.py", "pmy_failure_probe.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix == ".py":
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        selection="Exactly saved job764 selected_faces/selected_tets on v4 unsnapped one-round input",
        policy="One longest edge per marked face/tet; exact length then periodic class key; longest orbit first; disjoint complete stars",
        guard="Every child positive with q>=.1; unchanged 250000 background node/tet caps",
        rejection_order=["overlap", "quality", "budget"], quality_checked_even_when_overlap=True,
        estimator_scope="Same .002 normalized residual heuristic, tangent gradient on physical caps; not a distance certificate",
        witness_use="All choices and rejection evidence saved before reading eight held-out witness coordinates",
        candidate_generated=False, native_clipping_performed=False, solve_performed=False))
    data = json.loads(source.read_text())
    points, cells = np.array(data["points"]), np.array(data["tetrahedra"], dtype=np.int64)
    marks = np.load(a.selection/"selection.npz")
    result = dict(inputs=hashes, admitted=False, candidate_generated=False, native_clipping_performed=False,
                  admission_gates_unchanged=True, source=str(a.source), frozen_selection=str(a.selection))
    with Job(parts={"controls-and-inputs": 1, "edge-stars": 5, "background-audit": 2,
                    "criterion-observation": 2, "witness-observation": 2}) as job:
        with job.part("controls-and-inputs") as part:
            part.update(0, "Checking fixed input hashes and bisection controls", nodes=len(points), tets=len(cells))
            result["controls"] = controls()
            faces, counts, parents, boundary = preflight.inventory(cells)
            for key, actual in (("faces", faces), ("face_counts", counts), ("face_parents", parents)):
                np.testing.assert_array_equal(marks[key], actual)
            _, classes, _ = periodic_pairs(points, boundary, declared_bounds=([-.5]*3, [.5]*3))
            part.update(1, "Frozen selection indices and periodic classes verified", faces=int(marks["selected_faces"].sum()), tets=int(marks["selected_tets"].sum()))
        with job.part("edge-stars") as part:
            output = edge_pass(points, cells, faces, marks["selected_faces"], marks["selected_tets"], classes,
                lambda fraction, completed, total, counts, nodes, tets: part.update(fraction,
                    "Testing complete periodic edge stars", completed=completed, orbits=total,
                    accepted=counts["accepted"], quality_rejected=counts["quality"], overlap=counts["overlap"], nodes=nodes, tets=tets))
            refined, tets, owners, records, failed, decisions, coverage, arrays = output
            w.save(a.run/"edge-decisions.json", records)
            w.save(a.run/"quality-failure-proof.json", failed)
            np.savez_compressed(a.run/"edge-plan.npz", **arrays)
            np.savez_compressed(a.run/"background-edge-round.npz", points=refined, tetrahedra=tets, original_tet=owners)
            result.update(decisions=decisions, coverage=coverage,
                          overlapping_quality_failures=sum(not row["quality_ok"] and row["overlap"] for row in records),
                          all_quality_infeasible_orbits=sum(not row["quality_ok"] for row in records))
            result["frozen_decision_hashes"] = {name: w.digest(a.run/name) for name in (
                "edge-decisions.json", "quality-failure-proof.json", "edge-plan.npz", "background-edge-round.npz")}
            part.update(1, "Saved frozen decisions and all quality rejection proofs", accepted=decisions["accepted"], nodes=len(refined), tets=len(tets))
        with job.part("background-audit") as part:
            before, old_det, _ = preflight.metrics(points, cells)
            part.update(.4, "Original background gates verified", minimum_q=before["minimum_quality"])
            after, new_det, _ = preflight.metrics(refined, tets)
            error = float(np.max(abs(np.bincount(owners, weights=new_det)-old_det)/old_det))
            assert error < 1e-10
            np.testing.assert_array_equal(refined[:len(points)], points)
            result["background"] = dict(before=before, after=after, maximum_parent_volume_error=error,
                original_points_bitwise_preserved=True, passed=bool(after["quality_floor_met"] and after["within_caps"]))
            part.update(1, "Audited quality, exact periodic caps and fixed-domain subdivision", minimum_q=after["minimum_quality"], nodes=len(refined), tets=len(tets))
        with job.part("criterion-observation") as part:
            original, _ = criterion(case, points, cells, neck.raw_gradient)
            for kind in ("faces", "tets"):
                assert original[kind]["combined"] == prior["selection"][kind]["combined"]
            part.update(.5, "Replayed unchanged criterion on original background", faces=original["faces"]["combined"], tets=original["tets"]["combined"])
            updated, observations = criterion(case, refined, tets, neck.raw_gradient)
            np.savez_compressed(a.run/"criterion-observation.npz", **observations)
            result["criterion"] = dict(before=original, after=updated,
                unresolved_original_parent_tets=int(len(np.unique(owners[observations["selected_tets"]]))),
                scope="Counts and fractions refer to different subdivisions; residuals are first-order selection heuristics")
            part.update(1, "Observed remaining field error without changing selections", faces=updated["faces"]["combined"], tets=updated["tets"]["combined"])
        with job.part("witness-observation") as part:
            witnesses = json.loads((a.witnesses/"result.json").read_text())["witnesses"]
            assert len(witnesses) == 8
            original = residual.AffineBackground(points, cells, np.asarray(case.function(*points.T)))
            updated = residual.AffineBackground(refined, tets, np.asarray(case.function(*refined.T)))
            rows = []
            for index, witness in enumerate(witnesses):
                query = np.asarray(witness["probe"])
                before = residual.affine_boundary_distance(original, query, case.levels)
                saved = witness["neighborhood"]["backgrounds"]["refined"]["boundary"]
                np.testing.assert_allclose(before["nearest_distance"], saved["nearest_distance"], atol=1e-12, rtol=0)
                after = residual.affine_boundary_distance(updated, query, case.levels)
                containing = original.interpolation(query, witness["wall_level"])["containing_tets"]
                rows.append(dict(query_index=witness["query_index"], probe=witness["probe"],
                    before=before, after=after, original_containing_tets=containing,
                    containing_tet_split=bool((arrays["original_tet_edge"][containing] >= 0).any()),
                    updated_center=updated.interpolation(query, witness["wall_level"])))
                part.update((index+1)/len(witnesses), "Measured held-out affine boundary distances", completed=index+1, witnesses=len(witnesses), distance=after["nearest_distance"])
            result["witnesses"] = rows
            result["witness_summary"] = dict(count=len(rows),
                within_distance_gate=sum(row["after"]["nearest_is_global_for_affine_representation"] and row["after"]["nearest_distance"] <= .002 for row in rows),
                improved=sum(row["after"]["nearest_distance"] is not None and row["after"]["nearest_distance"] < row["before"]["nearest_distance"]-1e-12 for row in rows),
                worst_after=max(row["after"]["distance_lower_bound"] for row in rows))
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    assert all(w.digest(a.run/name) == digest for name, digest in result["frozen_decision_hashes"].items())
    result.update(complete=True, input_hashes_verified=True,
        interpretation="One background-only edge-star pass. A valid background is not an admitted PMY mesh. No native clipping, Tet10 count, timing win, physical topology certificate or FEA claim.")
    w.save(a.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-native-background-v4/one-round")
    parser.add_argument("--selection", type=Path, default=w.HERE/"runs/macbook-pmy-interpolation-control-v1")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--witnesses", type=Path, default=w.HERE/"runs/macbook-pmy-residual-probe-v1")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "source", "selection", "fixtures", "witnesses"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-pmy_edge_sampling_control.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **kw: None), timeout=TIMEOUT)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    findings = json.loads((args.run/"result.json").read_text())
    result("PMY edge decisions", **findings["decisions"])
    after = findings["background"]["after"]
    result("PMY background preflight", nodes=after["nodes"], tets=after["tets"], minimum_q=after["minimum_quality"], passed=findings["background"]["passed"])
    result("PMY held-out affine witnesses", **findings["witness_summary"])


if __name__ == "__main__":
    main()
