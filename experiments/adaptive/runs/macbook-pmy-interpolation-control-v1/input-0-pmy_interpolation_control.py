"""One field-only interpolation criterion and background resource preflight.

Compare three-way class sampling with normalized residuals at affine wall-chord
midpoints and wall-polygon vertex averages. These are selection heuristics, not
physical distance certificates. Freeze selection before reading saved witnesses.
Construct a background only if exact combinatorial counts fit the unchanged caps.
No native clipping, production policy change, mesh admission or FEA is performed.
"""

import argparse
from itertools import combinations, permutations
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
TOLERANCE = .002
TIMEOUT = 240


def wall_estimates(case, points, entities, gradient, *, faces):
    import numpy as np

    coordinates = points[entities]
    values = np.asarray(case.function(*points.T))[entities]
    maxima = np.zeros(len(entities))
    unresolved = np.zeros(len(entities), dtype=bool)
    active = np.zeros(len(entities), dtype=bool)
    records = []
    for level in case.levels:
        eligible = (values.min(1) <= level) & (values.max(1) >= level)
        indices = np.flatnonzero(eligible)
        xyz, field = coordinates[indices], values[indices]
        total = np.zeros((len(indices), 3))
        count = np.zeros(len(indices), dtype=np.int64)
        for i in range(entities.shape[1]):
            match = field[:, i] == level
            total[match] += xyz[match, i]
            count += match
        for i, j in combinations(range(entities.shape[1]), 2):
            crossing = (field[:, i]-level)*(field[:, j]-level) < 0
            fraction = (level-field[crossing, i])/(field[crossing, j]-field[crossing, i])
            total[crossing] += xyz[crossing, i]+fraction[:, None]*(xyz[crossing, j]-xyz[crossing, i])
            count += crossing
        valid = count >= (2 if faces else 3)
        indices, xyz, total, count = indices[valid], xyz[valid], total[valid], count[valid]
        probes = total/count[:, None]
        exact = np.asarray(case.function(*probes.T))
        gradients = np.asarray(gradient(probes))
        full_norm = np.linalg.norm(gradients, axis=1)
        locks = np.zeros_like(gradients, dtype=bool)
        if faces:
            for axis in range(3):
                locks[:, axis] = (xyz[:, :, axis] == -.5).all(1) | (xyz[:, :, axis] == .5).all(1)
            gradients[locks] = 0
        norm = np.linalg.norm(gradients, axis=1)
        residual = abs(exact-level)
        singular = (norm < 1e-12) & (residual > 1e-12)
        error = np.divide(residual, norm, out=np.zeros_like(residual), where=norm >= 1e-12)
        maxima[indices] = np.maximum(maxima[indices], error)
        unresolved[indices] |= singular
        active[indices] = True
        records.append(dict(entity=indices, level=np.full(len(indices), level), point=probes,
                            residual=residual, physical_gradient=full_norm,
                            constrained_gradient=norm, cap_locked=locks.any(1), estimate=error,
                            unresolved=singular))
    arrays = {key: np.concatenate([record[key] for record in records]) for key in records[0]}
    return dict(maximum=maxima, unresolved=unresolved, active=active,
                selected=(maxima > TOLERANCE) | unresolved, probes=arrays)


def local_center_growth(face_sequence):
    """Exact reference-tet count for one frozen face sequence and center insertion."""
    import numpy as np
    from pmy_neck_control import barycentric

    points = [np.array(p, dtype=float) for p in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))]
    descendants = [(0, 1, 2, 3)]
    for opposite in face_sequence:
        carrier = set(range(4))-{opposite}
        new = len(points)
        points.append(np.mean([points[i] for i in sorted(carrier)], axis=0))
        match = [i for i, tet in enumerate(descendants) if carrier.issubset(tet)]
        assert len(match) == 1
        index = match[0]
        tet = descendants[index]
        replacement = [tuple(new if i == node else i for i in tet) for node in sorted(carrier)]
        descendants = descendants[:index]+descendants[index+1:]+replacement
    points = np.array(points)
    b = barycentric(points, np.array(descendants), np.full(3, .25))
    hits = np.flatnonzero(b.min(1) >= -1e-9)
    assert len(hits)
    carrier = {descendants[hits[0]][j] for j in range(4) if b[hits[0], j] > 1e-9}
    assert 2 <= len(carrier) <= 4
    incident = sum(carrier.issubset(tet) for tet in descendants)
    return incident*(len(carrier)-1)


def exact_counts(points, cells, faces, parents, plan):
    import numpy as np

    sequence = [[] for _ in cells]
    for group in plan["groups"]:
        for face_index in group:
            carrier = set(map(int, faces[face_index]))
            for parent in parents[face_index]:
                if parent < 0:
                    continue
                opposite = [j for j, node in enumerate(cells[parent]) if node not in carrier]
                assert len(opposite) == 1
                sequence[parent].append(opposite[0])
    patterns = {}
    center_growth = 0
    for parent in np.flatnonzero(plan["selected_tets"]):
        key = tuple(sequence[parent])
        if key not in patterns:
            patterns[key] = local_center_growth(key)
        center_growth += patterns[key]
    faces_selected = int(plan["selected_faces"].sum())
    tetrahedra = len(cells)+2*sum(map(len, sequence))+center_growth
    nodes = len(points)+faces_selected+int(plan["selected_tets"].sum())
    return dict(nodes=nodes, tetrahedra=tetrahedra, original_nodes=len(points), original_tetrahedra=len(cells),
                selected_faces=faces_selected, selected_tets=int(plan["selected_tets"].sum()),
                face_tet_growth=2*sum(map(len, sequence)), center_tet_growth=center_growth,
                distinct_local_face_orders=len(patterns), within_caps=max(nodes, tetrahedra) <= w.CAP,
                scope="Exact stellar-subdivision combinatorics; not final clipped/curved mesh counts")


def controls():
    import numpy as np
    import pmy_barycenter_control as old
    import pmy_global_preclip_preflight as preflight
    from scorecard import Case

    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    tets = np.array([[0, 1, 2, 3]])
    faces = np.array(list(combinations(range(4), 3)))
    affine = Case("affine", lambda x, y, z: x, (-2., .1))
    gradient = lambda p: np.tile([1., 0., 0.], (len(p), 1))
    assert not wall_estimates(affine, points, faces, gradient, faces=True)["selected"].any()
    assert not wall_estimates(affine, points, tets, gradient, faces=False)["selected"].any()
    curved = Case("curved-same-class", lambda x, y, z: x*x+y*y+z*z-.1, (-2., 0.))
    for entities, is_face in ((faces, True), (tets, False)):
        centers = points[entities].mean(1)
        actual = curved.function(*centers.T)
        interpolated = curved.function(*points.T)[entities].mean(1)
        assert np.array_equal(old.classify(actual, curved.levels), old.classify(interpolated, curved.levels))
        check = wall_estimates(curved, points, entities, lambda p: 2*p, faces=is_face)
        assert check["selected"].any()
        scaled = Case("scaled", lambda x, y, z: 7*curved.function(x, y, z), (-14., 0.))
        alternative = wall_estimates(scaled, points, entities, lambda p: 14*p, faces=is_face)
        np.testing.assert_allclose(check["maximum"], alternative["maximum"], atol=1e-12, rtol=1e-12)
    cap_points = np.array([[0., 0., -.5], [.5, 0., -.5], [0., .5, -.5]])
    cap = Case("cap", lambda x, y, z: x*x+y*y+10*(z+.5)-.1, (-2., 0.))
    cap_gradient = lambda p: np.column_stack((2*p[:, 0], 2*p[:, 1], np.full(len(p), 10.)))
    cap_check = wall_estimates(cap, cap_points, np.array([[0, 1, 2]]), cap_gradient, faces=True)
    assert cap_check["probes"]["cap_locked"].all()
    assert (cap_check["probes"]["physical_gradient"] > 10.).all()
    np.testing.assert_allclose(cap_check["probes"]["constrained_gradient"], np.sqrt(.08), atol=1e-12)
    # Exercise all 65 ordered subsets of the four original faces. Geometry and
    # orientation may vary globally; the insertion counts are affine invariant.
    patterns = 0
    for length in range(5):
        for sequence in permutations(range(4), length):
            face_indices = list(sequence)
            face_array = np.array([sorted(set(range(4))-{i}) for i in range(4)])
            parent_array = np.column_stack((np.zeros(4, dtype=int), np.full(4, -1)))
            plan = dict(groups=[[i] for i in face_indices], selected_faces=np.isin(np.arange(4), face_indices),
                        selected_tets=np.array([True]), face_centers=points[face_array].mean(1), tet_centers=points[tets].mean(1))
            counts = exact_counts(points, tets, face_array, parent_array, plan)
            refined, refined_tets, _, _, _ = preflight.subdivide(points, tets, face_array, parent_array, plan)
            assert len(refined) == counts["nodes"] and len(refined_tets) == counts["tetrahedra"]
            patterns += 1
    assert patterns == 65

    return dict(affine_field_unselected=True, same_class_curved_error_detected=True,
                positive_field_rescaling_invariant=True, physical_cap_gradient_constrained=True,
                exact_subdivision_counts_controlled=True, ordered_face_patterns=patterns)


def worker(a):
    import numpy as np
    from microgen.shape import surface_functions
    import pmy_barycenter_control as old
    import pmy_global_preclip_preflight as preflight
    import pmy_neck_control as neck
    import tpms_family_fixtures as fixtures

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    source = a.source/"input.json"
    data = json.loads(source.read_text())
    prepared = json.loads((a.source.parent/"prepare-result.json").read_text())
    assert w.digest(source) == prepared["generated_input_hashes"]["one-round"]
    source_manifest = json.loads((a.source.parent/"manifest.json").read_text())
    for path in (Path(surface_functions.__file__), Path(neck.__file__)):
        assert w.digest(path) == source_manifest["inputs"][str(path.resolve())]
    case = fixtures.case(NAME, a.fixtures)
    frozen = json.loads((a.fixtures/"result.json").read_text())
    known = {Path(path).resolve(): digest for path, digest in frozen["inputs"].items()}
    specimen = list((a.fixtures/NAME).iterdir())
    assert all(w.digest(path) == known[path.resolve()] for path in specimen)
    files = [Path(__file__), Path(old.__file__), Path(preflight.__file__), Path(neck.__file__), Path(fixtures.__file__),
             Path(surface_functions.__file__), source, a.source.parent/"prepare-result.json", a.source.parent/"manifest.json",
             a.fixtures/"result.json", *specimen, a.witnesses/"result.json"]
    files += [w.HERE/name for name in ("microgen_workflow.py", "scorecard.py", "periodic_trace.py", "mmg_quality_localization.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix == ".py":
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        criterion="Keep class mismatches; add affine wall-chord midpoint / wall-polygon vertex-average normalized residual above .002",
        gradient="Physical gradient except constrained tangent gradient on physical cap faces",
        insertion_locations="Original face/tet barycenters; all original vertices fixed",
        estimator_scope="First-order selection heuristic, not a nearest-distance certificate",
        witness_use="Observer after selection is frozen; hashes are recorded beforehand",
        candidate_generated=False, native_clipping_performed=False, solve_performed=False))
    points, cells = np.array(data["points"]), np.array(data["tetrahedra"], dtype=np.int64)
    result = dict(inputs=hashes, admission_gates_unchanged=True, native_clipping_performed=False,
                  candidate_generated=False, admitted=False, source=str(a.source))
    with Job(parts={"controls": 1, "field-selection": 3, "resource-counts": 1, "background": 4, "witness-observation": 1}) as job:
        with job.part("controls") as part:
            result["controls"] = controls()
            part.update(1, "Affine, field-scaling and subdivision controls passed", controls=len(result["controls"]))
        with job.part("field-selection") as part:
            faces, incidence, parents, _ = preflight.inventory(cells)
            previous = preflight.selection(case, points, cells, faces, incidence)
            part.update(.25, "Evaluated current classification-only rule", faces=len(faces), tets=len(cells))
            face_estimates = wall_estimates(case, points, faces, neck.raw_gradient, faces=True)
            part.update(.6, "Measured affine wall-chord residuals", active_faces=int(face_estimates["active"].sum()))
            tet_estimates = wall_estimates(case, points, cells, neck.raw_gradient, faces=False)
            requested = previous["selected_faces"] | face_estimates["selected"]
            selected_faces, centers, groups, pairs, forced = preflight.paired_selection(
                points, faces, incidence, requested, previous["face_centers"])
            selected_tets = previous["selected_tets"] | tet_estimates["selected"]
            plan = dict(previous, selected_faces=selected_faces, selected_tets=selected_tets,
                        face_centers=centers, groups=groups, pairs=pairs, forced=forced)
            def counts(prior, estimates, selected):
                return dict(class_selected=int(prior.sum()), geometry_selected=int(estimates["selected"].sum()),
                    geometry_only=int((estimates["selected"] & ~prior).sum()), combined=int(selected.sum()),
                    wall_active=int(estimates["active"].sum()), unresolved_gradient=int(estimates["unresolved"].sum()),
                    maximum_normalized_residual=float(estimates["maximum"].max()))
            result["selection"] = dict(faces=counts(previous["selected_faces"], face_estimates, selected_faces),
                tets=counts(previous["selected_tets"], tet_estimates, selected_tets),
                forced_periodic_members=forced, paired_boundary_classes=len(pairs), frozen_before_witnesses=True)
            arrays = dict(faces=faces, face_counts=incidence, face_parents=parents,
                old_faces=previous["selected_faces"], old_tets=previous["selected_tets"],
                selected_faces=selected_faces, selected_tets=selected_tets,
                face_estimate=face_estimates["maximum"], tet_estimate=tet_estimates["maximum"])
            for prefix, record in (("face", face_estimates), ("tet", tet_estimates)):
                arrays.update({f"{prefix}_probe_{key}": value for key, value in record["probes"].items()})
            np.savez_compressed(a.run/"selection.npz", **arrays)
            result["selection_sha256"] = w.digest(a.run/"selection.npz")
            part.update(1, "Frozen geometry-aware selections with paired caps", faces=int(selected_faces.sum()), tets=int(selected_tets.sum()))
        job.result("field-selection", extra_faces=result["selection"]["faces"]["geometry_only"], extra_tets=result["selection"]["tets"]["geometry_only"])
        with job.part("resource-counts") as part:
            predicted = exact_counts(points, cells, faces, parents, plan)
            result["predicted_background"] = predicted
            part.update(1, "Computed exact proposed stellar-subdivision counts", nodes=predicted["nodes"], tets=predicted["tetrahedra"], within_caps=predicted["within_caps"])
        with job.part("background") as part:
            if predicted["within_caps"]:
                refined, tets, owners, _, stats = preflight.subdivide(points, cells, faces, parents, plan,
                    progress=lambda fraction, nodes, tets: part.update(.6*fraction, "Subdividing only the background", nodes=nodes, tets=tets))
                assert len(refined) == predicted["nodes"] and len(tets) == predicted["tetrahedra"]
                before, old_det, _ = preflight.metrics(points, cells)
                after, new_det, _ = preflight.metrics(refined, tets)
                error = float(np.max(abs(np.bincount(owners, weights=new_det)-old_det)/old_det))
                assert error < 1e-10
                result["background"] = dict(before=before, after=after, subdivision=stats,
                    maximum_parent_volume_error=error, passed=bool(after["quality_floor_met"] and after["within_caps"]))
                np.savez_compressed(a.run/"background-geometry-round.npz", points=refined, tetrahedra=tets, original_tet=owners)
                result["background_sha256"] = w.digest(a.run/"background-geometry-round.npz")
                part.update(1, "Checked background quality, volume and complete periodic caps", nodes=len(refined), tets=len(tets), minimum_q=after["minimum_quality"], below_q=after["below_quality_gate"])
            else:
                result["background"] = dict(passed=False, constructed=False, reason="Exact proposed background count exceeds unchanged caps")
                part.update(1, "Rejected background construction before exceeding caps", nodes=predicted["nodes"], tets=predicted["tetrahedra"])
        with job.part("witness-observation") as part:
            # No witness content is opened until every field-only choice is frozen.
            witnesses = json.loads((a.witnesses/"result.json").read_text())
            rows = []
            face_lookup = {tuple(face): i for i, face in enumerate(faces)}
            for i, witness in enumerate(witnesses["witnesses"]):
                containing = witness["neighborhood"]["backgrounds"]["refined"]["center"]["containing_tets"]
                surrounding = sorted({face_lookup[tuple(sorted(face))] for tet in cells[containing] for face in combinations(tet, 3)})
                rows.append(dict(query_index=witness["query_index"], containing_tets=containing,
                    class_rule_tet=bool(previous["selected_tets"][containing].any()),
                    class_rule_face=bool(previous["selected_faces"][surrounding].any()),
                    geometry_rule_tet=bool(tet_estimates["selected"][containing].any()),
                    geometry_rule_face=bool(face_estimates["selected"][surrounding].any()),
                    tet_estimates=tet_estimates["maximum"][containing].tolist(),
                    face_estimates=face_estimates["maximum"][surrounding].tolist()))
                part.update((i+1)/len(witnesses["witnesses"]), "Observed frozen witness neighborhoods after selection", completed=i+1, witnesses=len(witnesses["witnesses"]))
            result["witness_observation"] = rows
            job.result("witness-observation", class_covered=sum(r["class_rule_tet"] or r["class_rule_face"] for r in rows),
                       geometry_covered=sum(r["geometry_rule_tet"] or r["geometry_rule_face"] for r in rows))
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    assert w.digest(a.run/"selection.npz") == result["selection_sha256"]
    result.update(complete=True, input_hashes_verified=True,
        interpretation="Background-only preflight of one geometric sampling rule. Counts do not bound final Tet10 size; normalized residuals do not certify boundary distance. No native meshing or FEA claim.")
    w.save(a.run/"result.json", result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-native-background-v4/one-round")
    p.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    p.add_argument("--witnesses", type=Path, default=w.HERE/"runs/macbook-pmy-residual-probe-v1")
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    for name in ("run", "source", "fixtures", "witnesses"):
        setattr(a, name, getattr(a, name).resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run/"input-pmy_interpolation_control.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        a.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=TIMEOUT)
    w.save(a.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    findings = json.loads((a.run/"result.json").read_text())
    result("PMY interpolation selection", extra_faces=findings["selection"]["faces"]["geometry_only"],
           extra_tets=findings["selection"]["tets"]["geometry_only"])
    result("PMY background preflight", nodes=findings["predicted_background"]["nodes"],
           tets=findings["predicted_background"]["tetrahedra"], passed=findings["background"]["passed"])
    result("PMY witness coverage", class_covered=sum(r["class_rule_tet"] or r["class_rule_face"] for r in findings["witness_observation"]),
           geometry_covered=sum(r["geometry_rule_tet"] or r["geometry_rule_face"] for r in findings["witness_observation"]))


if __name__ == "__main__":
    main()
