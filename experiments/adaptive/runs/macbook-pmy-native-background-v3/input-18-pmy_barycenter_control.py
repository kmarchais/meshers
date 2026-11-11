"""One round of generic preclip classification checks on the frozen PMY patch.

Compare actual field values with the original affine interpolant at every face
and tet center. Insert only mismatches, skipping outer faces. Decisions use the
field and background alone; saved neck witnesses are read afterward as observers.
This is a local diagnostic, with unchanged admission gates and no admitted mesh.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w
import pmy_neck_control as neck


NAME = "tpms-pmy-sheet"
RUNS = w.HERE/"runs"
NODE_LIMIT = 2048
TET_LIMIT = 20000


def classify(values, levels):
    import numpy as np

    values = np.asarray(values)
    assert np.isfinite(values).all()
    return np.where(values < levels[0], -1, np.where(values > levels[1], 1, 0))


def selection(case, points, cells):
    """The complete one-round decision, made before any point is inserted."""
    import numpy as np

    faces, incidence = np.unique(np.sort(cells[:, np.asarray(neck.FACES)].reshape(-1, 3), axis=1),
                                 axis=0, return_counts=True)
    assert incidence.max() <= 2
    values = np.asarray(case.function(*points.T))
    records = []
    for kind, entities in (("face", faces), ("tet", cells)):
        centers = points[entities].mean(1)
        actual = np.asarray(case.function(*centers.T))
        linear = values[entities].mean(1)
        observed_class, affine_class = classify(actual, case.levels), classify(linear, case.levels)
        for index in np.flatnonzero(observed_class != affine_class):
            outer = kind == "face" and incidence[index] == 1
            records.append(dict(kind=kind, index=int(index), vertices=entities[index].tolist(),
                                point=centers[index].tolist(), actual_value=float(actual[index]),
                                affine_value=float(linear[index]), actual_class=int(observed_class[index]),
                                affine_class=int(affine_class[index]), selected=not outer,
                                reason="outer-facet-preservation" if outer else "classification-mismatch"))
    return dict(face_samples=len(faces), tet_samples=len(cells), records=records,
                class_labels={"-1": "below lower wall", "0": "inside band", "1": "above upper wall"},
                selected_faces=sum(r["selected"] and r["kind"] == "face" for r in records),
                selected_tets=sum(r["selected"] and r["kind"] == "tet" for r in records),
                skipped_outer_faces=sum(not r["selected"] for r in records),
                selection_scope="Original face and tet centers only; one frozen round; no witness coordinates")


def insert(points, cells, record):
    """Insert a face center in its full star, or a tet center's current carrier."""
    import numpy as np

    query = np.array(record["point"])
    distance = np.linalg.norm(points-query, axis=1)
    if distance.min() < 1e-12:
        return points, cells, dict(node=int(distance.argmin()), added=False, incident_tets=0)
    if record["kind"] == "face":
        carrier = np.array(record["vertices"])
        incident = np.flatnonzero(np.isin(cells, carrier).sum(1) == 3)
        assert len(incident) == 2, "An original internal face must still have two sides"
    else:
        hits, bary = neck.containing(points, cells, query, 1e-9)
        assert len(hits), "Original tet center escaped the unchanged background"
        carrier = cells[hits[0]][bary[0] > 1e-9]
        assert 2 <= len(carrier) <= 4
        incident = np.flatnonzero(np.isin(cells, carrier).sum(1) == len(carrier))
    new_id = len(points)
    planned = [np.where(tet == node, new_id, tet) for tet in cells[incident] for node in carrier]
    points = np.vstack((points, query))
    cells = np.vstack((np.delete(cells, incident, axis=0), neck.oriented(points, planned)))
    assert len(points) <= NODE_LIMIT and len(cells) <= TET_LIMIT, "Declared one-round local bound exceeded"
    return points, cells, dict(node=new_id, added=True, incident_tets=len(incident), carrier_size=len(carrier))


def refine(case, points, cells, plan, progress=None):
    import numpy as np
    from bistellar_probe import oriented_boundary

    original_points = points.copy()
    original_boundary = oriented_boundary(cells)
    before = neck.measure(points, cells)
    records = [r for r in plan["records"] if r["selected"]]
    operations = []
    for i, record in enumerate(records):
        points, cells, observation = insert(points, cells, record)
        actual = float(case.function(*points[observation["node"]]))
        np.testing.assert_allclose(actual, record["actual_value"], atol=1e-13, rtol=0)
        operations.append(dict(kind=record["kind"], original_index=record["index"], **observation))
        if progress is not None:
            progress((i+1)/len(records), len(points), len(cells))
    np.testing.assert_array_equal(points[:len(original_points)], original_points)
    assert oriented_boundary(cells) == original_boundary
    after = neck.measure(points, cells)
    assert after["oriented_faces_consistent"] and after["nonmanifold_faces"] == 0
    assert after["minimum_determinant"] > 0
    np.testing.assert_allclose(after["volume"], before["volume"], atol=1e-14, rtol=1e-12)
    return points, cells, operations


def controls():
    import numpy as np
    from scorecard import Case

    np.testing.assert_array_equal(classify([-.4, -.3, 0, .3, .4], (-.3, .3)), [-1, 0, 0, 0, 1])
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [0., 0., -1.]])
    cells = np.array([[0, 1, 2, 3], [0, 2, 1, 4]])
    affine = Case("affine", lambda x, y, z: x, (-.2, .8))
    assert not selection(affine, points, cells)["records"]
    face_case = Case("internal-face", lambda x, y, z: .01-(x-1/3)**2-(y-1/3)**2-z*z, (-2., 0.))
    plan = selection(face_case, points, cells)
    assert plan["selected_faces"] == 1 and plan["selected_tets"] == 0
    p, t, operations = refine(face_case, points, cells, plan)
    assert len(p) == 6 and len(t) == 6 and len(operations) == 1
    outer = selection(face_case, points[:4], cells[:1])
    assert outer["skipped_outer_faces"] == 1 and outer["selected_faces"] == 0
    for sign in (-1, 1):
        def hidden(x, y, z):
            return sign*(.02-(x-.25)**2-(y-.25)**2-(z-.25)**2)
        levels = (-2., 0.) if sign == 1 else (0., 2.)
        case = Case("hidden", hidden, levels)
        plan = selection(case, points[:4], cells[:1])
        assert plan["selected_tets"] == 1
        p, t, _ = refine(case, points[:4], cells[:1], plan)
        clipped, clipped_tets = neck.clip(case, p, t)
        assert not len(neck.containing(clipped, clipped_tets, np.full(3, .25))[0])
    return dict(all_three_band_classes=True, affine_field_unchanged=True,
                complete_internal_face_star=True, outer_mismatch_skipped=True,
                upper_hidden_void_detected=True, lower_hidden_void_detected=True)


def observe_witnesses(source, case, arms):
    """Observers are loaded only after the refinement decision and artifacts exist."""
    import numpy as np

    historical = json.loads((source/"result.json").read_text())
    point = np.array(historical["critical"]["point"])
    ends = historical["void_path_endpoints"]
    probes = np.vstack([np.array(end["point"])+np.linspace(0, 1, 129)[:, None]*(point-np.array(end["point"])) for end in ends])
    assert (classify(case.function(*probes.T), case.levels) == 1).all()
    results = {}
    for name, (points, cells) in arms.items():
        occupied = [i for i, p in enumerate(probes) if len(neck.containing(points, cells, p)[0])]
        results[name] = dict(samples=len(probes), material_samples=len(occupied), sample_indices=occupied,
                             critical_point_in_material=bool(len(neck.containing(points, cells, point)[0])))
    assert results["baseline"] == dict(samples=258, **historical["void_path_checks"]["before"])
    return results


def worker(a):
    import numpy as np
    from microgen.shape import surface_functions
    import tpms_family_fixtures as fixtures

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    archived = a.source/"background-before.npz"
    protocol = json.loads((a.fixtures/"protocol.json").read_text())
    frozen = json.loads((a.fixtures/"result.json").read_text())
    known = {Path(p).resolve(): h for p, h in frozen["inputs"].items()}
    fixture = a.fixtures/NAME
    assert all(w.digest(p) == known[p.resolve()] for p in fixture.iterdir())
    expected = next(h for p, h in protocol["source_hashes"].items() if Path(p).name == "surface_functions.py")
    assert w.digest(Path(surface_functions.__file__)) == expected
    files = [Path(__file__), Path(w.__file__), Path(neck.__file__), Path(fixtures.__file__),
             Path(surface_functions.__file__), archived, a.source/"result.json",
             a.fixtures/"protocol.json", a.fixtures/"result.json", *fixture.iterdir()]
    files += [w.HERE/n for n in ("scorecard.py", "bistellar_probe.py", "mmg_quality_localization.py")]
    hashes = {str(p.resolve()): w.digest(p) for p in files}
    with np.load(archived) as mesh:
        original_points, original_cells = mesh["points"], mesh["tetrahedra"]
    assert original_points.shape == (127, 3) and original_cells.shape == (396, 4)
    case = fixtures.case(NAME, a.fixtures)
    result = dict(inputs=hashes, rounds=1, physical_specimen_unchanged=True, original_nodes_unchanged=True,
                  outer_facets_unchanged=True, admitted=False, admission_gates_unchanged=True,
                  local_bounds=dict(nodes=NODE_LIMIT, tets=TET_LIMIT),
                  scope="One field-only barycenter-sign refinement on an archived local background. No global candidate.")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=180, memory_bytes=w.MEMORY,
           numerical_threads=w.NUMERICAL_THREAD_POLICY, local_bounds=result["local_bounds"],
           selection="All original face/tet centers; compare actual and affine three-way band classes; skip outer faces"))
    with Job(parts={"controls": 1, "selection": 1, "subdivision": 2, "clip-and-observe": 3}) as job:
        with job.part("controls") as part:
            result["controls"] = controls()
            part.update(1, "Classification and conformity controls passed", controls=len(result["controls"]))
        job.result("controls", passed=True)
        with job.part("selection") as part:
            plan = selection(case, original_points, original_cells)
            result["selection"] = plan
            w.save(a.run/"selection.json", plan)
            part.update(1, "Frozen one-round sample selection", faces=plan["selected_faces"], tets=plan["selected_tets"], skipped_outer=plan["skipped_outer_faces"])
        job.result("selection", faces=plan["selected_faces"], tets=plan["selected_tets"], skipped_outer=plan["skipped_outer_faces"])
        with job.part("subdivision") as part:
            points, cells, operations = refine(case, original_points.copy(), original_cells.copy(), plan,
                progress=lambda fraction, nodes, tets: part.update(fraction, "Conforming barycenter insertion", nodes=nodes, tets=tets))
            result["operations"] = operations
            result["background"] = dict(baseline=neck.measure(original_points, original_cells), one_round=neck.measure(points, cells))
            np.savez_compressed(a.run/"background-one-round.npz", points=points, tetrahedra=cells)
            part.update(1, "Preserved original vertices, outer facets and volume", nodes=len(points), tets=len(cells), minimum_q=result["background"]["one_round"]["minimum_quality"])
        job.result("subdivision", nodes=len(points), tets=len(cells), minimum_q=result["background"]["one_round"]["minimum_quality"])
        with job.part("clip-and-observe") as part:
            arms = {}
            result["clipped"] = {}
            for i, (name, p, t) in enumerate((("baseline", original_points, original_cells), ("one_round", points, cells))):
                cp, ct = neck.clip(case, p, t)
                arms[name] = (cp, ct)
                result["clipped"][name] = neck.measure(cp, ct)
                np.savez_compressed(a.run/("clipped-"+name+".npz"), points=cp, tetrahedra=ct)
                part.update((i+1)/3, "Clipped the unchanged field", nodes=len(cp), tets=len(ct), below_q=result["clipped"][name]["below_quality_gate"])
            result["void_path_observer"] = observe_witnesses(a.source, case, arms)
            original = result["void_path_observer"]["baseline"]["material_samples"]
            remaining = result["void_path_observer"]["one_round"]["material_samples"]
            result["false_material_samples_removed"] = original-remaining
            part.update(1, "Observed fixed path after all refinement choices", remaining=remaining, removed=original-remaining)
        job.result("clip-and-observe", remaining_false_samples=remaining, removed_false_samples=original-remaining,
                   below_q=result["clipped"]["one_round"]["below_quality_gate"])
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    result.update(complete=True, input_hashes_verified=True,
                  limits="Finite barycenter tests can miss features between samples. Subdivision does not guarantee quality after clipping. Path witnesses observe this one local control; they do not choose insertion points or establish global topology.")
    w.save(a.run/"result.json", result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--source", type=Path, default=RUNS/"macbook-pmy-neck-control-v1")
    p.add_argument("--fixtures", type=Path, default=RUNS/"macbook-tpms-family-fixtures-v1")
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    for name in ("run", "source", "fixtures"):
        setattr(a, name, getattr(a, name).resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run/"input-pmy_barycenter_control.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        a.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=180)
    w.save(a.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"


if __name__ == "__main__":
    main()
