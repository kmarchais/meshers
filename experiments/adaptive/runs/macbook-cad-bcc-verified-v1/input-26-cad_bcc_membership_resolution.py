"""Audit the original BCC cap ON ambiguities using only frozen records.

This does not reclassify CAD, modify tolerances, sample a new point, or establish
global CAD/field equivalence. Exact rational inequalities refer to the recorded
cylinder coordinates; CAD boolean semantics remain independently documented.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w
import cad_lattice_fixture as fixture

SOURCE = w.HERE/"runs/macbook-cad-bcc-coverage-v1"
DIAGNOSIS = w.HERE/"runs/macbook-cad-bcc-classifier-v1"


def cylinder_margins(points, records):
    import numpy as np
    offset = points[:, None, :]-np.array([row["center"] for row in records])[None, :, :]
    local = np.einsum("nki,kij->nkj", offset, np.array([row["rotation"] for row in records]))
    radial = np.array([row["radius"] for row in records])[None, :]-np.linalg.norm(local[:, :, 1:], axis=2)
    axial = np.array([row["height"] for row in records])[None, :]/2-abs(local[:, :, 0])
    return radial, axial


def exact_membership(point, record):
    """Strict polynomial membership tests, exact on the represented numbers."""
    q = [F(float(x))-F(float(c)) for x, c in zip(point, record["center"], strict=True)]
    rotation = [[F(float(value)) for value in row] for row in record["rotation"]]
    local = [sum(q[i]*rotation[i][j] for i in range(3)) for j in range(3)]
    axial = (F(float(record["height"]))/2)**2-local[0]**2
    radial = F(float(record["radius"]))**2-local[1]**2-local[2]**2
    a, b = [[F(float(value)) for value in row] for row in record["endpoints"]]
    d = [y-x for x, y in zip(a, b, strict=True)]
    v = [F(float(x))-y for x, y in zip(point, a, strict=True)]
    d2, dot = sum(x*x for x in d), sum(x*y for x, y in zip(v, d, strict=True))
    anchored_radial = F(float(record["radius"]))**2-sum(x*x for x in v)+dot*dot/d2
    return dict(raw_matrix_axial_squared_clearance=str(axial), raw_matrix_radial_squared_clearance=str(radial),
        anchored_segment_lower=str(dot), anchored_segment_upper=str(d2-dot),
        anchored_radial_squared_clearance=str(anchored_radial),
        raw_matrix_strict_interior=axial > 0 and radial > 0,
        anchored_strict_interior=dot > 0 and d2-dot > 0 and anchored_radial > 0)


def worker(args):
    import numpy as np
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    source, diagnosis = args.source, args.diagnosis
    names = {source: ("result.json", "primitive-receipt.json", "raw-primitive-receipt.json", "membership-probes.npz",
                     "membership-probes-input.npz", "actual-bcc.step"),
        diagnosis: ("result.json", "manifest.json", "fresh-inventory.json", "fresh-distances.json",
                    "fresh-classifications.json", "step-classifications.json", "step-inventory.json")}
    paths = [directory/name for directory, files in names.items() for name in files]
    paths.extend((Path(__file__).resolve(), Path(w.__file__).resolve(), Path(fixture.__file__).resolve()))
    hashes = {str(path): w.digest(path) for path in paths}
    for path in paths:
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/("input-"+path.name))
    old = json.loads((source/"result.json").read_text())
    evidence = json.loads((diagnosis/"result.json").read_text())
    assert evidence["complete"] and evidence["input_hashes_verified"]
    assert all(row["passed"] for row in evidence["controls"].values())
    assert w.digest(source/"actual-bcc.step") == old["cad"]["step_sha256"]
    assert not old["fixture_verified"] and not old["native_attempted"]
    w.save(args.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=60, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, new_cad_operations=False, new_samples=False,
        scope="Read-only fixture-verification proposal; original classifier failures remain archived"))
    with Job(parts={"load": 1, "ambiguities": 2, "archive": 1}) as job:
        with job.part("load") as part:
            with np.load(source/"membership-probes.npz") as saved:
                points, states, normalized, raw = [saved[key] for key in ("points", "cad_states", "normalized_field", "raw_field")]
            recipe = json.loads((source/"primitive-receipt.json").read_text())
            records = recipe["ordered_images"]
            recomputed = fixture.field(recipe["normalized_primitives"])(*points.T)
            recomputed = np.maximum(recomputed, abs(points).max(1)-.5)
            np.testing.assert_array_equal(recomputed, normalized)
            raw_recomputed = fixture.field(records, preserve_rotations=True)(*points.T)
            np.testing.assert_array_equal(np.maximum(raw_recomputed, abs(points).max(1)-.5), raw)
            radial, axial = cylinder_margins(points, records)
            best = np.minimum(radial, axial).argmax(1)
            mismatch_ids = {row["index"] for row in old["membership"]["mismatches"]}
            distances = {row["index"]: row for row in json.loads((diagnosis/"fresh-distances.json").read_text())}
            inventory = json.loads((diagnosis/"fresh-inventory.json").read_text())
            step = json.loads((diagnosis/"step-classifications.json").read_text())
            fresh = json.loads((diagnosis/"fresh-classifications.json").read_text())
            assert inventory["valid"] and inventory["solid_count"] == 1 and inventory["shells_closed"]
            assert not inventory["central_support_planes"]
            part.update(1, "Loaded every unchanged membership query", queries=len(points), ambiguities=len(mismatch_ids))
        rows, unresolved = [], []
        with job.part("ambiguities") as part:
            for index, point in enumerate(points):
                original = str(states[index])
                expected = "inside" if normalized[index] < 0 else "outside" if normalized[index] > 0 else "boundary"
                row = dict(index=index, point=point.tolist(), original_cad_state=original,
                    normalized_field=float(normalized[index]), raw_field=float(raw[index]),
                    fresh_classifier=fresh[index], step_corroboration=step[index], proposed_state=original,
                    resolution="unchanged-original-agreement")
                if index in mismatch_ids:
                    measured, record, k = distances[index], records[int(best[index])], int(best[index])
                    exact = exact_membership(point, record)
                    box = max(abs(F(float(x))) for x in point)-F(1, 2)
                    box_distance = float(np.linalg.norm(np.maximum(abs(point)-.5, 0))) if box > 0 else float(-box)
                    vertex = measured["nearest_vertex"][0]
                    vertex_point = np.asarray(inventory["vertices"][vertex["vertex"]]["point"])
                    vertex_distance = float(np.linalg.norm(point-vertex_point))
                    edge = measured["nearest_edge"][0]
                    actual_feet = measured["boundary_brep"]["solutions"]
                    foot_radial, foot_axial = cylinder_margins(np.array([item["point"] for item in actual_feet]), records)
                    foot_margin = np.minimum(foot_radial, foot_axial).max(1)
                    guard = fixture.ARITHMETIC*max(1., record["height"], float(np.linalg.norm(point-np.asarray(record["center"]))))
                    footprint = any(np.max(abs(np.asarray(item["point"]))) <= .5+guard and
                        np.min(abs(abs(np.asarray(item["point"]))-.5)) <= guard and
                        abs(item["euclidean_distance"]-box_distance) <= 1e-12 for item in actual_feet)
                    checks = dict(original_and_fresh_on=original == "boundary" and fresh[index]["reused"] == "boundary" and fresh[index]["fresh"] == "boundary",
                        unambiguous_box_side=box != 0, raw_cylinder_exact=exact["raw_matrix_strict_interior"],
                        anchored_cylinder_exact=exact["anchored_strict_interior"],
                        raw_margin_exceeds_recorded_arithmetic_scale=min(radial[index, k], axial[index, k]) > guard,
                        explicit_corner_vertex_within_own_tolerance=bool((abs(vertex_point) == .5).all() and vertex_distance <= vertex["vertex_tolerance"]),
                        known_cap_edge_point_within_own_tolerance=bool(np.count_nonzero(abs(abs(np.asarray(edge["foot"]))-.5) <= guard) >= 2 and edge["distance"] <= edge["edge_tolerance"]),
                        cad_distance_agrees_box_distance=measured["boundary_brep"]["done"] and abs(measured["boundary_brep"]["value"]-box_distance) <= 1e-12,
                        explicit_physical_cap_foot=footprint, cap_feet_in_cylinder_union=bool((foot_margin >= -guard).all()),
                        step_corroborates=step[index]["reused"] == expected and step[index]["fresh"] == expected,
                        analytic_side_agrees_field=expected == ("outside" if box > 0 else "inside"))
                    checks = {key: bool(value) for key, value in checks.items()}
                    accepted = all(checks.values())
                    row.update(resolution="cap-tolerance-ambiguity-resolved-by-recipe" if accepted else "unsupported-ambiguity",
                        proposed_state=expected if accepted else original, checks=checks, exact_recorded_membership=exact,
                        source_id=record["source_id"], periodic_image=record["image"], original_image_index=k,
                        raw_radial_margin=float(radial[index, k]), raw_axial_margin=float(axial[index, k]), arithmetic_scale=guard,
                        exact_box_signed_margin=str(box), analytic_domain_boundary_distance=box_distance,
                        saved_cad_boundary_distance=measured["boundary_brep"]["value"], physical_cap_feet=actual_feet,
                        cap_foot_union_margins=foot_margin.tolist(), edge_witness=edge, vertex_witness=vertex,
                        edge_extrema_failures_retained=len(measured["edge_extrema_failures"]))
                    if not accepted:
                        unresolved.append(index)
                else:
                    agrees = (original == "inside" and normalized[index] < 0) or (original == "outside" and normalized[index] > 0) or (original == "boundary" and abs(normalized[index]) <= 1e-8)
                    assert agrees, f"Unlisted original disagreement {index}"
                rows.append(row)
                if index % 128 == 0:
                    part.update(index/len(points), "Auditing fixed cap ambiguity evidence", queries=index, unresolved=len(unresolved))
            part.update(1, "Every original query retained", queries=len(points), unresolved=len(unresolved))
        with job.part("archive") as part:
            w.save(args.run/"queries.json", rows)
            counts = {name: sum(row["proposed_state"] == name for row in rows) for name in ("inside", "outside", "boundary")}
            result = dict(inputs=hashes, complete=True, original_queries=len(points), original_mismatch_count=len(mismatch_ids),
                resolved_ambiguities=len(mismatch_ids)-len(unresolved), unresolved_ambiguity_ids=unresolved,
                resolution_conditions_met=not unresolved, proposed_state_counts=counts, original_classifier_pass=False,
                fixture_verified=False, exact_cad_surface_sampling_pending=True, native_attempted=False,
                center=evidence["comparisons"]["fresh"]["center_distance"],
                policy="Retain original ON as tolerance ambiguity; resolve only cap probes with strict original-cylinder and box predicates, an explicit physical cap foot, local stored tolerance witnesses, and STEP corroboration. No CAD mutation or gate change.",
                limitations=["Finite query resolution does not prove global CAD/field equivalence.",
                    "Rational proofs apply to recorded matrix and endpoint-anchored cylinders, not every OCCT Boolean operation.",
                    "STEP alters topology tolerances and slightly changes numerical volume; it is corroboration, never replacement geometry.",
                    "Failed edge extrema remain failures; recorded edge feet are witnesses, not exhaustive closest-edge bounds.",
                    "Fresh original-CAD face/edge sampling and field consistency must pass before native meshing."])
            assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
            result["input_hashes_verified"] = True
            w.save(args.run/"result.json", result)
            part.update(1, "Archived resolution proposal without promoting fixture", resolved=result["resolved_ambiguities"])
            job.result("resolution", conditions_met=result["resolution_conditions_met"], resolved=result["resolved_ambiguities"], unresolved=len(unresolved))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--diagnosis", type=Path, default=DIAGNOSIS)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.source, args.diagnosis, args.run = args.source.resolve(), args.diagnosis.resolve(), args.run.resolve()
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import report, result
    report(0, "Auditing frozen BCC membership evidence", eta_s=60, memory_ceiling_gib=8)
    receipt = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=60)
    w.save(args.run/"worker.json", receipt)
    assert receipt["returncode"] == 0, "Read worker.log and worker.json"
    row = json.loads((args.run/"result.json").read_text())
    result("bcc-membership-resolution", conditions_met=row["resolution_conditions_met"],
        resolved=row["resolved_ambiguities"], unresolved=len(row["unresolved_ambiguity_ids"]), fixture_verified=False)


if __name__ == "__main__":
    main()
