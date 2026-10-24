"""Replay every saved preferred-only refinement and its original-domain maps."""

import json
from pathlib import Path

import numpy as np
import microgen_workflow as workflow


def verify_preflight(run, result):
    workflow.imports()
    from nested_refinement import interpolation_check
    from periodic_crease_audit import features
    from hp_shared_domain import compose, periodic_jumps

    source = Path(result["source"])
    original = dict(np.load(source / "mesh.npz"))
    original_fields = np.load(source / "fields.npz")["displacement"]
    current, fields = original, original_fields
    roots = np.arange(len(original["tetrahedra"]))
    mappings = np.repeat(np.eye(4)[None], len(roots), axis=0)
    offset = 0
    if result.get("continued_from"):
        continued = Path(result["continued_from"])
        prior = json.loads((continued / "result.json").read_text())
        assert prior["qualified"] and (prior["stop_reason"] == "round-limit" or (prior.get("projection_reaudit") and prior["stop_reason"] == "auditor-updated-geometry-passed"))
        assert Path(prior["source"]) == source and not prior["equilibrium_solve_performed"]
        current = dict(np.load(continued / "mesh.npz"))
        transfer = np.load(continued / "prolongation.npz")
        fields, roots, mappings = (transfer[key] for key in ("displacement", "parents", "parent_barycentric_maps"))
        interpolation_check(original["points"], original["tetrahedra"], original_fields,
            {**current, "displacement": fields, "parents": roots, "parent_barycentric_maps": mappings})
        offset = prior["total_rounds"] if "total_rounds" in prior else prior["completed_rounds"]
    pin = int(np.linalg.norm(original["points"], axis=1).argmin())
    assert pin == result["original_pin"]
    reports = []
    for row in result["rounds"]:
        assert row["level"] == offset+len(reports)+1
        directory = run / f"level-{row['level']}"
        assert workflow.digest(directory / "mesh.npz") == row["mesh_sha256"]
        assert workflow.digest(directory / "prolongation.npz") == row["prolongation_sha256"]
        mesh, transfer = dict(np.load(directory / "mesh.npz")), dict(np.load(directory / "prolongation.npz"))
        selected = np.bincount(transfer["parents"], minlength=len(current["tetrahedra"])) > 1
        feature = features(current)
        assert not (selected & ~feature["groups"]["curve_vertex_incident"]).any()
        assert int(selected.sum()) == row["bisected_curve_vertex_elements"]
        assert row["noncurve_bisected_elements"] == 0
        half = transfer["parent_barycentric_maps"] == .5
        child, replaced = np.nonzero(half.sum(axis=2) == 2)
        assert len(child) == 2 * int(selected.sum())
        endpoints = np.nonzero(half[child, replaced])[1].reshape(-1, 2)
        edges = current["tetrahedra"][transfer["parents"][child, None], endpoints]
        assert np.isin(edges, feature["curve_endpoints"]).any(axis=1).all()
        immediate = interpolation_check(current["points"], current["tetrahedra"], fields, {**mesh, **transfer})
        roots, mappings = compose(roots, mappings, transfer["parents"], transfer["parent_barycentric_maps"])
        original_check = interpolation_check(original["points"], original["tetrahedra"], original_fields,
            {**mesh, "displacement": transfer["displacement"], "parents": roots, "parent_barycentric_maps": mappings})
        assert len(mesh["points"]) == row["nodes"] <= 250000
        assert len(mesh["tetrahedra"]) == row["elements"] <= workflow.CAP
        assert len(mesh["points"])-len(current["points"]) == row["added_nodes"]
        geometry = workflow.geometry(run, "smooth-lattice", mesh["points"], mesh["tetrahedra"], mesh["surface"],
            box_projection_fallback=result.get("box_projection_fallback", False))
        assert geometry["geometry_admitted"] == row["geometry"]["geometry_admitted"]
        assert geometry["projection_failed"] == row["geometry"]["projection_failed"]
        if result.get("box_projection_fallback"):
            verify_box_witnesses(mesh, "smooth-lattice", row["geometry"])
        for key in ("sampled_two_sided_max", "relative_volume_error", "q_min_corner", "minimum_sampled_jacobian"):
            np.testing.assert_allclose(geometry[key], row["geometry"][key], atol=1e-13)
        periodic_jumps(mesh["points"], mesh["surface"], transfer["displacement"])
        np.testing.assert_array_equal(mesh["points"][pin], original["points"][pin])
        assert int(np.linalg.norm(mesh["points"], axis=1).argmin()) == pin
        assert row["selection"]["preferred_only"]
        assert row["selection"]["accepted_edge_orbits"] == row["selection"]["accepted_preferred_edge_orbits"]
        reports.append(dict(level=row["level"], nodes=len(mesh["points"]),
            bisected_curve_parents=int(selected.sum()), immediate=immediate, original=original_check,
            geometry_admitted=geometry["geometry_admitted"], failed_projections=geometry["projection_failed"]))
        current, fields = mesh, transfer["displacement"]
    final, transfer = np.load(run / "mesh.npz"), np.load(run / "prolongation.npz")
    for key in ("points", "tetrahedra", "surface"):
        np.testing.assert_array_equal(final[key], current[key])
    for key, expected in (("parents", roots), ("parent_barycentric_maps", mappings), ("displacement", fields)):
        np.testing.assert_array_equal(transfer[key], expected)
    assert not result["equilibrium_solve_performed"] and not result["native_candidate_promoted"] and not result["reference_promoted"]
    assert result["qualified"] == reports[-1]["geometry_admitted"]
    assert len(result["rounds"]) == result["completed_rounds"] <= result["maximum_rounds"] == 8
    assert (result["total_rounds"] if "total_rounds" in result else result["completed_rounds"]) == offset+len(reports)
    return dict(run=run.name, rounds=reports, original_pin=pin, stop=result["stop_reason"])


def verify_projection_witnesses(run, result):
    _, _, _, _, _, scorecard = workflow.imports()
    source = Path(result["source"])
    row = json.loads((source / "result.json").read_text())
    mesh = np.load(source / "mesh.npz")
    assert len(result["witnesses"]) == result["failed_projection_count"] == row["geometry"]["projection_failed"] > 0
    case, _, _ = workflow.specimen(result["case"])
    feasible = []
    for witness in result["witnesses"]:
        face = mesh["surface"][witness["face"]]
        weights = workflow.tri_values(scorecard.BARY)[witness["probe_index"]]
        probe = weights @ mesh["points"][face]
        np.testing.assert_allclose(probe, witness["probe"], atol=1e-15)
        np.testing.assert_array_equal(workflow.on_box_plane(abs(probe), .5), witness["locked"])
        old = np.asarray(witness["original_foot"])
        residual = abs(case.function(*old[None].T)[0]-witness["target_level"])/max(np.linalg.norm(scorecard.gradient(case.function, old[None])[0]), 1e-30)
        assert residual > 1e-9 or abs(old).max() > .5+1e-9
        np.testing.assert_allclose(np.linalg.norm(old-probe), witness["original_distance"], atol=1e-14)
        candidate = np.asarray(witness["constrained_foot"])
        value = abs(float(case.function(*candidate[None].T)[0]-witness["target_level"]))
        valid = bool(value < 1e-9 and abs(candidate).max() <= .5+1e-12)
        assert valid == witness["constrained_control"]["feasible_target"]
        np.testing.assert_allclose(np.linalg.norm(candidate-probe), witness["constrained_control"]["distance"], atol=1e-14)
        feasible.append(valid)
    assert all(feasible) == result["all_constrained_candidates_feasible"]
    assert result["source_rejection_preserved"] and not result["qualified"]
    assert not result["mesh_modified"] and not result["samples_regenerated"] and not result["stiffness_solve_performed"]
    return dict(run=run.name, failed_samples=len(feasible), feasible_bounded_targets=sum(feasible), source_rejection_preserved=True)


def verify_box_witnesses(mesh, name, geometry):
    _, _, _, _, _, scorecard = workflow.imports()
    triangle = mesh["points"][mesh["surface"][:, :3]]
    caps = np.any(np.all(workflow.on_box_plane(triangle, -.5), axis=1) | np.all(workflow.on_box_plane(triangle, .5), axis=1), axis=1)
    walls = mesh["surface"][~caps]
    weights = workflow.tri_values(scorecard.BARY)
    case, _, _ = workflow.specimen(name)
    levels = np.array([level for level in case.levels if np.isfinite(level)])
    recovered = 0
    attempts = geometry["box_projection_witnesses"]
    assert len(attempts) == geometry["box_projection_attempts"] == geometry["unbounded_projection_failed"]
    assert len({witness["sample_index"] for witness in attempts}) == len(attempts)
    for witness in attempts:
        face, local = divmod(witness["sample_index"], len(weights))
        probe = weights[local] @ mesh["points"][walls[face]]
        locked = workflow.on_box_plane(abs(probe), .5)
        level = levels[np.argmin(abs(case.function(*probe[None].T)[0]-levels))]
        old = np.asarray(witness["original_foot"])
        old_residual = abs(case.function(*old[None].T)[0]-level)/max(np.linalg.norm(scorecard.gradient(case.function, old[None])[0]), 1e-30)
        assert old_residual > 1e-9 or abs(old).max() > .5+1e-9
        candidate = np.asarray(witness["candidate"])
        residual = abs(case.function(*candidate[None].T)[0]-level)/max(np.linalg.norm(scorecard.gradient(case.function, candidate[None])[0]), 1e-30)
        valid = bool(np.isfinite(candidate).all() and residual <= 1e-9 and abs(candidate).max() <= .5+1e-12
            and np.all(abs(candidate[locked]-np.sign(probe[locked])*.5) <= 1e-12))
        assert valid == witness["feasible_target"]
        recovered += int(valid)
        np.testing.assert_allclose(np.linalg.norm(candidate-probe), witness["distance"], atol=1e-14)
    assert recovered == geometry["box_projection_recovered"]
    assert len(attempts)-recovered == geometry["projection_failed"]
    return recovered


def verify_projection_reaudit(run, result):
    workflow.imports()
    source = Path(result["audited_source"])
    old = json.loads((source / "result.json").read_text())
    assert not old["qualified"] and old["geometry"]["projection_failed"] > 0
    for name in ("mesh.npz", "prolongation.npz", result["case"]+"-samples.npz", result["case"]+"-specimen.json"):
        assert workflow.digest(run / name) == workflow.digest(source / name)
    mesh = np.load(run / "mesh.npz")
    legacy = workflow.geometry(run, result["case"], mesh["points"], mesh["tetrahedra"], mesh["surface"])
    corrected = workflow.geometry(run, result["case"], mesh["points"], mesh["tetrahedra"], mesh["surface"], box_projection_fallback=True)
    assert not legacy["geometry_admitted"] and legacy["projection_failed"] == old["geometry"]["projection_failed"]
    assert corrected["geometry_admitted"] == result["qualified"] and corrected["projection_failed"] == result["geometry"]["projection_failed"]
    assert result["geometry"]["projection_method"] == "newton-box-slsqp-v1"
    assert corrected["box_projection_attempts"] == legacy["projection_failed"]
    for key in ("sampled_two_sided_max", "relative_volume_error", "q_min_corner", "minimum_sampled_jacobian", "forward_samples", "reverse_samples"):
        np.testing.assert_allclose(corrected[key], result["geometry"][key], atol=1e-13)
        np.testing.assert_allclose(legacy[key], result["legacy_geometry"][key], atol=1e-13)
    recovered = verify_box_witnesses(mesh, result["case"], result["geometry"])
    assert result["source_rejection_preserved"] and not result["equilibrium_solve_performed"]
    return dict(run=run.name, legacy_failed=legacy["projection_failed"], recovered=recovered,
        qualified=result["qualified"], identical_mesh_fields_and_samples=True)
