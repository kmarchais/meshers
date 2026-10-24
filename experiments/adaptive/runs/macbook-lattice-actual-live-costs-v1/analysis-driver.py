"""Fresh actual lattice route with explicitly charged repair and star refinement."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

import numpy as np

import microgen_workflow as workflow


def prepare(source, output):
    from microgen_tet10_control import straight_elevation
    from nested_preflight import curve_vertex_star_edges
    from nested_refinement import refine, interpolation_check
    from periodic_crease_audit import features
    from homogenization import strain_tensor
    from hp_shared_domain import periodic_jumps
    from interior_relocation_trial import validate_repair_source

    repair = json.loads((source / "result.json").read_text())
    generation = Path(repair["source"])
    parent = json.loads((generation / "result.json").read_text())
    assert repair["case"] == "smooth-lattice" and repair["qualified"]
    assert repair["repair"]["remaining_bad_elements"] == 0 and repair["repair"]["boundary_coordinates_unchanged"]
    assert parent["stage_audit_enabled"] and not parent.get("saved_seed_replay", False)
    validate_repair_source(parent)
    assert parent["mmg_calls"] == repair["parent_mmg_calls"]
    assert workflow.digest(workflow.HERE / "microgen_tuned.py") == json.loads((generation / "manifest.json").read_text())["script_sha256"]
    paths = [Path(__file__)]
    paths.extend(workflow.HERE / name for name in ("microgen_tet10_control.py", "nested_preflight.py",
        "nested_refinement.py", "periodic_crease_audit.py", "hp_shared_domain.py", "test_straight_elevation.py",
        "test_nested_refinement.py", "microgen_tuned.py", "interior_relocation_trial.py"))
    paths.extend(directory / name for directory in (source, generation) for name in ("manifest.json", "result.json", "mesh.npz"))
    inputs = {}
    for directory in (source, generation):
        for name, expected in json.loads((directory / "manifest.json").read_text())["inputs"].items():
            path = Path(name)
            assert workflow.digest(path) == expected, str(path)
            inputs[str(path)] = expected
    for name in ("vtk-seed.json", "mmg-stage-1.json", "mmg-stage-2.json"):
        path = generation / name
        stage = json.loads(path.read_text())
        assert 0 < stage["elements"] <= workflow.CAP and stage["nodes"] <= 250000
        paths.append(path)
    for suffix in ("-specimen.json", "-samples.npz"):
        path = source / ("smooth-lattice" + suffix)
        assert workflow.digest(path) == inputs[str(generation / path.name)]
        paths.append(path)
        shutil.copyfile(path, output.parent / path.name)
    import microgen.remesh
    paths.extend((Path(microgen.remesh.__file__),
        Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))))
    inputs.update({str(path): workflow.digest(path) for path in paths})
    workflow.save(output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, output.parent / "analysis-driver.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_straight_elevation.py",
        "experiments/adaptive/test_nested_refinement.py"], check=True, timeout=60)
    start = time.perf_counter()
    original = np.load(source / "mesh.npz")
    elevation = straight_elevation(original["points"], original["tetrahedra"], original["surface"])
    if not elevation["within_node_cap"]:
        workflow.save(output, dict(qualified=False, status="straight-lift-node-cap", nodes=elevation["nodes"], source=str(source)))
        return
    points, cells, surface = (elevation[name] for name in ("points", "tetrahedra", "surface"))
    from periodic import periodic_pairs
    from meshers import Mesh
    matching, _, _ = periodic_pairs(points, surface)
    converted = Mesh(points, cells, surface, np.zeros(len(surface), dtype=np.uint8), matching,
        {"source": "Fresh actual VTK51/two-pass MMG, explicit fixed-boundary repair and straight elevation"}).to_fedoo()
    np.testing.assert_array_equal(converted.nodes, points)
    np.testing.assert_array_equal(converted.elements, cells)
    assert converted.elm_type == "tet10"
    fields = np.array([points @ strain_tensor(load).T for load in np.eye(6)])
    feature = features(dict(points=points, tetrahedra=cells, surface=surface))
    lift_s = time.perf_counter() - start
    phase = time.perf_counter()
    refined = refine(points, cells, surface, fields,
        preferred_edges=curve_vertex_star_edges(cells, feature["curve_endpoints"]))
    reproduction = interpolation_check(points, cells, fields, refined)
    refinement_s = time.perf_counter() - phase
    np.savez_compressed(output.parent / "elevated-mesh.npz", points=points, tetrahedra=cells, surface=surface)
    np.savez_compressed(output.parent / "mesh.npz", **{key: refined[key] for key in ("points", "tetrahedra", "surface")})
    np.savez_compressed(output.parent / "prolongation.npz", **{key: refined[key]
        for key in ("displacement", "parents", "parent_barycentric_maps")})
    geometry = workflow.geometry(output.parent, "smooth-lattice", refined["points"], refined["tetrahedra"], refined["surface"])
    periodic = periodic_jumps(refined["points"], refined["surface"], refined["displacement"])
    volume_change = geometry["volume"] / repair["geometry"]["volume"] - 1
    assert abs(volume_change) < 1e-10
    prepared_s = time.perf_counter() - start
    workflow.save(output, dict(case="smooth-lattice", route=repair["route"]+"-straight-tet10-curve-star",
        source=str(source), generation_source=str(generation), nodes=len(refined["points"]), elements=len(refined["tetrahedra"]),
        order=2, status="fresh-actual-preflight-passed" if geometry["geometry_admitted"] else "geometry-rejected",
        qualified=geometry["geometry_admitted"], geometry=geometry, periodic=periodic, reproduction=reproduction,
        refinement=refined["diagnostics"], relative_volume_change=volume_change,
        lifted_nodes=len(points), lift_and_conversion_s=lift_s, refinement_s=refinement_s, total_preparation_s=prepared_s,
        generation_and_audits_s=parent["total_workflow_s"], repair_and_audits_s=repair["total_saved_mesh_study_s"],
        charged_preparation_workflow_s=parent["total_workflow_s"]+repair["total_saved_mesh_study_s"]+prepared_s,
        fresh_actual_generation=True, saved_seed_replay=False, timing_comparable=False,
        extra_steps=["SLSQP interior repair", "manual straight Tet10/direct fedoo", "curve-vertex-incident star refinement"],
        control_fields="Six affine macrostrain fields for interpolation and periodic controls; no parent FE solve or equilibrium-field replay",
        field_solve_performed=False, reference_promoted=False,
        scope="Fresh Mac actual route through geometry preflight. Diagnostic intermediate audits are charged, so no production timing ratio; no stiffness or continuum claim."))
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


def solve_prepared(prepared, output):
    from hp_coupled_preflight import inherited_inputs
    from amg_reference import AuditedAMG
    from stress_recovery import solve, vector_first_reactions

    row = json.loads((prepared / "result.json").read_text())
    assert row["status"] == "fresh-actual-preflight-passed" and row["qualified"]
    assert row["fresh_actual_generation"] and not row["saved_seed_replay"] and not row["field_solve_performed"]
    paths = [Path(__file__), prepared / "result.json", prepared / "mesh.npz"]
    paths.extend(workflow.HERE / name for name in ("stress_recovery.py", "amg_reference.py", "hp_coupled_preflight.py"))
    inputs = {**inherited_inputs(prepared), **{str(path): workflow.digest(path) for path in paths}}
    for suffix in ("-specimen.json", "-samples.npz"):
        path = prepared / ("smooth-lattice" + suffix)
        inputs[str(path)] = workflow.digest(path)
        shutil.copyfile(path, output.parent / path.name)
    workflow.save(output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, output.parent / "analysis-driver.py")
    shutil.copyfile(prepared / "mesh.npz", output.parent / "mesh.npz")
    mesh = np.load(output.parent / "mesh.npz")
    assert len(mesh["points"]) <= 250000 and len(mesh["tetrahedra"]) <= workflow.CAP
    row.update(status="solving", prepared=str(prepared), completed_loads=0, field_solve_performed=True,
        solver="cached-elastic-block-amg-cg", solver_quadrature_points=15, saved_field_or_hierarchy_replayed=False,
        per_load_timeout_s=240, solve_stage_timeout_s=1200)
    workflow.save(output, row)
    def checkpoint(completed, partial, records):
        row.update(completed_loads=completed, partial_stiffness_columns=partial.tolist(), solver_records=records)
        workflow.save(output, row)
        signal.alarm(240)
    def expired(*_):
        raise TimeoutError("A fresh load exceeded the 240-second continuation gate")
    previous_handler = signal.signal(signal.SIGALRM, expired)
    phase = time.perf_counter()
    signal.alarm(240)
    try:
        fields, tensor, records = solve(mesh["points"], mesh["tetrahedra"], checkpoint,
            n_elm_gp=15, solver=AuditedAMG(elastic_blocks=True), reaction_reader=vector_first_reactions,
            release_reduced_cache=True)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_handler)
    np.savez_compressed(output.parent / "fields.npz", displacement=fields)
    solve_s = time.perf_counter() - phase
    row.pop("partial_stiffness_columns", None)
    row.update(status="solved-needs-reference", completed_loads=6, stiffness=tensor.tolist(), solver_records=records,
        solve_s=solve_s, charged_workflow_through_solve_s=row["charged_preparation_workflow_s"]+solve_s,
        field_archive_sha256=workflow.digest(output.parent / "fields.npz"),
        scope="Fresh Mac actual generation, explicit repair, straight elevation and star refinement followed by six fresh loads. Diagnostic intermediate audits are charged; no comparative timing, continuum or reference-promotion claim.")
    workflow.save(output, row)
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


def assess(source, audit, output):
    from hp_coupled_preflight import inherited_inputs
    from cubic_refinement import MODEL_REPORTS
    from periodic import tensor_difference

    started = time.perf_counter()
    row = json.loads((source / "result.json").read_text())
    assert row["qualified"] and row["completed_loads"] == 6 and row["fresh_actual_generation"]
    assert not row["saved_field_or_hierarchy_replayed"]
    assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
    audit_file = audit / (source.name + ".json")
    consistency_file = audit / (source.name + "-consistency.json")
    quadrature, consistency = (json.loads(path.read_text()) for path in (audit_file, consistency_file))
    measured = json.loads((audit / "manifest.json").read_text())["inputs"]
    assert all(workflow.digest(source / name) == measured[str(source / name)] for name in ("result.json", "mesh.npz", "fields.npz"))
    assert max(consistency[key] for key in ("reaction_relative_error", "hill_mandel_relative_error", "integrated_stiffness_skew_relative")) < 1e-7
    assert consistency["periodic_jump_max"] < 1e-9
    higher = quadrature["rules"][-1]
    assert higher["points"] == 125
    assert higher["energy_vs_solved_tensor"]["worst_relative_strain_energy_difference"] < 1e-4
    assert higher["energy_vs_previous_rule"]["worst_relative_strain_energy_difference"] < 1e-4
    paths = [Path(__file__), source / "result.json", source / "mesh.npz", source / "fields.npz", source / "worker.json",
        audit_file, consistency_file, audit / "manifest.json", audit / (source.name + "-worker.json")]
    paths.extend(workflow.HERE / "runs" / name for name in MODEL_REPORTS)
    inputs = {**inherited_inputs(source), **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, output.parent / "analysis-driver.py")
    tensor, models = np.asarray(row["stiffness"]), []
    for name in MODEL_REPORTS:
        path = workflow.HERE / "runs" / name
        report = json.loads(path.read_text())
        assert report["case"] == "smooth-lattice" and report["stiffness_target"] == .01
        for index, model in enumerate(report["models"]):
            models.append(dict(source=str(path), index=index, difference=tensor_difference(tensor, np.asarray(model["stiffness"]))))
    assert len(models) == 24
    worst = max(model["difference"]["worst_relative_strain_energy_difference"] for model in models)
    audit_s = json.loads((audit / (source.name + "-worker.json")).read_text())["wall_s"]
    assess_s = time.perf_counter()-started
    workflow.save(output, dict(fresh_actual_workflow_assessment=True, source=str(source), audit=str(audit),
        geometry_admitted=True, fresh_field_audit_passed=True, retained_model_count=24, models=models,
        stiffness_target=.01, worst_retained_model_disagreement=worst, passes_target_all_models=worst <= .01,
        charged_components_s=dict(generation_and_audits=row["generation_and_audits_s"], repair_and_audits=row["repair_and_audits_s"],
            elevation_refinement_and_admission=row["total_preparation_s"], solve=row["solve_s"], field_audit=audit_s, assessment=assess_s),
        charged_science_audited_workflow_s=row["charged_workflow_through_solve_s"]+audit_s+assess_s,
        fresh_actual_generation=True, saved_seed_replay=False, timing_comparable=False, reference_promoted=False,
        scope="One fresh actual lattice workflow with explicit repair, straight Tet10/direct fedoo and star refinement. Named stage timers charge intermediate diagnostics and serialization within their measured boundaries. Synthetic controls, untimed provenance checks, queue and inspection gaps are excluded. No repeated or matched-pair timing, default-pipeline, continuum or local-stress claim."))
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


def cost_receipt(assessment, localization, output):
    from hp_coupled_preflight import inherited_inputs

    report = json.loads((assessment / "result.json").read_text())
    solved = json.loads((Path(report["source"]) / "result.json").read_text())
    generation = Path(solved["generation_source"])
    diagnostic = localization / (generation.name+".json")
    worker = localization / (generation.name+"-worker.json")
    assert Path(json.loads(diagnostic.read_text())["source"]) == generation
    record = json.loads(worker.read_text())
    assert record["returncode"] == 0 and record["wall_s"] > 0
    paths = [Path(__file__), assessment / "result.json", diagnostic, worker, localization / "manifest.json"]
    inputs = {**inherited_inputs(assessment), **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, output.parent / "analysis-driver.py")
    costs = {**report["charged_components_s"], "quality_localization": record["wall_s"]}
    workflow.save(output, dict(fresh_actual_workflow_costs=True, assessment=str(assessment), localization=str(localization),
        generation_source=str(generation), charged_components_s=costs,
        charged_science_audited_workflow_s=sum(costs.values()),
        worst_retained_model_disagreement=report["worst_retained_model_disagreement"],
        passes_target_all_models=report["passes_target_all_models"], timing_comparable=False,
        scope="Fresh actual workflow named-stage costs plus its quality-localization diagnostic. Synthetic controls, untimed provenance checks, bookkeeping, queue and inspection gaps are excluded. No timing ratio or repeated comparison."))
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


def main():
    parser = argparse.ArgumentParser()
    stage = parser.add_mutually_exclusive_group(required=True)
    stage.add_argument("--source", type=Path)
    stage.add_argument("--prepared", type=Path)
    stage.add_argument("--assess", type=Path)
    stage.add_argument("--costs", type=Path)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--localization", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    if args.source:
        prepare(args.source, args.output)
    elif args.prepared:
        solve_prepared(args.prepared, args.output)
    elif args.assess:
        assert args.audit
        assess(args.assess, args.audit, args.output)
    else:
        assert args.localization
        cost_receipt(args.costs, args.localization, args.output)


if __name__ == "__main__":
    main()
