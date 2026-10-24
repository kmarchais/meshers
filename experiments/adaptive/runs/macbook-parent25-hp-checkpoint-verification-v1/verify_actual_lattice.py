"""Replay fresh actual lattice geometry, field admission and charged stages."""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh

import microgen_workflow as workflow


def verify_prepared(run, result):
    workflow.imports()
    from microgen_tet10_control import straight_elevation
    from nested_refinement import interpolation_check
    from homogenization import strain_tensor
    from hp_shared_domain import periodic_jumps
    from interior_relocation_trial import validate_repair_source

    assert result["fresh_actual_generation"] and not result["saved_seed_replay"]
    assert not result["field_solve_performed"]
    repair_path, generation = Path(result["source"]), Path(result["generation_source"])
    repair = json.loads((repair_path / "result.json").read_text())
    parent = json.loads((generation / "result.json").read_text())
    validate_repair_source(parent)
    assert Path(repair["source"]) == generation and parent["stage_audit_enabled"]
    assert workflow.digest(generation / "driver.py") == json.loads((generation / "manifest.json").read_text())["script_sha256"]
    stages = []
    for name in ("vtk-seed", "mmg-stage-1", "mmg-stage-2"):
        mesh = np.load(generation / f"{name}.npz")
        stage = json.loads((generation / f"{name}.json").read_text())
        assert len(mesh["points"]) == stage["nodes"] <= 250000
        assert len(mesh["tetrahedra"]) == stage["elements"] <= workflow.CAP
        geometry = workflow.geometry(generation, "smooth-lattice", mesh["points"], mesh["tetrahedra"], mesh["surface"])
        for key in ("sampled_two_sided_max", "relative_volume_error", "q_min_corner", "minimum_sampled_jacobian"):
            np.testing.assert_allclose(geometry[key], stage["geometry"][key], atol=1e-13)
        stages.append({"name": name, "nodes": stage["nodes"], "tets": stage["elements"]})
    original = np.load(repair_path / "mesh.npz")
    expected = straight_elevation(original["points"], original["tetrahedra"], original["surface"])
    elevated, final = np.load(run / "elevated-mesh.npz"), dict(np.load(run / "mesh.npz"))
    assert expected["within_node_cap"] and expected["nodes"] == result["lifted_nodes"]
    for key in ("points", "tetrahedra", "surface"):
        np.testing.assert_array_equal(expected[key], elevated[key])
    transfer = dict(np.load(run / "prolongation.npz"))
    fields = np.array([elevated["points"] @ strain_tensor(load).T for load in np.eye(6)])
    reproduced = interpolation_check(elevated["points"], elevated["tetrahedra"], fields, {**final, **transfer})
    periodic_jumps(final["points"], final["surface"], transfer["displacement"])
    geometry = workflow.geometry(run, "smooth-lattice", final["points"], final["tetrahedra"], final["surface"])
    assert geometry["geometry_admitted"] and len(final["points"]) == result["nodes"] <= 250000
    assert len(final["tetrahedra"]) == result["elements"] <= workflow.CAP
    for key in ("sampled_two_sided_max", "relative_volume_error", "q_min_corner", "minimum_sampled_jacobian"):
        np.testing.assert_allclose(geometry[key], result["geometry"][key], atol=1e-13)
    np.testing.assert_allclose(geometry["volume"]/repair["geometry"]["volume"]-1, result["relative_volume_change"], atol=1e-13)
    np.testing.assert_allclose(result["generation_and_audits_s"], parent["total_workflow_s"], atol=1e-12)
    np.testing.assert_allclose(result["repair_and_audits_s"], repair["total_saved_mesh_study_s"], atol=1e-12)
    np.testing.assert_allclose(result["charged_preparation_workflow_s"],
        parent["total_workflow_s"]+repair["total_saved_mesh_study_s"]+result["total_preparation_s"], atol=1e-12)
    assert 0 < result["total_preparation_s"] <= json.loads((run / "worker.json").read_text())["wall_s"]
    return dict(run=run.name, nodes=len(final["points"]), tets=len(final["tetrahedra"]),
        stages=stages, reproduction=reproduced, charged_preparation_s=result["charged_preparation_workflow_s"])


def verify_assessment(run, result):
    workflow.imports()
    from cubic_refinement import MODEL_REPORTS

    source, audit = Path(result["source"]), Path(result["audit"])
    solved = json.loads((source / "result.json").read_text())
    prepared = Path(solved["prepared"])
    assert workflow.digest(source / "mesh.npz") == workflow.digest(prepared / "mesh.npz")
    assert workflow.digest(source / "fields.npz") == solved["field_archive_sha256"]
    assert solved["completed_loads"] == len(solved["solver_records"]) == 6
    assert solved["fresh_actual_generation"] and not solved["saved_seed_replay"] and not solved["saved_field_or_hierarchy_replayed"]
    assert max(row["relative_residual"] for row in solved["solver_records"]) < 2e-9
    consistency = json.loads((audit / (source.name+"-consistency.json")).read_text())
    assert max(consistency[key] for key in ("reaction_relative_error", "hill_mandel_relative_error", "integrated_stiffness_skew_relative")) < 1e-7
    assert consistency["periodic_jump_max"] < 1e-9
    higher = json.loads((audit / (source.name+".json")).read_text())["rules"][-1]
    assert higher["points"] == 125
    assert max(higher[key]["worst_relative_strain_energy_difference"] for key in ("energy_vs_solved_tensor", "energy_vs_previous_rule")) < 1e-4
    tensor, errors = np.asarray(solved["stiffness"]), []
    expected = {(name, i) for name in MODEL_REPORTS for i in range(8)}
    actual = {(Path(model["source"]).name, model["index"]) for model in result["models"]}
    assert actual == expected and len(result["models"]) == result["retained_model_count"] == 24
    for model in result["models"]:
        reference = np.asarray(json.loads(Path(model["source"]).read_text())["models"][model["index"]]["stiffness"])
        values = eigvalsh((tensor+tensor.T-reference-reference.T)/2, (reference+reference.T)/2)
        np.testing.assert_allclose(values, model["difference"]["relative_energy_difference_eigenvalues"], atol=1e-12)
        errors.append(float(abs(values).max()))
    np.testing.assert_allclose(max(errors), result["worst_retained_model_disagreement"], atol=1e-12)
    assert result["passes_target_all_models"] == (max(errors) <= .01)
    costs = result["charged_components_s"]
    for key, field in (("generation_and_audits", "generation_and_audits_s"), ("repair_and_audits", "repair_and_audits_s"),
        ("elevation_refinement_and_admission", "total_preparation_s"), ("solve", "solve_s")):
        np.testing.assert_allclose(costs[key], solved[field], atol=1e-12)
    audit_s = json.loads((audit / (source.name+"-worker.json")).read_text())["wall_s"]
    np.testing.assert_allclose(costs["field_audit"], audit_s, atol=1e-12)
    assert 0 < costs["assessment"] <= json.loads((run / "worker.json").read_text())["wall_s"]
    np.testing.assert_allclose(sum(costs.values()), result["charged_science_audited_workflow_s"], atol=1e-12)
    assert not result["timing_comparable"] and not result["reference_promoted"]
    return dict(run=run.name, model_count=24, worst_disagreement=max(errors),
        passes_target_all_models=max(errors) <= .01, charged_workflow_s=result["charged_science_audited_workflow_s"])


def verify_costs(run, result):
    report = json.loads((Path(result["assessment"]) / "result.json").read_text())
    generation = Path(result["generation_source"])
    localization = Path(result["localization"])
    diagnostic = json.loads((localization / (generation.name+".json")).read_text())
    assert Path(diagnostic["source"]) == generation
    record = json.loads((localization / (generation.name+"-worker.json")).read_text())
    assert record["returncode"] == 0
    expected = {**report["charged_components_s"], "quality_localization": record["wall_s"]}
    assert expected.keys() == result["charged_components_s"].keys()
    for key, value in expected.items():
        assert value > 0
        np.testing.assert_allclose(value, result["charged_components_s"][key], atol=1e-12)
    np.testing.assert_allclose(sum(expected.values()), result["charged_science_audited_workflow_s"], atol=1e-12)
    assert result["worst_retained_model_disagreement"] == report["worst_retained_model_disagreement"]
    assert result["passes_target_all_models"] == report["passes_target_all_models"] and not result["timing_comparable"]
    return dict(run=run.name, charged_workflow_including_quality_localization_s=sum(expected.values()))
