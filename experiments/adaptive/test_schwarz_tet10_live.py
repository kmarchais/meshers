"""Keep the live workflow scope and failed accuracy visible in summaries."""

import pytest

from schwarz_tet10_live import comparison_summary


def inputs():
    common = {"qualified": True, "geometry": {"geometry_admitted": True},
        "case": "schwarz-solid", "elements": 100, "periodic": {}}
    parent = dict(common, first_pass_only=True, mmg_calls=[{"nofem": True}],
        resolution=30, first_pass_hmax_override=.1, generation_s=3., total_workflow_s=10.)
    solved = dict(common, order=2, nodes=200, solver="cached-elastic-block-amg-cg",
        solver_records=[{"relative_residual": 1e-12}] * 6, solve_s=8.)
    elevated = dict(solved, quadrature_sensitivity_passed=True, quadrature_sensitivity_limit=1e-4,
        solver_quadrature_points=4, lift_and_conversion_s=1., total_saved_mesh_study_s=20.)
    native = dict(solved, solver_quadrature_points=15, generation_s=2., total_workflow_s=15.)
    comparisons = {arm: {"stiffness_target": .001, "passes_target_under_all_declared_models": True}
        for arm in ("microgen", "meshers")}
    return parent, elevated, native, comparisons


def test_live_summary_charges_parent_and_distinguishes_phase_cost():
    report = comparison_summary(*inputs())
    assert report["science_audited_workflow_s"] == {"microgen": 30., "meshers": 15.}
    assert report["generation_lift_solve_s"] == {"microgen": 12., "meshers": 10.}
    assert report["single_pair_workflow_ratio_microgen_over_meshers"] == 2.
    assert report["repetitions"] == 1


def test_failed_accuracy_never_produces_timing_victory():
    parent, elevated, native, comparisons = inputs()
    comparisons["microgen"]["passes_target_under_all_declared_models"] = False
    report = comparison_summary(parent, elevated, native, comparisons)
    assert not report["both_pass_declared_reference_models"]
    assert report["single_pair_workflow_ratio_microgen_over_meshers"] is None
    assert report["single_pair_phase_ratio_microgen_over_meshers"] is None


def test_saved_or_weakened_control_is_rejected():
    parent, elevated, native, comparisons = inputs()
    elevated["quadrature_sensitivity_limit"] = .01
    with pytest.raises(AssertionError):
        comparison_summary(parent, elevated, native, comparisons)


def test_explicit_repair_is_charged_and_must_pass_geometry():
    parent, elevated, native, comparisons = inputs()
    repair = dict(parent, repair={"boundary_coordinates_unchanged": True,
        "remaining_bad_elements": 0, "signed_volume_difference": 0.},
        total_saved_mesh_study_s=4., repair_s=1.)
    parent["qualified"] = False
    report = comparison_summary(parent, elevated, native, comparisons, repair)
    assert report["science_audited_workflow_s"]["microgen"] == 34.
    assert report["generation_lift_solve_s"]["microgen"] == 13.
    assert report["explicit_interior_repair"]
    repair["geometry"] = {"geometry_admitted": False}
    with pytest.raises(AssertionError):
        comparison_summary(parent, elevated, native, comparisons, repair)
