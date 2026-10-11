from copy import deepcopy

from microgen_benchmark.qualification import admission, summarize


def config():
    return {
        "geometry_tolerance": 0.002,
        "volume_tolerance": 0.01,
        "minimum_quality": 0.1,
        "fea_target": 0.05,
        "error_targets": [0.01, 0.02, 0.05],
        "case_ids": ["gyroid"],
        "workflows": ["meshers-balanced"],
        "cells": [16, 32],
        "repeats": 2,
    }


def test_geometry_can_pass_while_fea_fails():
    row = {
        "sampled_geometry_distance": 0.001,
        "relative_volume_error": 0.003,
        "sampled_reverse_distance": 0.0015,
        "q_min": 0.2,
        "worst_energy_error": 0.2,
    }
    result = admission(row, config())
    assert result["geometry_eligible"] and result["quality_eligible"]
    assert not result["fea_eligible"]
    row["worst_energy_error"] = 0.02
    assert admission(row, config())["fea_eligible"]
    row["q_min"] = 0.01
    assert not admission(row, config())["fea_eligible"]
    assert admission(row, config())["geometry_eligible"]


def test_reverse_distance_and_reference_uncertainty_block_admission():
    row = {
        "sampled_geometry_distance": 0,
        "relative_volume_error": 0.009,
        "sampled_reverse_distance": 0,
    }
    assert not admission(row, config(), 0.002)["geometry_eligible"]
    row["relative_volume_error"], row["sampled_reverse_distance"] = 0, 0.003
    assert not admission(row, config())["geometry_eligible"]
    row.pop("sampled_reverse_distance")
    assert not admission(row, config())["geometry_eligible"]


def test_target_requires_all_repeats_and_reports_fastest_and_smallest():
    rows = [
        {
            "case": "gyroid",
            "workflow": "meshers-balanced",
            "cells": n,
            "repeat": r,
            "status": "ok",
            "geometry_eligible": True,
            "quality_eligible": True,
            "mesh_s": time,
            "tetrahedra": count,
            "worst_energy_error": 0.02,
        }
        for n, count, time in [(16, 1000, 2), (32, 2000, 1)]
        for r in range(2)
    ]
    doc = {"provenance": {"config": config()}, "rows": rows}
    target = next(t for t in summarize(doc)["targets"] if t["energy_target"] == 0.02)
    assert target["fastest"]["cells"] == 32
    assert target["smallest"]["cells"] == 16
    duplicated = deepcopy(doc)
    duplicated["rows"][1]["repeat"] = 0
    target = next(t for t in summarize(duplicated)["targets"] if t["energy_target"] == 0.02)
    assert target["smallest"]["cells"] == 32
    changed = deepcopy(doc)
    changed["rows"][0]["status"] = "timeout"
    target = next(t for t in summarize(changed)["targets"] if t["energy_target"] == 0.02)
    assert target["smallest"]["cells"] == 32
    changed["rows"] = changed["rows"][:1]
    assert all(t["status"] == "not reached" for t in summarize(changed)["targets"])
