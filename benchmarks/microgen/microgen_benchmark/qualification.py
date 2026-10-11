"""Application targets and observed cost to reach them, without interpolation."""

from math import log
from statistics import median


def admission(result, config, reference_uncertainty=0.0):
    """Keep geometry, element shape and solution accuracy as separate decisions."""
    reasons = []
    distance = result.get("sampled_geometry_distance")
    reverse = result.get("sampled_reverse_distance")
    if distance is None:
        reasons.append("Boundary distance unavailable: near-zero field gradient")
    elif distance > config["geometry_tolerance"]:
        reasons.append(f"Sampled boundary distance {distance:.6g}L exceeds tolerance")
    if reverse is not None and reverse > config["geometry_tolerance"]:
        reasons.append(f"Sampled target-to-mesh distance {reverse:.6g}L exceeds tolerance")
    upper = result["relative_volume_error"] + reference_uncertainty
    if upper > config["volume_tolerance"]:
        reasons.append(f"Volume error plus sampling uncertainty {upper:.2%} exceeds tolerance")
    geometry_ok = not reasons
    quality_ok = result.get("q_min", 0) >= config["minimum_quality"]
    energy = result.get("worst_energy_error")
    fea_ok = geometry_ok and quality_ok and energy is not None and energy <= config["fea_target"]
    return {
        "geometry_eligible": geometry_ok,
        "quality_eligible": quality_ok,
        "fea_eligible": fea_ok,
        "comparison_eligible": geometry_ok,
        "admission_reasons": reasons,
        "relative_volume_error_with_uncertainty": upper,
    }


def summarize(document):
    """Require every requested repeat before a configuration can reach a target."""
    config = document["provenance"]["config"]
    configurations = []
    for case in config["case_ids"]:
        for workflow in config["workflows"]:
            for cells in config["cells"]:
                rows = [
                    r
                    for r in document["rows"]
                    if r["case"] == case
                    and r["workflow"] == workflow
                    and r["cells"] == cells
                    and r["repeat"] >= 0
                ]
                complete = len(rows) == config["repeats"] and all(r["status"] == "ok" for r in rows)
                geometry = complete and all(
                    r.get("geometry_eligible", r.get("comparison_eligible", False)) for r in rows
                )
                energy = [
                    r["worst_energy_error"] for r in rows if r.get("worst_energy_error") is not None
                ]
                times = [r["mesh_s"] for r in rows if r.get("mesh_s") is not None]
                configurations.append(
                    {
                        "case": case,
                        "workflow": workflow,
                        "cells": cells,
                        "complete": complete,
                        "geometry_eligible": geometry,
                        "quality_eligible": complete
                        and all(r.get("quality_eligible", False) for r in rows),
                        "median_mesh_s": median(times) if times else None,
                        "max_tetrahedra": max((r.get("tetrahedra", 0) for r in rows), default=0),
                        "max_triangles": max((r.get("triangles", 0) for r in rows), default=0),
                        "worst_energy_error": max(energy)
                        if len(energy) == config["repeats"]
                        else None,
                    }
                )
    targets = []
    for case in config["case_ids"]:
        for workflow in config["workflows"]:
            candidates = [
                g
                for g in configurations
                if g["case"] == case and g["workflow"] == workflow and g["geometry_eligible"]
            ]
            for target in [None, *config["error_targets"]]:
                good = (
                    candidates
                    if target is None
                    else [
                        g
                        for g in candidates
                        if g["quality_eligible"]
                        and g["worst_energy_error"] is not None
                        and g["worst_energy_error"] <= target
                    ]
                )
                fastest = min(good, key=lambda g: g["median_mesh_s"], default=None)
                smallest = min(
                    good, key=lambda g: g["max_triangles"] or g["max_tetrahedra"], default=None
                )
                targets.append(
                    {
                        "case": case,
                        "workflow": workflow,
                        "energy_target": target,
                        "fastest": fastest,
                        "smallest": smallest,
                        "status": "reached" if good else "not reached",
                    }
                )
    convergence = []
    for case in config["case_ids"]:
        for workflow in config["workflows"]:
            sequence = sorted(
                [
                    g
                    for g in configurations
                    if g["case"] == case
                    and g["workflow"] == workflow
                    and g["geometry_eligible"]
                    and g["worst_energy_error"] is not None
                ],
                key=lambda g: g["cells"],
            )
            pairs = []
            for low, high in zip(sequence, sequence[1:], strict=False):
                ratio = high["max_tetrahedra"] / low["max_tetrahedra"]
                order = (
                    3 * log(low["worst_energy_error"] / high["worst_energy_error"]) / log(ratio)
                    if ratio > 1 and high["worst_energy_error"] > 0
                    else None
                )
                pairs.append(
                    {
                        "cells": [low["cells"], high["cells"]],
                        "effective_h_order": order,
                        "error_decreases": high["worst_energy_error"] < low["worst_energy_error"],
                    }
                )
            convergence.append(
                {
                    "case": case,
                    "workflow": workflow,
                    "levels": len(sequence),
                    "pairs": pairs,
                    "status": "decreasing"
                    if len(sequence) >= 3 and all(p["error_decreases"] for p in pairs)
                    else "insufficient levels"
                    if len(sequence) < 3
                    else "nonmonotonic",
                }
            )
    return {
        "configurations": configurations,
        "targets": targets,
        "convergence": convergence,
        "scope": "Observed configurations only. All repeats must pass; no interpolation or extrapolation. Geometry targets do not require FEA or element shape.",
    }
