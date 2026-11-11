"""Reproducible arithmetic audit of frozen Gyroid volume evidence; no quadrature.

Reference intervals are empirical sensitivity ranges, not certified enclosures.
Mesh volume integration roundoff is a separate, unmeasured numerical limitation.
"""

import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
import shutil

import microgen_workflow as w


RUNS = w.HERE/"runs"
VALIDATION = RUNS/"macbook-gyroid-volume-validation-v1"
CLOSURE = RUNS/"macbook-gyroid-volume-closure-v1"


def comparison(before, after, lower, upper):
    # Exact decimal representations of binary64 inputs prevent the decision
    # threshold from being rounded to another binary64 midpoint before testing.
    with localcontext() as context:
        context.prec = 100
        b, a, lo, hi = map(Decimal.from_float, (before, after, lower, upper))
        midpoint = (b+a)/2
        if b == a:
            verdict = "unchanged"
        elif hi < midpoint:
            verdict = "improves" if a < b else "worsens"
        elif lo > midpoint:
            verdict = "improves" if a > b else "worsens"
        else:
            verdict = "unresolved within empirical interval"
        side = "below both mesh volumes" if hi <= min(a, b) else (
            "above both mesh volumes" if lo >= max(a, b) else None)
        probes = [lo, hi, min(hi, max(lo, a)), min(hi, max(lo, b))]
        gains = [abs(b-value)-abs(a-value) for value in probes]
        clearance = midpoint-hi if midpoint > hi else lo-midpoint if midpoint < lo else Decimal(0)
        return dict(verdict=verdict, midpoint_decimal=str(midpoint),
            midpoint_clearance=float(clearance),
            absolute_volume_error_improvement_range=[float(min(gains)), float(max(gains))],
            reference_interval_same_side=side,
            same_side_reference_independent_gain=float(gains[0]) if side else None,
            exact_midpoint_test="For after<before, improvement iff reference<(before+after)/2; for after>before, reverse the inequality. Comparisons use exact decimal forms of the supplied binary64 endpoints and mesh volumes.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    args.run = args.run.resolve()
    args.run.mkdir(parents=True, exist_ok=False)
    previous = json.loads((VALIDATION/"result.json").read_text())
    closure = json.loads((CLOSURE/"result.json").read_text())
    paths = [Path(__file__), VALIDATION/"result.json", VALIDATION/"manifest.json",
             CLOSURE/"result.json", CLOSURE/"manifest.json"]
    records = []
    for row in closure["volumes"]:
        old = next(item for item in previous["rows"] if item["case"] == row["case"])
        estimates = dict(gauss64=row["gauss64_volume"], gauss128=row["gauss128"]["volume"],
                         adaptive_quadpack=row["adaptive_volume"])
        padding = max(old["adaptive"]["combined_reported_error_heuristic"], row["difference_128_64"])
        lower = min(estimates.values())-padding
        upper = max(estimates.values())+padding
        scores = row["score_comparison"]
        meshes = {item["stage"]: item for item in scores["rows"]}
        for item in meshes.values():
            path = Path(item["source"])
            assert w.digest(path) == item["source_sha256"]
            paths.append(path)
        comparisons = []
        for old_comparison in scores["comparisons"]:
            before, after = old_comparison["before"], old_comparison["after"]
            outcome = comparison(meshes[before]["volume"], meshes[after]["volume"], lower, upper)
            comparisons.append(dict(before=before, after=after,
                before_mesh_volume=meshes[before]["volume"], after_mesh_volume=meshes[after]["volume"],
                frozen_reference_absolute_error_improvement=old_comparison["old_absolute_error_improvement"],
                gauss128_absolute_error_improvement=old_comparison["gauss128_absolute_error_improvement"],
                **outcome))
        records.append(dict(case=row["case"], estimates=estimates,
            empirical_reference_interval=[lower, upper], padding=padding,
            padding_rule="Hull of Gauss64, Gauss128 and QUADPACK volume estimates, expanded by max(previous QUADPACK combined error heuristic, abs(Gauss128-Gauss64)). Not a rigorous coverage statement.",
            frozen_reference_volume=scores["frozen_reference_volume"],
            frozen_scramble_spread=scores["frozen_scramble_spread"],
            saved_mesh_scores=meshes, comparisons=comparisons))
    warning_rows = []
    for row in closure["warnings"]:
        warning_rows.append(dict(case=row["case"], z=row["z"],
            quadpack_volume=row["quadpack_replay"]["volume"],
            tanh_sinh_volume=row["tanh_sinh"]["volume"],
            gauss128_volume=row["gauss128_section"],
            quadpack_reported_error=row["quadpack_replay"]["reported_error"],
            tanh_sinh_reported_error=row["tanh_sinh"]["reported_error_sum"],
            observed_quadpack_minus_tanh=row["quadpack_replay"]["volume"]-row["tanh_sinh"]["volume"],
            observed_gauss128_minus_tanh=row["gauss128_section"]-row["tanh_sinh"]["volume"],
            observed_method_spread=row["observed_three_method_spread"],
            tanh_sinh_failed_intervals=row["tanh_sinh"]["failed_intervals"],
            transition_distance_z=row["transition"]["z_distance"],
            stationary_contact=row["transition"]["nearest_stationary_level_contact"]))
    assert len(records) == 4 and sum(len(row["comparisons"]) for row in records) == 12
    assert len(warning_rows) == 15
    hashes = {str(path.resolve()): w.digest(path) for path in paths}
    result = dict(complete=True, cases=records, warning_sections=warning_rows, inputs=hashes,
        scope="Only arithmetic over frozen measurements. Empirical reference intervals are not certified bounds or confidence intervals. Warning-section differences are observed errors relative to other numerical methods, not exact errors. Old references and admission stay unchanged.",
        mesh_volume_roundoff="Saved mesh volume values are reused without a new quadrature-roundoff audit. Reference convergence and mesh-volume numerical roundoff are separate limitations.",
        quadrature_performed=False, native_reference_promoted=False, admission_changed=False,
        robust_comparisons=sum(item["verdict"] in ("improves", "worsens", "unchanged")
            for row in records for item in row["comparisons"]))
    w.save(args.run/"result.json", result)
    shutil.copyfile(__file__, args.run/"input-gyroid_volume_rankings.py")
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    print(f"Saved {len(records)} cases, 12 midpoint comparisons and 15 warning sections to {args.run/'result.json'}")


if __name__ == "__main__":
    main()
