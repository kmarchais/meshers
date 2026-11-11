"""Independent adaptive quadrature and read-only Gyroid volume-score comparison.

Uses QUADPACK in physical coordinates, scalar line fractions and direct clipping
root isolation. The frozen Gauss experiment supplies widths and evidence only;
none of its quadrature nodes, endpoint mappings or transition panels are used.
No reference promotion, admission change, mesh motion or FEA solve occurs here.
"""

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

import microgen_workflow as w
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
import gyroid_analytic_volume as primary


RUNS = w.HERE / "runs"
PRIMARY = RUNS / "macbook-gyroid-analytic-volume-v1"
CASES = primary.CASES
TAU = 2 * math.pi
INNER_ABS = 2e-12
OUTER_ABS = 1e-10
REL = 1e-12
LIMIT = 150


def clipped_asin(value):
    return math.asin(max(-1., min(1., value)))


def line_fraction(angle, cosine_z, sine_z, half, *, negative=False):
    center = cosine_z * math.sin(angle)
    radius = math.hypot(math.cos(angle), sine_z)
    lower, upper = (-math.inf, 0.) if negative else (-half, half)
    if radius <= 1e-15:
        return float(lower <= center <= upper)
    return (clipped_asin((upper-center)/radius)
            - clipped_asin((lower-center)/radius)) / math.pi


def unique_angles(values):
    result = []
    for value in sorted(values):
        if not result or value-result[-1] > 2e-13:
            result.append(float(value))
    return result


def clipping_roots(cosine_z, sine_z, half, *, negative=False):
    """Bracket C +/- R = level directly, split at every branch critical point.

    d(C +/- R)/dY = cos(Y) * (cos(Z) -/+ sin(Y)/R).
    Squaring the second factor gives the critical sine magnitude below.
    Both signs are retained; the extra stationary candidates are harmless.
    This avoids the primary experiment's quadratic clipping-root formula.
    """
    k = cosine_z
    critical_sine = abs(k) * math.sqrt((2-k*k)/(1+k*k))
    base = math.asin(min(1., critical_sine))
    critical = unique_angles([-math.pi, -math.pi/2, 0., math.pi/2, math.pi,
                              base, math.pi-base, -base, -math.pi+base])
    roots, residuals, stationary_near_zeros = [], [], 0
    levels = (0.,) if negative else (-half, half)
    for sign in (-1., 1.):
        for level in levels:
            def function(angle):
                return (k*math.sin(angle)
                        + sign*math.hypot(math.cos(angle), sine_z) - level)
            values = [function(angle) for angle in critical]
            for angle, value in zip(critical, values, strict=True):
                if abs(value) <= 2e-13:
                    roots.append(angle)
                    residuals.append(abs(value))
                    stationary_near_zeros += 1
            for i, (left, right) in enumerate(zip(critical[:-1], critical[1:], strict=True)):
                if values[i]*values[i+1] < 0:
                    root = brentq(function, left, right, xtol=1e-14,
                                  rtol=4*np.finfo(float).eps)
                    roots.append(root)
                    residuals.append(abs(function(root)))
    points = unique_angles([*critical, *roots])
    return points[1:-1], dict(maximum_root_residual=max(residuals, default=0.),
        stationary_near_zeros=stationary_near_zeros, clipping_roots=len(roots))


def unpack_quad(answer):
    value, error, info = answer[:3]
    return float(value), float(error), dict(evaluations=int(info["neval"]),
        subintervals=int(info["last"]), message=str(answer[3]) if len(answer) > 3 else None)


def section(z, half, *, negative=False):
    cosine_z, sine_z = math.cos(TAU*z), math.sin(TAU*z)
    points, roots = clipping_roots(cosine_z, sine_z, half, negative=negative)
    value, error, info = unpack_quad(quad(
        lambda angle: line_fraction(angle, cosine_z, sine_z, half, negative=negative),
        -math.pi, math.pi, points=points, epsabs=TAU*INNER_ABS,
        epsrel=REL, limit=LIMIT, full_output=1))
    return value/TAU, error/TAU, dict(**info, **roots)


def controls():
    rows = []
    for half in (.1, .55, .9):
        value = line_fraction(0., 1., 0., half)
        expected = 2*math.asin(half)/math.pi
        assert abs(value-expected) < 1e-14
    assert line_fraction(math.pi/2, 1., 0., .55) == 0.
    assert line_fraction(math.pi/2, 1., 0., 1.1) == 1.
    # This control exercises the adaptive integrator and clipping roots, not
    # merely two ways to write the scalar fraction. Negative-half Gyroid area
    # is exactly 1/2 by (X,Y)->(X+pi,-Y), which changes the field's sign.
    for z in (-.41, -.13, .07, .32):
        value, error, info = section(z, .55, negative=True)
        assert abs(value-.5) < 2e-11 and info["message"] is None
        empty, _, empty_info = section(z, 0.)
        full, _, full_info = section(z, 1.6)
        assert abs(empty) < 1e-14 and abs(full-1.) < 2e-12
        assert empty_info["message"] is None and full_info["message"] is None
        rows.append(dict(z=z, negative_half_volume=value, quadpack_error=error))
    return dict(pure_sine=True, zero_amplitude_plane=True,
                zero_and_full_band=True, half_solid_sections=rows)


def adaptive_volume(width, knots, part, checkpoint):
    # Deliberately independent of the primary transition partition. Keep every
    # density interpolation knot exactly, plus a uniform physical z partition.
    panels = np.unique(np.concatenate((np.linspace(-.5, .5, 17), knots)))
    assert all(knot in panels for knot in knots)
    rows, cache = [], {}
    messages = Counter()
    examples = []
    total_evaluations = 0
    inner_max_error = 0.
    root_max_residual = 0.
    near_zero_count = 0
    last_report = time.monotonic()
    started = time.perf_counter()
    for index, (left, right) in enumerate(zip(panels[:-1], panels[1:], strict=True)):
        def integrand(z):
            nonlocal total_evaluations, inner_max_error, root_max_residual
            nonlocal near_zero_count, last_report
            if z not in cache:
                value, error, info = section(z, float(width(z)))
                cache[z] = value
                total_evaluations += info["evaluations"]
                inner_max_error = max(inner_max_error, error)
                root_max_residual = max(root_max_residual, info["maximum_root_residual"])
                near_zero_count += info["stationary_near_zeros"]
                if info["message"] is not None:
                    messages[info["message"]] += 1
                    if len(examples) < 12:
                        examples.append(dict(z=z, quadpack_error=error, **info))
            if time.monotonic()-last_report >= 1.:
                part.update(index/(len(panels)-1), "adaptive physical-coordinate quadrature",
                    sections=len(cache), inner_error=inner_max_error,
                    elapsed_s=time.perf_counter()-started)
                last_report = time.monotonic()
            return cache[z]
        value, error, info = unpack_quad(quad(integrand, float(left), float(right),
            epsabs=OUTER_ABS*float(right-left), epsrel=REL, limit=LIMIT, full_output=1))
        rows.append(dict(left=float(left), right=float(right), volume=value,
                         quadpack_error=error, **info))
        part.update((index+1)/(len(panels)-1), "completed physical z panel",
            sections=len(cache), inner_error=inner_max_error,
            elapsed_s=time.perf_counter()-started)
        checkpoint(dict(complete=False, outer_panels_completed=len(rows),
            outer_panel_count=len(panels)-1, panels=rows, sections=len(cache),
            inner_warning_counts=dict(messages)))
    outer_error = math.fsum(row["quadpack_error"] for row in rows)
    return dict(volume=math.fsum(row["volume"] for row in rows),
        seconds=time.perf_counter()-started, sections=len(cache),
        scalar_fraction_evaluations=total_evaluations, outer_panels=rows,
        frozen_knots_retained=len(knots), inner_max_reported_error=inner_max_error,
        outer_reported_error_sum=outer_error,
        combined_reported_error_heuristic=outer_error+inner_max_error,
        inner_warning_counts=dict(messages), inner_warning_examples=examples,
        outer_warning_count=sum(row["message"] is not None for row in rows),
        maximum_clipping_root_residual=root_max_residual,
        stationary_near_zero_count=near_zero_count,
        root_coverage="Direct branch roots bracketed between all analytic stationary candidates. Near-zero critical values retained. Floating-point root isolation is not interval certified.",
        uncertainty_scope="QUADPACK reports estimates, not bounds. The sum of outer estimates plus maximum observed inner estimate is a heuristic, not a confidence interval.")


def score_sources(name):
    native = RUNS / "macbook-native-face-fit-nine-v1"
    global_run = "macbook-volume-bias-pilot-v1" if name == "gyroid" else "macbook-volume-bias-seven-v1"
    return dict(baseline=native/(name+"-repeat-0-baseline")/"result.json",
                native_fit=native/(name+"-repeat-0-fit")/"result.json",
                global_correction=RUNS/global_run/name/"result.json",
                signed_envelope=RUNS/"macbook-volume-signed-envelope-six-v1"/name/"result.json")


def scores(name, reference, empirical_discrepancy, heuristic):
    sources = score_sources(name)
    records = {stage: json.loads(path.read_text()) for stage, path in sources.items()}
    original = json.loads((PRIMARY/(name+".json")).read_text())
    frozen_volume = original["frozen_reference"]["volume"]
    gauss_volume = original["estimate"]
    rows = []
    for stage, record in records.items():
        geometry = record["geometry"] if stage in ("baseline", "native_fit") else record["after"]["geometry"]
        volume = geometry["volume"]
        assert abs(geometry["relative_volume_error"]-(volume/frozen_volume-1)) < 2e-15
        rows.append(dict(stage=stage, source=str(sources[stage]), source_sha256=w.digest(sources[stage]),
            mesh_sha256=record["mesh_sha256"], volume=volume,
            old_relative_volume_error=geometry["relative_volume_error"],
            old_absolute_relative_volume_error=abs(geometry["relative_volume_error"]),
            old_absolute_volume_error=abs(volume-frozen_volume),
            adaptive_relative_volume_error=volume/reference-1,
            adaptive_absolute_volume_error=abs(volume-reference),
            gauss64_absolute_volume_error=abs(volume-gauss_volume),
            volume_float_ulp=math.ulp(volume)))
    by_stage = {row["stage"]: row for row in rows}
    comparisons = []
    for before, after in (("baseline", "native_fit"), ("native_fit", "global_correction"),
                          ("native_fit", "signed_envelope")):
        b, a = by_stage[before], by_stage[after]
        improvement = b["adaptive_absolute_volume_error"]-a["adaptive_absolute_volume_error"]
        gauss_improvement = b["gauss64_absolute_volume_error"]-a["gauss64_absolute_volume_error"]
        # Moving a reference by d can change a before/after absolute-error
        # difference by at most 2*d. This is sensitivity, not an error bound.
        observed_sensitivity = 2*empirical_discrepancy
        heuristic_sensitivity = 2*heuristic
        comparisons.append(dict(before=before, after=after,
            old_absolute_error_improvement=b["old_absolute_volume_error"]-a["old_absolute_volume_error"],
            adaptive_absolute_error_improvement=improvement,
            gauss64_absolute_error_improvement=gauss_improvement,
            observed_integration_discrepancy=empirical_discrepancy,
            twice_discrepancy_sensitivity=observed_sensitivity,
            absolute_improvement_over_twice_discrepancy=(abs(improvement)/observed_sensitivity
                if observed_sensitivity else None),
            sign_agrees_between_quadratures=improvement*gauss_improvement > 0,
            exceeds_observed_sensitivity_by_10x=abs(improvement) > 10*observed_sensitivity,
            exceeds_reported_heuristic_sensitivity_by_10x=abs(improvement) > 10*heuristic_sensitivity))
    native_volume = records["native_fit"]["geometry"]["volume"]
    repeated_before = [records[stage]["before"]["geometry"]["volume"]
                       for stage in ("global_correction", "signed_envelope")]
    return dict(frozen_reference_volume=frozen_volume,
        frozen_scramble_spread=original["frozen_reference"].get("volume_scramble_spread"),
        adaptive_reference_estimate=reference, primary_gauss64_estimate=gauss_volume,
        rows=rows, comparisons=comparisons,
        facewise=dict(available=False, reason="This mode was run only for spinodoid and Schwarz P/D, not these four Gyroid cases."),
        saved_mesh_volume_repeat_discrepancy=max(abs(value-native_volume) for value in repeated_before),
        mesh_volume_numerics="Saved Tet10 quadrature volumes are reused. Repetition discrepancy and float ulp are recorded separately; neither certifies quadrature roundoff. No new mesh-volume numerical audit is performed.",
        ranking_scope="Absolute volume-error improvements are compared with empirical integration discrepancy and reported heuristics. Neither is a rigorous error bound. Frozen scores and admission remain unchanged.")


def worker(args):
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    manifest = json.loads((PRIMARY/"manifest.json").read_text())
    # Verify the exact fixture definitions and source used by job771, including
    # its saved sample files. We do not rerun calibration or regenerate samples.
    for path, digest in manifest["inputs"].items():
        assert w.digest(Path(path)) == digest, f"Frozen input changed: {path}"
    paths = [Path(__file__), PRIMARY/"manifest.json", PRIMARY/"result.json"]
    for name in CASES:
        paths.extend([PRIMARY/(name+".json"), *score_sources(name).values()])
    hashes = {str(path.resolve()): w.digest(path) for path in paths}
    models = {name: primary.width_model(name) for name in CASES}
    rows = []
    with Job(parts={"controls": 1, **{name: 8 if "density" in name else 4 for name in CASES}}) as job:
        with job.part("controls") as part:
            part.update(0., "checking frozen inputs and independent analytic controls", force=True, cases=4)
            control = controls()
            part.update(1., "independent analytic controls passed", cases=4)
        w.save(args.run/"manifest.json", dict(inputs=hashes, verified_primary_inputs=manifest["inputs"],
            controls=control, scipy_version=scipy.__version__, inner_absolute_tolerance=INNER_ABS,
            outer_total_absolute_tolerance=OUTER_ABS, relative_tolerance=REL,
            quadpack_subinterval_limit=LIMIT, uniform_outer_panels=16,
            density_knots_preserved=True, timeout_s=180, memory_bytes=w.MEMORY,
            numerical_threads=w.NUMERICAL_THREAD_POLICY, reference_promotion=False,
            independent_method="QUADPACK adaptive Gauss-Kronrod in physical coordinates; direct clipping roots; uniform outer panels plus frozen density knots. No primary Gauss nodes, sin-squared map, or transition panels.",
            primary_xy_scope="Both coordinate reductions integrate the same fixed-z area by Fubini. The torus map (X,Y)->(pi/2-Y,pi/2-X) swaps the sin(Z),cos(Z) coefficients. Matched outer nodes and isometries can yield identical floating-point sums while retaining shared outer quadrature error. Their equality is not an independent outer convergence estimate."))
        for name in CASES:
            width, knots, _ = models[name]
            source = json.loads((PRIMARY/(name+".json")).read_text())
            with job.part(name) as part:
                part.update(0., "starting independent adaptive integration", force=True, sections=0)
                value = adaptive_volume(width, knots, part,
                    lambda partial: w.save(args.run/(name+"-partial.json"), partial))
            discrepancy = abs(value["volume"]-source["estimate"])
            empirical = max(discrepancy, source["empirical_stability"]["last_refinement_difference"])
            row = dict(case=name, adaptive=value, primary_gauss64_estimate=source["estimate"],
                adaptive_gauss_discrepancy=discrepancy, empirical_observed_scale=empirical,
                score_comparison=scores(name, value["volume"], empirical,
                                       max(empirical, value["combined_reported_error_heuristic"])),
                reference_promoted=False, gauss128_performed=False)
            w.save(args.run/(name+".json"), row)
            rows.append(row)
            w.save(args.run/"result.json", dict(complete=len(rows)==len(CASES), rows=rows,
                controls=control, inputs=hashes, reference_promoted=False,
                admission_recomputed=False, meshes_modified=False))
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    assert all(w.digest(Path(path)) == digest for path, digest in manifest["inputs"].items())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    args.run = args.run.resolve()
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-gyroid_volume_validation.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--worker", "--run", str(args.run)],
        args.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=180)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    saved = json.loads((args.run/"result.json").read_text())
    for row in saved["rows"]:
        value = row["adaptive"]
        result(row["case"]+" adaptive volume validation", volume=value["volume"],
            gauss_difference=row["adaptive_gauss_discrepancy"],
            reported_error=value["combined_reported_error_heuristic"],
            warnings=sum(value["inner_warning_counts"].values())+value["outer_warning_count"],
            seconds=value["seconds"], jsonl=str(args.run/"reported-results.jsonl"))
        for comparison in row["score_comparison"]["comparisons"]:
            result(row["case"]+" "+comparison["before"]+" to "+comparison["after"],
                old_abs_error_improvement=comparison["old_absolute_error_improvement"],
                new_abs_error_improvement=comparison["adaptive_absolute_error_improvement"],
                empirical_10x_margin=comparison["exceeds_observed_sensitivity_by_10x"],
                jsonl=str(args.run/"reported-results.jsonl"))


if __name__ == "__main__":
    main()
