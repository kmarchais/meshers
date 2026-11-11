"""One fixed Gauss128 closure and tanh-sinh audit of saved warning sections.

The 771 and 777 experiments stay immutable. This experiment compares numerical
methods, including failed tolerance estimates; it neither certifies intervals
nor promotes a reference or changes an admission decision.
"""

import argparse
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
from scipy.integrate import tanhsinh
import gyroid_analytic_volume as primary
import gyroid_volume_validation as validation


PREVIOUS = w.HERE / "runs/macbook-gyroid-volume-validation-v1"
TAU = 2*math.pi
TANH_ABS = 2e-13
TANH_REL = 2e-13
MAXLEVEL = 8


def tanh_section(z, half, *, negative=False):
    cosine, sine = math.cos(TAU*z), math.sin(TAU*z)
    inside, root_info = validation.clipping_roots(cosine, sine, half, negative=negative)
    edges = np.array([-math.pi, *inside, math.pi])

    def integrand(angle):
        center = cosine*np.sin(angle)
        radius = np.hypot(np.cos(angle), sine)
        lower, upper = (-np.inf, 0.) if negative else (-half, half)
        denominator = np.maximum(radius, 1e-15)
        low = np.clip((lower-center)/denominator, -1., 1.)
        high = np.clip((upper-center)/denominator, -1., 1.)
        fraction = (np.arcsin(high)-np.arcsin(low))/math.pi
        return np.where(radius > 1e-15, fraction,
                        ((center >= lower) & (center <= upper)).astype(float))

    answer = tanhsinh(integrand, edges[:-1], edges[1:],
        atol=TAU*TANH_ABS/(len(edges)-1), rtol=TANH_REL, maxlevel=MAXLEVEL)
    intervals = []
    for i in range(len(edges)-1):
        error = float(answer.error[i])
        intervals.append(dict(left=float(edges[i]), right=float(edges[i+1]),
            integral=float(answer.integral[i])/TAU,
            reported_error=error/TAU if math.isfinite(error) else None,
            status=int(answer.status[i]), success=bool(answer.success[i]),
            maxlevel=int(answer.maxlevel[i]), evaluations=int(answer.nfev[i])))
    finite_errors = all(row["reported_error"] is not None for row in intervals)
    return dict(volume=math.fsum(row["integral"] for row in intervals),
        reported_error_sum=math.fsum(row["reported_error"] for row in intervals) if finite_errors else None,
        all_intervals_converged=all(row["success"] for row in intervals),
        failed_intervals=sum(not row["success"] for row in intervals),
        intervals=intervals, root_info=root_info)


def controls():
    # Array-valued interval endpoints are the production call pattern. Check
    # orientation and summation independently with an exact polynomial integral.
    edges = np.array([-2., -.7, 0., .1, 1.])
    polynomial = tanhsinh(lambda x: x*x, edges[:-1], edges[1:],
                         atol=1e-13, rtol=1e-13, maxlevel=MAXLEVEL)
    expected = (edges[1:]**3-edges[:-1]**3)/3
    np.testing.assert_allclose(polynomial.integral, expected, atol=2e-13, rtol=0)
    assert np.all(polynomial.success)
    rows = []
    for z in (-.41, .07):
        half = tanh_section(z, .55, negative=True)
        full = tanh_section(z, 1.6)
        empty = tanh_section(z, 0.)
        assert abs(half["volume"]-.5) < 2e-12
        assert abs(full["volume"]-1.) < 2e-12 and abs(empty["volume"]) < 1e-14
        assert all(result["all_intervals_converged"] for result in (half, full, empty))
        rows.append(dict(z=z, half=half["volume"], full=full["volume"], empty=empty["volume"]))
    return dict(vectorized_interval_polynomial=True, half_full_empty_bands=rows)


def transition_diagnostic(name, z, half):
    old = json.loads((primary.RUNS/"macbook-gyroid-analytic-volume-v1"/(name+".json")).read_text())
    nearby = min(old["root_partition"]["roots"], key=lambda row: abs(row["z"]-z))
    cosine, sine = math.cos(TAU*z), math.sin(TAU*z)
    critical_sine = abs(cosine)*math.sqrt((2-cosine*cosine)/(1+cosine*cosine))
    base = math.asin(min(1., critical_sine))
    critical = validation.unique_angles([-math.pi, -math.pi/2, 0., math.pi/2, math.pi,
                                         base, math.pi-base, -base, -math.pi+base])
    stationary = []
    for angle in critical:
        cy, sy = math.cos(angle), math.sin(angle)
        radius = math.hypot(cy, sine)
        if radius <= 1e-15:
            continue
        for branch in (-1., 1.):
            derivative = cy*(cosine-branch*sy/radius)
            if abs(derivative) > 1e-11:
                continue
            second = -cosine*sy+branch*((sy*sy-cy*cy)/radius-sy*sy*cy*cy/radius**3)
            for level in (-half, half):
                gap = cosine*sy+branch*radius-level
                stationary.append(dict(angle=angle, branch=branch, level=level,
                    signed_gap=gap, second_derivative=second,
                    local_angular_scale=math.sqrt(2*abs(gap/second)) if abs(second) > 1e-15 else None))
    closest = min(stationary, key=lambda row: abs(row["signed_gap"]))
    return dict(nearest_saved_transition=nearby, z_distance=abs(nearby["z"]-z),
        nearest_stationary_level_contact=closest,
        interpretation="A small stationary branch-to-level gap produces nearby clipping-root coalescence or a narrow rounded near-contact. Its angular scale is a local quadratic estimate, not a root-separation certificate.")


def closure_scores(name, reference, discrepancy, heuristic):
    report = validation.scores(name, reference, discrepancy, heuristic)
    report["reference_kind"] = "Gauss128 closure estimate; no reference promotion"
    report["gauss128_reference_estimate"] = report.pop("adaptive_reference_estimate")
    for row in [*report["rows"], *report["comparisons"]]:
        for key in list(row):
            if key.startswith("adaptive_"):
                row["gauss128_"+key.removeprefix("adaptive_")] = row.pop(key)
    return report


def worker(args):
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    previous = json.loads((PREVIOUS/"result.json").read_text())
    manifest = json.loads((PREVIOUS/"manifest.json").read_text())
    for field in ("inputs", "verified_primary_inputs"):
        assert all(w.digest(Path(path)) == digest for path, digest in manifest[field].items())
    paths = [Path(__file__), PREVIOUS/"result.json", PREVIOUS/"manifest.json"]
    hashes = {str(path.resolve()): w.digest(path) for path in paths}
    warnings = [(row["case"], example) for row in previous["rows"]
                for example in row["adaptive"]["inner_warning_examples"]]
    total_reported = sum(sum(row["adaptive"]["inner_warning_counts"].values()) for row in previous["rows"])
    assert len(warnings) == total_reported == 15, "Audit all saved warnings, not a favorable subset"
    models = {name: primary.width_model(name) for name in primary.CASES}
    volume_rows, warning_rows = [], []
    parts = {"controls": 1, **{name+"-128": 5 if "density" in name else 1 for name in primary.CASES},
             "warning-sections": 5}
    with Job(parts=parts) as job:
        with job.part("controls") as part:
            part.update(0., "checking tanh-sinh interval controls and frozen inputs", force=True, warnings=15)
            control = controls()
            part.update(1., "tanh-sinh controls passed", warnings=15)
        w.save(args.run/"manifest.json", dict(inputs=hashes, verified_777_inputs=manifest["inputs"],
            verified_771_inputs=manifest["verified_primary_inputs"], controls=control,
            scipy_version=scipy.__version__, gauss_order=128, tanh_sinh_maxlevel=MAXLEVEL,
            tanh_sinh_absolute_tolerance=TANH_ABS, tanh_sinh_relative_tolerance=TANH_REL,
            warning_sections=15, timeout_s=180, memory_bytes=w.MEMORY,
            numerical_threads=w.NUMERICAL_THREAD_POLICY, reference_promoted=False,
            scope="One fixed order128 closure plus all15 saved warning sections. Independent tanh-sinh quadrature uses physical angles and direct clipping panels, without the Gauss sin-squared mapping. No precision or tolerance sweep."))
        for name in primary.CASES:
            width, knots, _ = models[name]
            panels, _ = primary.outer_breaks(name, width, knots)
            with job.part(name+"-128") as part:
                part.update(0., "fixed order128 analytic-line volume closure", force=True, order=128)
                started = time.perf_counter()
                value = primary.integrate(width, panels, 128, "x", lambda fraction, count:
                    part.update(fraction, "fixed order128 analytic-line volume closure", evaluations=count))
            source = next(row for row in previous["rows"] if row["case"] == name)
            scale = max(abs(value["volume"]-source["primary_gauss64_estimate"]),
                        abs(value["volume"]-source["adaptive"]["volume"]))
            row = dict(case=name, gauss128=value, seconds=time.perf_counter()-started,
                gauss64_volume=source["primary_gauss64_estimate"],
                adaptive_volume=source["adaptive"]["volume"],
                difference_128_64=abs(value["volume"]-source["primary_gauss64_estimate"]),
                difference_128_adaptive=abs(value["volume"]-source["adaptive"]["volume"]),
                empirical_observed_scale=scale,
                score_comparison=closure_scores(name, value["volume"], scale,
                    max(scale, source["adaptive"]["combined_reported_error_heuristic"])))
            volume_rows.append(row)
            w.save(args.run/(name+".json"), row)
        nodes, weights = primary.gauss(128)
        with job.part("warning-sections") as part:
            part.update(0., "auditing all saved QUADPACK warning sections", force=True, warnings=15)
            for index, (name, saved) in enumerate(warnings):
                z = saved["z"]
                half = float(models[name][0](z))
                started = time.perf_counter()
                qvalue, qerror, qinfo = validation.section(z, half)
                assert abs(qerror-saved["quadpack_error"]) < 1e-20
                assert qinfo["message"] == saved["message"]
                independent = tanh_section(z, half)
                gauss, evaluations = primary.section(z, half, "x", nodes, weights)
                spread = max(qvalue, independent["volume"], gauss)-min(qvalue, independent["volume"], gauss)
                row = dict(case=name, z=z, half_width=half, saved_warning=saved,
                    quadpack_replay=dict(volume=qvalue, reported_error=qerror, **qinfo),
                    tanh_sinh=independent, gauss128_section=gauss, gauss128_evaluations=evaluations,
                    quadpack_tanh_difference=abs(qvalue-independent["volume"]),
                    gauss_tanh_difference=abs(gauss-independent["volume"]),
                    observed_three_method_spread=spread,
                    reported_error_over_observed_spread=qerror/spread if spread else None,
                    transition=transition_diagnostic(name, z, half), seconds=time.perf_counter()-started)
                warning_rows.append(row)
                w.save(args.run/"warning-sections.json", dict(complete=len(warning_rows)==15, rows=warning_rows))
                part.update((index+1)/15, "compared saved warning section", sections=index+1,
                    observed_spread=spread, tanh_failures=independent["failed_intervals"])
        result = dict(complete=True, volumes=volume_rows, warnings=warning_rows, controls=control,
            largest_warning_section_spread=max(row["observed_three_method_spread"] for row in warning_rows),
            tanh_failed_intervals=sum(row["tanh_sinh"]["failed_intervals"] for row in warning_rows),
            reference_promoted=False, admission_recomputed=False, inputs_unchanged=True,
            scope="Observed agreement tests finite numerical methods. Error estimates and cross-method differences are empirical, not certified bounds. Warning-section agreement does not bound every unsaved section or all outer quadrature errors.")
        for field in ("inputs", "verified_primary_inputs"):
            assert all(w.digest(Path(path)) == digest for path, digest in manifest[field].items())
        assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
        w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    args.run = args.run.resolve()
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-gyroid_volume_closure.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--worker", "--run", str(args.run)],
        args.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=180)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    saved = json.loads((args.run/"result.json").read_text())
    for row in saved["volumes"]:
        result(row["case"]+" fixed volume closure", volume=row["gauss128"]["volume"],
            delta64=row["difference_128_64"], adaptive_delta=row["difference_128_adaptive"],
            seconds=row["seconds"], jsonl=str(args.run/"reported-results.jsonl"))
    result("saved warning sections", sections=len(saved["warnings"]),
        largest_observed_spread=saved["largest_warning_section_spread"],
        failed_tanh_intervals=saved["tanh_failed_intervals"],
        jsonl=str(args.run/"reported-results.jsonl"))


if __name__ == "__main__":
    main()
