"""Deterministic volume estimates after exact one-coordinate Gyroid reduction.

The x and y conditional line fractions are analytic. Remaining quadrature is
numerical: convergence and cross-route differences are empirical evidence, not
rigorous error bounds. Frozen specimen references and density widths are inputs,
never replaced or recalibrated by this experiment.
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
from scipy.optimize import brentq
import comparison_fixtures as fixtures
import linear_density_fixture


CASES = ("gyroid", "graded-gyroid-linear", "graded-gyroid-periodic",
         "graded-gyroid-density-linear")
RUNS = w.HERE/"runs"
DENSITY = RUNS/"macbook-linear-density-fixtures-v1"
ORDERS = (16, 32, 64)
TAU = 2*np.pi


def width_model(name):
    if name == "gyroid":
        return lambda z: np.asarray(z)*0+.55, np.array([-.5, .5]), "constant"
    if name == "graded-gyroid-linear":
        return lambda z: .55*(1+.5*np.asarray(z)), np.array([-.5, .5]), "linear"
    if name == "graded-gyroid-periodic":
        return lambda z: .55*(1+.25*np.cos(TAU*np.asarray(z))), np.array([-.5, .5]), "periodic"
    table = json.loads((DENSITY/"linear-density-calibration.json").read_text())
    z, half = np.asarray(table["z"]), np.asarray(table["half_width"])
    assert len(z) == 257 and np.all(np.diff(z) > 0)
    assert z[0] == -.5 and z[-1] == .5 and np.all((half > 0) & (half < 1))
    evaluate = lambda position: np.interp(position, z, half)
    np.testing.assert_array_equal(evaluate(z), half)
    return evaluate, z, "frozen piecewise linear"


def fraction_x(y, z, half, *, negative=False):
    """Exact x-line measure from sin(X)cos(Y)+cos(X)sin(Z)+sin(Y)cos(Z)."""
    radius = np.sqrt(np.cos(TAU*y)**2+np.sin(TAU*z)**2)
    center = np.sin(TAU*y)*np.cos(TAU*z)
    lower, upper = (-np.inf, 0.) if negative else (-half, half)
    denominator = np.maximum(radius, 1e-15)
    lo = np.clip((lower-center)/denominator, -1., 1.)
    hi = np.clip((upper-center)/denominator, -1., 1.)
    value = (np.arcsin(hi)-np.arcsin(lo))/np.pi
    return np.where(radius > 1e-15, value, ((center >= lower) & (center <= upper)).astype(float))


def fraction_y(x, z, half, *, negative=False):
    """Independent y-line measure, using cos-phase CDF rather than asin CDF."""
    radius = np.sqrt(np.cos(TAU*z)**2+np.sin(TAU*x)**2)
    center = np.sin(TAU*z)*np.cos(TAU*x)
    lower, upper = (-np.inf, 0.) if negative else (-half, half)
    denominator = np.maximum(radius, 1e-15)
    lo = np.clip((lower-center)/denominator, -1., 1.)
    hi = np.clip((upper-center)/denominator, -1., 1.)
    value = (np.arccos(lo)-np.arccos(hi))/np.pi
    return np.where(radius > 1e-15, value, ((center >= lower) & (center <= upper)).astype(float))


def unique_sorted(values, tolerance=2e-13):
    result = []
    for value in sorted(values):
        if not result or value-result[-1] > tolerance:
            result.append(float(value))
    return np.array(result)


def inner_breaks(z, half, axis, *, negative=False):
    # At a clipping transition, (L-k*t)^2 = 2-k^2-t^2.
    # t=sin(Y), k=cos(Z) for x reduction; t=cos(X), k=sin(Z) for y.
    k = np.cos(TAU*z) if axis == "x" else np.sin(TAU*z)
    levels = (0.,) if negative else (-half, half)
    angles = [-np.pi, -np.pi/2, 0., np.pi/2, np.pi]
    for level in levels:
        discriminant = (1+k*k)*(2-k*k)-level*level
        if discriminant < 0:
            continue
        for sign in (-1., 1.):
            t = (level*k+sign*np.sqrt(max(discriminant, 0.)))/(1+k*k)
            if not -1-1e-13 <= t <= 1+1e-13:
                continue
            t = np.clip(t, -1., 1.)
            if axis == "x":
                angle = float(np.arcsin(t))
                angles.extend((angle, np.pi-angle if angle >= 0 else -np.pi-angle))
            else:
                angle = float(np.arccos(t))
                angles.extend((-angle, angle))
    return unique_sorted(angles)


def stationary_points(cosine, sine, slope, a, b):
    """All derivative zeros of c*cos(2πz)+s*sin(2πz)-slope*z on [a,b]."""
    radius = TAU*math.hypot(cosine, sine)
    if radius == 0 or abs(slope) > radius:
        return []
    phase = math.atan2(sine, -cosine)
    angle = math.asin(max(-1., min(1., slope/radius)))
    result = []
    for base in (angle-phase, np.pi-angle-phase):
        for turn in range(-3, 4):
            point = (base+TAU*turn)/TAU
            if a < point < b:
                result.append(point)
    return result


def outer_breaks(name, width, knots):
    points = list(knots)+[-.5, -.25, 0., .25, .5]
    roots, near_tangencies = [], []
    segments = []
    # Clipping changes inner-panel topology when cos(Z) ± sin(Z) = ±h(z).
    # On every density segment the width is exactly affine. Split by every
    # analytic derivative zero before bracketing, so repeated roots are checked.
    for a, b in zip(knots[:-1], knots[1:], strict=True):
        slope = float((width(b)-width(a))/(b-a))
        for cosine in (-1., 1.):
            for sine in (-1., 1.):
                c, s, m = cosine, sine, slope
                if name == "graded-gyroid-periodic":
                    c, m = cosine-.55*.25, 0.
                critical = unique_sorted([a, b, *stationary_points(c, s, m, a, b)])
                function = lambda z: cosine*np.cos(TAU*z)+sine*np.sin(TAU*z)-float(width(z))
                values = np.array([function(z) for z in critical])
                for i, (left, right) in enumerate(zip(critical[:-1], critical[1:], strict=True)):
                    segments.append((float(left), float(right)))
                    if values[i]*values[i+1] < 0:
                        root = brentq(function, left, right, xtol=1e-14, rtol=4*np.finfo(float).eps)
                        points.append(root)
                        roots.append(dict(z=root, cosine=cosine, sine=sine,
                                          residual=abs(function(root)), kind="sign crossing"))
                for z, value in zip(critical, values, strict=True):
                    if abs(value) <= 1e-12:
                        points.append(z)
                        row = dict(z=float(z), cosine=cosine, sine=sine,
                                   residual=float(abs(value)), kind="near-zero endpoint or stationary point")
                        roots.append(row)
                        if a < z < b:
                            near_tangencies.append(row)
    # Keep every original lookup knot bitwise, including a possible root lying
    # only a few ulps away. Tiny adjacent panels are harmless; merging them
    # could otherwise shift a frozen interpolation breakpoint.
    partitions = np.unique(np.asarray(points, dtype=float))
    assert partitions[0] == -.5 and partitions[-1] == .5
    assert all(knot in partitions for knot in knots)
    return partitions, dict(frozen_knot_count=len(knots), outer_panels=len(partitions)-1,
                            monotone_brackets=len(segments), roots=roots,
                            near_tangent_candidates=near_tangencies,
                            coverage="Analytic derivative zeros split each affine-width or harmonic-width interval. Near-zero stationary points are retained, including possible repeated roots. Floating-point root isolation is not interval-certified.")


def gauss(order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return (nodes+1)/2, weights/2


def mapped_nodes(edges, nodes, weights):
    lengths = np.diff(edges)
    theta = np.pi*nodes/2
    positions = edges[:-1, None]+lengths[:, None]*np.sin(theta)[None, :]**2
    scaled_weights = lengths[:, None]*(np.pi/2*np.sin(np.pi*nodes)*weights)[None, :]
    return positions, scaled_weights


def section(z, half, axis, nodes, weights, *, negative=False):
    edges = inner_breaks(z, half, axis, negative=negative)
    angle, quadrature_weights = mapped_nodes(edges, nodes, weights)
    fraction = fraction_x(angle/TAU, z, half, negative=negative) if axis == "x" else fraction_y(angle/TAU, z, half, negative=negative)
    return float(np.sum(fraction*quadrature_weights)/TAU), fraction.size


def integrate(width, panels, order, axis, update):
    nodes, weights = gauss(order)
    reverse = axis == "y"
    edges = -panels[::-1] if reverse else panels
    positions, integration_weights = mapped_nodes(edges, nodes, weights)
    np.testing.assert_allclose(integration_weights.sum(), 1., atol=2e-14, rtol=0)
    values, evaluations = [], 0
    positions, integration_weights = positions.ravel(), integration_weights.ravel()
    for i, (coordinate, weight) in enumerate(zip(positions, integration_weights, strict=True)):
        z = -coordinate if reverse else coordinate
        value, count = section(z, float(width(z)), axis, nodes, weights)
        values.append(weight*value)
        evaluations += count
        update((i+1)/len(positions), evaluations)
    return dict(volume=math.fsum(values), order=order, line_axis=axis,
                reversed_outer_coordinate=reverse, line_fraction_evaluations=evaluations,
                outer_nodes=len(positions), outer_panels=len(panels)-1)


def midpoint(width, n, update):
    grid = (np.arange(n)+.5)/n-.5
    z = grid[None, :]
    half = width(z)
    total = 0.
    for start in range(0, n, 32):
        y = grid[start:start+32, None]
        total += float(fraction_x(y, z, half).sum())
        update(min(start+32, n)/n, min(start+32, n)*n)
    return dict(volume=total/(n*n), points_per_axis=n, rule="uniform midpoint, unsplit")


def controls():
    nodes, weights = gauss(32)
    # At z=0 and y=0 the x-line is exactly sin(X). At y=1/4,z=0
    # it is the constant plane value1, including the zero-amplitude branch.
    for half in (.1, .55, .9):
        expected = 2*np.arcsin(half)/np.pi
        np.testing.assert_allclose(fraction_x(0., 0., half), expected, atol=1e-14, rtol=0)
        np.testing.assert_allclose(fraction_y(0., 0., half), expected, atol=1e-14, rtol=0)
    assert fraction_x(.25, 0., .55) == 0 and fraction_x(.25, 0., 1.1) == 1
    for z in (-.41, -.13, .07, .32):
        for axis in ("x", "y"):
            negative, _ = section(z, .55, axis, nodes, weights, negative=True)
            np.testing.assert_allclose(negative, .5, atol=2e-12, rtol=0)
            empty, _ = section(z, 0., axis, nodes, weights)
            full, _ = section(z, 1.6, axis, nodes, weights)
            np.testing.assert_allclose(empty, 0., atol=1e-14, rtol=0)
            np.testing.assert_allclose(full, 1., atol=2e-14, rtol=0)
    return dict(exact_sine_band=True, zero_amplitude_plane=True,
                half_solid_symmetry=True, zero_width_and_full_band=True,
                full_band_bound="Cauchy-Schwarz gives |Gyroid|<=sqrt(S*(3-S))<=1.5.")


def identity(name, width):
    source = RUNS/"macbook-native-face-fit-nine-v1"/(name+"-repeat-0-fit")
    root = DENSITY if name == "graded-gyroid-density-linear" else RUNS/"macbook-comparison-fixtures-v1"
    target = fixtures.case(name, root)
    assert target.levels == (-1., 1.)
    with np.load(source/(name+"-samples.npz")) as samples:
        probes = np.concatenate([samples[key][::max(1, len(samples[key])//128)] for key in ("walls", "caps")])
    x, y, z = probes.T
    raw_x = np.cos(TAU*y)*np.sin(TAU*x)+np.sin(TAU*z)*np.cos(TAU*x)+np.sin(TAU*y)*np.cos(TAU*z)
    raw_y = np.cos(TAU*z)*np.sin(TAU*y)+np.sin(TAU*x)*np.cos(TAU*y)+np.sin(TAU*z)*np.cos(TAU*x)
    measured = target.function(x, y, z)
    errors = [float(abs(raw/width(z)-measured).max()) for raw in (raw_x, raw_y)]
    assert max(errors) < 1e-12
    if name == "graded-gyroid-density-linear":
        calibrated = linear_density_fixture.width(DENSITY)
        np.testing.assert_array_equal(width(z), calibrated(x, y, z))
    elif name.startswith("graded-"):
        np.testing.assert_array_equal(width(z), fixtures.half_width(name, x, y, z))
    return dict(saved_probes=len(probes), x_reduction_field_identity_error=errors[0],
                y_reduction_field_identity_error=errors[1], specimen_unchanged=True)


def worker(a):
    w.imports()
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    calibration = DENSITY/"linear-density-calibration.json"
    frozen_density = json.loads((DENSITY/"result.json").read_text())
    assert w.digest(calibration) == frozen_density["hashes"][calibration.name]
    files = [Path(__file__), Path(w.__file__), Path(fixtures.__file__),
             Path(linear_density_fixture.__file__), w.HERE/"scorecard.py", calibration, DENSITY/"result.json"]
    for name in CASES:
        source = RUNS/"macbook-native-face-fit-nine-v1"/(name+"-repeat-0-fit")
        files += [source/(name+"-samples.npz"), source/(name+"-specimen.json")]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for path in files:
        if path.suffix in (".py", ".json"):
            shutil.copyfile(path, a.run/("input-"+path.name))
    parts = {"controls": 1}
    for name in CASES:
        for order in ORDERS:
            parts[f"{name}-x-{order}"] = (order/16)**2
        parts[f"{name}-y-64"] = 16
        parts[f"{name}-midpoint"] = 2
    estimates = []
    with Job(parts=parts) as job:
        with job.part("controls") as part:
            part.update(0., "checking analytic reductions and frozen inputs", force=True, cases=4)
            control = controls()
            models = {name: width_model(name) for name in CASES}
            identities = {name: identity(name, models[name][0]) for name in CASES}
            splits = {name: outer_breaks(name, models[name][0], models[name][1]) for name in CASES}
            part.update(1., "analytic controls and root partitions ready", cases=4)
            job.result("analytic controls", passed=True)
        w.save(a.run/"manifest.json", dict(inputs=hashes, controls=control, identities=identities,
               root_partitions={name: splits[name][1] for name in CASES},
               quadrature_orders=list(ORDERS), independent_y_order=64, midpoint_order=1024,
               timeout_s=180, memory_bytes=w.MEMORY, numerical_threads=w.NUMERICAL_THREAD_POLICY,
               reference_promotion=False, density_target_not_used_as_truth=True))
        for name in CASES:
            width, knots, kind = models[name]
            panels, root_info = splits[name]
            reference_path = RUNS/"macbook-native-face-fit-nine-v1"/(name+"-repeat-0-fit")/(name+"-specimen.json")
            saved_reference = json.loads(reference_path.read_text())
            row = dict(case=name, width_kind=kind, frozen_knot_count=len(knots),
                       frozen_reference=saved_reference, reductions=[], root_partition=root_info)
            for order in ORDERS:
                with job.part(f"{name}-x-{order}") as part:
                    part.update(0., "integrating exact x-line fractions", force=True, order=order, evaluations=0)
                    started = time.perf_counter()
                    value = integrate(width, panels, order, "x", lambda fraction, count:
                                      part.update(fraction, "integrating exact x-line fractions", order=order, evaluations=count))
                    value["seconds"] = time.perf_counter()-started
                    row["reductions"].append(value)
                    job.result(f"{name} x{order}", volume=value["volume"], seconds=value["seconds"])
                    w.save(a.run/(name+".json"), row)
            with job.part(f"{name}-y-64") as part:
                part.update(0., "independent y reduction with reversed z", force=True, order=64, evaluations=0)
                started = time.perf_counter()
                independent = integrate(width, panels, 64, "y", lambda fraction, count:
                                        part.update(fraction, "independent y reduction with reversed z", order=64, evaluations=count))
                independent["seconds"] = time.perf_counter()-started
                row["independent_y"] = independent
                job.result(f"{name} y64", volume=independent["volume"], seconds=independent["seconds"])
            with job.part(f"{name}-midpoint") as part:
                part.update(0., "independent unsplit midpoint quadrature", force=True, n=1024, evaluations=0)
                row["midpoint"] = midpoint(width, 1024, lambda fraction, count:
                                           part.update(fraction, "independent unsplit midpoint quadrature", n=1024, evaluations=count))
                job.result(f"{name} midpoint", volume=row["midpoint"]["volume"])
            primary = [r["volume"] for r in row["reductions"]]
            last_delta = abs(primary[-1]-primary[-2])
            route_delta = abs(primary[-1]-independent["volume"])
            row.update(estimate=primary[-1],
                       empirical_stability=dict(first_refinement_difference=abs(primary[1]-primary[0]),
                           last_refinement_difference=last_delta, x_y_difference=route_delta,
                           conservative_observed_scale=max(last_delta, route_delta),
                           midpoint_difference=abs(primary[-1]-row["midpoint"]["volume"]),
                           scope="Observed differences only; not a rigorous error bound or confidence interval"),
                       difference_from_frozen_reference=primary[-1]-saved_reference["volume"],
                       source_reference_unchanged=True, reference_promoted=False)
            w.save(a.run/(name+".json"), row)
            estimates.append(row)
            w.save(a.run/"result.json", dict(complete=len(estimates)==len(CASES), rows=estimates,
                   inputs=hashes, controls=control, identities=identities, reference_promoted=False))
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    args.run = args.run.resolve()
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-gyroid_analytic_volume.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--worker", "--run", str(args.run)],
        args.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=180)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    saved = json.loads((args.run/"result.json").read_text())
    for row in saved["rows"]:
        result(row["case"]+" independent volume estimate", volume=row["estimate"],
               last_delta=row["empirical_stability"]["last_refinement_difference"],
               x_y_delta=row["empirical_stability"]["x_y_difference"],
               jsonl=str(args.run/"reported-results.jsonl"))


if __name__ == "__main__":
    main()
