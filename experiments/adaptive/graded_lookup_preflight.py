"""Check a native representation of the frozen density lookup, without meshing.

This uses existing comparisons/selects, not a new compiler opcode or a fitted
replacement field. It is a feasibility control, not a generation-speed result.
"""
import argparse
import hashlib
import json
from pathlib import Path
import signal
import time

import numpy as np
from meshers.expression import TraceError, compile_field

from comparison_fixtures import raw_gyroid


def native_width(table):
    knots = np.asarray(table["z"])
    half = np.asarray(table["half_width"])
    assert len(knots) == len(half) and np.all(np.diff(knots) > 0)

    def evaluate(x, y, z):
        def branch(lo, hi):
            if hi - lo == 1:
                slope = float((half[hi] - half[lo]) / (knots[hi] - knots[lo]))
                return float(half[lo]) + slope * (z - float(knots[lo]))
            middle = (lo + hi) // 2
            return np.where(z < float(knots[middle]), branch(lo, middle), branch(middle, hi))

        value = branch(0, len(knots) - 1)
        return np.where(z <= float(knots[0]), float(half[0]),
                        np.where(z >= float(knots[-1]), float(half[-1]), value))

    return evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    signal.alarm(45)
    table_path = Path(__file__).resolve().parents[2] / "handoff/macbook/linear-density-calibration-20261026.json"
    table = json.loads(table_path.read_text())
    knots, half = np.asarray(table["z"]), np.asarray(table["half_width"])

    def original(x, y, z):
        return raw_gyroid(x, y, z) / np.interp(z, knots, half)

    try:
        compile_field(original, gradient=False)
    except TraceError as error:
        original_error = str(error)
    else:
        raise AssertionError("The recorded unsupported interpolation no longer reproduces")

    width = native_width(table)

    def candidate(x, y, z):
        return raw_gyroid(x, y, z) / width(x, y, z)

    started = time.perf_counter()
    compiled_width = compile_field(width, gradient=False)
    compiled = compile_field(candidate, gradient=False)
    compile_seconds = time.perf_counter() - started
    assert not compiled.has_gradient
    # Every knot, adjacent representable values, endpoints, clamped exterior,
    # and independent random points. No frozen sample file is regenerated.
    rng = np.random.default_rng(20261026)
    z = np.r_[knots, np.nextafter(knots, -np.inf), np.nextafter(knots, np.inf), -1., 1.]
    edge_points = rng.uniform(-.5, .5, (len(z), 3))
    edge_points[:, 2] = z
    points = np.concatenate((edge_points, rng.uniform(-.5, .5, (8192, 3))))
    expected_width = np.interp(points[:, 2], knots, half)
    actual_width = compiled_width.values(points)
    expected_field = original(*points.T)
    actual_field = compiled.values(points)
    np.testing.assert_allclose(actual_width, expected_width, atol=1e-14, rtol=1e-14)
    np.testing.assert_allclose(actual_field, expected_field, atol=1e-12, rtol=1e-12)
    result = dict(passed=True, original_trace_error=original_error,
                  table_sha256=hashlib.sha256(table_path.read_bytes()).hexdigest(),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  probes=len(points), knots=len(knots), native_graph_nodes=compiled.node_count,
                  compile_seconds=compile_seconds,
                  maximum_width_difference=float(abs(actual_width - expected_width).max()),
                  maximum_field_difference=float(abs(actual_field - expected_field).max()),
                  native_analytic_gradient_enabled=False,
                  mesh_generated=False, speedup_measured=False,
                  scope="Existing native expression feasibility only. Both sides use the same frozen lookup; no mesh equivalence, admission or full generation speedup follows.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
