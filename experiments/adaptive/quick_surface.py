"""Indicative single-run timings of surface producers (not a timing study)."""

import sys
import time

import meshers
import scorecard

case = scorecard.CASES[sys.argv[1] if len(sys.argv) > 1 else "gyroid"]


def surface(balanced, fast=False, cells=24, iterations=5):
    mesh = meshers.generate_surface(
        case.function, bounds=(-0.5, 0.5) * 3, cells=cells, band=case.band,
        geometry_tolerance=1e-3, adaptive=True, iterations=iterations, balanced=balanced,
        fast=fast,
    )
    return mesh.points, mesh.triangles


for label, make in [
    ("fe-40", lambda: scorecard.flying_edges(case, 40)),
    ("fe-48", lambda: scorecard.flying_edges(case, 48)),
    ("adaptive", lambda: surface(False)),
    ("balanced", lambda: surface(True)),
    ("fast", lambda: surface(False, True, iterations=2)),
    ("fast-bal", lambda: surface(True, True, iterations=2)),
]:
    try:
        make()
    except Exception as error:  # noqa: BLE001 - report and continue
        print(f"{label:10s} failed: {error}")
        continue
    start = time.perf_counter()
    points, faces = make()
    seconds = time.perf_counter() - start
    m = scorecard.score(case, points, faces)
    print(f"{label:10s} {seconds:6.3f}s tri={m['triangles']:6d} H={m['hausdorff']:.2e} "
          f"dV={m['volume_error']:+.3%} amin={m['angle_min_deg']:.1f}", flush=True)
