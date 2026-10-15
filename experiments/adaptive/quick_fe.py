"""Indicative Flying Edges timings across grid sizes (not a timing study)."""

import sys
import time

import scorecard

case = scorecard.CASES[sys.argv[1]]
for n in map(int, sys.argv[2:]):
    scorecard.flying_edges(case, n)
    start = time.perf_counter()
    points, faces = scorecard.flying_edges(case, n)
    seconds = time.perf_counter() - start
    m = scorecard.score(case, points, faces)
    print(f"fe-{n:<4d} {seconds:6.3f}s tri={m['triangles']:6d} H={m['hausdorff']:.2e} "
          f"dV={m['volume_error']:+.3%}", flush=True)
