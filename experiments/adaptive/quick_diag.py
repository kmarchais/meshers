"""Print surface diagnostics timings for the fast mode (indicative)."""

import time

import meshers
import scorecard

case = scorecard.CASES["gyroid"]
for _ in range(3):
    start = time.perf_counter()
    mesh = meshers.generate_surface(
        case.function, bounds=(-0.5, 0.5) * 3, cells=24, band=case.band,
        geometry_tolerance=1e-3, adaptive=True, iterations=2, balanced=True, fast=True,
    )
    total = time.perf_counter() - start
    d = {k: round(v, 4) for k, v in mesh.diagnostics.items() if "second" in k}
    d.update({k: round(v, 4) for k, v in mesh.diagnostics["adaptive"]["seconds"].items() if v})
    print(f"{total:.3f}", d)
