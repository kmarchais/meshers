"""Forward/reverse distance breakdown for the fast modes (deterministic)."""

import sys

import meshers
import scorecard

for name in sys.argv[1:]:
    case = scorecard.CASES[name]
    for balanced in (False, True):
        mesh = meshers.generate_surface(
            case.function, bounds=(-0.5, 0.5) * 3, cells=24, band=case.band,
            geometry_tolerance=1e-3, adaptive=True, iterations=2, balanced=balanced, fast=True,
        )
        m = scorecard.score(case, mesh.points, mesh.triangles)
        print(name, balanced, f"fwd={m['forward_max']:.3e} rev={m['reverse_max']:.3e} "
              f"internal={mesh.diagnostics['maximum_sampled_surface_error']:.3e}")
