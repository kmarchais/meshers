"""Issue #4 BCC capsule lattice with Microgen's polynomial smooth union.

Hard `min` unions have concave creases (and stress singularities). Microgen's
`implicit_ops.smooth_union` uses the cubic polynomial smooth minimum below.
Run: uv run --no-sync python experiments/adaptive/smooth_lattice.py [blend ...]
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fea  # noqa: E402, F401  (frozen build path)
import meshers  # noqa: E402
import scorecard  # noqa: E402


def lattice(blend):
    def smin(a, b):
        h = np.maximum(blend - np.abs(a - b), 0.0) / blend
        return np.minimum(a, b) - h * h * h * blend / 6.0

    def field(x, y, z):
        d = None
        for c in scorecard.CORNERS:
            t = np.clip((x * c[0] + y * c[1] + z * c[2]) / (c @ c), 0, 1)
            r = np.sqrt((x - t * c[0]) ** 2 + (y - t * c[1]) ** 2 + (z - t * c[2]) ** 2)
            d = r if d is None else smin(d, r)
        return d - 0.12

    return field


for blend in map(float, sys.argv[1:] or ["0.03", "0.08"]):
    field = lattice(blend)
    case = scorecard.Case(f"smooth-lattice-{blend}", field, (-np.inf, 0.0))
    walls = scorecard.wall_samples(case)
    for cells in (8, 12, 16):
        for label, extra in (
            ("lattice", {"optimize_passes": 4}),
            ("adapt", {"optimize_passes": 1, "adapt": True, "max_length": 0.25}),
        ):
            try:
                mesh = meshers.generate(
                    field,
                    bounds=(-0.5, 0.5) * 3,
                    cells=cells,
                    minimum_quality=0.1,
                    geometry_tolerance=0.01,
                    max_tetrahedra=250_000,
                    threads=1,
                    **extra,
                )
            except Exception as error:  # noqa: BLE001 - rejections are results
                print(f"blend={blend} cells={cells:2d} {label:7s} REJECTED {str(error)[:60]}")
                continue
            r = scorecard.score_volume(case, mesh, walls)
            print(
                f"blend={blend} cells={cells:2d} {label:7s} ok tets={r['tetrahedra']:6d} "
                f"H={r['hausdorff']:.2e} dV={r['tet_volume_error']:+.4f} q_min={r['q_min']:.3f}"
            )
