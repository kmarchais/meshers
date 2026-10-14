"""Reproduce issues #3 and #4 (cells 8, tolerance 0.01, minimum quality 0.1).

The lattice pipeline is run with 0, 4 and 12 optimization passes, then the size-driven
adapter from the same 8-cell seed. Successful meshes get the independent two-sided
distance and volume checks of scorecard.py.
Run: uv run --no-sync python experiments/adaptive/issues.py
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "microgen-benchmark/benchmarks/microgen"))

import fea  # noqa: E402, F401  (frozen build path)
import meshers  # noqa: E402
import scorecard  # noqa: E402
from microgen_benchmark.cases import CASES  # noqa: E402

NAMES = ["gyroid", "schwarz-p", "neovius", "graded-gyroid", "thin-gyroid", "implicit-lattice"]


def wrap(case):
    levels = (-np.inf, case.band[1]) if case.implicit else tuple(case.band)
    return scorecard.Case(f"issue-{case.id}", case.field, levels)


def main():
    cases = {c.id: c for c in CASES}
    for name in NAMES:
        case = cases[name]
        local = wrap(case)
        walls = scorecard.wall_samples(local)
        runs = [(f"lattice passes={p}", {"optimize_passes": p}) for p in (0, 4, 12)]
        runs.append(("adapt", {"optimize_passes": 1, "adapt": True, "max_length": 0.25}))
        for label, extra in runs:
            try:
                mesh = meshers.generate(
                    case.field,
                    bounds=(-0.5, 0.5) * 3,
                    cells=8,
                    band=None if case.implicit else case.band,
                    minimum_quality=0.1,
                    geometry_tolerance=0.01,
                    max_tetrahedra=250_000,
                    threads=1,
                    **extra,
                )
            except Exception as error:  # noqa: BLE001 - rejections are results
                print(f"{name:17s} {label:18s} REJECTED {str(error)[:80]}")
                continue
            r = scorecard.score_volume(local, mesh, walls)
            print(
                f"{name:17s} {label:18s} ok tets={r['tetrahedra']:6d} H={r['hausdorff']:.2e} "
                f"dV={r['tet_volume_error']:+.4f} q_min={r['q_min']:.3f}"
            )


if __name__ == "__main__":
    main()
