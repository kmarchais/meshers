"""Issue #4 spinodoid in the PR #2 environment, which provides Microgen's Spinodoid.

Run from the benchmark worktree with this branch's frozen build first on the path:
uv run --locked python <this file> <build directory>
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, sys.argv[1])
sys.path.insert(1, str(Path(__file__).resolve().parent))

import meshers  # noqa: E402
import numpy as np  # noqa: E402
import scorecard  # noqa: E402
from microgen_benchmark.cases import CASES  # noqa: E402

case = {c.id: c for c in CASES}["spinodoid"]
local = scorecard.Case("issue-spinodoid", case.field, (-np.inf, 0.0))
walls = scorecard.wall_samples(local)
print("meshers from", meshers.__file__)
for cells, extra in (
    (8, {"optimize_passes": 4}),
    (8, {"optimize_passes": 1, "adapt": True, "max_length": 0.25}),
    (16, {"optimize_passes": 4}),
    (16, {"optimize_passes": 1, "adapt": True, "max_length": 0.25}),
):
    label = f"cells={cells} {'adapt' if extra.get('adapt') else 'lattice'}"
    start = time.perf_counter()
    try:
        mesh = meshers.generate(
            case.field,
            bounds=(-0.5, 0.5) * 3,
            cells=cells,
            minimum_quality=0.1,
            geometry_tolerance=0.01,
            max_tetrahedra=250_000,
            threads=1,
            **extra,
        )
    except Exception as error:  # noqa: BLE001 - rejections are results
        print(f"{label:20s} REJECTED {str(error)[:90]}")
        continue
    r = scorecard.score_volume(local, mesh, walls)
    print(
        f"{label:20s} ok tets={r['tetrahedra']:6d} H={r['hausdorff']:.2e} "
        f"dV={r['tet_volume_error']:+.4f} q_min={r['q_min']:.3f} "
        f"({time.perf_counter() - start:.1f}s incl. callbacks)"
    )
