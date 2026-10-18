"""Choose one activation scale from the admitted control's actual torus edges."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow

np, *_ = workflow.imports()
from periodic import periodic_pairs

source = ROOT / "experiments/adaptive/runs/g3-lattice-sizing-minbound-v1/g3-lattice-sizing-minbound-v1-uniform"
output = Path(__file__).with_name("result.json")
assert not output.exists()
row = json.loads((source / "result.json").read_text())
assert row["qualified"] and row["geometry"]["geometry_admitted"]
mesh = np.load(source / "mesh.npz")
_, classes, _ = periodic_pairs(mesh["points"], mesh["surface"])
tetrahedra = classes[mesh["tetrahedra"]]
edges = np.concatenate([tetrahedra[:, pair] for pair in ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))])
edges = np.unique(np.sort(edges, axis=1), axis=0)
delta = mesh["points"][edges[:, 1]] - mesh["points"][edges[:, 0]]
delta = (delta + 0.5) % 1 - 0.5
lengths = np.linalg.norm(delta, axis=1)
assert np.isfinite(lengths).all() and (lengths > 0).all()
workflow.save(output, {"inputs": {str(source / name): workflow.digest(source / name)
    for name in ("result.json", "mesh.npz")}, "driver_sha256": workflow.digest(Path(__file__)),
    "unique_torus_edges": len(edges), "edge_length_quantiles": np.quantile(lengths, [0, 0.1, 0.5, 0.9, 1]).tolist(),
    "selected_nominal_length": float(np.median(lengths)),
    "selection_rule": "One activation trial at the existing admitted mesh's median unique torus-edge length; no result-based scale sweep",
    "scope": "Observed geometric scale selection only, not a claim of optimal refinement or matched-accuracy improvement"})
