"""Check actual curve-edge sizes against the live pilot's requested metric."""

import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow
from periodic_crease_audit import features

started = time.perf_counter()
run = Path(__file__).parent
assert not (run / "result.json").exists()
source = ROOT / "experiments/adaptive/runs/g3-lattice-curve-grading-solve-v1/curve"
mesh = np.load(source / "mesh.npz")
row = json.loads((source / "result.json").read_text())
assert row["qualified"] and row["status"] == "solved-needs-reference"
cloud = np.load(source / "curve-samples.npz")["points"]
tree = cKDTree(cloud, boxsize=1)
data = features(mesh)
edges = np.asarray(data["curve_edges"])
locations = mesh["points"][edges]
distance, _ = tree.query(((locations + .5) % 1).reshape(-1, 3))
targets = np.minimum(.0025 + .8 * distance.reshape(-1, 2), .1)
lengths = np.linalg.norm(locations[:, 1] - locations[:, 0], axis=1)
ratios = lengths / targets.mean(axis=1)
cutoff = 4 / 3
rust = ROOT / "crates/meshers-core/src/tet_adapt.rs"
assert "const SPLIT: f64 = 4. / 3.;" in rust.read_text()
paths = [Path(__file__), source / "result.json", source / "mesh.npz",
         source / "curve-samples.npz", ROOT / "experiments/adaptive/periodic_crease_audit.py", rust]
inputs = {str(path): workflow.digest(path) for path in paths}
worst = []
for index in np.argsort(ratios)[-10:][::-1]:
    worst.append({"edge": edges[index].tolist(), "points": locations[index].tolist(),
        "length": float(lengths[index]), "endpoint_targets": targets[index].tolist(),
        "length_target_ratio": float(ratios[index])})
result = {"inputs": inputs, "curve_edges": len(edges),
    "edge_length_quantiles": np.quantile(lengths, [0, .1, .5, .9, 1]).tolist(),
    "length_target_ratio_quantiles": np.quantile(ratios, [0, .1, .5, .9, 1]).tolist(),
    "native_split_ratio": cutoff, "edges_above_split_ratio": int(np.count_nonzero(ratios > cutoff + 1e-8)),
    "fraction_above_split_ratio": float(np.mean(ratios > cutoff + 1e-8)), "worst_edges": worst,
    "seconds": time.perf_counter() - started,
    "scope": "Endpoint size-field realization only; not every native split must be admissible and no quality gate is relaxed"}
workflow.save(run / "result.json", result)
print(json.dumps({key: value for key, value in result.items() if key not in ("inputs", "worst_edges")}))
