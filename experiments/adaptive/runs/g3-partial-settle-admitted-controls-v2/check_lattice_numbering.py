"""Check whether lattice replay differences are node numbering rather than geometry."""

import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow

run = Path(__file__).parent
sources = [ROOT / "experiments/adaptive/runs/g3-lattice-adapt24-verified-v1/mesh.npz",
           run / "g3-lattice-adapt24-verified-v1/mesh.npz"]
inputs = {str(path): workflow.digest(path) for path in [Path(__file__), *sources]}
old, new = [np.load(path) for path in sources]
distance, mapping = cKDTree(old["points"]).query(new["points"])
assert len(np.unique(mapping)) == len(old["points"]) == len(new["points"])
assert distance.max() < 1e-12

def canonical(array):
    return np.array(sorted(map(tuple, np.sort(array, axis=1))))

matches = {key: np.array_equal(canonical(mapping[new[key]]), canonical(old[key]))
           for key in ("tetrahedra", "surface")}
assert all(matches.values())
output = run / "lattice-numbering.json"
assert not output.exists()
workflow.save(output, {"inputs": inputs, "maximum_matched_node_distance": float(distance.max()),
    "bijective_node_match": True, "canonical_connectivity_matches": matches,
    "scope": "Point coordinates and cell node sets up to numbering; no new equilibrium solve"})
assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())
print(json.dumps({"maximum_matched_node_distance": float(distance.max()), **matches}))
