"""Recover unclipped targets from the immutable live coarse indicator."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow

np, *_ = workflow.imports()
from periodic import periodic_pairs

source = ROOT / "experiments/adaptive/runs/g3-lattice-sizing-pilot-v1"
output = Path(__file__).with_name("result.json")
assert not output.exists()
data = np.load(source / "coarse-recovery.npz")
points, tetrahedra = data["points"], data["tetrahedra"]
faces = np.sort(tetrahedra[:, [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]].reshape(-1, 3), axis=1)
unique, counts = np.unique(faces, axis=0, return_counts=True)
_, classes, _ = periodic_pairs(points, unique[counts == 1])
volumes = np.linalg.det(points[tetrahedra[:, 1:]] - points[tetrahedra[:, :1]]) / 6
lengths = np.cbrt(6 * volumes)
element_targets = np.maximum(data["estimated_mass"] / (volumes * lengths**2), 1e-30)**(-1 / 5)
nodal = np.full(len(points), np.inf)
np.minimum.at(nodal, classes[tetrahedra.ravel()], np.repeat(element_targets, 4))
nodal = nodal[classes]
row = json.loads((source / "recovery.json").read_text())
raw = row["nominal_length"] * nodal / np.median(nodal)
np.testing.assert_allclose(np.clip(raw, 0.03, 0.1), data["target_length"], rtol=1e-13, atol=1e-14)
final = np.load(source / "recovery-mesh.npz")
final_volumes = np.linalg.det(final["points"][final["tetrahedra"][:, 1:]] -
    final["points"][final["tetrahedra"][:, :1]]) / 6
workflow.save(output, {"inputs": {str(source / name): workflow.digest(source / name)
    for name in ("coarse-recovery.npz", "recovery.json", "recovery-mesh.npz")},
    "driver_sha256": workflow.digest(Path(__file__)),
    "raw_target_quantiles": np.quantile(raw, [0, 0.01, 0.05, 0.1, 0.5, 0.9, 1]).tolist(),
    "coarse_element_length_quantiles": np.quantile(lengths, [0, 0.1, 0.5, 0.9, 1]).tolist(),
    "fraction_nodes_clipped_at_floor": float(np.mean(raw < 0.03)),
    "final_effective_uniform_length": float(np.cbrt(6 * final_volumes.sum() / len(final_volumes))),
    "scope": "Target-scale and clipping diagnostic only; not proof of which native sizing constraint controls each edge"})
