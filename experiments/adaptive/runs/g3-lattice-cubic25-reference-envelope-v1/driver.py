"""Read admitted cubic25 and score it against all unchanged model tensors."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow
from cubic_reference import cubic_node_count

np, *_ = workflow.imports()
from periodic import tensor_difference

run = Path(__file__).parent
assert not (run / "result.json").exists()
base = ROOT / "experiments/adaptive/runs"
source = base / "g3-lattice-cubic25-reference-blocked-symmetry-v2/result.json"
assessment = base / "g3-lattice-cubic25-assessment-v1/result.json"
manifest = base / "g3-lattice-cubic25-assessment-v1/manifest.json"
parent = base / "g3-lattice-cubic-parent23-preflight-v1"
paths = [base / name for name in (
    "g3-lattice-curve-grading-reference-comparison-v1.json",
    "g3-lattice-reference-sensitivity-v2.json",
    "g3-lattice-microgen-reference-comparison-v1.json")]
audit = json.loads(assessment.read_text())
row = json.loads(source.read_text())
assert audit["fresh_process_audit_passed"] and audit["retained_model_count"] == 24
assert Path(audit["source"]).resolve() == source.parent
assert row["qualified"] and row["completed_loads"] == 6 and row["quadrature_sensitivity_passed"]
frozen = {Path(path).resolve(): value for path, value in json.loads(manifest.read_text())["inputs"].items()}
for path in [source, *paths]:
    assert workflow.digest(path) == frozen[path]
models = []
for path in paths:
    report = json.loads(path.read_text())
    assert report["case"] == row["case"]
    for index, model in enumerate(report["models"]):
        models.append({"source": str(path), "index": index,
            "difference_normalized_by_model": tensor_difference(np.asarray(row["stiffness"]), np.asarray(model["stiffness"]))})
assert len(models) == 24
differences = [model["difference_normalized_by_model"]["worst_relative_strain_energy_difference"] for model in models]
parent_row = json.loads((parent / "result.json").read_text())
assert parent_row["qualified"] and parent_row["geometry"]["geometry_admitted"]
count = cubic_node_count(np.load(parent / "mesh.npz")["tetrahedra"])
inputs = [Path(__file__), source, assessment, manifest, *paths,
          parent / "result.json", parent / "mesh.npz", ROOT / "experiments/adaptive/cubic_reference.py"]
workflow.save(run / "result.json", {"inputs": {str(path): workflow.digest(path) for path in inputs},
    "source": str(source), "models": models, "model_count": len(models),
    "model_normalized_disagreement_range": [min(differences), max(differences)],
    "passes_one_percent_all_retained_models": max(differences) <= .01,
    "parent23_cubic_node_count": count, "parent23_within_existing_node_cap": count <= 250000,
    "scope": "Model-normalized disagreement, not an independent continuum bound; parent23 node-budget preflight, not a solve"})
print("Cubic25 model-normalized range:", min(differences), max(differences))
print("Parent23 exact cubic nodes:", count)
