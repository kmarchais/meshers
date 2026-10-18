"""Audit and compare the cap pilot without admitting its interrupted timing."""

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow

runs = ROOT / "experiments/adaptive/runs"
output = Path(__file__).with_name("result.json")
assert not output.exists()
source = runs / "g3-lattice-cap-size015-band005-solve-v1"
audit = runs / "g3-lattice-cap-size015-band005-final-audit-v1"
inputs = {str(Path(__file__)): workflow.digest(Path(__file__))}


def read(path):
    inputs[str(path)] = workflow.digest(path)
    return json.loads(path.read_text())


row = read(source / "result.json")
assert row["qualified"] and row["geometry"]["geometry_admitted"]
assert row["status"] == "solved-needs-reference" and len(row["solver_records"]) == 6
assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
manifest = read(audit / "manifest.json")
audited = {Path(path).resolve(): digest for path, digest in manifest["inputs"].items()}
for name in ("result.json", "mesh.npz", "fields.npz"):
    path = source / name
    inputs[str(path)] = workflow.digest(path)
    assert inputs[str(path)] == audited[path.resolve()]
consistency = read(audit / (source.name + "-consistency.json"))
assert max(consistency["reaction_relative_error"], consistency["hill_mandel_relative_error"],
    consistency["integrated_stiffness_skew_relative"]) < 1e-7
quadrature = read(audit / (source.name + ".json"))
assert quadrature["rules"][-1]["points"] == 125
assert quadrature["rules"][-1]["energy_vs_solved_tensor"]["worst_relative_strain_energy_difference"] < 1e-4
timing = read(runs / "g3-lattice-cap-size015-band005-timing-exclusion-v1.json")
assert not timing["timing_admitted"] and timing["source_result_sha256"] == inputs[str(source / "result.json")]
tensor = np.asarray(row["stiffness"])
errors = []
for name in ("g3-lattice-reference-sensitivity-v2.json",
    "g3-lattice-microgen-reference-comparison-v1.json", "g3-lattice-finer-reference-comparison-v1.json"):
    report = read(runs / name)
    for index, model in enumerate(report["models"]):
        reference = np.asarray(model["stiffness"])
        values = eigh(tensor - reference, reference, eigvals_only=True)
        errors.append({"source": name, "index": index, "length_proxy": model["length_proxy"],
            "assumed_order": model["assumed_order"], "relative_energy_eigenvalues": values.tolist(),
            "worst_relative_strain_energy_difference": float(abs(values).max())})
assert len(errors) == 24
values = [error["worst_relative_strain_energy_difference"] for error in errors]
workflow.save(output, {"inputs": inputs, "source": str(source), "nodes": row["nodes"],
    "elements": row["elements"], "models": errors, "retained_model_count": len(errors),
    "model_difference_range": [min(values), max(values)],
    "passes_one_percent_all_retained_models": max(values) <= .01,
    "timing_admitted": False,
    "scope": "Audited negative sizing pilot; all models retained; suspended-host timing excluded; no continuum certificate"})
