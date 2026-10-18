"""Compare audited lattice tensors without removing retained model families."""

import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "experiments/adaptive"))
import microgen_workflow as workflow


def witness(candidate, reference):
    values, vectors = eigh(candidate - reference, reference)
    index = int(np.argmax(abs(values)))
    strain = vectors[:, index] / np.linalg.norm(vectors[:, index])
    measured = float(strain @ candidate @ strain / (strain @ reference @ strain) - 1)
    assert abs(measured - values[index]) < 1e-12
    return {"worst_relative_strain_energy_difference": abs(measured),
        "signed_relative_energy_difference": measured,
        "engineering_strain_in_saved_tensor_basis": strain.tolist()}


def main():
    output = Path(__file__).with_name("result.json")
    assert not output.exists()
    for scale in (.96, 1, 1.04):
        assert abs(witness(np.eye(6) * scale, np.eye(6))[
            "worst_relative_strain_energy_difference"] - abs(scale - 1)) < 1e-14
    runs = ROOT / "experiments/adaptive/runs"
    inputs = {str(Path(__file__)): workflow.digest(Path(__file__))}

    def read(path):
        inputs[str(path)] = workflow.digest(path)
        return json.loads(path.read_text())

    models = []
    for name in ("g3-lattice-reference-sensitivity-v2.json",
        "g3-lattice-microgen-reference-comparison-v1.json",
        "g3-lattice-finer-reference-comparison-v1.json"):
        source = runs / name
        report = read(source)
        assert report["case"] == "smooth-lattice"
        for index, model in enumerate(report["models"]):
            tensor = np.asarray(model["stiffness"])
            assert tensor.shape == (6, 6) and np.isfinite(tensor).all()
            assert np.linalg.eigvalsh(tensor).min() > 0
            models.append((tensor, {"source": str(source), "index": index,
                "length_proxy": model["length_proxy"], "assumed_order": model["assumed_order"]}))
    assert len(models) == 24
    candidates = []
    audit = runs / "g3-lattice-finer-reference-field-audit-v1"
    sources = [(runs / name, audit) for name in (
        "g3-lattice-adapt24-solve-v1", "g3-lattice-adapt24-tol003-v1",
        "g3-lattice-adapt24-tol0028-v1", "g3-lattice-curved24-symmetric-v1",
        "g3-lattice-curved32-reference-v1", "g3-lattice-linear32-control-v1",
        "g3-lattice-tuned51-timing-v1")]
    sources.extend((runs / name, runs / audit_name) for name, audit_name in (
        ("g3-lattice-uniform34-quadrature4-v1", "g3-lattice-uniform34-field-audit-v1"),
        ("g3-lattice-microgen-straight-tet10-hausd001-amg-v2", "g3-lattice-microgen-tet10-field-audit-v1")))
    for source, audit in sources:
        row = read(source / "result.json")
        assert row["qualified"] and row["geometry"]["geometry_admitted"]
        assert row["case"] == "smooth-lattice" and row["elements"] <= workflow.CAP
        assert row["status"] in ("solved", "solved-needs-reference")
        assert len(row["solver_records"]) == 6
        assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
        if row.get("order") == 2 and row.get("solver_quadrature_points", 15) == 4:
            assert row["quadrature_sensitivity_passed"] and row["quadrature_sensitivity_limit"] <= 1e-4
        manifest = read(audit / "manifest.json")
        audited = {Path(path).resolve(): digest for path, digest in manifest["inputs"].items()}
        for name in ("result.json", "mesh.npz", "fields.npz"):
            path = source / name
            inputs[str(path)] = workflow.digest(path)
            assert inputs[str(path)] == audited[path.resolve()]
        checked = read(audit / (source.name + ".json"))
        assert Path(checked["source"]).resolve() == source.resolve()
        assert max(checked["reaction_relative_error"], checked["hill_mandel_relative_error"],
            checked["integrated_stiffness_skew_relative"]) < 1e-7
        tensor = np.asarray(row["stiffness"])
        errors = [{**metadata, **witness(tensor, reference)} for reference, metadata in models]
        worst = max(errors, key=lambda error: error["worst_relative_strain_energy_difference"])
        candidates.append({"source": str(source), "nodes": row["nodes"],
            "minimum_model_disagreement": min(error["worst_relative_strain_energy_difference"] for error in errors),
            "worst_model": worst,
            "passes_one_percent_all_retained_models": worst["worst_relative_strain_energy_difference"] <= .01})
    workflow.save(output, {"inputs": inputs, "retained_model_count": len(models),
        "candidates": candidates, "stiffness_target": .01,
        "scope": "Audited empirical model disagreements only; no model removed, continuum bound or live timing comparison"})


if __name__ == "__main__":
    main()
