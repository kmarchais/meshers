"""Summarize field-audited adaptive tensor refinement and declared model sensitivity."""

import argparse
import json
from pathlib import Path

import microgen_workflow as workflow


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--coarse-tolerance", type=float, required=True)
    parser.add_argument("--fine-tolerance", type=float, required=True)
    parser.add_argument("--field-audit", type=Path, required=True)
    parser.add_argument("--candidate", action="append", type=Path, default=[])
    parser.add_argument("--target", type=float, default=0.01)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    np, *_ = workflow.imports()
    from periodic import tensor_difference

    audit_manifest = args.field_audit / "manifest.json"
    audited = {Path(path).resolve(): value for path, value in
               json.loads(audit_manifest.read_text())["inputs"].items()}
    inputs = {str(audit_manifest): workflow.digest(audit_manifest)}
    def read_record(path):
        row = json.loads(path.read_text())
        assert row["qualified"] and row["geometry"]["geometry_admitted"]
        assert len(row["solver_records"]) == 6
        assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
        audit_path = args.field_audit / (path.parent.name + ".json")
        audit = json.loads(audit_path.read_text())
        assert Path(audit["source"]).resolve() == path.parent.resolve()
        assert max(audit["reaction_relative_error"], audit["hill_mandel_relative_error"],
                   audit["integrated_stiffness_skew_relative"]) < 1e-7
        for source in (path, path.parent / "mesh.npz", path.parent / "fields.npz"):
            value = workflow.digest(source)
            assert audited[source.resolve()] == value
            inputs[str(source)] = value
        inputs[str(audit_path)] = workflow.digest(audit_path)
        tensor = np.asarray(row["stiffness"])
        assert tensor.shape == (6, 6) and np.isfinite(tensor).all()
        assert np.allclose(tensor, tensor.T, atol=1e-10, rtol=0)
        assert np.linalg.eigvalsh(tensor).min() > 0
        return row, tensor

    coarse_row, coarse = read_record(args.coarse)
    fine_row, fine = read_record(args.fine)
    assert coarse_row["case"] == fine_row["case"]
    assert 0 < args.fine_tolerance < args.coarse_tolerance
    for row, tolerance in ((coarse_row, args.coarse_tolerance), (fine_row, args.fine_tolerance)):
        assert row.get("adapt_tolerance", tolerance) == tolerance
    ratios = {"target_length": (args.fine_tolerance / args.coarse_tolerance)**0.5,
              "element_density": (coarse_row["elements"] / fine_row["elements"])**(1 / 3)}
    models = []
    for proxy, ratio in ratios.items():
        assert 0 < ratio < 1
        for order in (2, 3, 4):
            extrapolated = (fine - ratio**order * coarse) / (1 - ratio**order)
            assert np.isfinite(extrapolated).all() and np.linalg.eigvalsh(extrapolated).min() > 0
            models.append({"length_proxy": proxy, "length_ratio": ratio, "assumed_order": order,
                "stiffness": extrapolated.tolist(), "fine_to_model": tensor_difference(fine, extrapolated)})
    candidates = []
    for path in args.candidate:
        row, tensor = read_record(path)
        assert row["case"] == fine_row["case"]
        differences = [tensor_difference(tensor, np.asarray(model["stiffness"]))[
            "worst_relative_strain_energy_difference"] for model in models]
        candidates.append({"source": str(path), "nodes": row["nodes"], "elements": row["elements"],
            "total_workflow_s": row.get("total_workflow_s", row.get("total_audit_workflow_s")),
            "difference_from_fine": tensor_difference(tensor, fine),
            "declared_model_difference_range": [min(differences), max(differences)],
            "passes_target_under_all_declared_models": max(differences) <= args.target})
    workflow.save(args.output, {"case": fine_row["case"], "inputs": inputs,
        "script_sha256": workflow.digest(Path(__file__)), "primary_reference_stiffness": fine.tolist(),
        "primary_reference_source": str(args.fine), "fine_pair_difference": tensor_difference(coarse, fine),
        "models": models, "candidates": candidates, "stiffness_target": args.target,
        "scope": "Empirical adaptive refinement sensitivity; assumptions are not certified continuum error bounds"})


if __name__ == "__main__":
    main()
