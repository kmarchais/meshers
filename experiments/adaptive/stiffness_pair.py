"""Compare admitted tensors or a previously justified reference summary."""

import argparse
import json
from pathlib import Path

import microgen_workflow as workflow


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=Path, required=True)
    normalization = parser.add_mutually_exclusive_group(required=True)
    normalization.add_argument("--second", type=Path)
    normalization.add_argument("--reference-summary", type=Path)
    parser.add_argument("--orders", nargs="+", type=int, default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.orders or (args.reference_summary and all(order > 0 for order in args.orders))
    if args.output.exists():
        raise FileExistsError(args.output)
    np, *_ = workflow.imports()
    from periodic import tensor_difference
    tensors, inputs = [], {}
    for path in ((args.first, args.second) if args.second else (args.first,)):
        row = json.loads(path.read_text())
        assert row["qualified"] and row["geometry"]["geometry_admitted"]
        assert row["status"] in ("solved", "solved-needs-reference")
        if row.get("order") == 2 and row.get("solver_quadrature_points", 15) == 4:
            assert row.get("quadrature_sensitivity_passed")
            assert row["quadrature_sensitivity_limit"] <= 1e-4
        assert len(row["solver_records"]) == 6
        assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
        tensor = np.asarray(row["stiffness"])
        assert tensor.shape == (6, 6) and np.isfinite(tensor).all()
        assert np.linalg.eigvalsh(tensor).min() > 0
        tensors.append(tensor)
        inputs[str(path)] = workflow.digest(path)
    if args.reference_summary:
        summary = json.loads(args.reference_summary.read_text())
        reference = summary["references"][row["case"]]
        assert reference["qualified_inputs"]
        tensor = np.asarray(reference["stiffness"])
        assert tensor.shape == (6, 6) and np.isfinite(tensor).all()
        assert np.linalg.eigvalsh(tensor).min() > 0
        tensors.append(tensor)
        inputs[str(args.reference_summary)] = workflow.digest(args.reference_summary)
    report = {"inputs": inputs, "script_sha256": workflow.digest(Path(__file__)),
        "scope": ("Previously justified reference comparison; not continuum certification" if args.reference_summary else
            "Pairwise admitted mesh tensor difference; no justified fine reference inferred"),
        "second_tensor_is_normalization": True, "difference": tensor_difference(*tensors)}
    if args.orders:
        levels = {record["cells"]: record for record in summary["rows"]
            if record["case"] == row["case"] and record["route"] == "meshers-curved"}
        coarse_row, fine_row = levels[16], levels[20]
        for record in (coarse_row, fine_row):
            assert record["qualified"] and record["geometry"]["geometry_admitted"]
            assert len(record["solver_records"]) == 6
            assert max(solver["relative_residual"] for solver in record["solver_records"]) < 2e-9
        coarse, fine = (np.asarray(record["stiffness"]) for record in (coarse_row, fine_row))
        ratios = {"uniform_seed": 16 / 20,
            "element_density": (coarse_row["elements"] / fine_row["elements"])**(1 / 3)}
        models = []
        for proxy, ratio in ratios.items():
            assert 0 < ratio < 1
            for assumed_order in args.orders:
                extrapolated = (fine - ratio**assumed_order * coarse) / (1 - ratio**assumed_order)
                models.append({"length_proxy": proxy, "length_ratio": ratio, "assumed_order": assumed_order,
                    "difference": tensor_difference(tensors[0], extrapolated)})
        differences = [model["difference"]["worst_relative_strain_energy_difference"] for model in models]
        report.update(declared_order_models=models, declared_order_difference_range=[min(differences), max(differences)],
            passes_one_percent_under_all_declared_models=max(differences) <= 0.01)
    workflow.save(args.output, report)


if __name__ == "__main__":
    main()
