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
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    np, *_ = workflow.imports()
    from periodic import tensor_difference
    tensors, inputs = [], {}
    for path in ((args.first, args.second) if args.second else (args.first,)):
        row = json.loads(path.read_text())
        assert row["qualified"] and row["geometry"]["geometry_admitted"]
        assert len(row["solver_records"]) == 6
        assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
        tensor = np.asarray(row["stiffness"])
        assert tensor.shape == (6, 6) and np.isfinite(tensor).all()
        assert np.linalg.eigvalsh(tensor).min() > 0
        tensors.append(tensor)
        inputs[str(path)] = workflow.digest(path)
    if args.reference_summary:
        reference = json.loads(args.reference_summary.read_text())["references"][row["case"]]
        assert reference["qualified_inputs"]
        tensor = np.asarray(reference["stiffness"])
        assert tensor.shape == (6, 6) and np.isfinite(tensor).all()
        assert np.linalg.eigvalsh(tensor).min() > 0
        tensors.append(tensor)
        inputs[str(args.reference_summary)] = workflow.digest(args.reference_summary)
    workflow.save(args.output, {"inputs": inputs, "script_sha256": workflow.digest(Path(__file__)),
        "scope": ("Previously justified reference comparison; not continuum certification" if args.reference_summary else
            "Pairwise admitted mesh tensor difference; no justified fine reference inferred"),
        "second_tensor_is_normalization": True, "difference": tensor_difference(*tensors)})


if __name__ == "__main__":
    main()
