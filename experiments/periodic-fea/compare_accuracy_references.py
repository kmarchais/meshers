"""Compare the same candidates against two independent volume topologies."""

import argparse
import sys
from pathlib import Path

import accuracy_metrics
from run import digest, write_json

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report


def compare(references, candidates, output, geometry_gate):
    # Read-only inputs and atomic fixed-name outputs make an interrupted run restartable.
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, reference in enumerate(references):
        accuracy_metrics.report = lambda fraction, message, **kwargs: report(
            (index + fraction) / len(references), f"Reference {index + 1}/{len(references)}: {message}", **kwargs)
        rows.append(accuracy_metrics.compare(reference, candidates, output / f"reference-{index}.json",
                                             geometry_gate=geometry_gate))
    write_json(output / "comparisons.json", {"state": "complete", "source_sha256": digest(Path(__file__)),
                                             "common_geometry_gate": geometry_gate, "comparisons": rows})
    report(1, "Accuracy comparisons against both reference topologies complete", eta_s=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, nargs=2, required=True)
    parser.add_argument("--candidates", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--common-geometry-gate", type=float, default=.002)
    args = parser.parse_args()
    compare(args.references, args.candidates, args.output, args.common_geometry_gate)
