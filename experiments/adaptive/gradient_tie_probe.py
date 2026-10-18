"""Check compiled gradients at smooth-union branch ties against the fixed field."""

import argparse
import os
from pathlib import Path
import resource
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def symmetric_lattice(coordinates_x, coordinates_y, coordinates_z):
    import numpy as np
    import scorecard

    distance = None
    for corner in scorecard.CORNERS:
        parameter = np.clip((coordinates_x * corner[0] + coordinates_y * corner[1]
                             + coordinates_z * corner[2]) / 0.75, 0, 1)
        radius = np.sqrt((coordinates_x - parameter * corner[0])**2
            + (coordinates_y - parameter * corner[1])**2 + (coordinates_z - parameter * corner[2])**2)
        if distance is None:
            distance = radius
        else:
            difference = np.abs(distance - radius)
            blend = np.maximum(0.08 - difference, 0.0) / 0.08
            distance = (distance + radius - difference) / 2 - blend**3 * 0.08 / 6
    return distance - 0.12


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    resource.setrlimit(resource.RLIMIT_AS, (workflow.MEMORY, workflow.MEMORY))
    with Job(parts=["gradient-ties"]) as job:
        with job.part("gradient-ties") as part:
            part.update(0, "Checking smooth-union derivatives at branch ties", probes=6)
            np, _, meshers, _, _, scorecard = workflow.imports()
            case, _, _ = workflow.specimen("smooth-lattice")
            original = meshers.compile_field(case.function)
            symmetric = meshers.compile_field(symmetric_lattice)
            queries = np.array([[0., 0.17, -0.17], [0., -0.17448645899797485, 0.1744864589979749],
                [0.17, 0., -0.17], [-0.17, 0.17, 0.], [0., -0.17, -0.17], [0.17, 0., 0.17]])
            numeric = scorecard.gradient(case.function, queries)
            original_gradient = np.array([original.gradient(*point) for point in queries])
            symmetric_gradient = np.array([symmetric.gradient(*point) for point in queries])
            random = np.random.default_rng(20261018).uniform(-0.5, 0.5, (4096, 3))
            expected, actual = case.function(*random.T), symmetric.values(random)
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
            assert np.array_equal(actual < 0, expected < 0)
            original_error = float(abs(original_gradient - numeric).max())
            symmetric_error = float(abs(symmetric_gradient - numeric).max())
            assert original_error > 1e-3 and symmetric_error < 1e-7
            workflow.save(args.output, {"script_sha256": workflow.digest(Path(__file__)),
                "specimen_source_sha256": workflow.digest(Path(workflow.__file__)),
                "queries": queries.tolist(), "finite_difference_gradient": numeric.tolist(),
                "original_compiled_gradient": original_gradient.tolist(),
                "symmetric_compiled_gradient": symmetric_gradient.tolist(),
                "original_max_gradient_error": original_error, "symmetric_max_gradient_error": symmetric_error,
                "maximum_field_value_difference": float(abs(actual - expected).max()),
                "field_value_verification_points": len(random),
                "identity": "min(a,b) = (a+b-abs(a-b))/2; physical smooth-union field unchanged"})
            job.result("gradient-ties", original_error=original_error, symmetric_error=symmetric_error)


if __name__ == "__main__":
    main()
