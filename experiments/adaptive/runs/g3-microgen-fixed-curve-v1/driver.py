"""Check whether rejected Schwarz probes lie on unchanged VTK cap edges."""

import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job
import microgen_workflow as workflow


def main():
    import numpy as np
    from scipy.spatial import cKDTree
    from cap_quality_bound import caps

    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--diagnostic", type=Path, required=True)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (workflow.MEMORY, workflow.MEMORY))
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    inputs = {str(path): workflow.digest(path) for path in (args.seed, args.mesh, args.diagnostic)}
    workflow.save(args.run / "manifest.json", {"inputs": inputs,
        "script_sha256": workflow.digest(Path(__file__)), "memory_bytes": workflow.MEMORY,
        "threads": 1, "distance_gate": 0.002})
    with Job(parts=["retained-curves"]) as job:
        with job.part("retained-curves") as part:
            part.update(0, "Matching rejected cap-curve edges to the original VTK seed", meshes=2)
            seed, mesh = np.load(args.seed), np.load(args.mesh)
            diagnostic = json.loads(args.diagnostic.read_text())
            assert diagnostic["case"] == "schwarz-solid"
            seed_caps, _ = caps(seed)
            seed_edges = {tuple(sorted((int(face[first]), int(face[second]))))
                for face in seed_caps for first, second in ((0, 1), (1, 2), (2, 0))}
            tree = cKDTree(seed["points"])
            examples = []
            for probe in diagnostic["geometry"]["worst_boundary_samples"]:
                coordinates = np.asarray(probe["coordinates"])
                point = np.asarray(probe["probe"])
                for axis in np.flatnonzero(probe["locked"]):
                    edge_nodes = np.flatnonzero(abs(coordinates[:, axis] - point[axis]) < 1e-12)
                    if len(edge_nodes) != 2:
                        continue
                    distances, seed_nodes = tree.query(coordinates[edge_nodes])
                    retained = bool(distances.max() < 1e-10
                        and tuple(sorted(map(int, seed_nodes))) in seed_edges)
                    lower = np.maximum(point - 0.002, -0.5)
                    upper = np.minimum(point + 0.002, 0.5)
                    maximum_absolute_coordinate = np.maximum(abs(lower), abs(upper))
                    minimum_field_in_coordinate_box = float(np.cos(
                        2 * np.pi * maximum_absolute_coordinate).sum())
                    examples.append({"probe": point.tolist(), "sampled_distance": probe["distance"],
                        "final_edge_nodes": mesh["surface"][probe["face"]][edge_nodes].tolist(),
                        "seed_edge_nodes": seed_nodes.tolist(), "seed_edge_retained": retained,
                        "maximum_endpoint_matching_error": float(distances.max()),
                        "analytic_field_lower_bound_in_gate_coordinate_box": minimum_field_in_coordinate_box,
                        "fixed_edge_exceeds_gate": retained and minimum_field_in_coordinate_box > 0})
            assert examples and any(row["fixed_edge_exceeds_gate"] for row in examples)
            workflow.save(args.run / "result.json", {"case": "schwarz-solid", "examples": examples,
                "scope": "Unchanged VTK cap-edge geometry obstruction; no stiffness or local stress certification"})
            assert all(workflow.digest(Path(path)) == value for path, value in inputs.items())
            job.result("retained-curves", fixed_obstructions=sum(
                row["fixed_edge_exceeds_gate"] for row in examples), probes=len(examples))


if __name__ == "__main__":
    main()
