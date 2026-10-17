"""Recheck saved G3 artifacts without changing their reports or running FEA."""

import argparse
import json
import os
from pathlib import Path
import resource
import sys

ROOT = Path(__file__).resolve().parents[2]
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))

from pueue_progress import Job
import numpy as np
from scipy.linalg import eigvalsh

import microgen_workflow as experiment


def volume(points, tetrahedra):
    coordinates = points[tetrahedra]
    if tetrahedra.shape[1] == 4:
        determinants = np.linalg.det(coordinates[:, 1:] - coordinates[:, :1])
        assert np.all(determinants > 0)
        return float(determinants.sum() / 6)
    gradient = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    barycentric = np.vstack([[0.25]*4, np.full((4, 4), 1/6) + np.eye(4) / 3])
    weights = [-2/15] + [3/40]*4
    result = 0.0
    for values, weight in zip(barycentric, weights, strict=True):
        derivatives = np.vstack([(4 * values - 1)[:, None] * gradient,
                                 [4 * (values[first] * gradient[second] + values[second] * gradient[first])
                                  for first, second in ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))]])
        determinants = np.linalg.det(np.einsum("na,tnb->tab", derivatives, coordinates))
        assert np.all(determinants > 0)
        result += weight * determinants.sum()
    return float(result)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (experiment.MEMORY, experiment.MEMORY))
    experiment.imports()
    run = args.run.resolve()
    summary = json.loads((run / "summary.json").read_text())
    manifest = json.loads((run / "manifest.json").read_text())
    original_hashes = {path.name: experiment.digest(path) for path in run.iterdir() if path.is_file()}
    records = []
    labels = [f"{row['case']}-{row['route']}-{row['cells']}" for row in summary["rows"]]
    with Job(parts=["provenance"] + labels) as job:
        with job.part("provenance") as part:
            part.update(0, "Checking frozen source and inherited hashes", artifacts=len(original_hashes))
            assert experiment.digest(run / "driver.py") == manifest["script_sha256"]
            assert experiment.digest(Path(experiment.__file__)) == manifest["script_sha256"]
            for name, provenance in json.loads((run / "inherited.json").read_text()).items():
                assert experiment.digest(Path(provenance["source"])) == provenance["sha256"]
                assert experiment.digest(run / name) == provenance["sha256"]
            for path in (run / "microgen-source").iterdir():
                assert experiment.digest(path) == manifest["source_sha256"][path.name]
        for label, row in zip(labels, summary["rows"], strict=True):
            with job.part(label) as part:
                part.update(0, "Rechecking " + label, nodes=row.get("nodes", 0))
                record = {"name": label, "status": row["status"]}
                mesh_path = run / (label + ".npz")
                if mesh_path.exists():
                    data = np.load(mesh_path)
                    points, tetrahedra = data["points"], data["tetrahedra"]
                    assert len(points) == row["nodes"] and len(tetrahedra) == row["elements"]
                    independent_volume = volume(points, tetrahedra)
                    difference = abs(independent_volume - row["geometry"]["volume"])
                    assert difference < 1e-11, (label, difference)
                    record["independent_volume_difference"] = difference
                    if row.get("qualified"):
                        from periodic import periodic_pairs
                        periodic_pairs(points + 0.5, data["surface"][:, :3])
                    if row["route"] == "meshers-curved":
                        _, _, _, _, _, scorecard = experiment.imports()
                        samples = np.load(run / f"{row['case']}-samples.npz")
                        queries = np.concatenate([samples["walls"], samples["caps"]])
                        part.update(0.5, "Checking finer quadratic boundary tessellation", nodes=len(points))
                        vertices, faces = experiment.tessellate(points, data["surface"], 16)
                        reverse = scorecard.distance_to_mesh(vertices, faces, queries)
                        maximum = float(reverse.max())
                        record["reverse_max_subdivision16"] = maximum
                        record["forward_and_reverse16_max"] = max(maximum, row["geometry"]["sampled_forward_max"])
                        assert record["forward_and_reverse16_max"] < 0.002
                        del vertices, faces, reverse
                if "stiffness" in row:
                    stiffness = np.array(row["stiffness"])
                    assert len(row["solver_records"]) == 6
                    assert max(solve["relative_residual"] for solve in row["solver_records"]) < 2e-9
                    assert eigvalsh(stiffness).min() > 0
                    if row["case"] in summary["references"]:
                        reference = np.array(summary["references"][row["case"]]["stiffness"])
                        spectrum = eigvalsh(stiffness - reference, reference)
                        stored = row["mixed_load_reference_disagreement"]
                        assert np.allclose(spectrum, stored["relative_energy_difference_eigenvalues"], atol=1e-13, rtol=1e-12)
                        record["mixed_load_reference_disagreement"] = float(np.max(abs(spectrum)))
                records.append(record)
                job.result(label, mesh_checked=int(mesh_path.exists()), tensor_checked=int("stiffness" in row))
    assert original_hashes == {path.name: experiment.digest(path) for path in run.iterdir() if path.is_file()}
    experiment.save(args.output, {"run": str(run), "summary_sha256": experiment.digest(run / "summary.json"),
                                 "artifacts_unchanged": True, "records": records,
                                 "scope": "New G3 artifacts only. Historical saved fields were not transferred and were not revalidated."})


if __name__ == "__main__":
    main()
