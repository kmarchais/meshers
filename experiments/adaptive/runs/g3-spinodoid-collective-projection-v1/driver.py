"""Test simultaneous wall-root moves on a saved rejected linear mesh."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    np, *_, scorecard = workflow.imports()
    from periodic import periodic_pairs

    archive = np.load(args.source / "mesh.npz")
    points, tets, surface = (archive[key] for key in ("points", "tetrahedra", "surface"))
    row = json.loads((args.source / "result.json").read_text())
    assert row["case"] == "spinodoid" and tets.shape[1] == 4 and len(tets) <= workflow.CAP
    _, classes, _ = periodic_pairs(points, surface)
    triangles = points[surface]
    cap = np.any(np.all(workflow.on_box_plane(triangles, -.5), axis=1) |
        np.all(workflow.on_box_plane(triangles, .5), axis=1), axis=1)
    wall = np.unique(surface[~cap])
    groups, index = np.unique(classes[wall], return_index=True)
    representatives = wall[index]
    origins = points[representatives]
    frozen = json.loads((args.source / "fourier-coefficients.json").read_text())
    coefficients = np.array(frozen["coefficients"])
    modes, amplitude, phase = coefficients[:, :3], coefficients[:, 3], coefficients[:, 4]

    def evaluate(coordinates):
        angle = 2 * np.pi * ((coordinates + .5) @ modes.T) + phase
        return frozen["threshold"] - np.cos(angle) @ amplitude

    def gradient(coordinates):
        angle = 2 * np.pi * ((coordinates + .5) @ modes.T) + phase
        return 2 * np.pi * (np.sin(angle) * amplitude) @ modes

    case, _, _ = workflow.specimen("spinodoid")
    queries = np.random.default_rng(20261018).uniform(-.5, .5, (4096, 3))
    np.testing.assert_allclose(evaluate(queries), case.function(*queries.T), atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(gradient(queries[:32]), scorecard.gradient(case.function, queries[:32]),
        atol=2e-7, rtol=2e-7)
    target = origins.copy()
    locked = workflow.on_box_plane(abs(origins), .5)
    for _ in range(40):
        derivative = gradient(target)
        derivative[locked] = 0
        step = evaluate(target)[:, None] * derivative / np.maximum((derivative**2).sum(1), 1e-30)[:, None]
        target -= step
        if np.linalg.norm(step, axis=1).max() < 1e-12:
            break
    residual = abs(evaluate(target)) / np.maximum(np.linalg.norm(gradient(target), axis=1), 1e-30)
    valid = (residual < 1e-9) & np.all(abs(target) <= .5 + 1e-12, axis=1)
    np.testing.assert_array_equal(target[locked], origins[locked])
    displacement = np.zeros_like(points)
    displacement[groups[valid]] = target[valid] - origins[valid]
    displacement = displacement[classes]
    candidates = []
    for fraction in (1.0, .5, .25):
        moved = points + fraction * displacement
        np.testing.assert_array_equal(moved[workflow.on_box_plane(abs(points), .5)],
            points[workflow.on_box_plane(abs(points), .5)])
        _, _, periodic = periodic_pairs(moved, surface)
        metrics = workflow.geometry(args.run, "spinodoid", moved, tets, surface)
        label = str(fraction).replace(".", "p")
        np.savez_compressed(args.run / (label + "-mesh.npz"), points=moved, tetrahedra=tets, surface=surface)
        candidates.append({"fraction": fraction, "geometry": metrics, "periodic": periodic})
    workflow.save(args.run / "result.json", {"source": str(args.source),
        "wall_classes": len(groups), "failed_root_classes": int((~valid).sum()),
        "maximum_full_move": float(np.linalg.norm(displacement, axis=1).max()),
        "candidates": candidates,
        "scope": "Saved-mesh simultaneous projection diagnostic, not a live workflow or continuum certificate"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [args.source / name for name in ("mesh.npz", "result.json", "fourier-coefficients.json",
        "spinodoid-samples.npz", "spinodoid-specimen.json")]
    paths += [Path(__file__), Path(workflow.__file__)]
    hashes = {str(path): workflow.digest(path) for path in paths}
    for name in ("spinodoid-samples.npz", "spinodoid-specimen.json"):
        shutil.copyfile(args.source / name, args.run / name)
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "timeout_s": 240, "threads": 1, "fractions": [1.0, .5, .25]})
    with Job(parts=["simultaneous-root-moves"]) as job:
        with job.part("simultaneous-root-moves") as part:
            command = [sys.executable, __file__, "--run", str(args.run), "--source", str(args.source), "--worker"]
            record = workflow.bounded(command, args.run / "worker.log", part, timeout=240)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("simultaneous-root-moves", seconds=record["wall_s"],
                admitted=sum(candidate["geometry"]["geometry_admitted"] for candidate in row["candidates"]))
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
