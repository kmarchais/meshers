"""Test bounded periodic volume offsets without changing native meshing code."""

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


def worker(source, run):
    np, *_ = workflow.imports()
    from periodic import periodic_pairs
    from scipy.optimize import brentq

    row = json.loads((source / "result.json").read_text())
    assert row["case"] == "smooth-lattice" and row["qualified"]
    archive = np.load(source / "mesh.npz")
    points, tetrahedra, surface = (archive[name] for name in ("points", "tetrahedra", "surface"))
    assert tetrahedra.shape[1] == 4 and len(tetrahedra) <= workflow.CAP
    _, classes, _ = periodic_pairs(points, surface)
    triangles = points[surface]
    cap = np.zeros(len(surface), dtype=bool)
    for axis in range(3):
        for value in (-0.5, 0.5):
            cap |= np.all(abs(triangles[:, :, axis] - value) < 1e-10, axis=1)
    wall = np.unique(surface[~cap])
    case, _, _ = workflow.specimen(row["case"])
    locations = points[wall]
    gradient = np.empty_like(locations)
    step = 1e-6
    for axis in range(3):
        shift = np.eye(3)[axis] * step
        gradient[:, axis] = (case.function(*(locations + shift).T) -
            case.function(*(locations - shift).T)) / (2 * step)
    length = np.linalg.norm(gradient, axis=1)
    assert (length > 1e-10).all()
    direction = np.zeros_like(points)
    direction[wall] = gradient / length[:, None]
    locked = abs(abs(points) - 0.5) < 1e-10
    direction[locked] = 0
    grouped = np.zeros_like(direction)
    np.add.at(grouped, classes, direction)
    counts = np.bincount(classes, minlength=len(points))
    direction = grouped[classes] / counts[classes, None]
    assert abs(direction[locked]).max(initial=0) == 0

    def volume(offset):
        moved = points + offset * direction
        return float(np.linalg.det(moved[tetrahedra[:, 1:]] - moved[tetrahedra[:, :1]]).sum() / 6)

    baseline = volume(0)
    assert abs(baseline / row["geometry"]["volume"] - 1) < 1e-10
    reference = json.loads((run / "smooth-lattice-specimen.json").read_text())
    target = reference["volume"]
    bracket = (-0.002, 0.002)
    assert volume(bracket[0]) < baseline < volume(bracket[1])
    candidates = []
    for fraction in (0.5, 1.0):
        requested = baseline + fraction * (target - baseline)
        offset = brentq(lambda value: volume(value) - requested, *bracket, xtol=1e-14)
        moved = points + offset * direction
        assert np.array_equal(moved[locked], points[locked])
        assert np.max(abs(moved)) <= 0.5 + 1e-12
        _, _, periodic = periodic_pairs(moved, surface)
        metrics = workflow.geometry(run, row["case"], moved, tetrahedra, surface)
        assert abs(metrics["volume"] / requested - 1) < 1e-10
        label = "half" if fraction == 0.5 else "full"
        np.savez_compressed(run / (label + "-mesh.npz"), points=moved, tetrahedra=tetrahedra, surface=surface)
        candidates.append({"label": label, "fraction": fraction, "normal_offset": offset,
            "maximum_displacement": float(np.linalg.norm(moved - points, axis=1).max()),
            "geometry": metrics, "periodic": periodic})
    workflow.save(run / "result.json", {"source": str(source), "baseline_volume": baseline,
        "volume_target": target, "target_scramble_spread": reference["volume_scramble_spread"],
        "candidates": candidates,
        "scope": "Saved Tet4 geometry perturbation only; sampled volume target is uncertain; no stiffness or live workflow claim"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.source, args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    original = workflow.HERE / "runs/g3-microgen-defaults-v4"
    files = [args.source / name for name in ("result.json", "mesh.npz")]
    for name in ("smooth-lattice-specimen.json", "smooth-lattice-samples.npz"):
        files.append(original / name)
        shutil.copyfile(original / name, args.run / name)
    files += [Path(__file__), Path(workflow.__file__)]
    hashes = {str(path): workflow.digest(path) for path in files}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 180, "offset_bracket": [-0.002, 0.002],
        "volume_correction_fractions": [0.5, 1.0]})
    with Job(parts=["volume-offsets"]) as job:
        with job.part("volume-offsets") as part:
            record = workflow.bounded([sys.executable, __file__, "--worker", "--source", str(args.source),
                "--run", str(args.run)], args.run / "worker.log", part, timeout=180)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("volume-offsets", seconds=record["wall_s"],
                admitted=sum(candidate["geometry"]["geometry_admitted"] for candidate in row["candidates"]),
                peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
