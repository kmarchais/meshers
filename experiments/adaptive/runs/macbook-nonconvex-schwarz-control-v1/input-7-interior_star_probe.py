"""One fixed-boundary interior-star feasibility check on genuine MMG output."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np
from scipy.optimize import minimize

import microgen_workflow as workflow
from cap_quality_bound import bounds, self_check
from job_parts import worker_part
from mmg_quality_localization import quality

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

SOURCE = workflow.HERE / "runs/g3-gyroid-tuned43-v3"


def relocate_star(points, cells, node, target=.1):
    selected = np.flatnonzero((cells == node).any(axis=1))
    assert len(selected) and target > 0
    star = cells[selected]
    positions = np.argmax(star == node, axis=1)
    fixed = np.array([cell[np.arange(4) != position] for cell, position in zip(star, positions, strict=True)])
    ceiling, _ = bounds(points[fixed])
    initial = quality(points, star)
    constraint_target = target + 1e-9
    row = {"node": int(node), "incident_tets": len(star), "quality_target": target,
        "initial_minimum_quality": float(initial.min()), "fixed_opposite_face_upper_bound": float(ceiling.min()),
        "target_attained": False, "status": "fixed-face-obstructed", "numerical_quality_guard": 1e-9}
    if ceiling.min() < target:
        return points[node].copy(), row
    origin = points[node].copy()
    scale = float(np.linalg.norm(points[np.unique(star)] - origin, axis=1).max())
    assert scale > 0
    coordinates = (points[star] - origin) / scale
    initial_squared = sum(np.sum((coordinates[:, first] - coordinates[:, second])**2, axis=1)
        for first in range(4) for second in range(first))
    def constraints(displacement):
        moved = coordinates.copy()
        moved[np.arange(len(star)), positions] = displacement
        determinants = np.linalg.det(moved[:, 1:] - moved[:, :1])
        squared = sum(np.sum((moved[:, first] - moved[:, second])**2, axis=1)
            for first in range(4) for second in range(first))
        return determinants / initial_squared**1.5 - constraint_target / (12 * np.sqrt(3)) * (squared / initial_squared)**1.5
    started = time.monotonic()
    result = minimize(lambda displacement: np.dot(displacement, displacement), np.zeros(3), method="SLSQP",
        constraints={"type": "ineq", "fun": constraints}, options={"maxiter": 300, "ftol": 1e-13})
    candidate = origin + scale * result.x
    moved = points.copy()
    moved[node] = candidate
    final = quality(moved, star)
    determinants = np.linalg.det(moved[star[:, 1:]] - moved[star[:, :1]])
    attained = bool(determinants.min() > 0 and final.min() >= target)
    row.update(status="feasible-numerically" if attained else "not-found-not-impossibility-proof",
        target_attained=attained, solver_success=bool(result.success), solver_message=str(result.message),
        iterations=int(result.nit), seconds=time.monotonic() - started,
        final_minimum_quality=float(final.min()), maximum_move=float(np.linalg.norm(candidate - origin)),
        minimum_normalized_constraint=float(constraints(result.x).min()),
        candidate_coordinates=candidate.tolist(), minimum_signed_determinant=float(determinants.min()))
    return candidate, row


def worker(run):
    self_check()
    row = json.loads((SOURCE / "result.json").read_text())
    assert row["case"] == "gyroid" and row["route"].startswith("microgen-vtk") and len(row["mmg_calls"]) == 2
    assert row["geometry"]["sampled_two_sided_max"] <= .002 and abs(row["geometry"]["relative_volume_error"]) <= .01
    mesh = np.load(SOURCE / "mesh.npz")
    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    assert len(cells) <= workflow.CAP and cells.shape[1] == 4
    measured = quality(points, cells)
    assert abs(measured.min() - row["geometry"]["q_min_corner"]) < 1e-12
    boundary = np.zeros(len(points), dtype=bool)
    boundary[np.unique(surface)] = True
    worst = int(measured.argmin())
    free = cells[worst][~boundary[cells[worst]]]
    assert len(free) == 1, "Reconsider the diagnostic if the recorded one-free-vertex configuration changed"
    node = int(free[0])
    candidate, result = relocate_star(points, cells, node)
    moved = points.copy()
    moved[node] = candidate
    np.testing.assert_array_equal(moved[boundary], points[boundary])
    old_determinants = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])
    new_determinants = np.linalg.det(moved[cells[:, 1:]] - moved[cells[:, :1]])
    signed_volume_difference = float((new_determinants.sum() - old_determinants.sum()) / 6)
    assert abs(signed_volume_difference) < 1e-12
    workflow.save(run / "result.json", {"case": "gyroid", "source": str(SOURCE), "worst_element": worst,
        "nodes": len(points), "elements": len(cells), "star": result,
        "all_boundary_coordinates_unchanged": True, "connectivity_unchanged": True,
        "signed_volume_difference": signed_volume_difference,
        "minimum_signed_determinant_after_trial": float(new_determinants.min()),
        "global_minimum_quality_after_trial": float(quality(moved, cells).min()),
        "qualified": False, "stiffness_solve_performed": False,
        "scope": "One interior-vertex feasibility diagnostic on frozen actual MMG output; not a repaired workflow, Microgen default, Meshers seed, or admitted comparator"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), SOURCE / "result.json", SOURCE / "mesh.npz"]
    paths.extend(workflow.HERE / name for name in ("job_parts.py", "test_job_parts.py", "test_interior_star_probe.py",
        "microgen_workflow.py", "cap_quality_bound.py", "mmg_quality_localization.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "maximum_vertices_moved": 1, "quality_target": .1,
        "one_slsqp_attempt": True, "surface_and_specimen_unchanged": True})
    with Job(parts=["controls", "star-feasibility"]) as job:
        for stage in ("controls", "star-feasibility"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_job_parts.py", "experiments/adaptive/test_interior_star_probe.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker"])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                values = {"seconds": record["wall_s"], "peak_gib": record["peak_rss_bytes"] / 1024**3}
                if stage == "star-feasibility":
                    result = json.loads((args.run / "result.json").read_text())["star"]
                    values.update(attained=int(result["target_attained"]), quality=result.get("final_minimum_quality"),
                        upper_bound=result["fixed_opposite_face_upper_bound"])
                job.result(stage, **values)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
