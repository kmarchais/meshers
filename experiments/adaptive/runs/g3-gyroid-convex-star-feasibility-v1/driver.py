"""Convex fixed-neighbor quality feasibility and checked tangent bounds."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import warnings

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np
from scipy.optimize import linprog, minimize

import microgen_workflow as workflow
from job_parts import worker_part
from mmg_quality_localization import quality

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

SOURCE = workflow.HERE / "runs/g3-gyroid-firstpass-fem43-interior-relocation-v1"
BOUND_MARGIN = 1e-8


class StarConstraints:
    def __init__(self, points, cells, node, target):
        self.cells = cells[(cells == node).any(axis=1)]
        assert len(self.cells) and target > 0
        self.positions = np.argmax(self.cells == node, axis=1)
        self.origin = points[node].copy()
        self.scale = float(np.linalg.norm(points[np.unique(self.cells)] - self.origin, axis=1).max())
        assert self.scale > 0
        self.coordinates = (points[self.cells] - self.origin) / self.scale
        self.squared = sum(np.sum((self.coordinates[:, first] - self.coordinates[:, second])**2, axis=1)
            for first in range(4) for second in range(first))
        self.coefficient = target / (12 * np.sqrt(3))
        self.lower = (-.5 - self.origin) / self.scale
        self.upper = (.5 - self.origin) / self.scale

    def evaluate(self, position):
        coordinates = self.coordinates.copy()
        coordinates[np.arange(len(self.cells)), self.positions] = position
        edges = coordinates[:, 1:] - coordinates[:, :1]
        determinants = np.linalg.det(edges)
        determinant_gradients = np.zeros_like(coordinates)
        determinant_gradients[:, 1] = np.cross(edges[:, 1], edges[:, 2])
        determinant_gradients[:, 2] = np.cross(edges[:, 2], edges[:, 0])
        determinant_gradients[:, 3] = np.cross(edges[:, 0], edges[:, 1])
        determinant_gradients[:, 0] = -determinant_gradients[:, 1:].sum(axis=1)
        determinant_gradient = determinant_gradients[np.arange(len(self.cells)), self.positions]
        squared = sum(np.sum((coordinates[:, first] - coordinates[:, second])**2, axis=1)
            for first in range(4) for second in range(first))
        squared_gradient = 8 * position - 2 * coordinates.sum(axis=1)
        values = determinants / self.squared**1.5 - self.coefficient * (squared / self.squared)**1.5
        gradients = determinant_gradient / self.squared[:, None]**1.5
        gradients -= 1.5 * self.coefficient * np.sqrt(squared / self.squared)[:, None] * squared_gradient / self.squared[:, None]
        return values, gradients


def checked_upper_bound(model, position, weights):
    weights = np.asarray(weights, dtype=float)
    assert weights.shape == (len(model.cells),) and np.isfinite(weights).all()
    assert weights.min() >= 0 and abs(weights.sum() - 1) < 1e-12
    values, gradients = model.evaluate(position)
    gradient = weights @ gradients
    support = np.maximum(gradient * model.lower, gradient * model.upper).sum()
    return float(weights @ values - gradient @ position + support)


def diagnose(points, cells, node, target=.1):
    model = StarConstraints(points, cells, node, target)
    initial, _ = model.evaluate(np.zeros(3))
    def constraints(variables):
        return model.evaluate(variables[:3])[0] + variables[3]
    def derivatives(variables):
        return np.column_stack((model.evaluate(variables[:3])[1], np.ones(len(model.cells))))
    start = np.r_[np.zeros(3), max(0, -initial.min()) + 1e-6]
    result = minimize(lambda variables: variables[3], start, jac=lambda variables: np.array([0., 0., 0., 1.]),
        method="SLSQP", bounds=[*zip(model.lower, model.upper, strict=True), (0, None)],
        constraints={"type": "ineq", "fun": constraints, "jac": derivatives},
        options={"maxiter": 1000, "ftol": 1e-13})
    position = result.x[:3]
    values, gradients = model.evaluate(position)
    intercepts = values - gradients @ position
    count = len(model.cells)
    inequalities = np.zeros((6, count + 3))
    for axis in range(3):
        inequalities[2 * axis, :count] = gradients[:, axis] * model.lower[axis]
        inequalities[2 * axis + 1, :count] = gradients[:, axis] * model.upper[axis]
        inequalities[2 * axis:2 * axis + 2, count + axis] = -1
    with warnings.catch_warnings(record=True) as caught:
        linear = linprog(np.r_[intercepts, np.ones(3)], A_ub=inequalities, b_ub=np.zeros(6),
            A_eq=np.r_[np.ones(count), np.zeros(3)][None], b_eq=[1.],
            bounds=[(0, None)] * count + [(None, None)] * 3, method="highs-ds", options={"threads": 1})
    assert linear.success, linear.message
    weights = np.maximum(linear.x[:count], 0)
    weights /= weights.sum()
    upper = checked_upper_bound(model, position, weights)
    moved = points.copy()
    moved[node] = model.origin + model.scale * position
    signed = np.linalg.det(moved[model.cells[:, 1:]] - moved[model.cells[:, :1]])
    candidate_quality = float(quality(moved, model.cells).min())
    attained = bool(signed.min() > 0 and candidate_quality >= target)
    obstructed = upper < -BOUND_MARGIN
    assert not (attained and obstructed)
    return {"node": int(node), "incident_tets": count, "quality_target": target,
        "initial_minimum_quality": float(quality(points, model.cells).min()),
        "candidate_minimum_quality": candidate_quality, "candidate_minimum_signed_determinant": float(signed.min()),
        "candidate_coordinates": moved[node].tolist(), "positive_candidate_at_target": attained,
        "slack_optimizer_success": bool(result.success), "slack_optimizer_message": str(result.message),
        "slack_iterations": int(result.nit), "optimal_slack_reported": float(result.x[3]),
        "minimum_slack_constraint": float(constraints(result.x).min()),
        "checked_tangent_upper_bound": upper, "required_negative_margin": BOUND_MARGIN,
        "fixed_neighbor_obstruction": obstructed, "weights": weights.tolist(),
        "normalized_expansion_point": position.tolist(), "normalized_box_lower": model.lower.tolist(),
        "normalized_box_upper": model.upper.tolist(), "physical_origin": model.origin.tolist(), "scale": model.scale,
        "linear_solver_success": bool(linear.success), "linear_solver_warnings": [str(item.message) for item in caught],
        "scope": "Floating-point checked concave-tangent bound for this single moving vertex and fixed neighbors inside the comparison box; not interval certification, coupled-relocation impossibility, or continuum accuracy"}


def worker(run):
    source = json.loads((SOURCE / "result.json").read_text())
    assert not source["qualified"] and source["repair"]["remaining_bad_elements"] == 49
    mesh = np.load(SOURCE / "mesh.npz")
    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    assert len(cells) <= workflow.CAP
    measured = quality(points, cells)
    worst = int(measured.argmin())
    boundary = np.zeros(len(points), dtype=bool)
    boundary[np.unique(surface)] = True
    free = cells[worst][~boundary[cells[worst]]]
    assert len(free) == 1
    report = diagnose(points, cells, int(free[0]))
    model = StarConstraints(points, cells, int(free[0]), .1)
    independently_checked = checked_upper_bound(model, np.asarray(report["normalized_expansion_point"]), report["weights"])
    assert abs(independently_checked - report["checked_tangent_upper_bound"]) < 1e-14
    workflow.save(run / "result.json", {"case": "gyroid", "source": str(SOURCE), "worst_element": worst,
        "diagnostic": report, "global_geometry_unchanged": True, "source_mesh_unchanged": True,
        "qualified": False, "stiffness_solve_performed": False,
        "scope": "One frozen failed star, no vertex moves applied to source mesh, no new remeshing or stiffness workflow"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), SOURCE / "result.json", SOURCE / "mesh.npz"]
    paths.extend(workflow.HERE / name for name in ("test_star_feasibility.py", "job_parts.py", "test_job_parts.py",
        "microgen_workflow.py", "mmg_quality_localization.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "highs_threads": 1, "worker_timeout_s": 120, "quality_target": .1,
        "negative_margin": BOUND_MARGIN, "one_star_only": True,
        "packages": {name: importlib.metadata.version(name) for name in ("numpy", "scipy")}})
    with Job(parts=["controls", "frozen-star"]) as job:
        for stage in ("controls", "frozen-star"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_star_feasibility.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker"])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
