"""Replay failed forward samples and test feasible box-constrained target feet."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from scipy.optimize import minimize

import microgen_workflow as workflow


def constrained_foot(function, probe, level, locked, initial):
    scale = max(float(np.linalg.norm(initial-probe)), 1e-3)
    lower, upper = np.full(3, -.5), np.full(3, .5)
    lower[locked] = upper[locked] = np.sign(probe[locked])*.5
    start = np.clip(initial, lower, upper)
    def residual(displacement):
        point = probe+scale*displacement
        return float(function(*point[None].T)[0]-level)/scale
    result = minimize(lambda displacement: float(displacement @ displacement), (start-probe)/scale,
        method="SLSQP", bounds=list(zip((lower-probe)/scale, (upper-probe)/scale)),
        constraints={"type": "eq", "fun": residual}, options={"maxiter": 300, "ftol": 1e-13})
    foot = probe+scale*result.x
    remaining = abs(float(function(*foot[None].T)[0]-level))
    feasible = bool(remaining < 1e-9 and abs(foot).max() <= .5+1e-12)
    return foot, dict(solver_success=bool(result.success), message=str(result.message), iterations=int(result.nit),
        absolute_field_residual=remaining, feasible_target=feasible, distance=float(np.linalg.norm(foot-probe)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    _, _, _, _, _, scorecard = workflow.imports()
    source = json.loads((args.source / "result.json").read_text())
    paths = [Path(__file__), args.source / "mesh.npz", args.source / "result.json", Path(workflow.__file__), Path(scorecard.__file__),
        workflow.HERE / "test_projection_failure_probe.py"]
    paths.extend(args.source / (source["case"]+suffix) for suffix in ("-samples.npz", "-specimen.json"))
    inputs = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    shutil.copyfile(workflow.HERE / "test_projection_failure_probe.py", args.output.parent / "test_projection_failure_probe.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_projection_failure_probe.py"], check=True, timeout=30)
    mesh = np.load(args.source / "mesh.npz")
    points, surface = mesh["points"], mesh["surface"]
    assert surface.shape[1] == 6
    triangle = points[surface[:, :3]]
    caps = np.any(np.all(workflow.on_box_plane(triangle, -.5), axis=1) | np.all(workflow.on_box_plane(triangle, .5), axis=1), axis=1)
    wall_indices = np.flatnonzero(~caps)
    weights = workflow.tri_values(scorecard.BARY)
    probes = np.einsum("sn,fnd->fsd", weights, points[surface[~caps]]).reshape(-1, 3)
    locked = workflow.on_box_plane(abs(probes), .5)
    case, _, _ = workflow.specimen(source["case"])
    levels = np.array([value for value in case.levels if np.isfinite(value)])
    values = case.function(*probes.T)
    targets = levels[np.argmin(abs(values[:, None]-levels), axis=1)]
    foot = probes.copy()
    for _ in range(40):
        gradient = scorecard.gradient(case.function, foot)
        gradient[locked] = 0
        residual = case.function(*foot.T)-targets
        step = residual[:, None]*gradient/np.maximum((gradient**2).sum(axis=1), 1e-30)[:, None]
        foot -= step
        if np.linalg.norm(step, axis=1).max() < 1e-12:
            break
    residual = abs(case.function(*foot.T)-targets)/np.maximum(np.linalg.norm(scorecard.gradient(case.function, foot), axis=1), 1e-30)
    failed = (residual > 1e-9) | np.any(abs(foot) > .5+1e-9, axis=1)
    assert int(failed.sum()) == source["geometry"]["projection_failed"]
    rows = []
    for index in np.flatnonzero(failed):
        wall, sample = divmod(int(index), len(weights))
        candidate, control = constrained_foot(case.function, probes[index], targets[index], locked[index], foot[index])
        rows.append(dict(face=int(wall_indices[wall]), probe_index=sample, probe=probes[index].tolist(),
            original_foot=foot[index].tolist(), original_distance=float(np.linalg.norm(foot[index]-probes[index])),
            original_normalized_residual=float(residual[index]), original_box_overshoot=float(max(abs(foot[index]).max()-.5, 0)),
            locked=locked[index].tolist(), target_level=float(targets[index]),
            constrained_foot=candidate.tolist(), constrained_control=control))
    workflow.save(args.output, dict(projection_failure_localization=True, source=str(args.source), case=source["case"],
        failed_projection_count=int(failed.sum()), forward_sample_count=len(probes), witnesses=rows,
        all_constrained_candidates_feasible=all(row["constrained_control"]["feasible_target"] for row in rows),
        mesh_modified=False, samples_regenerated=False, source_rejection_preserved=True,
        qualified=False, stiffness_solve_performed=False,
        scope="Same failed forward probes, original projector replay and box-constrained feasible-target trials. Feasible feet are distance upper bounds, not certified closest points. The existing geometry rejection is not waived or rewritten."))
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
