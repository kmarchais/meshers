"""Geometry-only damped projection control; never applies or admits a mesh."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as workflow
import numpy as np


def project(function, points, gradient, maximum_step):
    origin = points.copy()
    locked = workflow.on_box_plane(abs(origin), .5)
    result = origin.copy()
    for _ in range(60):
        previous = result.copy()
        residual = function(*result.T)
        derivative = gradient(function, result)
        derivative[locked] = 0
        step = residual[:, None] * derivative / np.maximum((derivative**2).sum(1), 1e-30)[:, None]
        length = np.linalg.norm(step, axis=1)
        step *= np.minimum(1., maximum_step / np.maximum(length, 1e-30))[:, None]
        step[length < 1e-12] = 0
        accepted = length < 1e-12
        for _ in range(20):
            trial = result - step
            pending = ~accepted & (abs(trial).max(axis=1) <= .5)
            good = np.zeros(len(result), dtype=bool)
            good[pending] = abs(function(*trial[pending].T)) < abs(residual[pending])
            result[good] = trial[good]
            accepted |= good
            step[~accepted] *= .5
            if accepted.all():
                break
        if np.linalg.norm(result - previous, axis=1).max() < 1e-12:
            break
    residual_distance = abs(function(*result.T)) / np.maximum(np.linalg.norm(gradient(function, result), axis=1), 1e-30)
    np.testing.assert_array_equal(result[locked], origin[locked])
    return result, residual_distance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.source, args.output = args.source.resolve(), args.output.resolve()
    assert not args.output.exists()
    _, _, _, _, pv, scorecard = workflow.imports()
    from periodic import periodic_pairs
    paths = [Path(__file__), Path(workflow.__file__), Path(scorecard.__file__),
        args.source / "mesh.npz", args.source / "result.json", args.source / "spinodoid-samples.npz",
        workflow.HERE / "runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json"]
    inputs = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    plane = lambda x, y, z: x - .1
    test = np.array([[.2, .5, .1], [-.1, -.5, -.1]])
    projected, error = project(plane, test, scorecard.gradient, 1 / 144)
    np.testing.assert_allclose(projected[:, 0], .1, atol=1e-10)
    assert error.max() < 1e-10
    mesh = np.load(args.source / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    _, classes, _ = periodic_pairs(points, surface)
    triangles = points[surface]
    cap = np.any(np.all(workflow.on_box_plane(triangles, -.5), axis=1)
        | np.all(workflow.on_box_plane(triangles, .5), axis=1), axis=1)
    wall_nodes = np.unique(surface[~cap])
    _, representative = np.unique(classes[wall_nodes], return_index=True)
    roots = wall_nodes[representative]
    canonical = points[roots].copy()
    canonical[workflow.on_box_plane(canonical, .5)] = -.5
    case, _, _ = workflow.specimen("spinodoid")
    projected, residual = project(case.function, canonical, scorecard.gradient, 1 / 144)
    moved = points.copy()
    # Each representative's displacement is shared by its exact periodic class.
    order = np.argsort(classes[roots])
    group = np.searchsorted(classes[roots][order], classes[wall_nodes])
    moved[wall_nodes] += (projected - canonical)[order[group]]
    _, _, periodic = periodic_pairs(moved, surface)
    packed = np.column_stack((np.full(len(surface), 3), surface)).ravel()
    samples = np.load(args.source / "spinodoid-samples.npz")
    queries = np.concatenate([samples["walls"], samples["caps"]])
    _, feet = pv.PolyData(moved, packed).find_closest_cell(queries, return_closest_point=True)
    distances = np.linalg.norm(queries - feet, axis=1)
    determinant = np.linalg.det(moved[cells[:, 1:]] - moved[cells[:, :1]])
    workflow.save(args.output, {"source": str(args.source), "case": "spinodoid",
        "wall_nodes": len(wall_nodes), "periodic_representatives": len(roots),
        "maximum_step_per_iteration": 1 / 144, "maximum_iterations": 60, "maximum_backtracks": 20,
        "known_plane_control_passed": True, "unresolved_representatives": int((residual > 1e-9).sum()),
        "maximum_field_residual_distance": float(residual.max()),
        "maximum_displacement": float(np.linalg.norm(moved - points, axis=1).max()),
        "reverse_maximum": float(distances.max()), "reverse_samples_above_gate": int((distances > .002).sum()),
        "minimum_corner_determinant": float(determinant.min()), "nonpositive_tets": int((determinant <= 0).sum()),
        "signed_volume_relative_change": float(determinant.sum() / np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]]).sum() - 1),
        "periodic": periodic, "source_mesh_unchanged": True, "qualified": False,
        "mmg_calls": [], "stiffness_solve_performed": False,
        "scope": "Hypothetical boundary-only projection with fixed connectivity; no volume safeguard or applied mesh; reverse distance alone cannot admit this trial."})
    assert all(workflow.digest(Path(path)) == expected for path, expected in inputs.items())
    print(f"Projected geometry only: reverse {distances.max():.9g}, nonpositive tets {(determinant <= 0).sum()}, unresolved {(residual > 1e-9).sum()}")


if __name__ == "__main__":
    main()
