"""One periodic cap-interior vertex feasibility test, without applying moves."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import warnings

import microgen_workflow as workflow
from cap_corner_probe import polygon_corners
from cap_quality_bound import bounds, caps
from star_feasibility import StarConstraints, checked_upper_bound
from mmg_quality_localization import quality

import numpy as np
from scipy.optimize import linprog, minimize


class PlaneStar:
    def __init__(self, points, cells, faces, nodes, axis, target):
        self.axes = [other for other in range(3) if other != axis]
        self.nodes = np.asarray(nodes)
        self.volumes = []
        self.guards = []
        for node in nodes:
            incident = cells[(cells == node).any(axis=1)]
            guard = min(1e-12, float(quality(points, incident).min()) * .01)
            assert guard > 0
            self.volumes.append(StarConstraints(points, incident, node, guard))
            self.guards.append(guard)
        self.scale = max(model.scale for model in self.volumes)
        self.lower = np.max((-.5 - points[self.nodes][:, self.axes]) / self.scale, axis=0)
        self.upper = np.min((.5 - points[self.nodes][:, self.axes]) / self.scale, axis=0)
        self.triangles, opposite = [], []
        for node in nodes:
            selected = faces[(faces == node).any(axis=1)].copy()
            for face in selected:
                vertices = points[face][:, self.axes]
                signed = np.linalg.det(vertices[1:] - vertices[:1])
                assert signed != 0
                if signed < 0:
                    face[1], face[2] = face[2], face[1]
                face = np.roll(face, -int(np.flatnonzero(face == node).item()))
                self.triangles.append(face)
                opposite.append((points[face[1:]][:, self.axes] - points[node, self.axes]) / self.scale)
        self.triangles = np.asarray(self.triangles)
        self.opposite = np.asarray(opposite)
        self.target = target
        self.cells = [*self.triangles, *(cell for model in self.volumes for cell in model.cells)]

    def evaluate(self, position):
        first, second = self.opposite[:, 0], self.opposite[:, 1]
        a, b = first - position, second - position
        determinant = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
        squared = (a * a).sum(1) + (b * b).sum(1) + ((first - second)**2).sum(1)
        derivative = np.column_stack((first[:, 1] - second[:, 1], second[:, 0] - first[:, 0]))
        values = [2 * np.sqrt(3) * determinant - self.target * squared]
        gradients = [2 * np.sqrt(3) * derivative - self.target * (4 * position - 2 * (first + second))]
        for model in self.volumes:
            query = np.zeros(3)
            query[self.axes] = position * self.scale / model.scale
            current, derivative = model.evaluate(query)
            values.append(current)
            gradients.append(derivative[:, self.axes] * self.scale / model.scale)
        return np.concatenate(values), np.vstack(gradients)


def diagnose(points, cells, faces, nodes, axis, target=.100001):
    model = PlaneStar(points, cells, faces, nodes, axis, target)
    initial, _ = model.evaluate(np.zeros(2))
    result = minimize(lambda x: x[2], [0., 0., max(0., -float(initial.min())) + 1e-6],
        jac=lambda x: np.array([0., 0., 1.]), method="SLSQP",
        bounds=[*zip(model.lower, model.upper, strict=True), (0, None)],
        constraints={"type": "ineq", "fun": lambda x: model.evaluate(x[:2])[0] + x[2],
            "jac": lambda x: np.column_stack((model.evaluate(x[:2])[1], np.ones(len(model.cells))))},
        options={"maxiter": 1000, "ftol": 1e-13})
    position = result.x[:2]
    values, gradients = model.evaluate(position)
    intercepts = values - gradients @ position
    count = len(values)
    inequalities = np.zeros((4, count + 2))
    for dimension in range(2):
        inequalities[2 * dimension, :count] = gradients[:, dimension] * model.lower[dimension]
        inequalities[2 * dimension + 1, :count] = gradients[:, dimension] * model.upper[dimension]
        inequalities[2 * dimension:2 * dimension + 2, count + dimension] = -1
    with warnings.catch_warnings(record=True) as caught:
        linear = linprog(np.r_[intercepts, np.ones(2)], A_ub=inequalities, b_ub=np.zeros(4),
            A_eq=np.r_[np.ones(count), np.zeros(2)][None], b_eq=[1.],
            bounds=[(0, None)] * count + [(None, None)] * 2, method="highs-ds",
            options={"threads": 1, "time_limit": 30})
    assert linear.success, linear.message
    weights = np.maximum(linear.x[:count], 0)
    weights /= weights.sum()
    upper = checked_upper_bound(model, position, weights)
    moved = points.copy()
    for node in nodes:
        moved[node, model.axes] += position * model.scale
    triangle_coordinates = moved[model.triangles][:, :, model.axes]
    signed = np.linalg.det(triangle_coordinates[:, 1:] - triangle_coordinates[:, :1])
    cap_quality = float(bounds(moved[model.triangles])[0].min())
    star_cells = np.unique(np.vstack([current.cells for current in model.volumes]), axis=0)
    determinants = np.linalg.det(moved[star_cells[:, 1:]] - moved[star_cells[:, :1]])
    attained = bool(signed.min() > 0 and cap_quality >= target and determinants.min() > 0
        and all(quality(moved, current.cells).min() >= guard
            for current, guard in zip(model.volumes, model.guards, strict=True)))
    assert not (attained and upper < -1e-8)
    return {"nodes": list(map(int, nodes)), "axis": axis, "target_cap_quality": target,
        "volume_seed_guards": model.guards, "candidate_coordinates": moved[model.nodes].tolist(),
        "candidate_minimum_cap_quality": cap_quality, "candidate_minimum_volume_determinant": float(determinants.min()),
        "positive_candidate_at_cap_target": attained, "checked_tangent_upper_bound": upper,
        "fixed_neighbor_obstruction": upper < -1e-8, "weights": weights.tolist(),
        "normalized_expansion_point": position.tolist(), "scale": model.scale,
        "slack_optimizer_success": bool(result.success), "slack_optimizer_message": str(result.message),
        "linear_solver_warnings": [str(item.message) for item in caught],
        "incident_cap_triangles": model.triangles.tolist(), "incident_volume_tets": len(star_cells)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.source, args.output = args.source.resolve(), args.output.resolve()
    assert not args.output.exists()
    paths = [Path(__file__), args.source / "mesh.npz", args.source / "result.json"]
    paths.extend(workflow.HERE / name for name in ("star_feasibility.py", "cap_corner_probe.py",
        "cap_quality_bound.py", "mmg_quality_localization.py", "test_cap_plane_star.py", "microgen_workflow.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.output.parent / "analysis-inputs.json", hashes)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_cap_plane_star.py"],
        check=True, timeout=45)
    workflow.imports()
    from facet_incidence_audit import audit
    from periodic import periodic_pairs

    mesh = np.load(args.source / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert audit(points, cells, surface)["oriented_facet_consistency_passed"]
    _, classes, _ = periodic_pairs(points, surface)
    cap_faces, triangles = caps(mesh)
    worst = int(bounds(triangles)[0].argmin())
    face = cap_faces[worst]
    axis = next(axis for axis in range(3) if np.all(workflow.on_box_plane(abs(points[face, axis]), .5)))
    selected = np.all(workflow.on_box_plane(points[cap_faces, axis], points[face[0], axis]), axis=1)
    axes = [other for other in range(3) if other != axis]
    corners, ambiguous = polygon_corners(points[:, axes], cap_faces[selected])
    boundary = {row["node"] for row in corners} | set(ambiguous)
    free = [int(node) for node in face if node not in boundary]
    assert 1 <= len(free) <= 2
    diagnostics = []
    for node in free:
        nodes = np.flatnonzero(classes == classes[node])
        assert len(nodes) == 2 and np.all(workflow.on_box_plane(abs(points[nodes, axis]), .5))
        diagnostics.append(diagnose(points, cells, cap_faces, nodes, axis))
    workflow.save(args.output, {"source": str(args.source), "case": "spinodoid", "worst_cap_face": face.tolist(),
        "initial_worst_cap_ceiling": float(bounds(triangles)[0][worst]), "diagnostics": diagnostics,
        "source_mesh_unchanged": True, "qualified": False, "mmg_calls": [], "stiffness_solve_performed": False,
        "scope": "One worst cap face and each of at most two periodic free-vertex orbits, fixed neighbors and cap polygon; cap target plus seed orientation guards are not final tet-quality admission. No move is applied."})
    assert all(workflow.digest(Path(path)) == expected for path, expected in hashes.items())
    print("Completed fixed-neighbor cap-plane feasibility; no source edits or FE")


if __name__ == "__main__":
    main()
