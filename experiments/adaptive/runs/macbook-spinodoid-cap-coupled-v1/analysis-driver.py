"""Bounded independent planar motion of two periodic cap-interior orbits."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as workflow
from cap_quality_bound import bounds, caps
from mmg_quality_localization import quality

import numpy as np
from scipy.optimize import minimize


class CoupledPatch:
    def __init__(self, points, cells, faces, orbits, axis, target=.1):
        self.points, self.axis, self.target = points, axis, target
        self.axes = [other for other in range(3) if other != axis]
        self.orbits = [np.asarray(orbit, dtype=int) for orbit in orbits]
        assert len(orbits) == 2 and all(len(orbit) == 2 for orbit in orbits)
        self.nodes = np.concatenate(self.orbits)
        assert len(np.unique(self.nodes)) == 4 and target > 0
        owner = np.full(len(points), -1, dtype=int)
        for group, orbit in enumerate(self.orbits):
            owner[orbit] = group
            assert np.all(workflow.on_box_plane(abs(points[orbit, axis]), .5))
            np.testing.assert_allclose(points[orbit[0], self.axes], points[orbit[1], self.axes], atol=1e-10, rtol=0)
        self.cell_indices = np.flatnonzero((owner[cells] >= 0).any(axis=1))
        self.cells = cells[self.cell_indices].copy()
        self.tet_owner = owner[self.cells]
        self.scale = float(np.linalg.norm(points[self.cells[:, 1:]] - points[self.cells[:, :1]], axis=2).max())
        assert self.scale > 0
        self.tet_base = (points[self.cells] - points[self.cells[:, :1]]) / self.scale
        self.original_squared = sum(((self.tet_base[:, a] - self.tet_base[:, b])**2).sum(1)
            for a in range(4) for b in range(a))
        old_quality = quality(points, self.cells)
        assert (old_quality > 0).all()
        self.volume_guards = np.minimum(1e-12, old_quality * .01)
        self.coefficient = self.volume_guards / (12 * np.sqrt(3))
        selected = faces[(owner[faces] >= 0).any(axis=1)].copy()
        coordinates = points[selected][:, :, self.axes]
        signed = np.linalg.det(coordinates[:, 1:] - coordinates[:, :1])
        assert (signed != 0).all(), "Moving vertices must be interior to one cap plane"
        negative = signed < 0
        selected[negative, 1], selected[negative, 2] = selected[negative, 2].copy(), selected[negative, 1].copy()
        self.faces, self.tri_owner = selected, owner[selected]
        self.tri_base = (points[selected][:, :, self.axes] - points[selected[:, :1]][:, :, self.axes]) / self.scale
        self.lower = np.concatenate([np.max((-.5 - points[orbit][:, self.axes]) / self.scale, axis=0) for orbit in self.orbits])
        self.upper = np.concatenate([np.min((.5 - points[orbit][:, self.axes]) / self.scale, axis=0) for orbit in self.orbits])

    def evaluate(self, variables):
        displacement = np.vstack((np.asarray(variables).reshape(2, 2), np.zeros(2)))
        tri = self.tri_base + displacement[self.tri_owner]
        edge = tri[:, 1:] - tri[:, :1]
        determinant = edge[:, 0, 0] * edge[:, 1, 1] - edge[:, 0, 1] * edge[:, 1, 0]
        determinant_gradient = np.zeros_like(tri)
        determinant_gradient[:, 1] = np.column_stack((edge[:, 1, 1], -edge[:, 1, 0]))
        determinant_gradient[:, 2] = np.column_stack((-edge[:, 0, 1], edge[:, 0, 0]))
        determinant_gradient[:, 0] = -determinant_gradient[:, 1:].sum(axis=1)
        squared = sum(((tri[:, a] - tri[:, b])**2).sum(1) for a in range(3) for b in range(a))
        target = self.target + 1e-6  # Search margin; direct checks use the unchanged requested target.
        tri_values = 2 * np.sqrt(3) * determinant - target * squared
        tri_gradient = 2 * np.sqrt(3) * determinant_gradient - target * (6 * tri - 2 * tri.sum(axis=1)[:, None])
        shifts = np.zeros((3, 3))
        shifts[:, self.axes] = displacement
        tet = self.tet_base + shifts[self.tet_owner]
        edge = tet[:, 1:] - tet[:, :1]
        determinant = np.linalg.det(edge)
        determinant_gradient = np.zeros_like(tet)
        determinant_gradient[:, 1] = np.cross(edge[:, 1], edge[:, 2])
        determinant_gradient[:, 2] = np.cross(edge[:, 2], edge[:, 0])
        determinant_gradient[:, 3] = np.cross(edge[:, 0], edge[:, 1])
        determinant_gradient[:, 0] = -determinant_gradient[:, 1:].sum(axis=1)
        squared = sum(((tet[:, a] - tet[:, b])**2).sum(1) for a in range(4) for b in range(a))
        ratio = squared / self.original_squared
        tet_values = determinant / self.original_squared**1.5 - self.coefficient * ratio**1.5
        tet_gradient = determinant_gradient / self.original_squared[:, None, None]**1.5
        tet_gradient -= (1.5 * self.coefficient * np.sqrt(ratio) / self.original_squared)[:, None, None] * (
            8 * tet - 2 * tet.sum(axis=1)[:, None])
        tri_jacobian = np.column_stack([np.einsum("tn,tnk->tk", self.tri_owner == group, tri_gradient) for group in range(2)])
        tet_jacobian = np.column_stack([np.einsum("tn,tnk->tk", self.tet_owner == group, tet_gradient)[:, self.axes] for group in range(2)])
        return np.r_[tri_values, tet_values], np.vstack((tri_jacobian, tet_jacobian))

    def moved(self, variables):
        moved = self.points.copy()
        for orbit, displacement in zip(self.orbits, np.asarray(variables).reshape(2, 2), strict=True):
            moved[np.ix_(orbit, self.axes)] += displacement * self.scale
        return moved

    def check(self, variables):
        moved = self.moved(variables)
        triangles = moved[self.faces]
        planar = triangles[:, :, self.axes]
        signed_area = np.linalg.det(planar[:, 1:] - planar[:, :1])
        with np.errstate(invalid="ignore", divide="ignore"):
            cap_quality = bounds(triangles)[0]
            tet_quality = quality(moved, self.cells)
        determinants = np.linalg.det(moved[self.cells[:, 1:]] - moved[self.cells[:, :1]])
        finite = bool(np.isfinite(cap_quality).all() and np.isfinite(tet_quality).all())
        periodic_error = max(float(abs((moved[orbit] - self.points[orbit]) - (moved[orbit[0]] - self.points[orbit[0]])).max())
            for orbit in self.orbits)
        feasible = bool(finite and cap_quality.min() >= self.target and signed_area.min() > 0
            and determinants.min() > 0 and (tet_quality >= self.volume_guards).all() and periodic_error < 1e-12)
        return {"feasible": feasible, "metrics_finite": finite,
            "minimum_cap_quality": float(cap_quality.min()) if np.isfinite(cap_quality).all() else None,
            "minimum_volume_quality": float(tet_quality.min()) if np.isfinite(tet_quality).all() else None,
            "minimum_volume_determinant": float(determinants.min()), "minimum_signed_cap_determinant": float(signed_area.min()),
            "periodic_displacement_error": periodic_error, "candidate_coordinates": moved[self.nodes].tolist()}


def search(model, extra_start=None):
    smooth = []
    for group, orbit in enumerate(model.orbits):
        neighbors = np.unique(model.faces[(model.tri_owner == group).any(axis=1)])
        smooth.extend((model.points[neighbors][:, model.axes].mean(axis=0) - model.points[orbit[0], model.axes]) / model.scale)
    smooth = np.asarray(smooth)
    starts = [("original", np.zeros(4))]
    if extra_start is not None:
        starts.append(("previous-single-orbit", np.asarray(extra_start)))
    starts.extend((("half-smoothing", smooth / 2), ("smoothing", smooth)))
    rng = np.random.default_rng(20261024)
    while len(starts) < 8:
        starts.append((f"perturbation-{len(starts)}", smooth / 2 + rng.normal(0, .1, 4)))
    attempts, chosen = [], None
    for name, initial in starts:
        initial = np.clip(initial, model.lower, model.upper)
        values, _ = model.evaluate(initial)
        count = len(values)
        result = minimize(lambda x: x[-1], np.r_[initial, max(0., -float(values.min())) + 1e-6],
            jac=lambda x: np.r_[np.zeros(4), 1.], method="SLSQP",
            bounds=[*zip(model.lower, model.upper, strict=True), (0, None)],
            constraints={"type": "ineq", "fun": lambda x: model.evaluate(x[:4])[0] + x[-1],
                "jac": lambda x: np.column_stack((model.evaluate(x[:4])[1], np.ones(count)))},
            options={"maxiter": 500, "ftol": 1e-13})
        report = model.check(result.x[:4])
        attempts.append({"start": name, "initial_variables": initial.tolist(), "variables": result.x[:4].tolist(),
            "solver_success": bool(result.success), "solver_message": str(result.message),
            "iterations": int(result.nit), "reported_slack": float(result.x[-1]), **report})
        if report["feasible"]:
            chosen = result.x[:4]
            break
    return chosen, attempts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--witness", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.source, args.witness, args.output = args.source.resolve(), args.witness.resolve(), args.output.resolve()
    assert not args.output.exists()
    witness = json.loads(args.witness.read_text())
    assert witness["source_mesh_sha256"] == workflow.digest(args.source / "mesh.npz")
    prior = Path(witness["source"]) / "result.json"
    assert witness["source_result_sha256"] == workflow.digest(prior)
    paths = [Path(__file__), args.source / "mesh.npz", args.source / "result.json", args.witness, prior]
    paths.extend(workflow.HERE / name for name in ("microgen_workflow.py", "cap_corner_probe.py",
        "cap_quality_bound.py", "mmg_quality_localization.py", "facet_incidence_audit.py", "test_cap_coupled_patch.py"))
    inputs = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    shutil.copyfile(workflow.HERE / "test_cap_coupled_patch.py", args.output.parent / "test_cap_coupled_patch.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_cap_coupled_patch.py"],
        check=True, timeout=45)
    workflow.imports()
    from cap_corner_probe import polygon_corners
    from facet_incidence_audit import audit
    from periodic import periodic_pairs

    mesh = np.load(args.source / "mesh.npz")
    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    incidence = audit(points, cells, surface)
    assert incidence["oriented_facet_consistency_passed"] and not incidence["nonpositive_corner_determinants"]
    _, classes, _ = periodic_pairs(points, surface)
    certificate = witness["diagnostics"][0]
    roots = [certificate["nodes"][0], *certificate["other_free_cap_vertices_in_support"]]
    assert len(roots) == 2
    axis = json.loads(prior.read_text())["diagnostics"][0]["axis"]
    orbits = [np.flatnonzero(classes == classes[node]) for node in roots]
    cap_faces, _ = caps(mesh)
    for side in (-.5, .5):
        selected = cap_faces[np.all(workflow.on_box_plane(points[cap_faces, axis], side), axis=1)]
        rows, ambiguous = polygon_corners(points[:, [other for other in range(3) if other != axis]], selected)
        boundary = {row["node"] for row in rows} | set(ambiguous)
        assert not any(node in boundary for orbit in orbits for node in orbit)
    model = CoupledPatch(points, cells, cap_faces, orbits, axis)
    previous = json.loads(prior.read_text())["diagnostics"][0]
    extra = np.r_[np.asarray(previous["normalized_expansion_point"]) * previous["scale"] / model.scale, [0., 0.]]
    chosen, attempts = search(model, extra)
    periodic = None
    if chosen is not None:
        _, _, periodic = periodic_pairs(model.moved(chosen), surface)
    workflow.save(args.output, {"source": str(args.source), "case": json.loads((args.source / "result.json").read_text())["case"],
        "roots": roots, "orbits": [orbit.tolist() for orbit in orbits], "axis": axis, "unknowns": 4,
        "incident_cap_triangles": len(model.faces), "incident_volume_tets": len(model.cells),
        "target_cap_quality": model.target, "search_target_cap_quality": model.target + 1e-6,
        "maximum_starts": 8, "maximum_iterations_per_start": 500, "attempts": attempts,
        "candidate_found": chosen is not None, "candidate_periodic_check": periodic,
        "infeasibility_proved": False, "source_mesh_unchanged": True, "qualified": False,
        "mmg_calls": [], "stiffness_solve_performed": False,
        "scope": "Independent motion of two periodic cap-interior orbits with fixed outer curves; directly checked local candidate or inconclusive nonlinear search, never a global infeasibility certificate. No move is applied."})
    assert all(workflow.digest(Path(path)) == expected for path, expected in inputs.items())
    print("Completed four-variable cap patch preflight; no source edit, MMG, or FE")


if __name__ == "__main__":
    main()
