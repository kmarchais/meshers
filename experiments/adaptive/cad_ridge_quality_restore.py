"""One interior-quality restoration solve with the repaired curved wall fixed.

Only frozen field-rank proposal0 is replayed. Its interior corner stars expand
the cavity once. A single constrained solve minimizes sorted-quality slack;
the existing mesh admission criteria remain separate from this research guard.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

import microgen_workflow as w
import cad_ridge_cavity_control as cavity
import cad_intersection_curved_control as curved
import cad_union_comparison as base


SOURCE = w.HERE/"runs/macbook-cad-intersection-curved-v1/native-curved-tet10"
PROPOSALS = w.HERE/"runs/macbook-cad-ridge-general-v1"
ROUND_OFF = 1e-13
MAX_ITERATIONS = 120


def quality_gradient(coordinates):
    """Signed MMG corner quality and analytic derivatives for each tet vertex."""
    import numpy as np

    x = np.asarray(coordinates)
    edges = x[:, 1:]-x[:, :1]
    determinant = np.einsum("ij,ij->i", edges[:, 0], np.cross(edges[:, 1], edges[:, 2]))
    squared = sum(((x[:, a]-x[:, b])**2).sum(1) for a in range(4) for b in range(a))
    assert (squared > 0).all()
    det_gradient = np.empty_like(x)
    det_gradient[:, 1] = np.cross(edges[:, 1], edges[:, 2])
    det_gradient[:, 2] = np.cross(edges[:, 2], edges[:, 0])
    det_gradient[:, 3] = np.cross(edges[:, 0], edges[:, 1])
    det_gradient[:, 0] = -det_gradient[:, 1:].sum(1)
    squared_gradient = 8*x-2*x.sum(1)[:, None]
    factor = 12*np.sqrt(3)/squared**1.5
    quality = factor*determinant
    gradient = factor[:, None, None]*(det_gradient-
        (1.5*determinant/squared)[:, None, None]*squared_gradient)
    return quality, gradient, determinant


def controls():
    import numpy as np
    from nested_refinement import positive_quality

    coordinates = np.array([
        [[0., 0., 0.], [1., 0., 0.], [.5, .8, 0.], [.4, .3, .9]],
        [[-.2, .1, .3], [.6, .2, .2], [.1, .7, .1], [.05, .2, .5]],
    ])
    quality, gradient, determinant = quality_gradient(coordinates)
    assert (determinant > 0).all()
    np.testing.assert_allclose(quality, positive_quality(coordinates), atol=2e-15, rtol=0)
    epsilon = 1e-6
    finite = np.empty_like(gradient)
    for node in range(4):
        for axis in range(3):
            plus, minus = coordinates.copy(), coordinates.copy()
            plus[:, node, axis] += epsilon
            minus[:, node, axis] -= epsilon
            finite[:, node, axis] = (quality_gradient(plus)[0]-quality_gradient(minus)[0])/(2*epsilon)
    difference = float(abs(gradient-finite).max())
    assert difference < 2e-8, "Analytic quality gradient differs from central differences"
    np.testing.assert_allclose(gradient.sum(1), 0, atol=2e-14, rtol=0)
    np.testing.assert_allclose(np.einsum("tij,tij->t", gradient, coordinates), 0, atol=2e-14, rtol=0)
    return dict(signed_quality_matches_existing_positive_metric=True,
        analytic_gradient_maximum_absolute_difference=difference, finite_difference_step=epsilon,
        translation_gradient_sum_zero=True, scale_direction_gradient_zero=True,
        fixed_synthetic_tetrahedra=len(coordinates))


def inventory(points, cells, surface, proposal):
    import numpy as np

    original = np.asarray(proposal["parents"])
    boundary_nodes = set(map(int, surface.ravel()))
    interior = np.array([int(node) for node in np.unique(cells[original, :4]) if int(node) not in boundary_nodes])
    expanded = np.flatnonzero(np.isin(cells[:, :4], interior).any(1) | np.isin(np.arange(len(cells)), original))
    local_nodes = np.unique(cells[expanded])
    physical = cavity.boundary_keys(surface)
    interface = {key: value for key, value in cavity.boundary_keys(cavity.boundary(cells[expanded])).items()
                 if key not in physical}
    interface_nodes = np.unique([node for face in interface.values() for node in face])
    assert not np.isin(interior, interface_nodes).any()
    assert not np.isin(interior, surface).any()
    assert len(interior) == 14 and len(expanded) == 213, "Frozen proposal0 cavity changed"
    clearance = float(.5-abs(points[local_nodes]).max())
    assert clearance > 1e-10, "Interior restoration cavity reaches periodic caps"
    shortest = []
    for node in interior:
        incident = cells[np.any(cells[:, :4] == node, axis=1), :4]
        neighbors = np.setdiff1d(np.unique(incident), [node])
        shortest.append(float(np.linalg.norm(points[neighbors]-points[node], axis=1).min()))
    return dict(interior=interior, expanded=expanded, local_nodes=local_nodes,
        interface_nodes=interface_nodes, shortest=np.asarray(shortest), interface_faces=len(interface),
        cap_clearance=clearance, physical_nodes=np.unique(surface))


class Restoration:
    def __init__(self, repaired_points, original_points, cells, plan):
        import numpy as np

        self.free = plan["interior"]
        self.radii = .1*plan["shortest"]
        self.cells = cells[plan["expanded"], :4]
        self.coordinates = repaired_points[self.cells].copy()
        index = {int(node): position for position, node in enumerate(self.free)}
        self.owner = np.array([[index.get(int(node), -1) for node in tet] for tet in self.cells])
        self.target = np.sort(quality_gradient(original_points[self.cells])[0])
        self.dimension = 3*len(self.free)
        self.last_x = self.last = None

    def evaluate(self, variables):
        import numpy as np

        if self.last_x is not None and np.array_equal(variables, self.last_x):
            return self.last
        direction = np.asarray(variables[:self.dimension]).reshape(-1, 3)
        delta = np.vstack((direction*self.radii[:, None], np.zeros(3)))
        coordinates = self.coordinates+delta[self.owner]
        quality, gradient, determinant = quality_gradient(coordinates)
        order = np.argsort(quality, kind="stable")
        jacobian = np.zeros((len(quality), self.dimension+1))
        for node in range(len(self.free)):
            jacobian[:, 3*node:3*node+3] = (gradient*(self.owner == node)[:, :, None]).sum(1)*self.radii[node]
        jacobian = jacobian[order]
        jacobian[:, -1] = 1
        residual = quality[order]+variables[-1]-self.target
        radius_values = 1-(direction**2).sum(1)
        radius_jacobian = np.zeros((len(self.free), self.dimension+1))
        for node in range(len(self.free)):
            radius_jacobian[node, 3*node:3*node+3] = -2*direction[node]
        values = np.r_[residual, radius_values]
        jac = np.vstack((jacobian, radius_jacobian))
        self.last_x = np.asarray(variables).copy()
        self.last = (values, jac, quality, determinant, direction)
        return self.last

    def points(self, repaired_points, cells, variables):
        import numpy as np

        moved = repaired_points.copy()
        delta = np.asarray(variables[:self.dimension]).reshape(-1, 3)*self.radii[:, None]
        displacements = {int(node): value for node, value in zip(self.free, delta, strict=True)}
        for node, value in displacements.items():
            moved[node] += value
        for (a, b), midpoint in cavity.midpoint_map(cells).items():
            if a in displacements or b in displacements:
                moved[midpoint] += .5*(displacements.get(a, 0)+displacements.get(b, 0))
        return moved


def summarize_quality(values):
    import numpy as np

    return dict(zip(("min", "p01", "p10", "median", "p90", "p99", "max"),
                    map(float, np.quantile(values, [0, .01, .1, .5, .9, .99, 1]))))


def worker(args):
    np, *_ = w.imports()
    from scipy.optimize import minimize
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    from nested_refinement import positive_quality
    from quadratic_edge_removal_probe import measure
    from periodic_trace import periodic_pairs
    from comparison_mesh_audit import quadratic_conformity
    from facet_incidence_audit import audit
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    source, proposals = args.source, args.proposals
    manifest = json.loads((proposals/"manifest.json").read_text())
    assert w.digest(proposals/"frozen-proposals.json") == manifest["proposal_sha256"]
    proposal = json.loads((proposals/"frozen-proposals.json").read_text())["frozen"][0]
    assert proposal["proposal_index"] == 0
    trial = proposals/"proposal-00"
    trial_result = json.loads((trial/"result.json").read_text())
    assert trial_result["accepted_minimum"]
    assert w.digest(trial/"candidate-cavity.npz") == trial_result["candidate_sha256"]
    saved = json.loads((source/"result.json").read_text())
    assert saved["qualified"] and w.digest(source/"mesh.npz") == saved["mesh_sha256"]
    with np.load(source/"mesh.npz") as data:
        original_points, cells, surface = (data[key] for key in ("points", "tetrahedra", "surface"))
    repaired = original_points.copy()
    with np.load(trial/"candidate-cavity.npz") as data:
        np.testing.assert_array_equal(data["parents"], proposal["parents"])
        np.testing.assert_array_equal(data["cells"], cells[data["parents"]])
        np.testing.assert_array_equal(original_points[data["nodes"]], data["points_before"])
        repaired[data["nodes"]] = data["points_after"]
    paths = {Path(__file__), Path(w.__file__), Path(cavity.__file__), Path(curved.__file__), Path(base.__file__),
        Path(base.fixture.__file__), Path(curved.probe.__file__), source/"result.json", source/"mesh.npz", proposals/"manifest.json",
        proposals/"frozen-proposals.json", trial/"result.json", trial/"candidate-cavity.npz", trial/"local-metrics.npz"}
    paths.update(w.HERE/name for name in ("scorecard.py", "nested_refinement.py", "quadratic_edge_removal_probe.py",
        "quadrature_audit.py", "periodic_trace.py", "comparison_mesh_audit.py", "facet_incidence_audit.py", "box_projection.py"))
    paths.update(source/(base.NAME+suffix) for suffix in ("-samples.npz", "-specimen.json"))
    hashes = {str(path.resolve()): w.digest(path) for path in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(source/(base.NAME+suffix), args.run/(base.NAME+suffix))
    np.savez_compressed(args.run/"repaired-before-restoration.npz", points=repaired, tetrahedra=cells, surface=surface)
    result = dict(inputs=hashes, source=str(source), proposals=str(proposals), proposal_index=0,
        complete=False, ordinary_admitted=False, research_sorted_quality_feasible=False,
        restoration_accepted=False, numerical_threads=1, vtk_backend=vtkSMPTools.GetBackend(),
        production_change=False, finite_element_solve_performed=False, optimization_performed=False,
        generation_timing_claim=False,
        scope="One fixed-boundary interior restoration attempt on frozen field-rank proposal0; sorted dominance is an additional research guard")
    w.save(args.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=240, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        optimizer="SLSQP with analytic signed-MMG quality Jacobian and sorted constraints", starts=1,
        maximum_iterations=MAX_ITERATIONS, objective="minimize shared nonnegative sorted-quality slack",
        radius_shortest_incident_edge=.1, radius_source="Original 751 corner coordinates",
        sorted_comparison_roundoff=ROUND_OFF, ordinary_quality_gate=.1,
        boundary="Repaired physical Tri6 boundary and expanded cavity interface fixed bitwise"))
    with Job(parts={"controls": 1, "cavity": 1, "restore": 4, "independent-checks": 2, "global-admission": 4}) as job:
        with job.part("controls") as part:
            part.update(0, "Checking analytic quality gradients", controls=0)
            result["controls"] = controls()
            part.update(1, "Finite-difference gradient controls passed", gradient_error=result["controls"]["analytic_gradient_maximum_absolute_difference"])
        job.result("controls", gradient_error=result["controls"]["analytic_gradient_maximum_absolute_difference"])
        with job.part("cavity") as part:
            part.update(0, "Expanding complete stars of interior corners", original_tets=len(proposal["parents"]))
            plan = inventory(original_points, cells, surface, proposal)
            model = Restoration(repaired, original_points, cells, plan)
            result["cavity"] = dict(interior_corners=plan["interior"].tolist(), tetrahedron_indices=plan["expanded"].tolist(),
                tetrahedra=len(plan["expanded"]), corner_degrees_of_freedom=model.dimension,
                interface_faces=plan["interface_faces"], cap_clearance=plan["cap_clearance"],
                displacement_radii=model.radii.tolist())
            variables = np.zeros(model.dimension+1)
            initial_quality = model.evaluate(variables)[2]
            variables[-1] = max(0., float((model.target-np.sort(initial_quality)).max()))+1e-12
            result["initial_sorted_quality_deficit"] = float((model.target-np.sort(initial_quality)).max())
            result["initialization"] = dict(all_corner_displacements_zero=True, slack=float(variables[-1]))
            part.update(1, "Frozen expanded interface and displacement radii", tets=len(plan["expanded"]), dofs=model.dimension)
        job.result("cavity", tetrahedra=len(plan["expanded"]), dofs=model.dimension,
            sorted_deficit=result["initial_sorted_quality_deficit"])
        w.save(args.run/"result.json", result)
        with job.part("restore") as part:
            history = []
            objective_gradient = np.zeros(model.dimension+1)
            objective_gradient[-1] = 1
            started = time.perf_counter()
            part.update(0, "Starting one sorted-quality slack solve", iterations=0, slack=variables[-1])

            def callback(value):
                constraints, _, qualities, determinant, direction = model.evaluate(value)
                row = dict(iteration=len(history)+1, slack=float(value[-1]),
                    maximum_sorted_deficit=float((model.target-np.sort(qualities)).max()),
                    maximum_radius_fraction=float(np.linalg.norm(direction, axis=1).max()),
                    minimum_corner_determinant=float(determinant.min()),
                    minimum_constraint=float(constraints.min()))
                history.append(row)
                part.update(min(len(history)/MAX_ITERATIONS, .99), "Restoring sorted corner quality",
                    iterations=len(history), slack=row["slack"], sorted_deficit=row["maximum_sorted_deficit"])

            solution = minimize(lambda value: value[-1], variables, jac=lambda _: objective_gradient,
                method="SLSQP", bounds=[(-1., 1.)]*model.dimension+[(0., None)],
                constraints={"type": "ineq", "fun": lambda value: model.evaluate(value)[0],
                             "jac": lambda value: model.evaluate(value)[1]}, callback=callback,
                options={"maxiter": MAX_ITERATIONS, "ftol": 1e-12, "disp": False})
            result["optimizer"] = dict(success=bool(solution.success), status=int(solution.status), message=str(solution.message),
                iterations=int(solution.nit), evaluations=int(solution.nfev), slack=float(solution.x[-1]),
                elapsed_s=time.perf_counter()-started, history=history, variables=solution.x.tolist(),
                one_start=True, no_observer_objective=True)
            result["optimization_performed"] = True
            moved = model.points(repaired, cells, solution.x)
            np.savez_compressed(args.run/"candidate.npz", points=moved, tetrahedra=cells, surface=surface)
            result["candidate_sha256"] = w.digest(args.run/"candidate.npz")
            part.update(1, "Optimizer returned; independent checks follow", iterations=int(solution.nit), slack=float(solution.x[-1]))
        job.result("restore", optimizer_success=bool(solution.success), iterations=int(solution.nit), slack=float(solution.x[-1]))
        w.save(args.run/"result.json", result)
        with job.part("independent-checks") as part:
            part.update(0, "Checking fixed curved domain and complete quality vectors", tets=len(plan["expanded"]))
            np.testing.assert_array_equal(moved[plan["physical_nodes"]], repaired[plan["physical_nodes"]])
            np.testing.assert_array_equal(moved[plan["interface_nodes"]], repaired[plan["interface_nodes"]])
            outside = np.setdiff1d(np.arange(len(moved)), plan["local_nodes"])
            np.testing.assert_array_equal(moved[outside], repaired[outside])
            original_q = positive_quality(original_points[cells[:, :4]])
            repaired_q = positive_quality(repaired[cells[:, :4]])
            restored_q = positive_quality(moved[cells[:, :4]])
            global_delta = np.sort(restored_q)-np.sort(original_q)
            local_delta = np.sort(restored_q[plan["expanded"]])-np.sort(original_q[plan["expanded"]])
            strict = bool((global_delta >= -ROUND_OFF).all())
            if (local_delta >= -ROUND_OFF).all():
                assert strict, "Fixed outside cells did not preserve local sorted dominance"
            old_jac, old_vol = measure(repaired[cells[plan["expanded"]]])
            new_jac, new_vol = measure(moved[cells[plan["expanded"]]])
            volume_change = float(new_vol.sum()-old_vol.sum())
            volume_relative = abs(volume_change)/abs(float(old_vol.sum()))
            radius_fraction = np.linalg.norm((moved[model.free]-repaired[model.free])/model.radii[:, None], axis=1)
            radius_passed = bool(radius_fraction.max() <= 1+1e-10)
            topology = quadratic_conformity(moved, cells, surface)
            assert topology["passed"]
            facets = audit(moved, cells, surface)
            _, _, periodic = periodic_pairs(moved, surface, declared_bounds=([-.5]*3, [.5]*3))
            result["research_sorted_quality_feasible"] = strict
            result["checks"] = dict(physical_boundary_bitwise_unchanged=True, quadratic_interface_bitwise_unchanged=True,
                outside_cavity_bitwise_unchanged=True, counts_and_connectivity_unchanged=True,
                original_quality=summarize_quality(original_q), repaired_quality=summarize_quality(repaired_q),
                restored_quality=summarize_quality(restored_q), minimum_global_sorted_quality_change=float(global_delta.min()),
                minimum_local_sorted_quality_change=float(local_delta.min()), regressed_global_rank_entries=int((global_delta < -ROUND_OFF).sum()),
                maximum_radius_fraction=float(radius_fraction.max()), radius_passed=radius_passed,
                repaired_local_minimum_jacobian=float(old_jac.min()), restored_local_minimum_jacobian=float(new_jac.min()),
                local_volume_change=volume_change, local_volume_relative_difference=volume_relative,
                fixed_domain_volume_balance_passed=bool(volume_relative < 1e-10), periodic=periodic,
                quadratic_conformity=topology, facet_incidence=facets)
            np.savez_compressed(args.run/"quality-vectors.npz", original=original_q, repaired=repaired_q,
                restored=restored_q, global_sorted_difference=global_delta, local_sorted_difference=local_delta,
                restored_local_jacobians=new_jac, repaired_local_jacobians=old_jac)
            part.update(1, "Completed independent fixed-domain checks", sorted_feasible=strict,
                minimum_jacobian=float(new_jac.min()), volume_balance=volume_relative)
        job.result("independent-checks", sorted_feasible=strict, regressed_ranks=result["checks"]["regressed_global_rank_entries"],
            radius_passed=radius_passed, minimum_jacobian=float(new_jac.min()))
        w.save(args.run/"result.json", result)
        with job.part("global-admission") as part:
            result["ordinary_admission"] = curved.archive(args.run, "restored-candidate", (moved, cells, surface), part)
            result["ordinary_admitted"] = result["ordinary_admission"]["qualified"]
            result["restoration_accepted"] = bool(result["ordinary_admitted"] and strict and radius_passed and volume_relative < 1e-10)
        job.result("global-admission", ordinary_admitted=result["ordinary_admitted"],
            restoration_accepted=result["restoration_accepted"], optimizer_success=bool(solution.success))
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    result.update(complete=True, input_hashes_verified=True, original_meshes_and_samples_unchanged=True)
    w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--proposals", type=Path, default=PROPOSALS)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "source", "proposals"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import report, result
    report(0, "Starting bounded fixed-domain quality restoration", eta_s=240, memory_ceiling_gib=w.MEMORY/1024**3)
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=240)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    row = json.loads((args.run/"result.json").read_text())
    result("quality-restoration", optimizer_success=row["optimizer"]["success"], slack=row["optimizer"]["slack"],
        sorted_quality_feasible=row["research_sorted_quality_feasible"], ordinary_admitted=row["ordinary_admitted"],
        restoration_accepted=row["restoration_accepted"], minimum_sorted_change=row["checks"]["minimum_global_sorted_quality_change"],
        q_min=row["checks"]["restored_quality"]["min"])


if __name__ == "__main__":
    main()
