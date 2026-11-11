"""One quality-restoration continuation with a dual-selected interior star added.

The physical candidate starts at the archived 773 solution. Only corner 3419,
selected by the frozen stationarity certificate, is newly freed. This is one
continuation, not a new zero start or a search over cavities and parameters.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

import cad_ridge_quality_restore as restore

w = restore.w
SOURCE = w.HERE/"runs/macbook-cad-ridge-quality-restore-v1"
CERTIFICATE = w.HERE/"runs/macbook-cad-ridge-stationarity-v1"


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

    prior = json.loads((args.source/"result.json").read_text())
    certificate = json.loads((args.certificate/"result.json").read_text())
    assert prior["complete"] and prior["ordinary_admitted"] and certificate["complete"]
    assert certificate["input_hashes_verified"] and not certificate["next_control"]["direction_applied"]
    for path, digest in certificate["inputs"].items():
        assert w.digest(Path(path)) == digest, f"Changed stationarity input: {path}"
    assert str(args.source/"result.json") in certificate["inputs"]
    assert w.digest(args.source/"candidate.npz") == prior["candidate_sha256"]
    original = Path(prior["source"])
    proposals = Path(prior["proposals"])
    for path, digest in prior["inputs"].items():
        assert w.digest(Path(path)) == digest, f"Changed prior restoration input: {path}"
    with np.load(original/"mesh.npz") as archive:
        points, cells, surface = (archive[key] for key in ("points", "tetrahedra", "surface"))
    with np.load(args.source/"repaired-before-restoration.npz") as archive:
        repaired = archive["points"]
        np.testing.assert_array_equal(cells, archive["tetrahedra"])
        np.testing.assert_array_equal(surface, archive["surface"])
    with np.load(args.source/"candidate.npz") as archive:
        initial = archive["points"]
        np.testing.assert_array_equal(cells, archive["tetrahedra"])
        np.testing.assert_array_equal(surface, archive["surface"])
    proposal = json.loads((proposals/"frozen-proposals.json").read_text())["frozen"][0]
    plan = restore.inventory(points, cells, surface, proposal)
    selected = certificate["next_control"]["selected_corner"]
    assert selected == 3419, "Only the approved frozen dual-reaction corner is admitted"
    star = np.flatnonzero((cells[:, :4] == selected).any(1))
    neighbors = np.setdiff1d(np.unique(cells[star, :4]), [selected])
    shortest = float(np.linalg.norm(points[neighbors]-points[selected], axis=1).min())
    plan["interior"] = np.r_[plan["interior"], selected]
    plan["expanded"] = np.union1d(plan["expanded"], star)
    plan["shortest"] = np.r_[plan["shortest"], shortest]
    plan["local_nodes"] = np.unique(cells[plan["expanded"]])
    physical = restore.cavity.boundary_keys(surface)
    interface = {key: face for key, face in restore.cavity.boundary_keys(
        restore.cavity.boundary(cells[plan["expanded"]])).items() if key not in physical}
    plan["interface_nodes"] = np.unique([node for face in interface.values() for node in face])
    plan["interface_faces"] = len(interface)
    plan["cap_clearance"] = float(.5-abs(points[plan["local_nodes"]]).max())
    assert len(plan["interior"]) == 15 and len(plan["expanded"]) == 224
    assert plan["interior"].tolist() == certificate["next_control"]["free_corners"]
    assert not np.isin(plan["interior"], plan["interface_nodes"]).any()
    assert not np.isin(plan["interior"], surface).any()
    assert plan["cap_clearance"] > 1e-10
    model = restore.Restoration(repaired, points, cells, plan)
    variables = np.r_[prior["optimizer"]["variables"][:-1], [0., 0., 0.], prior["optimizer"]["slack"]]
    np.testing.assert_array_equal(model.points(repaired, cells, variables), initial)
    np.testing.assert_array_equal(initial[plan["physical_nodes"]], repaired[plan["physical_nodes"]])
    np.testing.assert_array_equal(initial[plan["interface_nodes"]], repaired[plan["interface_nodes"]])
    initial_deficit = float((model.target-np.sort(model.evaluate(variables)[2])).max())
    variables[-1] = max(0., initial_deficit)+1e-12
    paths = {Path(path) for path in prior["inputs"]} | {Path(path) for path in certificate["inputs"]}
    paths.update((Path(__file__).resolve(), args.source/"result.json", args.source/"candidate.npz",
        args.source/"repaired-before-restoration.npz", args.certificate/"result.json",
        args.certificate/"manifest.json", args.certificate/"derivative-certificate.npz"))
    hashes = {str(path): w.digest(path) for path in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(original/(restore.base.NAME+suffix), args.run/(restore.base.NAME+suffix))
    np.savez_compressed(args.run/"initial-candidate.npz", points=initial, tetrahedra=cells, surface=surface)
    w.save(args.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=240, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        selected_new_corner=selected, selection="Frozen dual reaction certificate, no geometric query selection",
        optimizer="One SLSQP continuation with analytic sorted-MMG constraints", maximum_iterations=restore.MAX_ITERATIONS,
        coordinate_start="Archived 773 solution, plus zero displacement of newly freed corner 3419",
        radius=.1, radius_source="Original 751 shortest incident corner edges",
        radius_center="Repaired pre-773 corner coordinates; bounds apply to total motion",
        sorted_comparison_roundoff=restore.ROUND_OFF, ordinary_quality_gate=.1,
        boundary="Repaired physical Tri6 boundary and new 224-tet cavity interface fixed bitwise"))
    result = dict(inputs=hashes, source=str(args.source), certificate=str(args.certificate), complete=False,
        optimization_performed=False, finite_element_solve_performed=False, ordinary_admitted=False,
        research_sorted_quality_feasible=False, restoration_accepted=False, numerical_threads=1,
        vtk_backend=vtkSMPTools.GetBackend(), production_change=False, generation_timing_claim=False,
        stationarity_diagnosis=dict(prior_slack=prior["optimizer"]["slack"],
            prior_stationarity_residual=certificate["stationarity_residual"], selected_new_corner=selected,
            expanded_directional_derivative=certificate["next_control"]["directional_quality_derivative_minimum"],
            scope="First-order diagnostic only, independent of this attempt's feasibility"),
        initial_sorted_quality_deficit=initial_deficit,
        initialization=dict(archived_773_coordinates_bitwise_reproduced=True, newly_free_corner_displacement_zero=True,
            slack=float(variables[-1]), radius_center_is_original_repaired_mesh=True),
        cavity=dict(interior_corners=plan["interior"].tolist(), tetrahedron_indices=plan["expanded"].tolist(),
            tetrahedra=len(plan["expanded"]), corner_degrees_of_freedom=model.dimension,
            interface_faces=len(interface), cap_clearance=plan["cap_clearance"], displacement_radii=model.radii.tolist()))
    w.save(args.run/"result.json", result)
    with Job(parts={"controls": 1, "restore": 4, "independent-checks": 2, "global-admission": 4}) as job:
        with job.part("controls") as part:
            part.update(0, "Checking gradients and frozen continuation", free_corners=len(model.free), tets=len(model.cells))
            result["controls"] = restore.controls()
            part.update(1, "Gradient and continuation controls passed", gradient_error=result["controls"]["analytic_gradient_maximum_absolute_difference"])
        job.result("controls", free_corners=len(model.free), tets=len(model.cells), initial_deficit=initial_deficit)
        with job.part("restore") as part:
            history = []
            objective_gradient = np.r_[np.zeros(model.dimension), 1.]
            started = time.perf_counter()
            part.update(0, "Continuing from archived 773 coordinates", iterations=0, slack=float(variables[-1]))

            def callback(value):
                constraints, _, quality, determinant, direction = model.evaluate(value)
                row = dict(iteration=len(history)+1, slack=float(value[-1]),
                    maximum_sorted_deficit=float((model.target-np.sort(quality)).max()),
                    maximum_radius_fraction=float(np.linalg.norm(direction, axis=1).max()),
                    minimum_corner_determinant=float(determinant.min()), minimum_constraint=float(constraints.min()))
                history.append(row)
                part.update(min(len(history)/restore.MAX_ITERATIONS, .99), "Restoring quality with one expanded star",
                    iterations=len(history), slack=row["slack"], sorted_deficit=row["maximum_sorted_deficit"])

            solution = minimize(lambda value: value[-1], variables, jac=lambda _: objective_gradient,
                method="SLSQP", bounds=[(-1., 1.)]*model.dimension+[(0., None)],
                constraints={"type": "ineq", "fun": lambda value: model.evaluate(value)[0],
                             "jac": lambda value: model.evaluate(value)[1]}, callback=callback,
                options={"maxiter": restore.MAX_ITERATIONS, "ftol": 1e-12, "disp": False})
            result["optimizer"] = dict(success=bool(solution.success), status=int(solution.status), message=str(solution.message),
                iterations=int(solution.nit), evaluations=int(solution.nfev), slack=float(solution.x[-1]),
                elapsed_s=time.perf_counter()-started, variables=solution.x.tolist(), history=history,
                single_archived_start=True, no_observer_objective=True)
            result["optimization_performed"] = True
            moved = model.points(repaired, cells, solution.x)
            np.savez_compressed(args.run/"candidate.npz", points=moved, tetrahedra=cells, surface=surface)
            result["candidate_sha256"] = w.digest(args.run/"candidate.npz")
            part.update(1, "Continuation returned; independent checks follow", iterations=int(solution.nit), slack=float(solution.x[-1]))
        job.result("restore", optimizer_success=bool(solution.success), initial_slack=initial_deficit, slack=float(solution.x[-1]))
        w.save(args.run/"result.json", result)
        with job.part("independent-checks") as part:
            part.update(0, "Checking fixed boundary and full quality vectors", tets=len(plan["expanded"]))
            np.testing.assert_array_equal(moved[plan["physical_nodes"]], initial[plan["physical_nodes"]])
            np.testing.assert_array_equal(moved[plan["interface_nodes"]], initial[plan["interface_nodes"]])
            outside = np.setdiff1d(np.arange(len(moved)), plan["local_nodes"])
            np.testing.assert_array_equal(moved[outside], initial[outside])
            original_q = positive_quality(points[cells[:, :4]])
            repaired_q = positive_quality(repaired[cells[:, :4]])
            initial_q = positive_quality(initial[cells[:, :4]])
            restored_q = positive_quality(moved[cells[:, :4]])
            global_delta = np.sort(restored_q)-np.sort(original_q)
            local_delta = np.sort(restored_q[plan["expanded"]])-np.sort(original_q[plan["expanded"]])
            strict = bool((global_delta >= -restore.ROUND_OFF).all())
            if (local_delta >= -restore.ROUND_OFF).all():
                assert strict
            old_jac, old_vol = measure(initial[cells[plan["expanded"]]])
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
                original_quality=restore.summarize_quality(original_q), repaired_quality=restore.summarize_quality(repaired_q),
                initial_773_quality=restore.summarize_quality(initial_q), restored_quality=restore.summarize_quality(restored_q),
                minimum_global_sorted_quality_change=float(global_delta.min()), minimum_local_sorted_quality_change=float(local_delta.min()),
                regressed_global_rank_entries=int((global_delta < -restore.ROUND_OFF).sum()),
                maximum_radius_fraction=float(radius_fraction.max()), radius_passed=radius_passed,
                initial_local_minimum_jacobian=float(old_jac.min()), restored_local_minimum_jacobian=float(new_jac.min()),
                local_volume_change=volume_change, local_volume_relative_difference=volume_relative,
                fixed_domain_volume_balance_passed=bool(volume_relative < 1e-10), periodic=periodic,
                quadratic_conformity=topology, facet_incidence=facets)
            np.savez_compressed(args.run/"quality-vectors.npz", original=original_q, repaired=repaired_q,
                initial_773=initial_q, restored=restored_q, global_sorted_difference=global_delta,
                local_sorted_difference=local_delta, restored_local_jacobians=new_jac, initial_local_jacobians=old_jac)
            part.update(1, "Fixed-domain checks complete", sorted_feasible=strict, volume_balance=volume_relative,
                minimum_jacobian=float(new_jac.min()))
        job.result("independent-checks", sorted_feasible=strict, radius_passed=radius_passed,
            regressed_ranks=result["checks"]["regressed_global_rank_entries"], minimum_jacobian=float(new_jac.min()))
        w.save(args.run/"result.json", result)
        with job.part("global-admission") as part:
            result["ordinary_admission"] = restore.curved.archive(args.run, "restored-candidate", (moved, cells, surface), part)
            result["ordinary_admitted"] = result["ordinary_admission"]["qualified"]
            result["restoration_accepted"] = bool(result["ordinary_admitted"] and strict and radius_passed and volume_relative < 1e-10)
        job.result("global-admission", ordinary_admitted=result["ordinary_admitted"], restoration_accepted=result["restoration_accepted"])
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    result.update(complete=True, input_hashes_verified=True, prior_results_and_samples_unchanged=True)
    w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "source", "certificate"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import report, result
    report(0, "Starting one expanded-cavity continuation", eta_s=240, memory_ceiling_gib=w.MEMORY/1024**3)
    receipt = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=240)
    w.save(args.run/"worker.json", receipt)
    assert receipt["returncode"] == 0, "Read worker.log and worker.json"
    row = json.loads((args.run/"result.json").read_text())
    result("quality-continuation", optimizer_success=row["optimizer"]["success"], initial_slack=row["initial_sorted_quality_deficit"],
        slack=row["optimizer"]["slack"], sorted_quality_feasible=row["research_sorted_quality_feasible"],
        ordinary_admitted=row["ordinary_admitted"], restoration_accepted=row["restoration_accepted"],
        minimum_sorted_change=row["checks"]["minimum_global_sorted_quality_change"])


if __name__ == "__main__":
    main()
