"""One bounded, dual-guided cavity-growth continuation of archived CAD repair.

At most four complete-star expansions and 512 tetrahedra. Every expansion
follows a nonnegative stationarity certificate and a derivative check. All
physical boundary coordinates and original displacement radii stay fixed.
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
SOURCE = w.HERE/"runs/macbook-cad-ridge-quality-expand-v1"
LIMITS = dict(growth_steps=4, cavity_tetrahedra=512, iterations=120, timeout_s=240,
    active_tolerance=1e-9, certificate_tolerance=1e-9, sorted_gap=1e-10,
    interior_reaction=1e-8, radius_margin=1e-7, directional_gain=1e-10)


def make_plan(data, free, parents):
    import numpy as np

    expanded = np.union1d(parents, np.flatnonzero(np.isin(data.cells[:, :4], free).any(1)))
    local_nodes = np.unique(data.cells[expanded])
    boundary = restore.cavity.boundary_keys(restore.cavity.boundary(data.cells[expanded]))
    interface = {key: face for key, face in boundary.items() if key not in data.physical_faces}
    interface_nodes = np.unique([node for face in interface.values() for node in face])
    assert not np.isin(free, interface_nodes).any() and not np.isin(free, data.surface).any()
    return dict(interior=np.asarray(free), expanded=expanded, local_nodes=local_nodes,
        interface_nodes=interface_nodes, interface_faces=len(interface),
        cap_clearance=float(.5-abs(data.original[local_nodes]).max()),
        shortest=np.array([data.shortest(int(node)) for node in free]))


def audit_state(data, model, plan, variables, folder):
    import numpy as np
    from nested_refinement import positive_quality
    from quadratic_edge_removal_probe import measure
    from comparison_mesh_audit import quadratic_conformity

    folder.mkdir()
    points = model.points(data.repaired, data.cells, variables)
    np.savez_compressed(folder/"candidate.npz", points=points, tetrahedra=data.cells, surface=data.surface)
    quality = positive_quality(points[data.cells[:, :4]])
    delta = np.sort(quality)-np.sort(data.original_quality)
    local_delta = np.sort(quality[plan["expanded"]])-model.target
    jacobian, volume = measure(points[data.cells[plan["expanded"]]])
    _, old_volume = measure(data.initial[data.cells[plan["expanded"]]])
    volume_difference = float(abs(volume.sum()-old_volume.sum())/abs(old_volume.sum()))
    radius = np.linalg.norm((points[model.free]-data.repaired[model.free])/model.radii[:, None], axis=1)
    outside = np.setdiff1d(np.arange(len(points)), plan["local_nodes"])
    try:
        topology = quadratic_conformity(points, data.cells, data.surface)
    except AssertionError as error:
        topology = dict(passed=False, error=str(error))
    checks = dict(physical_boundary_bitwise_fixed=bool(np.array_equal(points[data.physical_nodes], data.initial[data.physical_nodes])),
        expanded_interface_bitwise_fixed=bool(np.array_equal(points[plan["interface_nodes"]], data.initial[plan["interface_nodes"]])),
        outside_cavity_bitwise_fixed=bool(np.array_equal(points[outside], data.initial[outside])),
        quality_gate=bool(quality.min() >= .1), positive125=bool(jacobian.min() > 0),
        total_radius=bool(radius.max() <= 1+1e-10), volume_balance=volume_difference < 1e-10,
        quadratic_conformity=bool(topology["passed"]))
    row = dict(checks=checks, local_checks_passed=all(checks.values()), quality=restore.summarize_quality(quality),
        global_cdf_passed=bool((delta >= -restore.ROUND_OFF).all()), minimum_global_sorted_change=float(delta.min()),
        regressed_global_rank_entries=int((delta < -restore.ROUND_OFF).sum()),
        individually_regressed_tetrahedra=int((quality < data.original_quality-restore.ROUND_OFF).sum()),
        local_sorted_deficit=float(-local_delta.min()), maximum_radius_fraction=float(radius.max()),
        minimum_local_jacobian=float(jacobian.min()), local_relative_volume_difference=volume_difference,
        free_corners=model.free.tolist(), tetrahedron_indices=plan["expanded"].tolist(),
        tetrahedra=len(model.cells), interface_faces=plan["interface_faces"], cap_clearance=plan["cap_clearance"],
        variables=variables.tolist(), candidate_sha256=w.digest(folder/"candidate.npz"), topology=topology)
    np.savez_compressed(folder/"quality-vectors.npz", quality=quality, global_sorted_difference=delta,
        local_sorted_difference=local_delta, local_minimum_jacobians=jacobian)
    w.save(folder/"state.json", row)
    return points, row


def nonnegative_dual(matrix, objective):
    import numpy as np
    from scipy.optimize import nnls

    raw, raw_residual = nnls(matrix, objective, maxiter=5*max(matrix.shape))
    total = float(raw.sum())
    weights = raw/total if total > 0 else raw
    normalized_residual = float(np.linalg.norm(matrix@weights-objective))
    return raw, weights, dict(raw_dual_weights=raw.tolist(), raw_dual_sum=total,
        raw_dual_residual=float(raw_residual), normalized_dual_weights=weights.tolist(),
        normalized_dual_sum=float(weights.sum()), normalized_dual_residual=normalized_residual,
        passed=bool(total > 0 and normalized_residual <= LIMITS["certificate_tolerance"]))


def stationarity(data, model, plan, variables, points, folder):
    import numpy as np

    values = variables.copy()
    _, _, quality, _, direction = model.evaluate(values)
    values[-1] = max(0., float((model.target-np.sort(quality)).max()))
    residual, jacobian, quality, _, _ = model.evaluate(values)
    order = np.argsort(quality, kind="stable")
    gap = float(np.diff(quality[order]).min())
    radius = float(np.linalg.norm(direction, axis=1).max())
    row = dict(passed=False, required_slack=float(values[-1]), solver_slack=float(variables[-1]),
        sorted_gap=gap, maximum_radius_fraction=radius, selected_interior_corners=[])
    if gap <= LIMITS["sorted_gap"]:
        row["stop_reason"] = "unsupported-sorted-tie"
    elif radius >= 1-LIMITS["radius_margin"]:
        row["stop_reason"] = "active-radius-bound"
    else:
        active = np.flatnonzero(residual[:len(quality)] < LIMITS["active_tolerance"])
        matrix = jacobian[active].T
        objective = np.r_[np.zeros(model.dimension), 1.]
        try:
            raw, weights, fitted = nonnegative_dual(matrix, objective)
            row.update(active_ranks=active.tolist(), active_tetrahedra=plan["expanded"][order[active]].tolist(),
                active_residuals=residual[active].tolist(), **fitted)
            np.savez_compressed(folder/"stationarity.npz", active_jacobian=matrix.T, raw_weights=raw,
                normalized_weights=weights, active_ranks=active)
            if not fitted["passed"]:
                row["stop_reason"] = "nonnegative-stationarity-certificate-failed"
            else:
                row["passed"] = True
                selected_cells = model.cells[order[active]]
                _, gradients, _ = restore.quality_gradient(points[selected_cells])
                reactions = {}
                for weight, cell, gradient in zip(weights, selected_cells, gradients, strict=True):
                    for node, vector in zip(cell, gradient, strict=True):
                        reactions[int(node)] = reactions.get(int(node), np.zeros(3))+weight*vector
                records = []
                for node, vector in sorted(reactions.items()):
                    role = "free" if node in model.free else "wall" if node in data.physical_set else "fixed-interior"
                    scaled = float(np.linalg.norm(vector)*.1*data.shortest(node))
                    records.append(dict(node=node, role=role, reaction=vector.tolist(), scaled_norm=scaled))
                    if role == "fixed-interior" and scaled > LIMITS["interior_reaction"]:
                        row["selected_interior_corners"].append(node)
                row["corner_reactions"] = records
                if not row["selected_interior_corners"]:
                    row["stop_reason"] = "no-supported-interior-reaction"
        except (RuntimeError, ValueError) as error:
            row.update(passed=False, stop_reason="nonnegative-dual-fit-failed", error=str(error))
    w.save(folder/"certificate.json", row)
    return row


def improving_direction(model, variables, folder):
    import numpy as np

    _, _, quality, _, _ = model.evaluate(variables)
    variables[-1] = max(0., float((model.target-np.sort(quality)).max()))
    residual, jacobian, quality, _, _ = model.evaluate(variables)
    active = np.flatnonzero(residual[:len(quality)] < LIMITS["active_tolerance"])
    matrix = jacobian[active, :-1]
    vector = np.linalg.lstsq(matrix, np.ones(len(active)), rcond=None)[0]
    scale = float(np.linalg.norm(vector.reshape(-1, 3), axis=1).max())
    if scale > 0:
        vector /= scale
    gain = matrix@vector
    gap = float(np.diff(np.sort(quality)).min())
    row = dict(passed=bool(scale > 0 and gain.min() > LIMITS["directional_gain"] and gap > LIMITS["sorted_gap"]),
        active_ranks=active.tolist(), active_gradient_rank=int(np.linalg.matrix_rank(matrix)), sorted_gap=gap,
        minimum_quality_directional_derivative=float(gain.min()), maximum_quality_directional_derivative=float(gain.max()),
        direction_applied=False, initial_required_slack=float(variables[-1]))
    np.savez_compressed(folder/"growth-direction.npz", active_quality_jacobian=matrix, normalized_direction=vector)
    w.save(folder/"growth-direction.json", row)
    variables[-1] += 1e-12
    return row


def continue_once(model, variables, part, step):
    import numpy as np
    from scipy.optimize import minimize

    history = []
    started = time.perf_counter()
    objective = np.r_[np.zeros(model.dimension), 1.]

    def callback(value):
        residual, _, quality, _, direction = model.evaluate(value)
        row = dict(iteration=len(history)+1, slack=float(value[-1]),
            required_slack=float((model.target-np.sort(quality)).max()), minimum_constraint=float(residual.min()),
            maximum_radius_fraction=float(np.linalg.norm(direction, axis=1).max()))
        history.append(row)
        part.update((step-1+min(len(history)/LIMITS["iterations"], .99))/LIMITS["growth_steps"],
            "Continuing quality restoration in expanded cavity", growth=step, iterations=len(history),
            tets=len(model.cells), slack=row["slack"])

    solution = minimize(lambda value: value[-1], variables, jac=lambda _: objective, method="SLSQP",
        bounds=[(-1., 1.)]*model.dimension+[(0., None)],
        constraints={"type": "ineq", "fun": lambda value: model.evaluate(value)[0],
            "jac": lambda value: model.evaluate(value)[1]}, callback=callback,
        options={"maxiter": LIMITS["iterations"], "ftol": 1e-12, "disp": False})
    return solution.x, dict(success=bool(solution.success), status=int(solution.status), message=str(solution.message),
        iterations=int(solution.nit), evaluations=int(solution.nfev), slack=float(solution.x[-1]),
        elapsed_s=time.perf_counter()-started, history=history, single_continuation_start=True)


def worker(args):
    np, *_ = w.imports()
    import scipy
    from nested_refinement import positive_quality
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    prior = json.loads((args.source/"result.json").read_text())
    assert prior["complete"] and prior["ordinary_admitted"]
    for path, digest in prior["inputs"].items():
        assert w.digest(Path(path)) == digest, f"Changed frozen input: {path}"
    assert w.digest(args.source/"candidate.npz") == prior["candidate_sha256"]
    original_restore = Path(prior["source"])
    first = json.loads((original_restore/"result.json").read_text())
    original = Path(first["source"])
    with np.load(original/"mesh.npz") as archive:
        points, cells, surface = (archive[key] for key in ("points", "tetrahedra", "surface"))
    with np.load(original_restore/"repaired-before-restoration.npz") as archive:
        repaired = archive["points"]
    with np.load(args.source/"candidate.npz") as archive:
        initial = archive["points"]
        np.testing.assert_array_equal(cells, archive["tetrahedra"])
        np.testing.assert_array_equal(surface, archive["surface"])
    lengths = {}

    def shortest(node):
        if node not in lengths:
            star = cells[(cells[:, :4] == node).any(1), :4]
            neighbors = np.setdiff1d(np.unique(star), [node])
            lengths[node] = float(np.linalg.norm(points[neighbors]-points[node], axis=1).min())
        return lengths[node]

    data = SimpleNamespace(original=points, repaired=repaired, initial=initial, cells=cells, surface=surface,
        physical_nodes=np.unique(surface), physical_set=set(map(int, surface.ravel())),
        physical_faces=restore.cavity.boundary_keys(surface), original_quality=positive_quality(points[cells[:, :4]]),
        shortest=shortest)
    plan = make_plan(data, prior["cavity"]["interior_corners"], prior["cavity"]["tetrahedron_indices"])
    assert len(plan["expanded"]) == 224 and len(plan["interior"]) == 15
    model = restore.Restoration(repaired, points, cells, plan)
    variables = np.asarray(prior["optimizer"]["variables"])
    np.testing.assert_array_equal(model.points(repaired, cells, variables), initial)
    paths = {Path(path) for path in prior["inputs"]} | {Path(__file__).resolve(),
        args.source/"result.json", args.source/"candidate.npz"}
    hashes = {str(path): w.digest(path) for path in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(original/(restore.base.NAME+suffix), args.run/(restore.base.NAME+suffix))
    w.save(args.run/"manifest.json", dict(inputs=hashes, limits=LIMITS, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, native_mesh_caps=w.CAP,
        numpy_version=np.__version__, scipy_version=scipy.__version__,
        radius_fraction=.1, radius_lengths="Original 751 shortest corner edges", radius_center="Original repaired proposal0",
        start="Archived 778 solution", growth="All fixed interior corners with positive-dual scaled reaction above 1e-8",
        whole_growth_rejected_on_cap=True, multistart=False,
        scope="Research CDF guard is additional to ordinary mesh admission; local slack across different cavity sizes is not the comparison metric"))
    result = dict(inputs=hashes, complete=False, source=str(args.source), production_change=False,
        finite_element_solve_performed=False, generation_timing_claim=False, stages=[], ordinary_admitted=False,
        numerical_threads=1, vtk_backend=vtkSMPTools.GetBackend(), original_quality=restore.summarize_quality(data.original_quality),
        repaired_quality=restore.summarize_quality(positive_quality(repaired[cells[:, :4]])))
    with Job(parts={"controls": 1, "continuations": 4, "global-admission": 3}) as job:
        with job.part("controls") as part:
            part.update(0, "Checking fixed continuation and analytic gradients", tets=len(model.cells))
            result["controls"] = restore.controls()
            target = np.array([0., 1.])
            _, _, positive = nonnegative_dual(np.array([[-1., 1.], [1., 1.]]), target)
            _, _, negative = nonnegative_dual(np.array([[1., 2.], [1., 1.]]), target)
            assert positive["passed"] and not negative["passed"]
            result["controls"]["nonnegative_dual"] = dict(positive=positive, negative=negative,
                scope="Reject a fixture whose exact unconstrained weights are [2,-1]")
            part.update(1, "Fixed-input and gradient controls passed", gradient_error=result["controls"]["analytic_gradient_maximum_absolute_difference"])
        with job.part("continuations") as part:
            part.update(0, "Auditing archived 778 state", growth=0, tets=len(model.cells))
            optimizer = None
            for step in range(LIMITS["growth_steps"]+1):
                folder = args.run/f"stage-{step:02d}"
                moved, state = audit_state(data, model, plan, variables, folder)
                if optimizer is not None:
                    w.save(folder/"optimizer.json", optimizer)
                result["stages"].append(dict(stage=step, **state, optimizer=optimizer))
                job.result(f"stage-{step}", tets=len(model.cells), cdf_passed=state["global_cdf_passed"],
                    global_sorted_change=state["minimum_global_sorted_change"], local_checks=state["local_checks_passed"])
                w.save(args.run/"result.json", result)
                if not state["local_checks_passed"]:
                    reason = "candidate-local-check-failed"
                    break
                if state["global_cdf_passed"]:
                    reason = "full-global-cdf-guard-passed"
                    break
                certificate = stationarity(data, model, plan, variables, moved, folder)
                result["stages"][-1]["certificate"] = certificate
                if not certificate["passed"] or not certificate["selected_interior_corners"]:
                    reason = certificate["stop_reason"]
                    break
                if step == LIMITS["growth_steps"]:
                    reason = "growth-step-limit"
                    break
                added = certificate["selected_interior_corners"]
                next_plan = make_plan(data, np.r_[model.free, added], plan["expanded"])
                growth = dict(added_corners=added, proposed_tetrahedra=len(next_plan["expanded"]),
                    proposed_cap_clearance=next_plan["cap_clearance"], accepted=False)
                result["stages"][-1]["growth"] = growth
                if len(next_plan["expanded"]) > LIMITS["cavity_tetrahedra"] or next_plan["cap_clearance"] <= 1e-10:
                    reason = "whole-growth-resource-or-periodic-cap-rejection"
                    w.save(folder/"growth.json", growth)
                    break
                next_model = restore.Restoration(repaired, points, cells, next_plan)
                next_variables = np.r_[variables[:-1], np.zeros(3*len(added)), variables[-1]]
                np.testing.assert_array_equal(next_model.points(repaired, cells, next_variables), moved)
                direction = improving_direction(next_model, next_variables, folder)
                growth["direction"] = direction
                if not direction["passed"]:
                    reason = "expanded-first-order-direction-check-failed"
                    w.save(folder/"growth.json", growth)
                    break
                growth["accepted"] = True
                w.save(folder/"growth.json", growth)
                plan, model = next_plan, next_model
                variables, optimizer = continue_once(model, next_variables, part, step+1)
            result.update(stop_reason=reason, growth_steps=len(result["stages"])-1,
                final_state=state, research_sorted_quality_feasible=state["global_cdf_passed"])
            part.update(1, "Bounded cavity continuation stopped", growth=result["growth_steps"],
                global_sorted_change=state["minimum_global_sorted_change"], cdf_passed=state["global_cdf_passed"])
        w.save(args.run/"result.json", result)
        with job.part("global-admission") as part:
            result["ordinary_admission"] = restore.curved.archive(args.run, "final-candidate", (moved, cells, surface), part)
            result["ordinary_admitted"] = result["ordinary_admission"]["qualified"]
            result["restoration_accepted"] = bool(result["ordinary_admitted"] and state["local_checks_passed"] and state["global_cdf_passed"])
        job.result("global-admission", ordinary_admitted=result["ordinary_admitted"], restoration_accepted=result["restoration_accepted"])
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    result.update(complete=True, input_hashes_verified=True, frozen_inputs_unchanged=True)
    w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.run, args.source = args.run.resolve(), args.source.resolve()
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import report, result
    report(0, "Starting bounded dual-guided cavity continuation", eta_s=LIMITS["timeout_s"], memory_ceiling_gib=w.MEMORY/1024**3)
    receipt = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=LIMITS["timeout_s"])
    w.save(args.run/"worker.json", receipt)
    assert receipt["returncode"] == 0, "Read worker.log and worker.json"
    row = json.loads((args.run/"result.json").read_text())
    result("adaptive-cavity", growth_steps=row["growth_steps"], stop_reason=row["stop_reason"],
        global_sorted_change=row["final_state"]["minimum_global_sorted_change"],
        ordinary_admitted=row["ordinary_admitted"], cdf_passed=row["research_sorted_quality_feasible"],
        restoration_accepted=row["restoration_accepted"])


if __name__ == "__main__":
    main()
