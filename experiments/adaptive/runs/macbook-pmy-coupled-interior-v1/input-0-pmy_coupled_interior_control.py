"""One coupled interior solve per certificate-selected PMY component.

The four fixed cavities come from seven checked single-vertex obstructions.
All physical and cavity boundary vertices stay fixed. This coupled problem is
nonconvex: a failed solve is not an impossibility result. One original start,
no cavity growth, no radius sweep, no retries and no new admission claim.
"""

import argparse
import importlib.metadata
import json
from pathlib import Path
import shutil
import sys
import time
import traceback

import microgen_workflow as w


SOURCE = w.HERE / "runs/macbook-pmy-interior-quality-repair-v3"
CERTIFICATE = w.HERE / "runs/macbook-pmy-star-obstruction-v1"
TARGET = .100001
FTOL = 1e-8
MAX_ITERATIONS = 100
WORKER_TIMEOUT_S = 180
SOLVE_BUDGET_S = 150  # Leave time to save and independently audit accepted moves.


class Coupled:
    def __init__(self, points, cells, nodes):
        import numpy as np

        self.cells = cells
        self.nodes = np.asarray(nodes, dtype=int)
        self.coordinates = points[cells].copy()
        self.origins = points[self.nodes].copy()
        self.indices = [np.where(cells == node) for node in self.nodes]
        self.scales = np.array([np.median(np.linalg.norm(
            points[np.setdiff1d(np.unique(cells[rows]), node)]-points[node], axis=1))
            for node, (rows, _) in zip(self.nodes, self.indices, strict=True)])
        assert (self.scales > 0).all()

    def evaluate(self, variables):
        import numpy as np
        from cad_ridge_quality_restore import quality_gradient

        x = self.coordinates.copy()
        positions = self.origins + self.scales[:, None]*variables.reshape(-1, 3)
        for position, (rows, slots) in zip(positions, self.indices, strict=True):
            x[rows, slots] = position
        q, gradient, determinant = quality_gradient(x)
        jacobian = np.zeros((len(self.cells), len(variables)))
        for index, ((rows, slots), scale) in enumerate(zip(self.indices, self.scales, strict=True)):
            jacobian[rows, 3*index:3*index+3] = gradient[rows, slots]*scale
        return q-TARGET, jacobian, determinant

    def solve(self, callback):
        import numpy as np
        from scipy.optimize import minimize

        return minimize(lambda u: .5*np.dot(u, u), np.zeros(3*len(self.nodes)), jac=lambda u: u,
            method="SLSQP", constraints=[dict(type="ineq", fun=lambda u: self.evaluate(u)[0],
                jac=lambda u: self.evaluate(u)[1])], callback=callback,
            options=dict(maxiter=MAX_ITERATIONS, ftol=FTOL, disp=False))


def controls(path):
    import numpy as np
    from cad_ridge_quality_restore import controls as quality_controls
    from nested_refinement import positive_quality

    report = dict(complete=False, stage="existing-analytic-quality-controls")
    w.save(path, report)
    report.update(quality_controls())
    points = np.array([[0., 0., 0.], [1., 0., 0.], [.5, np.sqrt(3)/2, 0.],
                       [.5, np.sqrt(3)/6, np.sqrt(2/3)], [0., 0., 0.], [0., 0., 0.]])
    points[4] = points[:4].mean(0)
    points[5] = points[[4, 1, 2, 3]].mean(0)
    cells = np.array([[5, 1, 2, 3], [4, 5, 2, 3], [4, 1, 5, 3], [4, 1, 2, 5],
                      [0, 4, 2, 3], [0, 1, 4, 3], [0, 1, 2, 4]])
    model = Coupled(points, cells, [4, 5])
    probe = np.array([.003, -.002, .001, -.001, .002, -.003])
    values, jacobian, determinant = model.evaluate(probe)
    step = 1e-6
    finite = np.column_stack([(model.evaluate(probe+step*axis)[0]-model.evaluate(probe-step*axis)[0])/(2*step)
                             for axis in np.eye(len(probe))])
    moved = points.copy()
    moved[model.nodes] = model.origins + model.scales[:, None]*probe.reshape(-1, 3)
    independent = positive_quality(moved[cells])
    before_det = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
    report.update(stage="coupled-assembly", maximum_jacobian_difference=float(abs(jacobian-finite).max()),
        maximum_independent_quality_difference=float(abs(values+TARGET-independent).max()),
        minimum_probe_determinant=float(determinant.min()),
        volume_relative_change=float(abs(determinant.sum()-before_det.sum())/before_det.sum()),
        coupled_variables=len(probe), synthetic_tets=len(cells))
    w.save(path, report)
    assert report["maximum_jacobian_difference"] < 2e-8
    assert report["maximum_independent_quality_difference"] < 3e-15
    assert determinant.min() > 0 and report["volume_relative_change"] < 1e-13
    report.update(complete=True, stage="passed", nonlinear_control_solves=0)
    w.save(path, report)


def plan(points, cells, surface, certificate):
    import numpy as np
    from facet_incidence_audit import FACES

    choices, selected = [], set()
    for row in certificate["rows"]:
        anchor = min((r for r in row["anchors"] if r["numerical_separation"]),
            key=lambda r: max(r["checked_tangent_upper_bound"], r["independent_tangent_upper_bound"]))
        chosen = [row["node"], *anchor["fixed_interior_vertices_in_support"]]
        selected.update(chosen)
        choices.append(dict(failed_node=row["node"], anchor=anchor["anchor"],
            support_elements=anchor["support_global_elements"], movable_nodes=chosen,
            checked_upper_bound=anchor["checked_tangent_upper_bound"],
            independent_upper_bound=anchor["independent_tangent_upper_bound"]))
    assert sorted(selected) == [158, 1816, 1832, 1846, 5457, 5480, 8580, 9763, 11699, 12892, 12922]
    assert not selected.intersection(map(int, surface.ravel()))
    stars = {node: set(map(int, np.flatnonzero((cells == node).any(1)))) for node in selected}
    groups = []
    remaining = set(selected)
    while remaining:
        group = {min(remaining)}
        pending = list(group)
        while pending:
            node = pending.pop()
            adjacent = {other for other in remaining-group if stars[node] & stars[other]}
            group.update(adjacent)
            pending.extend(sorted(adjacent))
        remaining -= group
        groups.append(sorted(group))
    components = []
    for nodes in groups:
        ids = np.array(sorted(set().union(*(stars[node] for node in nodes))), dtype=int)
        faces = cells[ids][:, FACES].reshape(-1, 3)
        unique, counts = np.unique(np.sort(faces, axis=1), axis=0, return_counts=True)
        boundary = unique[counts == 1]
        assert (counts <= 2).all() and not np.isin(nodes, boundary).any()
        components.append(dict(nodes=nodes, incident_elements=ids.tolist(),
            elements=len(ids), cavity_boundary_triangles=len(boundary),
            all_movable_nodes_strictly_inside_fixed_cavity=True))
    assert [c["elements"] for c in components] == [144, 96, 90, 124]
    for index, component in enumerate(components):
        assert not any(set(component["incident_elements"]) & set(other["incident_elements"])
                       for other in components[index+1:])
    return dict(selection=choices, movable_nodes=sorted(selected), components=components,
        component_order="smallest node id", component_incident_unions_disjoint=True,
        cavity_growth=False, starts_per_component=1, max_iterations=MAX_ITERATIONS)


def worker(a):
    import numpy as np
    from facet_incidence_audit import audit
    from mmg_quality_localization import quality
    from periodic_trace import periodic_pairs
    from pmy_interior_quality_repair import quality_summary

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    with Job().part("control") as part:
        part.update(0, "checking coupled derivatives and frozen selection", force=True,
                    timeout_remaining_s=WORKER_TIMEOUT_S)
        controls(a.run / "controls.json")
        prior = json.loads((a.source / "result.json").read_text())
        certificate = json.loads((a.certificate / "result.json").read_text())
        assert prior["complete"] and not prior["admitted"]
        assert certificate["complete"] and certificate["checked_fixed_neighbor_obstructions"] == 7
        source = a.source / "repaired-linear.npz"
        assert w.digest(source) == prior["output_sha256"]
        cm = json.loads((a.certificate / "manifest.json").read_text())
        assert cm["inputs"][str(source)] == w.digest(source)
        with np.load(source) as archive:
            points, cells, surface = [archive[k].copy() for k in ("points", "tetrahedra", "surface")]
        original = points.copy()
        assert max(len(points), len(cells)) <= w.CAP
        before_topology = audit(points, cells, surface)
        assert before_topology["oriented_facet_consistency_passed"] and not before_topology["nonpositive_corner_determinants"]
        before_periodic = periodic_pairs(points, surface, declared_bounds=([-.5]*3, [.5]*3))
        initial_quality = quality(points, cells)
        assert int((initial_quality < .1).sum()) == 18
        selection = plan(points, cells, surface, certificate)
        w.save(a.run / "plan.json", selection)
        records = []
        started = time.monotonic()
        budget_exhausted = False
        for index, component in enumerate(selection["components"]):
            if time.monotonic()-started >= SOLVE_BUDGET_S:
                budget_exhausted = True
                break
            nodes = np.asarray(component["nodes"], dtype=int)
            ids = np.asarray(component["incident_elements"], dtype=int)
            model = Coupled(points, cells[ids], nodes)
            before = quality(points, cells[ids])
            row = dict(component=index, nodes=nodes.tolist(), elements=len(ids),
                minimum_before=float(before.min()), bad_before=int((before < .1).sum()), accepted=False)
            iterations = 0

            def callback(u):
                nonlocal iterations
                iterations += 1
                elapsed = time.monotonic()-started
                part.update(.1 + .75*(index+min(iterations/MAX_ITERATIONS, .99))/len(selection["components"]),
                    f"coupled component {index+1}/4, iteration {iterations}",
                    timeout_remaining_s=max(SOLVE_BUDGET_S-elapsed, 0),
                    minimum_constraint=float(model.evaluate(u)[0].min()),
                    accepted_components=sum(r["accepted"] for r in records))
                if elapsed >= SOLVE_BUDGET_S:
                    raise TimeoutError("Shared coupled solve budget exhausted")

            try:
                solved = model.solve(callback)
                row.update(optimizer_success=bool(solved.success), optimizer_status=int(solved.status),
                    optimizer_message=str(solved.message), iterations=int(solved.nit), evaluations=int(solved.nfev))
                if not np.isfinite(solved.x).all():
                    row.update(status="rejected-nonfinite-optimizer-point")
                else:
                    candidate = points.copy()
                    candidate[nodes] = model.origins+model.scales[:, None]*solved.x.reshape(-1, 3)
                    after = quality(candidate, cells[ids])
                    before_det = np.linalg.det(points[cells[ids, 1:]]-points[cells[ids, :1]])
                    after_det = np.linalg.det(candidate[cells[ids, 1:]]-candidate[cells[ids, :1]])
                    if not (np.isfinite(after).all() and np.isfinite(after_det).all()):
                        raise FloatingPointError("Nonfinite independent candidate quality or determinant")
                    relative_volume = float(abs(after_det.sum()-before_det.sum())/before_det.sum())
                    guards = dict(optimizer_success=bool(solved.success), finite=bool(np.isfinite(candidate[nodes]).all()),
                        all_cavity_quality_at_least_point_one=bool((after >= .1).all()),
                        positive_cavity_determinants=bool((after_det > 0).all()),
                        fixed_cavity_volume_preserved=relative_volume <= 1e-12,
                        no_previously_good_cavity_tet_fails=bool((after[before >= .1] >= .1).all()))
                    row.update(status="accepted" if all(guards.values()) else "rejected", guards=guards,
                        minimum_after=float(after.min()), bad_after=int((after < .1).sum()),
                        constraint_minimum=float(model.evaluate(solved.x)[0].min()),
                        proposed_positions=candidate[nodes].tolist(), original_positions=points[nodes].tolist(),
                        scaled_displacements=solved.x.reshape(-1, 3).tolist(), node_scales=model.scales.tolist(),
                        maximum_displacement=float(np.linalg.norm(candidate[nodes]-points[nodes], axis=1).max()),
                        signed_cavity_volume_relative_change=relative_volume,
                        good_neighbor_quality_decreases=int((after[before >= .1] < before[before >= .1]).sum()))
                    if all(guards.values()):
                        points[nodes] = candidate[nodes]
                        row["accepted"] = True
            except (AssertionError, FloatingPointError, ValueError, TimeoutError) as error:
                row.update(status="rejected-optimizer-exception", error_type=type(error).__name__,
                    error=str(error), iterations=iterations, traceback=traceback.format_exc())
                if isinstance(error, TimeoutError):
                    budget_exhausted = True
            records.append(row)
            w.save(a.run / "partial.json", dict(complete=False, records=records, budget_exhausted=budget_exhausted))
            part.update(.1 + .75*(index+1)/4, f"recorded coupled component {index+1}/4",
                accepted_components=sum(r["accepted"] for r in records),
                remaining_bad=int((quality(points, cells) < .1).sum()))
            if budget_exhausted:
                break
        solve_loop_s = time.monotonic()-started
        part.update(.9, "auditing accepted components against the frozen mesh", force=True)
        final_quality = quality(points, cells)
        final_topology = audit(points, cells, surface)
        final_periodic = periodic_pairs(points, surface, declared_bounds=([-.5]*3, [.5]*3))
        boundary = np.unique(surface)
        changed = np.flatnonzero((points.view(np.uint64) != original.view(np.uint64)).any(1))
        accepted_nodes = [node for r in records if r["accepted"] for node in r["nodes"]]
        initial_volume = float(np.linalg.det(original[cells[:, 1:]]-original[cells[:, :1]]).sum()/6)
        final_volume = float(np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]]).sum()/6)
        guards = dict(boundary_coordinates_bitwise_unchanged=bool(np.array_equal(
            points[boundary].view(np.uint64), original[boundary].view(np.uint64))),
            only_accepted_component_nodes_moved=bool(np.isin(changed, accepted_nodes).all()),
            positive_volume=final_topology["nonpositive_corner_determinants"] == 0,
            coherent_facets=final_topology["oriented_facet_consistency_passed"],
            fixed_polyhedral_volume=abs(final_volume-initial_volume)/initial_volume <= 1e-12,
            periodic_pairs_exactly_unchanged=all(np.array_equal(x, y) for x, y in zip(before_periodic[0], final_periodic[0], strict=True)),
            periodic_classes_exactly_unchanged=bool(np.array_equal(before_periodic[1], final_periodic[1])),
            no_new_quality_failures=bool((final_quality[initial_quality >= .1] >= .1).all()))
        w.save(a.run / "audit.json", dict(invariant_guards=guards, topology=final_topology,
            periodic=final_periodic[2], initial_volume=initial_volume, final_volume=final_volume))
        assert all(guards.values()), guards
        np.savez_compressed(a.run / "repaired-linear.npz", points=points, tetrahedra=cells, surface=surface)
        with np.load(a.run / "repaired-linear.npz") as saved, np.load(source) as frozen:
            assert np.array_equal(saved["tetrahedra"], frozen["tetrahedra"])
            assert np.array_equal(saved["surface"], frozen["surface"])
            assert saved["points"].shape == frozen["points"].shape
            assert np.array_equal(saved["points"].view(np.uint64), points.view(np.uint64))
        guards.update(connectivity_and_surface_arrays_unchanged=True, node_and_element_counts_unchanged=True)
        w.save(a.run / "result.json", dict(complete=True, admitted=False, source_sha256=w.digest(source),
            output_sha256=w.digest(a.run / "repaired-linear.npz"), nodes=len(points), elements=len(cells),
            initial_quality=quality_summary(initial_quality), final_quality=quality_summary(final_quality),
            initial_bad_elements=np.flatnonzero(initial_quality < .1).tolist(),
            final_bad_elements=np.flatnonzero(final_quality < .1).tolist(),
            final_linear_quality_threshold_passed=bool((final_quality >= .1).all()), records=records,
            selected_nodes=selection["movable_nodes"], component_count=4, visited_components=len(records),
            accepted_components=sum(r["accepted"] for r in records), moved_nodes=changed.tolist(),
            solve_loop_s=solve_loop_s, solve_budget_exhausted=budget_exhausted, invariant_guards=guards,
            initial_polyhedral_volume=initial_volume, final_polyhedral_volume=final_volume,
            initial_topology=before_topology, final_topology=final_topology,
            periodic_before=before_periodic[2], periodic_after=final_periodic[2],
            scope="One nonconvex coupled solve per fixed certificate-selected component. Failed components are rejected whole and do not prove impossibility. Saved Tet4 diagnostic only; no Tet10 admission, FEA, continuum or generation-time comparison."))
        part.update(1, "saved coupled control and independent audits", force=True,
            bad_tets=int((final_quality < .1).sum()), minimum_quality=float(final_quality.min()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--source", default=SOURCE, type=Path)
    parser.add_argument("--certificate", default=CERTIFICATE, type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = parser.parse_args()
    a.run, a.source, a.certificate = a.run.resolve(), a.source.resolve(), a.certificate.resolve()
    if a.worker:
        try:
            return worker(a)
        except BaseException as error:
            w.save(a.run / "failure.json", dict(type=type(error).__name__, error=str(error),
                traceback=traceback.format_exc(), complete=False, admitted=False))
            raise
    a.run.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__).resolve(), a.source / "result.json", a.source / "repaired-linear.npz",
             a.certificate / "result.json", a.certificate / "manifest.json"]
    files += [w.HERE / name for name in ("cad_ridge_quality_restore.py", "pmy_interior_quality_repair.py",
        "nested_refinement.py", "mmg_quality_localization.py", "facet_incidence_audit.py",
        "periodic_trace.py", "microgen_workflow.py", "job_parts.py")]
    hashes = {str(path): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".json"):
            shutil.copyfile(path, a.run / f"input-{index}-{path.name}")
    w.save(a.run / "manifest.json", dict(inputs=hashes, numerical_threads=w.NUMERICAL_THREAD_POLICY,
        memory_bytes=w.MEMORY, worker_timeout_s=WORKER_TIMEOUT_S, shared_solve_budget_s=SOLVE_BUDGET_S,
        starts_per_component=1, max_iterations_per_component=MAX_ITERATIONS, target_quality=TARGET,
        optimizer_ftol=FTOL, admission_quality=.1, fixed_boundary=True, fixed_connectivity=True,
        expected_movable_nodes=11, expected_component_tets=[144, 96, 90, 124],
        no_radius_sweep=True, no_cavity_growth=True, no_restarts=True, nonconvex_problem=True,
        whole_component_rejection=True, admitted=False,
        packages={name: importlib.metadata.version(name) for name in ("numpy", "scipy")}))
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    from job_parts import worker_part

    with Job(parts={"worker": 1, "control": 4}) as job:
        with worker_part(job, "worker") as part:
            record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--run", str(a.run),
                "--source", str(a.source), "--certificate", str(a.certificate), "--worker"],
                a.run / "worker.log", part, timeout=WORKER_TIMEOUT_S, memory_bytes=w.MEMORY)
        w.save(a.run / "worker.json", record)
        job.result("worker", returncode=record["returncode"], seconds=record["wall_s"],
                   peak_gib=record["peak_rss_bytes"]/1024**3)
        assert record["returncode"] == 0, "Read worker.log, controls.json, failure.json, partial.json and audit.json"
        assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
        row = json.loads((a.run / "result.json").read_text())
        row.update(input_hashes_unchanged=True, resources=record)
        w.save(a.run / "result.json", row)
        job.result("PMY coupled interior control", accepted_components=row["accepted_components"],
            remaining_bad=len(row["final_bad_elements"]), minimum_quality=row["final_quality"]["minimum"], admitted=False)


if __name__ == "__main__":
    main()
