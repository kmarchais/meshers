"""One fixed-boundary Tet4 quality repair of the rejected PMY snapshot.

For a single movable vertex, the signed determinant D is affine and the sum
S of squared edge lengths is convex quadratic. The MMG quality constraint
q >= t is D - t*S**1.5/(12*sqrt(3)) >= 0, a convex feasible set. Each complete
incident star gets one nearest-feasible-point solve. This diagnostic neither
generates Tet10 elements nor changes admission gates or native defaults.
"""

import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import sys
import time
import traceback

import microgen_workflow as w


SOURCE = w.HERE / "runs/macbook-pmy-linear-failure-v1/candidate/linear-final-linear.npz"
SOURCE_SHA256 = "06e370191be2c0dcdc43684407f65b70c63548300068c0047b71a13be6b08fae"
TARGET = .100000001
MAX_ITERATIONS = 100
OPTIMIZATION_TIMEOUT_S = 180
WORKER_TIMEOUT_S = 240


class Star:
    """Dimensionless fixed-threshold constraints for one full interior star."""

    def __init__(self, points, cells, node, target=TARGET):
        import numpy as np

        self.node = int(node)
        self.origin = points[node].copy()
        self.slots = np.argmax(cells == node, axis=1)
        assert ((cells == node).sum(1) == 1).all()
        neighbors = np.setdiff1d(np.unique(cells), node)
        distances = np.linalg.norm(points[neighbors] - self.origin, axis=1)
        self.scale = float(np.median(distances))
        assert self.scale > 0
        self.shortest_edge = float(distances.min())
        x = (points[cells] - self.origin) / self.scale
        edges = x[:, 1:] - x[:, :1]
        self.det0 = np.einsum("ij,ij->i", edges[:, 0], np.cross(edges[:, 1], edges[:, 2]))
        gradient = np.empty_like(x)
        gradient[:, 1] = np.cross(edges[:, 1], edges[:, 2])
        gradient[:, 2] = np.cross(edges[:, 2], edges[:, 0])
        gradient[:, 3] = np.cross(edges[:, 0], edges[:, 1])
        gradient[:, 0] = -gradient[:, 1:].sum(1)
        self.det_gradient = gradient[np.arange(len(cells)), self.slots]
        self.squared0 = sum(((x[:, i] - x[:, j]) ** 2).sum(1)
                            for i in range(4) for j in range(i))
        self.squared_gradient0 = (8*x - 2*x.sum(1)[:, None])[np.arange(len(cells)), self.slots]
        assert (self.det0 > 0).all() and (self.squared0 > 0).all()
        self.factor = 12*np.sqrt(3) / self.squared0**1.5
        self.target = target

    def constraint(self, displacement):
        import numpy as np

        squared = self.squared0 + self.squared_gradient0 @ displacement + 3*np.dot(displacement, displacement)
        return self.factor*(self.det0 + self.det_gradient @ displacement) - self.target*(squared/self.squared0)**1.5

    def jacobian(self, displacement):
        import numpy as np

        squared = self.squared0 + self.squared_gradient0 @ displacement + 3*np.dot(displacement, displacement)
        return (self.factor[:, None]*self.det_gradient
                - (1.5*self.target*np.sqrt(squared/self.squared0)/self.squared0)[:, None]
                * (self.squared_gradient0 + 6*displacement))

    def solve(self, callback=None):
        import numpy as np
        from scipy.optimize import minimize

        return minimize(lambda u: .5*np.dot(u, u), np.zeros(3), jac=lambda u: u,
            method="SLSQP", constraints=[dict(type="ineq", fun=self.constraint, jac=self.jacobian)],
            callback=callback, options=dict(maxiter=MAX_ITERATIONS, ftol=1e-12, disp=False))


def quality_summary(values):
    import numpy as np

    return dict(zip(("minimum", "p01", "p10", "median", "p90", "p99", "maximum"),
                    map(float, np.quantile(values, [0, .01, .1, .5, .9, .99, 1]))))


def controls():
    import numpy as np
    from cad_ridge_quality_restore import controls as quality_controls
    from mmg_quality_localization import quality

    report = quality_controls()
    points = np.array([[0., 0., 0.], [1., 0., 0.], [.5, np.sqrt(3)/2, 0.],
                       [.5, np.sqrt(3)/6, np.sqrt(2/3)], [0., 0., 0.]])
    points[4] = .001*points[0] + .333*points[1:4].sum(0)
    cells = np.array([[4, 1, 2, 3], [0, 4, 2, 3], [0, 1, 4, 3], [0, 1, 2, 4]])
    star = Star(points, cells, 4)
    epsilon = 1e-6
    probe = np.array([.007, -.009, .005])
    finite = np.column_stack([(star.constraint(probe + epsilon*axis)
                              - star.constraint(probe - epsilon*axis))/(2*epsilon) for axis in np.eye(3)])
    difference = float(abs(star.jacobian(probe) - finite).max())
    assert difference < 2e-8, "Fixed-threshold analytic Jacobian failed central differences"
    moved = points.copy()
    moved[4] = star.origin + star.scale*probe
    physical = moved[cells]
    det = np.linalg.det(physical[:, 1:] - physical[:, :1]) / star.scale**3
    squared = sum(((physical[:, i]-physical[:, j])**2).sum(1)
                  for i in range(4) for j in range(i)) / star.scale**2
    direct = star.factor*det - star.target*(squared/star.squared0)**1.5
    np.testing.assert_allclose(star.constraint(probe), direct, atol=3e-15, rtol=0)
    old = quality(points, cells)
    assert old.min() < .1
    solved = star.solve()
    moved[4] = star.origin + star.scale*solved.x
    new = quality(moved, cells)
    before_det = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])
    after_det = np.linalg.det(moved[cells[:, 1:]] - moved[cells[:, :1]])
    assert solved.success and new.min() >= .1 and (after_det > 0).all()
    assert abs(after_det.sum()-before_det.sum()) <= 1e-13*before_det.sum()
    # No tetrahedron can exceed quality 1. This deliberately infeasible solve
    # must never be accepted by the independent quality guard.
    impossible = Star(points, cells, 4, target=1.01)
    rejected = impossible.solve()
    assert float(impossible.constraint(rejected.x).min()) < 0
    report.update(convex_constraint_jacobian_maximum_absolute_difference=difference,
        affine_determinant_and_quadratic_edge_formula_match_direct_evaluation=True,
        fixed_star_repair_minimum_before=float(old.min()), fixed_star_repair_minimum_after=float(new.min()),
        fixed_star_repair_iterations=int(solved.nit), fixed_star_signed_volume_preserved=True,
        impossible_quality_above_one_rejected=True, impossible_optimizer_success=bool(rejected.success))
    return report


def repair(a):
    import numpy as np
    from facet_incidence_audit import audit
    from mmg_quality_localization import quality
    from periodic_trace import periodic_pairs

    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    with Job().part("repair") as part:
        part.update(0, "checking fixed-star controls and frozen source", force=True, eta_s=WORKER_TIMEOUT_S)
        w.save(a.run / "controls.json", controls())
        assert w.digest(a.source) == SOURCE_SHA256
        with np.load(a.source) as archive:
            points, cells, surface = [archive[k].copy() for k in ("points", "tetrahedra", "surface")]
        original = points.copy()
        assert cells.shape == (66637, 4) and points.shape == (21387, 3) and surface.shape[1] == 3
        assert max(len(points), len(cells)) <= w.CAP
        before_topology = audit(points, cells, surface)
        assert before_topology["oriented_facet_consistency_passed"]
        assert before_topology["nonpositive_corner_determinants"] == 0
        periodic_before = periodic_pairs(points, surface, declared_bounds=([-.5]*3, [.5]*3))
        initial_quality = quality(points, cells)
        current_quality = initial_quality.copy()
        initial_bad = np.flatnonzero(initial_quality < .1)
        assert len(initial_bad) == 22
        boundary = np.unique(surface)
        nodes = np.setdiff1d(np.unique(cells[initial_bad]), boundary)
        stars = {int(node): np.flatnonzero((cells == node).any(1)) for node in nodes}
        nodes = sorted(map(int, nodes), key=lambda node: (float(initial_quality[stars[node]].min()), node))
        assert nodes and len(nodes) <= 88
        w.save(a.run / "plan.json", dict(initial_bad_elements=initial_bad.tolist(),
            initial_bad_associated_interior_nodes=nodes, ordering="initial incident minimum quality, then node id",
            incident_element_counts={str(n): len(stars[n]) for n in nodes},
            one_pass=True, one_start=True, max_iterations=MAX_ITERATIONS, fixed_target=TARGET))
        records = []
        started = time.monotonic()
        exhausted = False
        accepted = 0
        for index, node in enumerate(nodes):
            elapsed = time.monotonic() - started
            if elapsed >= OPTIMIZATION_TIMEOUT_S:
                exhausted = True
                break
            ids = stars[node]
            before = current_quality[ids].copy()
            row = dict(node=node, incident_elements=ids.tolist(), minimum_before=float(before.min()),
                       below_threshold_before=int((before < .1).sum()), accepted=False)
            if (before >= .1).all():
                row.update(status="skipped-already-feasible")
            else:
                star = Star(points, cells[ids], node)
                iterations = 0

                def callback(displacement):
                    nonlocal iterations
                    iterations += 1
                    elapsed = time.monotonic() - started
                    part.update(.05 + .85*(index + min(iterations/MAX_ITERATIONS, .99))/len(nodes),
                        f"interior node {index+1}/{len(nodes)}, iteration {iterations}",
                        eta_s=max(OPTIMIZATION_TIMEOUT_S-elapsed, 0),
                        bad_tets=int((current_quality < .1).sum()), accepted_moves=accepted,
                        minimum_constraint=float(star.constraint(displacement).min()))
                    if elapsed >= OPTIMIZATION_TIMEOUT_S:
                        raise TimeoutError("Bounded 180 s repair budget exhausted")

                try:
                    solved = star.solve(callback)
                except TimeoutError as error:
                    row.update(status="time-budget-exhausted", error=str(error), iterations=iterations)
                    exhausted = True
                    records.append(row)
                    w.save(a.run / "partial.json", dict(records=records, accepted_moves=accepted,
                        time_budget_exhausted=True, complete=False))
                    break
                if not np.isfinite(solved.x).all():
                    row.update(status="rejected-nonfinite-optimizer-point", optimizer_status=int(solved.status),
                        optimizer_message=str(solved.message), iterations=int(solved.nit))
                    records.append(row)
                    w.save(a.run / "partial.json", dict(records=records, accepted_moves=accepted,
                        time_budget_exhausted=False, complete=False))
                    continue
                candidate = points.copy()
                candidate[node] = star.origin + star.scale*solved.x
                after = quality(candidate, cells[ids])
                det_before = np.linalg.det(points[cells[ids, 1:]] - points[cells[ids, :1]])
                det_after = np.linalg.det(candidate[cells[ids, 1:]] - candidate[cells[ids, :1]])
                volume_relative_change = float(abs(det_after.sum()-det_before.sum())/det_before.sum())
                guards = dict(optimizer_success=bool(solved.success), finite=bool(np.isfinite(candidate[node]).all()),
                    all_incident_quality_at_least_point_one=bool((after >= .1).all()),
                    positive_incident_determinants=bool((det_after > 0).all()),
                    fixed_star_volume_preserved=bool(volume_relative_change <= 1e-12),
                    no_previously_good_incident_tet_fails=bool((after[before >= .1] >= .1).all()))
                row.update(status="accepted" if all(guards.values()) else "rejected", guards=guards,
                    optimizer_status=int(solved.status), optimizer_message=str(solved.message),
                    iterations=int(solved.nit), evaluations=int(solved.nfev),
                    minimum_after=float(after.min()), below_threshold_after=int((after < .1).sum()),
                    constraint_minimum=float(star.constraint(solved.x).min()),
                    signed_star_volume_relative_change=volume_relative_change,
                    original_position=points[node].tolist(), proposed_position=candidate[node].tolist(),
                    displacement=float(np.linalg.norm(candidate[node]-points[node])),
                    median_incident_edge=star.scale, shortest_incident_edge=star.shortest_edge,
                    native_finite_difference_step_over_shortest_edge=(1/16*.002)/star.shortest_edge,
                    native_last_step_over_shortest_edge=(1/16*.1/8)/star.shortest_edge,
                    good_neighbor_quality_decreases=int((after[before >= .1] < before[before >= .1]).sum()))
                if all(guards.values()):
                    points[node] = candidate[node]
                    current_quality[ids] = after
                    row["accepted"] = True
                    accepted += 1
            records.append(row)
            w.save(a.run / "partial.json", dict(records=records, accepted_moves=accepted,
                time_budget_exhausted=exhausted, complete=False))
            part.update(.05 + .85*(index+1)/len(nodes), f"checked interior node {index+1}/{len(nodes)}",
                eta_s=max(OPTIMIZATION_TIMEOUT_S-(time.monotonic()-started), 0),
                bad_tets=int((current_quality < .1).sum()), accepted_moves=accepted,
                minimum_quality=float(current_quality.min()))
        repair_loop_s = time.monotonic()-started
        part.update(.92, "independent whole-mesh geometry, quality and periodic audits", force=True, eta_s=15)
        final_quality = quality(points, cells)
        final_topology = audit(points, cells, surface)
        periodic_after = periodic_pairs(points, surface, declared_bounds=([-.5]*3, [.5]*3))
        initial_volume = float(np.linalg.det(original[cells[:, 1:]]-original[cells[:, :1]]).sum()/6)
        final_volume = float(np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]]).sum()/6)
        changed = np.flatnonzero((original.view(np.uint64) != points.view(np.uint64)).any(1))
        invariant_guards = dict(boundary_coordinates_bitwise_unchanged=bool(np.array_equal(
            original[boundary].view(np.uint64), points[boundary].view(np.uint64))),
            only_initial_selected_interior_nodes_moved=bool(np.isin(changed, nodes).all()),
            positive_volume=final_topology["nonpositive_corner_determinants"] == 0,
            coherent_facets=final_topology["oriented_facet_consistency_passed"],
            fixed_polyhedral_volume=abs(final_volume-initial_volume)/initial_volume <= 1e-12,
            periodic_pairs_exactly_unchanged=all(np.array_equal(x, y) for x, y in zip(periodic_before[0], periodic_after[0], strict=True)),
            periodic_classes_exactly_unchanged=bool(np.array_equal(periodic_before[1], periodic_after[1])),
            no_new_quality_failures=bool((final_quality[initial_quality >= .1] >= .1).all()),
            quality_bookkeeping_matches_independent_recomputation=bool(np.allclose(final_quality, current_quality, atol=2e-15, rtol=0)))
        assert all(invariant_guards.values()), invariant_guards
        np.savez_compressed(a.run / "repaired-linear.npz", points=points, tetrahedra=cells, surface=surface)
        with np.load(a.run / "repaired-linear.npz") as saved:
            assert np.array_equal(saved["tetrahedra"], cells) and np.array_equal(saved["surface"], surface)
            assert np.array_equal(saved["points"].view(np.uint64), points.view(np.uint64))
            with np.load(a.source) as frozen:
                assert np.array_equal(saved["tetrahedra"], frozen["tetrahedra"])
                assert np.array_equal(saved["surface"], frozen["surface"])
                assert saved["points"].shape == frozen["points"].shape
        invariant_guards.update(connectivity_and_surface_arrays_unchanged=True, node_and_element_counts_unchanged=True)
        w.save(a.run / "result.json", dict(complete=True, admitted=False,
            source_sha256=w.digest(a.source), output_sha256=w.digest(a.run / "repaired-linear.npz"),
            nodes=len(points), elements=len(cells), initial_quality=quality_summary(initial_quality),
            final_quality=quality_summary(final_quality), initial_bad_elements=initial_bad.tolist(),
            final_bad_elements=np.flatnonzero(final_quality < .1).tolist(),
            final_linear_quality_threshold_passed=bool((final_quality >= .1).all()),
            selected_nodes=nodes, visited_nodes=len(records), accepted_moves=accepted,
            moved_nodes=changed.tolist(), time_budget_exhausted=exhausted,
            optimization_and_loop_s=repair_loop_s, invariant_guards=invariant_guards,
            initial_polyhedral_volume=initial_volume, final_polyhedral_volume=final_volume,
            initial_topology=before_topology, final_topology=final_topology,
            periodic_before=periodic_before[2], periodic_after=periodic_after[2], records=records,
            scope="Saved rejected Tet4 diagnostic with unchanged boundary and connectivity. No Tet10 admission, FEA, continuum, native meshing-time or competitor claim."))
        part.update(1, "saved interior repair and independent audits", force=True, eta_s=0,
                    bad_tets=int((final_quality < .1).sum()), minimum_quality=float(final_quality.min()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--source", default=SOURCE, type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = parser.parse_args()
    a.run, a.source = a.run.resolve(), a.source.resolve()
    if a.worker:
        try:
            return repair(a)
        except BaseException as error:
            w.save(a.run / "failure.json", dict(type=type(error).__name__, error=str(error),
                traceback=traceback.format_exc(), complete=False, admitted=False))
            raise
    assert w.digest(a.source) == SOURCE_SHA256
    a.run.mkdir(parents=True, exist_ok=False)
    source_result = json.loads((a.source.parent / "result.json").read_text())
    assert any(r["stage"] == "final-linear" and r["sha256"] == SOURCE_SHA256
               for r in source_result["native_linear_trace"])
    files = [Path(__file__).resolve(), a.source, a.source.parent / "result.json", a.source.parent / "input.json"]
    files += [w.HERE / name for name in ("cad_ridge_quality_restore.py", "cad_ridge_cavity_control.py",
        "cad_intersection_curved_control.py", "cad_union_comparison.py", "nested_refinement.py",
        "mmg_quality_localization.py", "facet_incidence_audit.py", "periodic_trace.py",
        "microgen_workflow.py", "job_parts.py")]
    hashes = {str(path): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".json"):
            shutil.copyfile(path, a.run / f"input-{index}-{path.name}")
    w.save(a.run / "manifest.json", dict(inputs=hashes, source=str(a.source), source_sha256=SOURCE_SHA256,
        platform=platform.system(), architecture=platform.machine(), python=sys.version,
        packages={p: importlib.metadata.version(p) for p in ("numpy", "scipy")},
        numerical_threads=w.NUMERICAL_THREAD_POLICY, memory_bytes=w.MEMORY,
        worker_timeout_s=WORKER_TIMEOUT_S, optimization_timeout_s=OPTIMIZATION_TIMEOUT_S,
        max_iterations_per_node=MAX_ITERATIONS, target_quality=TARGET, admission_quality=.1,
        method="One deterministic pass, nearest feasible point per complete incident star, analytic SLSQP constraints",
        fixed_boundary=True, fixed_connectivity=True, fixed_counts=True, saved_tet4_only=True,
        native_options=source_result["options"], admitted=False,
        timing_scope="Python diagnostic repair only; no fresh native meshing or competitor timing"))
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job
    from job_parts import worker_part

    with Job(parts={"worker": 1, "repair": 4}) as job:
        with worker_part(job, "worker") as part:
            record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--run", str(a.run),
                "--source", str(a.source), "--worker"], a.run / "worker.log", part,
                timeout=WORKER_TIMEOUT_S, memory_bytes=w.MEMORY)
        w.save(a.run / "worker.json", record)
        job.result("worker", seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"]/1024**3,
                   returncode=record["returncode"])
        assert record["returncode"] == 0, "Read worker.log, failure.json and partial.json"
        assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
        row = json.loads((a.run / "result.json").read_text())
        row.update(input_hashes_unchanged=True, resources=record)
        w.save(a.run / "result.json", row)
        job.result("PMY fixed-boundary interior repair", minimum_quality=row["final_quality"]["minimum"],
            remaining_bad=len(row["final_bad_elements"]), accepted_moves=row["accepted_moves"], admitted=False)


if __name__ == "__main__":
    main()
