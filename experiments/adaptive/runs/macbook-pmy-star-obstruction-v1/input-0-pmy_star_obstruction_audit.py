"""Checked tangent separation for seven archived PMY fixed-neighbor stars.

For concave quality constraints h_i, every feasible point obeys every tangent
halfspace. Nonnegative normalized weights give an upper bound on their weighted
sum. The gradient balance error is bounded over a geometry-derived finite ball:
q=12*sqrt(3)*D/S**1.5, D<=|grad D|*r, S>=3*r*r imply
r<=sqrt(4*|grad D|/q_target), about any fixed opposite-face centroid.

This is a floating-point bound checked with an independent direct evaluator.
It is not interval certification or an obstruction to coupled vertex motion.
No nonlinear solve, mesh displacement or mesh file is produced.
"""

import argparse
import importlib.metadata
import json
from pathlib import Path
import shutil
import sys
import traceback
import warnings

import microgen_workflow as w


SOURCE = w.HERE / "runs/macbook-pmy-interior-quality-repair-v3"
ORIGINAL = w.HERE / "runs/macbook-pmy-linear-failure-v1/candidate/linear-final-linear.npz"
TARGET = .1
NEGATIVE_MARGIN = 1e-10
TIMEOUT = 60


def tangent_audit(points, cells, node, position, target=TARGET):
    import numpy as np
    from scipy.optimize import linprog
    from pmy_interior_quality_repair import Star
    from star_feasibility import StarConstraints

    model = Star(points, cells, node, target=target)
    u = (np.asarray(position)-model.origin)/model.scale
    values, gradients = model.constraint(u), model.jacobian(u)
    # Four equalities yield a basic LP solution with at most four weights.
    with warnings.catch_warnings(record=True) as caught:
        solution = linprog(values, A_eq=np.vstack((gradients.T, np.ones(len(cells)))),
            b_eq=[0., 0., 0., 1.], bounds=(0, None), method="highs-ds",
            options=dict(threads=1, time_limit=2., primal_feasibility_tolerance=1e-9,
                         dual_feasibility_tolerance=1e-9))
    report = dict(target=target, anchor_position=np.asarray(position).tolist(),
        normalized_anchor=u.tolist(), physical_origin=model.origin.tolist(), scale=model.scale,
        linear_success=bool(solution.success), linear_status=int(solution.status),
        linear_message=str(solution.message), linear_iterations=int(solution.nit),
        linear_warnings=[str(item.message) for item in caught],
        numerical_separation=False, formal_interval_certificate=False)
    if not solution.success:
        return report
    raw = solution.x
    assert np.isfinite(raw).all() and raw.min() >= -1e-12 and raw.sum() > 0
    weights = np.maximum(raw, 0)
    weights /= weights.sum()
    residual = weights @ gradients
    weighted_value = float(weights @ values)
    centers = -model.squared_gradient0/6
    radii = np.sqrt(4*np.linalg.norm(model.det_gradient, axis=1)/target)
    upper_bounds = weighted_value + (centers-u) @ residual + np.linalg.norm(residual)*radii
    chosen = int(np.argmin(upper_bounds))
    radius = float(radii[chosen])
    center = centers[chosen]
    upper = float(upper_bounds[chosen])

    # Existing independent implementation reconstructs each determinant and
    # edge sum directly at the anchor. It uses a different coordinate scale.
    independent = StarConstraints(points, cells, node, target)
    iv, ig = independent.evaluate((np.asarray(position)-independent.origin)/independent.scale)
    iv = iv*12*np.sqrt(3)
    ig = ig*12*np.sqrt(3)*model.scale/independent.scale
    value_difference = float(abs(values-iv).max())
    gradient_difference = float(abs(gradients-ig).max())
    independent_residual = weights @ ig
    independent_value = float(weights @ iv)
    independent_upper = float(independent_value + independent_residual @ (center-u)
                              + np.linalg.norm(independent_residual)*radius)
    assert value_difference < 5e-12 and gradient_difference < 5e-11
    support = np.flatnonzero(weights > 0)
    report.update(weights=weights.tolist(), weight_sum=float(weights.sum()),
        raw_minimum_weight=float(raw.min()), positive_weight_support=support.tolist(),
        positive_weight_count=len(support), weighted_constraint_at_anchor=weighted_value,
        weighted_tangent_intercept=float(weighted_value-residual @ u),
        balanced_gradient_residual=residual.tolist(), balanced_gradient_residual_norm=float(np.linalg.norm(residual)),
        selected_radius_element=chosen, normalized_ball_center=center.tolist(), normalized_ball_radius=radius,
        physical_ball_center=(model.origin+model.scale*center).tolist(), physical_ball_radius=model.scale*radius,
        stationarity_residual_ball_contribution=float(np.linalg.norm(residual)*radius),
        center_shift_contribution=float(residual @ (center-u)), checked_tangent_upper_bound=upper,
        independent_tangent_upper_bound=independent_upper,
        independent_value_maximum_absolute_difference=value_difference,
        independent_gradient_maximum_absolute_difference=gradient_difference,
        numerical_separation=bool(max(upper, independent_upper) < -NEGATIVE_MARGIN),
        required_negative_margin=NEGATIVE_MARGIN,
        support_constraints_at_anchor=values[support].tolist(),
        support_gradients_at_anchor=gradients[support].tolist())
    return report


def controls(path):
    import numpy as np
    from mmg_quality_localization import quality
    from pmy_interior_quality_repair import Star

    points = np.array([[0., 0., 0.], [1., 0., 0.], [.5, np.sqrt(3)/2, 0.],
                       [.5, np.sqrt(3)/6, np.sqrt(2/3)], [0., 0., 0.]])
    points[4] = .001*points[0] + .333*points[1:4].sum(0)
    cells = np.array([[4, 1, 2, 3], [0, 4, 2, 3], [0, 1, 4, 3], [0, 1, 2, 4]])
    center = points[:4].mean(0)
    moved = points.copy()
    moved[4] = center
    q = quality(moved, cells)
    det = np.linalg.det(moved[cells[:, 1:]]-moved[cells[:, :1]])
    report = dict(complete=False, feasible_centroid_quality=q.tolist(), feasible_centroid_determinants=det.tolist())
    w.save(path, report)
    assert (q >= TARGET).all() and (det > 0).all()
    report["feasible_initial_anchor"] = tangent_audit(points, cells, 4, points[4])
    w.save(path, report)
    assert not report["feasible_initial_anchor"]["numerical_separation"]
    report["feasible_centroid_anchor"] = tangent_audit(points, cells, 4, center)
    w.save(path, report)
    assert report["feasible_centroid_anchor"]["linear_success"]
    assert not report["feasible_centroid_anchor"]["numerical_separation"]
    model = Star(points, cells, 4, target=TARGET)
    centers = -model.squared_gradient0/6
    radii = np.sqrt(4*np.linalg.norm(model.det_gradient, axis=1)/TARGET)
    u = (center-model.origin)/model.scale
    report["feasible_centroid_inside_all_derived_balls"] = bool((np.linalg.norm(u-centers, axis=1) <= radii).all())
    w.save(path, report)
    assert report["feasible_centroid_inside_all_derived_balls"]
    report["impossible_quality_above_one"] = tangent_audit(points, cells, 4, center, target=1.01)
    w.save(path, report)
    assert report["impossible_quality_above_one"]["numerical_separation"]
    report["complete"] = True
    w.save(path, report)


def worker(a):
    import numpy as np
    from facet_incidence_audit import audit
    from mmg_quality_localization import quality

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    with Job().part("audit") as part:
        part.update(0, "checking feasible and impossible synthetic stars", force=True, timeout_remaining_s=TIMEOUT)
        controls(a.run / "controls.json")
        saved = json.loads((a.source / "result.json").read_text())
        assert saved["complete"] and not saved["admitted"] and saved["accepted_moves"] == 2
        assert w.digest(a.original) == saved["source_sha256"]
        assert w.digest(a.source / "repaired-linear.npz") == saved["output_sha256"]
        with np.load(a.original) as archive:
            points, cells, surface = [archive[key].copy() for key in ("points", "tetrahedra", "surface")]
        assert max(len(points), len(cells)) <= w.CAP
        topology = audit(points, cells, surface)
        assert topology["oriented_facet_consistency_passed"] and not topology["nonpositive_corner_determinants"]
        boundary = set(map(int, np.unique(surface)))
        failed = [row for row in saved["records"] if row["status"] == "rejected"]
        assert len(failed) == 7
        rows = []
        for recorded in saved["records"]:
            node = recorded["node"]
            if recorded["status"] == "accepted":
                points[node] = recorded["proposed_position"]
                continue
            assert recorded["status"] == "rejected" and node not in boundary
            np.testing.assert_array_equal(points[node], recorded["original_position"])
            ids = np.flatnonzero((cells == node).any(1))
            np.testing.assert_array_equal(ids, recorded["incident_elements"])
            initial_quality = quality(points, cells[ids])
            assert abs(initial_quality.min()-recorded["minimum_before"]) < 2e-15
            row = dict(node=node, incident_elements=ids.tolist(), initial_minimum_quality=float(initial_quality.min()),
                archived_optimizer_status=recorded["optimizer_status"], archived_optimizer_minimum_quality=recorded["minimum_after"],
                anchors=[], all_other_vertices_fixed=True)
            for anchor, position in (("original", recorded["original_position"]), ("returned", recorded["proposed_position"])):
                result = tangent_audit(points, cells[ids], node, position)
                result["anchor"] = anchor
                if result["linear_success"]:
                    support = np.asarray(result["positive_weight_support"], dtype=int)
                    result["support_global_elements"] = ids[support].tolist()
                    result["support_cells"] = cells[ids[support]].tolist()
                    result["support_initial_quality"] = initial_quality[support].tolist()
                    other = set(map(int, np.unique(cells[ids[support]])))-{node}
                    result["fixed_boundary_vertices_in_support"] = sorted(other & boundary)
                    result["fixed_interior_vertices_in_support"] = sorted(other-boundary)
                row["anchors"].append(result)
                w.save(a.run / "partial.json", dict(complete=False, rows=rows, current=row))
                part.update(.1 + .85*(len(rows)+(1 if anchor == "original" else 2)/2)/len(failed),
                    f"checked star {len(rows)+1}/{len(failed)}, {anchor} tangent",
                    negative_bound=result.get("checked_tangent_upper_bound"),
                    separations=sum(any(t["numerical_separation"] for t in r["anchors"]) for r in rows))
            row["fixed_neighbor_obstruction"] = any(r["numerical_separation"] for r in row["anchors"])
            row["classification"] = ("checked-single-vertex-fixed-neighbor-obstruction"
                if row["fixed_neighbor_obstruction"] else "unresolved-by-two-archived-tangents")
            rows.append(row)
            w.save(a.run / "partial.json", dict(complete=False, rows=rows))
        # Replaying only archived accepted moves must reproduce the saved v3
        # result exactly. Failed-star certificates use their actual old state.
        with np.load(a.source / "repaired-linear.npz") as repaired:
            np.testing.assert_array_equal(points.view(np.uint64), repaired["points"].view(np.uint64))
            np.testing.assert_array_equal(cells, repaired["tetrahedra"])
            np.testing.assert_array_equal(surface, repaired["surface"])
        result = dict(complete=True, admitted=False, failed_stars=len(rows),
            checked_fixed_neighbor_obstructions=sum(row["fixed_neighbor_obstruction"] for row in rows),
            exact_admission_target=TARGET, optimization_clearance_used=False, rows=rows,
            replay_of_archived_accepted_moves_bitwise_exact=True, nonlinear_solves=0,
            mesh_written=False, source_file_writes_performed=False, formal_interval_certificate=False,
            scope="Floating-point concave-tangent separation with stationarity residual bounded over a geometry-derived ball and independently recomputed tangents. Only one interior vertex is free per star. This does not establish impossibility of coupled interior motion, topology changes or boundary changes; an unresolved LP is not a geometric obstruction.")
        w.save(a.run / "result.json", result)
        part.update(1, "saved fixed-neighbor obstruction evidence", force=True,
            separations=result["checked_fixed_neighbor_obstructions"], failed_stars=len(rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--original", type=Path, default=ORIGINAL)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = parser.parse_args()
    a.run, a.source, a.original = a.run.resolve(), a.source.resolve(), a.original.resolve()
    if a.worker:
        try:
            return worker(a)
        except BaseException as error:
            w.save(a.run / "failure.json", dict(type=type(error).__name__, error=str(error),
                traceback=traceback.format_exc(), complete=False))
            raise
    a.run.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__).resolve(), a.original, a.source / "result.json", a.source / "repaired-linear.npz"]
    files += [w.HERE / name for name in ("pmy_interior_quality_repair.py", "star_feasibility.py",
        "facet_incidence_audit.py", "mmg_quality_localization.py", "microgen_workflow.py", "job_parts.py")]
    hashes = {str(path): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".json"):
            shutil.copyfile(path, a.run / f"input-{index}-{path.name}")
    w.save(a.run / "manifest.json", dict(inputs=hashes, numerical_threads=w.NUMERICAL_THREAD_POLICY,
        highs_threads=1, memory_bytes=w.MEMORY, worker_timeout_s=TIMEOUT, per_lp_time_limit_s=2,
        exact_admission_target=TARGET, required_negative_margin=NEGATIVE_MARGIN,
        maximum_actual_stars=7, expansion_points=["original", "archived-returned"], maximum_actual_lps=14,
        nonlinear_solves=0, no_mesh_output=True, formal_interval_certificate=False,
        radius_derivation="For the fixed opposite-face centroid c, r=norm(u-c), S>=3r^2 and D<=norm(grad D)*r. Hence q>=t implies r<=sqrt(4*norm(grad D)/t). Weighted tangent residual is maximized on this containing ball.",
        packages={name: importlib.metadata.version(name) for name in ("numpy", "scipy")}))
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    from job_parts import worker_part

    with Job(parts={"worker": 1, "audit": 4}) as job:
        with worker_part(job, "worker") as part:
            record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--run", str(a.run),
                "--source", str(a.source), "--original", str(a.original), "--worker"],
                a.run / "worker.log", part, timeout=TIMEOUT, memory_bytes=w.MEMORY)
        w.save(a.run / "worker.json", record)
        job.result("worker", returncode=record["returncode"], seconds=record["wall_s"],
                   peak_gib=record["peak_rss_bytes"]/1024**3)
        assert record["returncode"] == 0, "Read worker.log, controls.json, failure.json and partial.json"
        assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
        row = json.loads((a.run / "result.json").read_text())
        row.update(input_hashes_unchanged=True, resources=record)
        w.save(a.run / "result.json", row)
        job.result("PMY fixed-neighbor obstruction audit", failed_stars=row["failed_stars"],
            checked_obstructions=row["checked_fixed_neighbor_obstructions"], nonlinear_solves=0)


if __name__ == "__main__":
    main()
