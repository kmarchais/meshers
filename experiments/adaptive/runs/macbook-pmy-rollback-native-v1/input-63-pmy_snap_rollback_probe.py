"""Read-only analytic PMY snapping and periodic validity-rollback diagnosis.

This is a NumPy replay of the scalar equations, not a bitwise native trace.
Native clipping, production guards, input meshes and fixtures remain unchanged.
Every rollback restores both original coordinates and original normalized field
values over complete periodic vertex orbits, followed by a global validity scan.
"""

import argparse
from itertools import product
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
TIMEOUT = 180


def replay(points, case):
    import numpy as np
    from pmy_neck_control import raw_gradient

    probes = np.array(list(product((0., .23, .61, 1.), repeat=3)))-.5
    center = sum(case.levels)/2
    scale = max((case.levels[1]-case.levels[0])/2, float(abs(case.function(*probes.T)-center).max()))
    threshold = (case.levels[1]-case.levels[0])/2/scale

    def value(p):
        return (case.function(*np.where(p == .5, -.5, p).T)-center)/scale

    initial = np.asarray(value(points))
    levels = np.copysign(threshold, initial)
    lock = abs(points) == .5
    gradient = raw_gradient(points)/scale
    gradient[lock] = 0
    norms = np.linalg.norm(gradient, axis=1)
    eligible = (norms >= 1e-12) & (abs(initial-levels) <= .2/16*norms)
    active = eligible.copy()
    trial = points.copy()
    for _ in range(8):
        gradient = raw_gradient(trial)/scale
        gradient[lock] = 0
        squared = (gradient**2).sum(1)
        active &= squared >= 1e-20
        trial[active] -= ((value(trial)[active]-levels[active])/squared[active])[:, None]*gradient[active]
    distance = np.linalg.norm(trial-points, axis=1)
    residual = abs(value(trial)-levels)
    accepted = eligible & (abs(trial) <= .5).all(1) & (distance <= .2/16) & (residual < 1e-11)
    moved, values = points.copy(), initial.copy()
    moved[accepted], values[accepted] = trial[accepted], levels[accepted]
    changed_points = np.any(moved != points, axis=1)
    changed_values = values != initial
    return moved, values, initial, threshold, dict(scale=scale, threshold=threshold,
        accepted=int(accepted.sum()), eligible=int(eligible.sum()), moved_points=int(changed_points.sum()),
        changed_values=int(changed_values.sum()), value_only_changes=int((~changed_points & changed_values).sum()),
        maximum_move=float(np.linalg.norm(moved-points, axis=1).max()),
        maximum_accepted_normalized_residual=float(residual[accepted].max(initial=0)),
        locked_coordinates_preserved=bool(np.array_equal(moved[lock], points[lock]))), accepted


def inspect(points, cells, values, threshold, original_det):
    import numpy as np
    from mmg_quality_localization import quality

    determinant = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
    active = ~((values[cells] > threshold+1e-12).all(1) | (values[cells] < -threshold-1e-12).all(1))
    local_h = np.cbrt(2*original_det)
    cutoff = local_h*local_h*local_h*.01
    inverted = determinant <= 0
    subvolume = abs(determinant) < cutoff
    bad = inverted | (active & subvolume)
    q = quality(points, cells)
    summary = dict(inverted=int(inverted.sum()), active_sub2_percent=int((active & subvolume).sum()),
        inactive_positive_sub2_percent=int((~active & subvolume & ~inverted).sum()),
        guard_failures=int(bad.sum()), active_tets=int(active.sum()),
        minimum_signed_determinant=float(determinant.min()),
        minimum_signed_ratio=float((determinant/original_det).min()),
        minimum_quality_magnitude=float(q.min()), positive_below_q=int(((q < .1) & ~inverted).sum()),
        signed_box_volume=float(determinant.sum()/6))
    return summary, dict(determinant=determinant, ratio=determinant/original_det, active=active,
                         inverted=inverted, subvolume=subvolume, bad=bad, quality_magnitude=q, cutoff=cutoff)


def rollback(original, cells, original_values, first_points, first_values, classes, threshold, progress):
    import numpy as np

    original_det = np.linalg.det(original[cells[:, 1:]]-original[cells[:, :1]])
    assert original_det.min() > 0
    points, values = first_points.copy(), first_values.copy()
    initial_changed = np.any(points != original, axis=1) | (values != original_values)
    rounds, all_reverted, previous_bad = [], np.zeros(len(points), dtype=bool), set()
    for iteration in range(int(initial_changed.sum())+1):
        state, arrays = inspect(points, cells, values, threshold, original_det)
        bad = np.flatnonzero(arrays["bad"])
        if not len(bad):
            return points, values, rounds, all_reverted, state, arrays
        changed = np.any(points != original, axis=1) | (values != original_values)
        touched = np.unique(cells[bad])
        orbit_ids = np.unique(classes[touched[changed[touched]]])
        restore = np.isin(classes, orbit_ids) & changed
        assert restore.any(), "Valid original has an unresolved no-progress guard failure: replay/rollback discrepancy"
        stars = np.flatnonzero(restore[cells].any(1))
        rows = dict(round=iteration+1, before=state, bad_tets=bad.tolist(),
            newly_bad_tets=sorted(set(map(int, bad))-previous_bad) if iteration else [],
            periodic_orbits=orbit_ids.tolist(), restored_vertices=np.flatnonzero(restore).tolist(),
            complete_affected_star_tets=stars.tolist(),
            restored_coordinate_vertices=int((restore & np.any(points != original, axis=1)).sum()),
            restored_value_only_vertices=int((restore & ~np.any(points != original, axis=1) & (values != original_values)).sum()))
        points[restore], values[restore] = original[restore], original_values[restore]
        all_reverted |= restore
        np.testing.assert_array_equal(points[all_reverted], original[all_reverted])
        np.testing.assert_array_equal(values[all_reverted], original_values[all_reverted])
        rows["remaining_changed_vertices"] = int((np.any(points != original, axis=1) | (values != original_values)).sum())
        rounds.append(rows)
        previous_bad = set(map(int, bad))
        progress(int(all_reverted.sum()), int(initial_changed.sum()), len(rounds), len(bad), len(stars))
    raise AssertionError("Monotone restoration exceeded its finite original-state bound")


def controls():
    import numpy as np

    original = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    tet = np.array([[0, 1, 2, 3]])
    det = np.array([1.])
    flat = original.copy()
    flat[3, 2] = .01
    active = np.full(4, .5)
    state, _ = inspect(flat, tet, active, 1., det)
    assert state["guard_failures"] == 1 and state["inverted"] == 0
    state, _ = inspect(flat, tet, np.full(4, 2.), 1., det)
    assert state["guard_failures"] == 0 and state["inactive_positive_sub2_percent"] == 1
    flat[3, 2] = -1.
    assert inspect(flat, tet, np.full(4, 2.), 1., det)[0]["inverted"] == 1
    periodic = np.vstack((original, original+[3., 0., 0.]))
    cells = np.vstack((tet, tet+4))
    moved = periodic.copy()
    moved[3, 2] = -1.
    initial = np.full(8, .5)
    values = initial.copy()
    values[[3, 7]] = 1.
    fixed, restored, rounds, reverted, state, _ = rollback(periodic, cells, initial, moved, values,
        np.tile(np.arange(4), 2), 1., lambda *args: None)
    np.testing.assert_array_equal(fixed, periodic)
    np.testing.assert_array_equal(restored, initial)
    assert reverted[[3, 7]].all() and rounds[0]["restored_value_only_vertices"] == 1
    assert state["guard_failures"] == 0
    # A has determinant P.x. B has determinant Q.x-P.x. Initially only A
    # fails: P=-1, Q=-.5. Restoring P=1 makes B fail; restoring Q=2 repairs it.
    chain = np.array([[0., 0., 0.], [0., 1., 0.], [0., 0., 1.], [1., 0., 0.], [2., 0., 0.]])
    cells = np.array([[0, 1, 2, 3], [3, 4, 1, 2]])
    moved = chain.copy()
    moved[[3, 4], 0] = [-1., -.5]
    initial = np.full(5, .5)
    fixed, restored, rounds, _, state, _ = rollback(chain, cells, initial, moved, initial.copy(),
        np.arange(5), 1., lambda *args: None)
    assert len(rounds) == 2 and rounds[0]["bad_tets"] == [0] and rounds[1]["newly_bad_tets"] == [1]
    np.testing.assert_array_equal(fixed, chain)
    np.testing.assert_array_equal(restored, initial)
    assert state["guard_failures"] == 0
    return dict(sign_inversion_rejected_even_if_inactive=True, active_sub2_percent_rejected=True,
        inactive_positive_sub2_percent_matches_native_skip=True, periodic_orbit_restoration=True,
        unchanged_coordinate_value_restored=True, two_round_neighbor_cascade_rescanned=True)


def worker(a):
    import numpy as np
    from microgen.shape import surface_functions
    import pmy_global_preclip_preflight as preflight
    import pmy_neck_control as neck
    import pmy_snap_topology_probe as topology
    import tpms_family_fixtures as fixtures
    from periodic_trace import periodic_pairs

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    source = a.source/"edge/input.json"
    native = json.loads((a.source/"edge/result.json").read_text())
    prepared = json.loads((a.source/"prepare-result.json").read_text())
    parent_manifest = json.loads((a.source/"manifest.json").read_text())
    assert native["status"] == "native-rejected" and "snapping inverted" in native["error"]
    assert w.digest(source) == prepared["generated_input_hashes"]["edge"]
    case = fixtures.case(NAME, a.fixtures)
    files = [Path(__file__), Path(neck.__file__), Path(topology.__file__), Path(preflight.__file__),
        Path(fixtures.__file__), Path(surface_functions.__file__), source, a.source/"result.json",
        a.source/"manifest.json", a.source/"prepare-result.json", a.source/"edge/result.json",
        a.source/"edge/native-output.json", *list((a.fixtures/NAME).iterdir()),
        w.ROOT/"crates/meshers-core/src/implicit.rs", w.ROOT/"crates/meshers-core/src/geometry.rs",
        w.ROOT/"crates/meshers-core/examples/pmy_background.rs"]
    files += [w.HERE/name for name in ("microgen_workflow.py", "periodic_trace.py", "mmg_quality_localization.py")]
    for path in (Path(neck.__file__), Path(surface_functions.__file__),
                 w.ROOT/"crates/meshers-core/src/implicit.rs",
                 w.ROOT/"crates/meshers-core/examples/pmy_background.rs",
                 *list((a.fixtures/NAME).iterdir())):
        assert w.digest(path) == parent_manifest["inputs"][str(path)]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs"):
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, native_error=native["error"], native_options=native["options"],
        guard="Nonpositive signed determinant on every supplied tet; abs(det)<local_h^3*.01 only on clipping-active tets; local_h=(2*original_det)^(1/3)",
        rollback="Restore original point and normalized field value for all changed periodic vertex orbits touching bad tets; rescan every tet until stable",
        scope="Analytic NumPy first-snap and proposed validity rollback, not actual native trace or candidate admission",
        production_changed=False, native_clipping_performed=False))
    data = json.loads(source.read_text())
    points, cells = np.asarray(data["points"]), np.asarray(data["tetrahedra"], dtype=np.int64)
    result = dict(inputs=hashes, native_error=native["error"], admitted=False, native_trace=False,
                  guard_unchanged=True, native_clipping_performed=False)
    with Job(parts={"controls": 1, "first-snap": 2, "periodic-rollback": 3, "final-audit": 2}) as job:
        with job.part("controls") as part:
            result["controls"] = controls()
            part.update(1, "Sign, local-volume, periodic-value and cascade controls passed", controls=len(result["controls"]))
        with job.part("first-snap") as part:
            original, original_det, boundary = preflight.metrics(points, cells)
            _, classes, _ = periodic_pairs(points, boundary, declared_bounds=([-.5]*3, [.5]*3))
            canonical = np.where(points == .5, -.5, points)
            rank = np.floor((canonical+.5)*1e11+.5).astype(np.int64)
            _, rank_classes = np.unique(rank, axis=0, return_inverse=True)
            # Native rollback uses rank(original_point), not the moved position.
            assert all(len(np.unique(classes[rank_classes == key])) == 1 for key in np.unique(rank_classes))
            assert len(np.unique(rank_classes)) == len(np.unique(classes))
            part.update(.25, "Verified original background and native periodic rank classes", nodes=len(points), tets=len(cells))
            moved, values, original_values, threshold, snapping, accepted = replay(points, case)
            old_replay, _ = topology.snapped(points, case)
            np.testing.assert_array_equal(moved, old_replay)
            first, first_arrays = inspect(moved, cells, values, threshold, original_det)
            result.update(original_background=original, snapping=snapping, first_snap=first)
            np.savez_compressed(a.run/"first-snap.npz", original_points=points, tetrahedra=cells,
                original_values=original_values, original_determinant=original_det, snapped_points=moved,
                snapped_values=values, accepted=accepted, periodic_classes=classes, native_rank=rank, **first_arrays)
            details = []
            for index in np.flatnonzero(first_arrays["bad"]):
                tet = cells[index]
                details.append(dict(tet=int(index), vertices=tet.tolist(), active=bool(first_arrays["active"][index]),
                    signed_ratio=float(first_arrays["ratio"][index]), original_points=points[tet].tolist(),
                    snapped_points=moved[tet].tolist(), original_values=original_values[tet].tolist(), snapped_values=values[tet].tolist()))
            w.save(a.run/"first-failure-proof.json", details)
            part.update(1, "Archived every first-snap guard failure", inverted=first["inverted"], bad_tets=first["guard_failures"])
        with job.part("periodic-rollback") as part:
            fixed, restored_values, rounds, reverted, final, final_arrays = rollback(
                points, cells, original_values, moved, values, rank_classes, threshold,
                lambda undone, total, round, bad, stars: part.update(undone/max(total, 1),
                    "Restored periodic orbits and rescanned complete background", restored=undone,
                    initially_changed=total, rounds=round, bad_tets=bad, affected_star_tets=stars))
            w.save(a.run/"rollback-rounds.json", rounds)
            result.update(rollback_rounds=len(rounds), rollback_history=rounds, final=final,
                restored_vertices=int(reverted.sum()), restored_periodic_orbits=int(len(np.unique(rank_classes[reverted]))),
                retained_snapped_coordinates=int(np.any(fixed != points, axis=1).sum()),
                retained_changed_values=int((restored_values != original_values).sum()),
                newly_failing_tets_after_restoration=sum(len(row["newly_bad_tets"]) for row in rounds))
            part.update(1, "Global rescan reached a stable guard-valid state", rounds=len(rounds), restored=int(reverted.sum()), guard_failures=final["guard_failures"])
        with job.part("final-audit") as part:
            np.testing.assert_array_equal(fixed[reverted], points[reverted])
            np.testing.assert_array_equal(restored_values[reverted], original_values[reverted])
            lock = abs(points) == .5
            np.testing.assert_array_equal(fixed[lock], points[lock])
            _, _, periodic = periodic_pairs(fixed, boundary, declared_bounds=([-.5]*3, [.5]*3))
            np.testing.assert_allclose(final_arrays["determinant"].sum()/6, 1., atol=1e-12, rtol=0)
            assert final["guard_failures"] == 0 and final["inverted"] == 0
            original_margin = np.linalg.norm(points[cells[:, 1:]]-points[cells[:, :1]], axis=2).max(1)**3
            negative = first_arrays["inverted"]
            result["failure_robustness"] = dict(
                minimum_inverted_absolute_ratio=float(abs(first_arrays["ratio"][negative]).min()) if negative.any() else None,
                minimum_inverted_det_over_local_edge_cubed=float((abs(first_arrays["determinant"][negative])/original_margin[negative]).min()) if negative.any() else None,
                interpretation="Magnitude diagnostics only; no interval or bitwise native-trace certificate")
            result["final_audit"] = dict(periodic=periodic, locked_coordinates_exact=True,
                restored_coordinates_exact=True, restored_normalized_field_values_exact=True,
                all_background_tets_rescanned=True, unchanged_tetrahedron_connectivity=True,
                stable_no_remaining_cascade=True, minimum_background_quality_is_not_native_admission=True)
            np.savez_compressed(a.run/"rollback-state.npz", points=fixed, tetrahedra=cells,
                values=restored_values, reverted=reverted, **final_arrays)
            part.update(1, "Audited restored values, fixed caps and final signed determinants", minimum_q=final["minimum_quality_magnitude"], below_q=final["positive_below_q"], minimum_ratio=final["minimum_signed_ratio"])
    result["artifacts"] = {name: w.digest(a.run/name) for name in (
        "first-snap.npz", "first-failure-proof.json", "rollback-rounds.json", "rollback-state.npz")}
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    result.update(complete=True, input_hashes_verified=True,
        next_native_control="Before clip_background returns, collect changed original-point ranks from guard-failing supplied tets and restore complete periodic orbits, including original values. Rescan before clipping. Preserve errors if no progress; retain existing topology rollback and unchanged volume/quality/geometry gates.",
        trace_requirement="A native control must record actual pre/post validity-rollback failures and restored counts. This analytic probe does not identify native first-snap coordinates bitwise.")
    w.save(a.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-edge-native-v1")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "source", "fixtures"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-pmy_snap_rollback_probe.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **kw: None), timeout=TIMEOUT)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    findings = json.loads((args.run/"result.json").read_text())
    result("PMY analytic first snap", **findings["first_snap"])
    result("PMY periodic validity rollback", rounds=findings["rollback_rounds"], restored=findings["restored_vertices"],
        newly_failing=findings["newly_failing_tets_after_restoration"], retained_moves=findings["retained_snapped_coordinates"])
    result("PMY rollback final background", **findings["final"])


if __name__ == "__main__":
    main()
