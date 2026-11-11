"""Compare absolute and relative determinant guards on frozen PMY backgrounds.

This is an analytic replay of the first scalar snap pass. It changes no native
guard, generates no mesh, and makes no admission claim. All original positive
determinants, signed post-snap determinants, and local quality values are retained.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


def diagnose(points, cells, case):
    import numpy as np
    from mmg_quality_localization import quality
    from pmy_snap_topology_probe import snapped

    moved, snapping = snapped(points, case)
    before = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
    after = np.linalg.det(moved[cells[:, 1:]]-moved[cells[:, :1]])
    ratio = after/before
    floor = (1/16)**3*.01
    changed = np.any(moved != points, axis=1)
    values = case.function(*np.where(moved == .5, -.5, moved).T)
    original_values = case.function(*np.where(points == .5, -.5, points).T)
    values[changed] = np.where(original_values[changed] >= 0, case.levels[1], case.levels[0])
    scale = snapping["scale"]
    used = ~((values[cells] > case.levels[1]+scale*1e-12).all(1)
             | (values[cells] < case.levels[0]-scale*1e-12).all(1))
    absolute_fail = abs(after) < floor
    relative_fail = ratio < .02
    negative = after <= 0
    q0, q1 = quality(points, cells), quality(moved, cells)
    moderate = absolute_fail & (ratio >= .02) & ~negative
    details = []
    for index in np.flatnonzero(used & (absolute_fail | negative))[:16]:
        tet = cells[index]
        details.append(dict(tet=int(index), vertices=tet.tolist(),
            original_determinant=float(before[index]), snapped_determinant=float(after[index]),
            determinant_ratio=float(ratio[index]), minimum_original_quality=float(q0[index]),
            snapped_quality_magnitude=float(q1[index]),
            absolute_guard_rejects=bool(absolute_fail[index]), relative_guard_rejects=bool(relative_fail[index]),
            positive_orientation=bool(not negative[index]), contains_added_samples=bool((tet >= 9009).any()),
            original_points=points[tet].tolist(), snapped_points=moved[tet].tolist(),
            displacement=np.linalg.norm(moved[tet]-points[tet], axis=1).tolist()))
    counts = lambda mask: dict(all_background=int(mask.sum()), used_by_clipping=int((used & mask).sum()))
    result = dict(nodes=len(points), tetrahedra=len(cells), snapping=snapping,
        absolute_cutoff=floor, original_minimum_determinant=float(before.min()),
        original_below_absolute_cutoff=int((before < floor).sum()),
        largest_cutoff_fraction_of_original=float((floor/before).max()),
        original_minimum_quality=float(q0.min()), snapped_minimum_signed_determinant=float(after.min()),
        minimum_signed_determinant_ratio=float(ratio.min()),
        original_bcc_det_ratio_cutoff=.02,
        absolute_guard_failure=counts(absolute_fail), relative_guard_failure=counts(relative_fail),
        nonpositive_orientation=counts(negative),
        positive_moderate_shrink_absolute_only=counts(moderate),
        absolute_only_with_quality_floor=counts(moderate & (q1 >= .1)),
        first_native_guard_failure=int(np.flatnonzero(used & absolute_fail)[0]) if (used & absolute_fail).any() else None,
        witnesses=details)
    arrays = dict(original_points=points, snapped_points=moved, tetrahedra=cells,
                  original_determinants=before, snapped_determinants=after,
                  original_quality=q0, snapped_quality_magnitude=q1, clipping_active=used)
    return result, arrays


def worker(a):
    import numpy as np
    import pmy_snap_topology_probe as topology
    import pmy_neck_control as neck
    import tpms_family_fixtures as fixtures

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    files = [Path(__file__), Path(topology.__file__), Path(neck.__file__), Path(fixtures.__file__),
             w.ROOT/"crates/meshers-core/src/implicit.rs", w.ROOT/"crates/meshers-core/examples/pmy_background.rs",
             a.source/"baseline/input.json", a.source/"one-round/input.json", a.source/"result.json",
             a.fixtures/"tpms-pmy-sheet/tpms-pmy-sheet-definition.json"]
    hashes = {str(p.resolve()): w.digest(p) for p in files}
    case = fixtures.case("tpms-pmy-sheet", a.fixtures)
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=180, memory_bytes=w.MEMORY,
           numerical_threads=w.NUMERICAL_THREAD_POLICY, scope=__doc__))
    results = {}
    with Job(parts={"baseline": 1, "one-round": 1}) as job:
        for name in ("baseline", "one-round"):
            with job.part(name) as part:
                source = json.loads((a.source/name/"input.json").read_text())
                points, cells = np.array(source["points"]), np.array(source["tetrahedra"], dtype=np.int64)
                row, arrays = diagnose(points, cells, case)
                np.savez_compressed(a.run/(name+"-determinants.npz"), **arrays)
                row["arrays_sha256"] = w.digest(a.run/(name+"-determinants.npz"))
                results[name] = row
                part.update(1, "Compared signed volume and quality before and after snapping",
                    negative=row["nonpositive_orientation"]["used_by_clipping"],
                    absolute_failures=row["absolute_guard_failure"]["used_by_clipping"],
                    relative_failures=row["relative_guard_failure"]["used_by_clipping"])
            job.result(name, pre_snap_rejections=row["original_below_absolute_cutoff"],
                       nonpositive=row["nonpositive_orientation"]["used_by_clipping"],
                       moderate_absolute_only=row["positive_moderate_shrink_absolute_only"]["used_by_clipping"])
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    w.save(a.run/"result.json", dict(complete=True, input_hashes_verified=True, arms=results,
        guard_unchanged=True, generated_candidate=False,
        scope="Analytic first-snap replay. Absolute criterion uses abs(det); signed orientation is independently reported. No native rollback, clipping, optimization or admission."))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-native-background-v3")
    p.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    for name in ("run", "source", "fixtures"):
        setattr(a, name, getattr(a, name).resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run/"input-pmy_snap_guard_probe.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        a.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=180)
    w.save(a.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    for name, row in json.loads((a.run/"result.json").read_text())["arms"].items():
        result(name, pre_snap_rejections=row["original_below_absolute_cutoff"],
               nonpositive=row["nonpositive_orientation"]["used_by_clipping"],
               moderate_absolute_only=row["positive_moderate_shrink_absolute_only"]["used_by_clipping"],
               absolute_failures=row["absolute_guard_failure"]["used_by_clipping"],
               relative_failures=row["relative_guard_failure"]["used_by_clipping"])


if __name__ == "__main__":
    main()
