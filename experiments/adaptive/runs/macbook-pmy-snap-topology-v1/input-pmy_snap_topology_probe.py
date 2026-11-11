"""Isolate post-snap cap triangulation in the PMY background control.

This replays only the bounded scalar snapping equations and Builder rank rule.
It does not generate a candidate or claim native bitwise snapping equivalence.
"""

import argparse
from itertools import product
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


def snapped(points, case):
    import numpy as np
    from pmy_neck_control import raw_gradient

    probes = np.array(list(product((0., .23, .61, 1.), repeat=3))) - .5
    center = sum(case.levels) / 2
    scale = max((case.levels[1]-case.levels[0])/2,
                float(abs(case.function(*probes.T)-center).max()))
    threshold = (case.levels[1]-case.levels[0])/2/scale
    def value(p):
        canonical = np.where(p == .5, -.5, p)
        return (case.function(*canonical.T)-center)/scale

    initial = value(points)
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
        g2 = (gradient**2).sum(1)
        active &= g2 >= 1e-20
        trial[active] -= ((value(trial)[active]-levels[active])/g2[active])[:, None]*gradient[active]
    accepted = eligible & (abs(trial) <= .5).all(1) & (np.linalg.norm(trial-points, axis=1) <= .2/16)
    accepted &= abs(value(trial)-levels) < 1e-11
    result = points.copy()
    result[accepted] = trial[accepted]
    return result, dict(scale=scale, initial_eligible=int(eligible.sum()), accepted=int(accepted.sum()),
                        maximum_move=float(np.linalg.norm(result-points, axis=1).max()),
                        locked_coordinates_preserved=bool(np.array_equal(result[lock], points[lock])))


def worker(a):
    import numpy as np
    import pmy_native_background_control as control
    import tpms_family_fixtures as fixtures

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    source = a.source/"baseline"/"input.json"
    data = json.loads(source.read_text())
    points = np.array(data["points"])
    original = np.array(data["tetrahedra"])
    case = fixtures.case(control.NAME, a.fixtures)
    files = [source, Path(__file__), Path(control.__file__), w.HERE/"pmy_neck_control.py",
             Path(fixtures.__file__), w.ROOT/"crates/meshers-core/src/implicit.rs",
             w.ROOT/"crates/meshers-core/src/lib.rs", w.ROOT/"crates/meshers-core/examples/pmy_background.rs"]
    hashes = {str(p.resolve()): w.digest(p) for p in files}
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=180, memory_bytes=w.MEMORY,
           numerical_threads=w.NUMERICAL_THREAD_POLICY, scope=__doc__))
    with Job(parts={"snap-replay": 1, "cap-topology": 1}) as job:
        with job.part("snap-replay") as part:
            raw = control.raw_bcc_tetrahedra(points)
            np.testing.assert_array_equal(raw, original)
            moved, snap = snapped(points, case)
            part.update(1, "Replayed scalar snapping without changing cap coordinates", moved=snap["accepted"])
        with job.part("cap-topology") as part:
            after = control.raw_bcc_tetrahedra(moved)
            old_set = {tuple(sorted(map(int, t))) for t in raw}
            new_set = {tuple(sorted(map(int, t))) for t in after}
            removed, added = sorted(old_set-new_set), sorted(new_set-old_set)
            assert len(removed) == len(added)
            groups = {}
            for kind, entries in (("removed", removed), ("added", added)):
                for tet in entries:
                    center = [i for i in tet if i >= 17**3]
                    assert len(center) == 1
                    face = [i for i in tet if i != center[0]]
                    caps = [(axis, side) for axis in range(3) for side in (-.5, .5)
                            if (points[face, axis] == side).all()]
                    assert len(caps) == 1, "Difference extends beyond a physical cap pyramid"
                    key = (center[0], *caps[0])
                    groups.setdefault(key, {"removed": [], "added": []})[kind].append(face)
            witnesses = []
            per_cap = {}
            for (center, axis, side), record in sorted(groups.items()):
                assert len(record["removed"]) == len(record["added"]) == 2
                old_diagonal = sorted(set(record["removed"][0]) & set(record["removed"][1]))
                new_diagonal = sorted(set(record["added"][0]) & set(record["added"][1]))
                assert old_diagonal != new_diagonal
                vertices = sorted(set(sum(record["removed"], [])))
                assert set(vertices) == set(sum(record["added"], []))
                per_cap[f"axis{axis}:{side}"] = per_cap.get(f"axis{axis}:{side}", 0)+1
                witnesses.append(dict(center=center, axis=axis, side=side, vertices=vertices,
                    original=points[vertices].tolist(), snapped=moved[vertices].tolist(),
                    old_diagonal=old_diagonal, new_diagonal=new_diagonal))
            part.update(1, "Localized every changed background tet to a cap diagonal", changed_caps=len(groups), removed_tets=len(removed))
        job.result("cap-topology", changed_caps=len(groups), removed_tets=len(removed), all_differences_on_caps=True)
    np.savez_compressed(a.run/"snap-replay.npz", original_points=points, snapped_points=moved,
                        frozen_tetrahedra=raw, dynamic_tetrahedra=after)
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    result = dict(complete=True, input_hashes_verified=True, snap=snap,
        original_background_nodes=len(points), original_background_tets=len(raw),
        removed_tets=len(removed), added_tets=len(added), changed_cap_pyramids=len(groups),
        all_geometric_tet_differences_are_cap_diagonal_changes=True, per_cap=per_cap,
        witnesses=witnesses[:12], artifact_sha256=w.digest(a.run/"snap-replay.npz"),
        scope="First scalar snapping pass, before native topology rollback. Analytic replay, not a bitwise native trace. No clipping or admission.",
        implication="Ordinary clip_lattice chooses cap diagonals after snapping; an externally supplied frozen tetrahedral background fixes them before snapping.",
        next_control="Compare frozen baseline versus its field-only refined background, label fixed triangulation explicitly, retain ordinary native baseline separately and all admission gates.")
    w.save(a.run/"result.json", result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-native-background-v2")
    p.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    for name in ("run", "source", "fixtures"):
        setattr(a, name, getattr(a, name).resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run/"input-pmy_snap_topology_probe.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        a.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=180)
    w.save(a.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"


if __name__ == "__main__":
    main()
