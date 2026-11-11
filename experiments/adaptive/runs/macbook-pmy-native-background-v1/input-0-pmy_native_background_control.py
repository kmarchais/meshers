"""Run frozen generic PMY preclip sampling through native mesh generation.

The parent submits this driver through qrun. Each preparation/generation/audit
worker has a 180-second, 8-GiB process-tree bound and standard Job progress.
The ordinary native BCC run and the same externally supplied BCC partition must
agree before the refined background is tested. No admission gates are relaxed.
"""

import argparse
from functools import partial
import json
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import traceback

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
RUNS = w.HERE / "runs"
TIMEOUT = 180


def reordered_background(original, refined):
    """Restore Builder's original node ordering without moving any coordinate."""
    import numpy as np

    p0, t0 = original
    p1, t1 = refined
    lookup = {tuple(p): i for i, p in enumerate(p0)}
    corners = [(-.5+x/16, -.5+y/16, -.5+z/16)
               for z in range(17) for y in range(17) for x in range(17)]
    centers = [(-.5+(x+.5)/16, -.5+(y+.5)/16, -.5+(z+.5)/16)
               for z in range(16) for y in range(16) for x in range(16)]
    order0 = np.array([lookup[p] for p in corners+centers])
    assert len(order0) == len(p0) == 9009 and len(np.unique(order0)) == len(p0)
    np.testing.assert_array_equal(p1[:len(p0)], p0)
    output = []
    for points, cells, order in ((p0, t0, order0),
                                (p1, t1, np.r_[order0, np.arange(len(p0), len(p1))])):
        inverse = np.empty(len(order), dtype=np.int64)
        inverse[order] = np.arange(len(order))
        output.append((points[order], inverse[cells]))
    return output


def prepare(a):
    import numpy as np
    import pmy_global_preclip_preflight as preflight
    import tpms_family_fixtures as fixtures

    frozen = json.loads((a.background / "result.json").read_text())
    assert frozen["preflight_passed"] and not frozen["clipping_performed"]
    assert w.digest(a.background / "selection.npz") == frozen["selection_sha256"]
    fixture = a.fixtures / NAME
    known = {Path(p).resolve(): h for p, h in frozen["inputs"].items()}
    assert all(w.digest(p) == known[p.resolve()] for p in fixture.iterdir())
    case = fixtures.case(NAME, a.fixtures)
    with np.load(a.background / "background-before.npz") as f:
        original = f["points"], f["tetrahedra"]
    with np.load(a.background / "background-one-round.npz") as f:
        refined = f["points"], f["tetrahedra"]
    faces, count, parents, _ = preflight.inventory(original[1])
    plan = preflight.selection(case, *original, faces, count)
    p, t, _, _, _ = preflight.subdivide(*original, faces, parents, plan)
    np.testing.assert_array_equal(p, refined[0])
    np.testing.assert_array_equal(t, refined[1])
    after, _, _ = preflight.metrics(*refined)
    assert after["quality_floor_met"] and after["within_caps"]
    assert after["minimum_quality"] == frozen["background"]["one_round"]["minimum_quality"]
    baseline, one_round = reordered_background(original, refined)
    probes = np.random.default_rng(20261111).uniform(-.5, .5, (1024, 3))
    for name, arrays in (("builtin", baseline), ("baseline", baseline), ("one-round", one_round)):
        points, cells = arrays
        directory = a.run / name
        directory.mkdir()
        for source in fixture.iterdir():
            shutil.copyfile(source, directory / source.name)
        data = dict(points=points.tolist(), tetrahedra=cells.tolist(),
                    levels=list(case.levels), probes=probes.tolist(),
                    builtin_background=name == "builtin")
        (directory / "input.json").write_text(json.dumps(data, separators=(",", ":")))
    result = dict(status="prepared", complete=True,
                  frozen_selection_and_subdivision_reproduced_exactly=True,
                  background=after, field_probe_count=len(probes),
                  source_array_hashes={name: w.digest(a.background / name)
                      for name in ("background-before.npz", "background-one-round.npz", "selection.npz")},
                  generated_input_hashes={name: w.digest(a.run / name / "input.json")
                      for name in ("builtin", "baseline", "one-round")})
    w.save(a.run / "prepare-result.json", result)


def worker(a):
    if a.worker == "prepare":
        return prepare(a)
    import numpy as np
    import pmy_neck_control as neck
    import tpms_family_fixtures as fixtures
    from comparison_mesh_audit import observe
    from periodic_trace import periodic_pairs

    directory = a.run / a.worker
    output = directory / "native-output.json"
    subprocess.run([str(a.binary), str(directory / "input.json"), str(output)], check=True)
    record = json.loads(output.read_text())
    points = record.pop("points", None)
    cells = record.pop("tetrahedra", None)
    surface = record.pop("surface", None)
    supplied_pairs = record.pop("periodic_pairs", None)
    record.pop("boundary_tags", None)
    record["native_output_sha256"] = w.digest(output)
    case = fixtures.case(NAME, a.fixtures)
    inputs = json.loads((directory / "input.json").read_text())
    probes = np.array(inputs["probes"])
    measured = np.array(record.pop("probe_values"))
    gradient = np.array(record.pop("probe_gradients"))
    expected = case.function(*probes.T)
    np.testing.assert_allclose(measured, expected, atol=3e-14, rtol=3e-14)
    np.testing.assert_allclose(gradient, neck.raw_gradient(probes), atol=1e-12, rtol=1e-12)
    record["field_identity"] = dict(probes=len(probes),
        maximum_value_difference=float(abs(measured-expected).max()),
        maximum_gradient_difference=float(abs(gradient-neck.raw_gradient(probes)).max()),
        frozen_levels=list(case.levels), passed=True)
    record.update(admitted=False, complete=False)
    if points is None:
        record["complete"] = True
        w.save(directory / "result.json", record)
        return
    points, cells, surface = np.array(points), np.array(cells, dtype=np.int64), np.array(surface, dtype=np.int64)
    mesh_path = directory / "mesh.npz"
    np.savez_compressed(mesh_path, points=points, tetrahedra=cells, surface=surface)
    record.update(nodes=len(points), elements=len(cells), mesh_sha256=w.digest(mesh_path), order=2)
    w.save(directory / "result.json", record)
    if max(len(points), len(cells)) > w.CAP:
        record.update(status="cap-rejected", complete=True)
        w.save(directory / "result.json", record)
        return
    try:
        w.imports()  # Register the original periodic observer's import directory.
        auditor = partial(periodic_pairs, declared_bounds=([-.5]*3, [.5]*3))
        record["geometry"], record["periodic"], record["quality"] = observe(
            directory, NAME, case, points, cells, surface, True, periodic_auditor=auditor)
        checked, _, _ = auditor(points, surface)
        assert all({tuple(map(int, p)) for p in supplied} == {tuple(map(int, p)) for p in independent}
                   for supplied, independent in zip(supplied_pairs, checked, strict=True))
        record["native_and_independent_pairs_agree"] = True
        record["admitted"] = record["geometry"]["geometry_admitted"]
        record["status"] = "admitted" if record["admitted"] else "admission-rejected"
        # The saved witness is only an observer, after all mesh decisions are over.
        witness = json.loads((a.observers / "result.json").read_text())["critical"]["point"]
        record["known_void_saddle_occupancy"] = neck.saved_occupancy(directory, np.array(witness))
    except Exception as error:
        record.update(status="audit-failed", error=str(error), traceback=traceback.format_exc())
    record["complete"] = True
    w.save(directory / "result.json", record)


def compare_baselines(run):
    import numpy as np

    results = [json.loads((run / name / "result.json").read_text()) for name in ("builtin", "baseline")]
    if any(not (run / name / "mesh.npz").exists() for name in ("builtin", "baseline")):
        return dict(passed=False, reason="A baseline did not return a mesh", statuses=[r["status"] for r in results])
    with np.load(run / "builtin" / "mesh.npz") as a, np.load(run / "baseline" / "mesh.npz") as b:
        exact = all(np.array_equal(a[k], b[k]) for k in ("points", "tetrahedra", "surface"))
        maximum = float(np.max(abs(a["points"]-b["points"]))) if a["points"].shape == b["points"].shape else None
    return dict(passed=exact, arrays_bitwise_identical=exact, maximum_same_index_coordinate_difference=maximum,
                reason="External BCC must reproduce ordinary BCC through the same native runner before refinement")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--binary", type=Path, default=w.ROOT / "target/release/examples/pmy_background")
    p.add_argument("--background", type=Path, default=RUNS / "macbook-pmy-global-preclip-preflight-v1")
    p.add_argument("--fixtures", type=Path, default=RUNS / "macbook-tpms-family-fixtures-v1")
    p.add_argument("--observers", type=Path, default=RUNS / "macbook-pmy-neck-control-v1")
    p.add_argument("--worker", choices=("prepare", "builtin", "baseline", "one-round"), help=argparse.SUPPRESS)
    a = p.parse_args()
    for key in ("run", "binary", "background", "fixtures", "observers"):
        setattr(a, key, getattr(a, key).resolve())
    if a.worker:
        return worker(a)
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job
    from job_parts import worker_part
    from microgen.shape import surface_functions

    a.run.mkdir(parents=True, exist_ok=False)
    files = [Path(__file__), a.binary, a.background / "result.json", a.background / "selection.npz",
             a.background / "background-before.npz", a.background / "background-one-round.npz",
             a.observers / "result.json", a.fixtures / "result.json", a.fixtures / "protocol.json",
             Path(surface_functions.__file__), *list((a.fixtures / NAME).iterdir()),
             w.ROOT / "crates/meshers-core/src/implicit.rs",
             w.ROOT / "crates/meshers-core/examples/pmy_background.rs"]
    files += [w.HERE / name for name in ("microgen_workflow.py", "job_parts.py", "pmy_global_preclip_preflight.py",
              "pmy_barycenter_control.py", "pmy_neck_control.py", "tpms_family_fixtures.py", "comparison_mesh_audit.py",
              "periodic_trace.py", "scorecard.py", "quadratic_edge_removal_probe.py", "facet_incidence_audit.py")]
    hashes = {str(path): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs"):
            shutil.copyfile(path, a.run / f"input-{index}-{path.name}")
    w.save(a.run / "manifest.json", dict(inputs=hashes, worker_timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        unchanged_gates=dict(distance=.002, absolute_relative_volume=.01, minimum_quality=.1,
            positive_curved_jacobians=True, full_periodic_caps=True, coherent_facets=True),
        builtin_options_geometry_tolerance=1.,
        timing_scope="Single diagnostic native call; excludes preclip selection, preparation and audits. No workflow speed claim.",
        native_path="Experimental externally supplied background; ordinary production default unchanged",
        physical_specimen_unchanged=True, frozen_background_replay=True, solve_performed=False))
    rows = []
    tasks = ("prepare", "builtin", "baseline", "one-round")
    def expired(*_):
        raise TimeoutError("Whole PMY control deadline exceeded")
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(len(tasks) * (TIMEOUT+5) + 30)
    with Job(parts={name: 1 for name in tasks}) as job:
        for name in tasks:
            with worker_part(job, name) as part:
                record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker", name],
                                   a.run / f"{name}-worker.log", part, timeout=TIMEOUT)
            w.save(a.run / f"{name}-worker.json", record)
            path = a.run / "prepare-result.json" if name == "prepare" else a.run / name / "result.json"
            row = json.loads(path.read_text()) if path.exists() else dict(status="no-result")
            rows.append(dict(arm=name, worker=record, result=row))
            job.result(name, status=row["status"], returncode=record["returncode"], admitted=row.get("admitted", False))
            w.save(a.run / "partial.json", dict(rows=rows, complete=False))
            if record["returncode"] != 0:
                break
            if name == "baseline":
                comparison = compare_baselines(a.run)
                w.save(a.run / "baseline-equivalence.json", comparison)
                if not comparison["passed"]:
                    break
    signal.alarm(0)
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    w.save(a.run / "result.json", dict(rows=rows, complete=True, input_hashes_verified=True,
        refined_attempted=any(r["arm"] == "one-round" for r in rows),
        admitted=any(r["arm"] == "one-round" and r["result"].get("admitted", False) for r in rows),
        scope="Controlled native clipping experiment, not a production sampling policy or stiffness result"))


if __name__ == "__main__":
    main()
