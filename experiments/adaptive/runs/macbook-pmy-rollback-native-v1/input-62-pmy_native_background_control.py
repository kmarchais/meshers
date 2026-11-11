"""Run frozen generic PMY preclip sampling through native mesh generation.

The parent submits this driver through qrun. Each preparation/generation/audit
worker has a 180-second, 8-GiB process-tree bound and standard Job progress.
The ordinary native BCC and externally supplied frozen BCC are separate controls.
By default their complete arrays must agree. An explicit fixed-triangulation
ablation requires saved evidence that post-snap cap diagonals explain the mismatch.
No admission gates are relaxed.
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


def raw_bcc_tetrahedra(points, n=16):
    """Replay clip_lattice's pre-clipping tuples, including periodic cap pivots."""
    import numpy as np

    def node(x, y, z):
        return x+(n+1)*(y+(n+1)*z)

    def rank(i):
        # Rust rounds nonnegative normalized ranks to nearest, ties upward.
        return tuple(int((0. if x == .5 else x+.5)*1e11+.5) for x in points[i])

    faces = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (3, 7, 6, 2), (0, 4, 7, 3), (1, 2, 6, 5))
    tets = []
    for z in range(n):
        for y in range(n):
            for x in range(n):
                corners = (node(x, y, z), node(x+1, y, z), node(x+1, y+1, z), node(x, y+1, z),
                           node(x, y, z+1), node(x+1, y, z+1), node(x+1, y+1, z+1), node(x, y+1, z+1))
                center = (n+1)**3+x+n*(y+n*z)
                for fi, template in enumerate(faces):
                    face = [corners[i] for i in template]
                    coordinate = (z, z, y, y, x, x)[fi]
                    positive = fi % 2 == 1
                    boundary = coordinate == (n-1 if positive else 0)
                    if boundary:
                        pivot = min(range(4), key=lambda j: rank(face[j]))
                        for k in (1, 2):
                            tets.append((center, face[pivot], face[(pivot+k) % 4], face[(pivot+k+1) % 4]))
                    elif positive:
                        neighbor = center+(n*n, n*n, n, n, 1, 1)[fi]
                        tets.extend((center, neighbor, face[j], face[(j+1) % 4]) for j in range(4))
    return np.array(tets, dtype=np.int64)


def reordered_background(original, refined, owners):
    """Restore native nodes and raw tet traversal, preserving every geometric tet."""
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
    mapped = []
    for points, cells, order in ((p0, t0, order0),
                                (p1, t1, np.r_[order0, np.arange(len(p0), len(p1))])):
        inverse = np.empty(len(order), dtype=np.int64)
        inverse[order] = np.arange(len(order))
        mapped.append((points[order], inverse[cells]))
    raw = raw_bcc_tetrahedra(mapped[0][0])
    parent_by_key = {tuple(sorted(tet)): i for i, tet in enumerate(mapped[0][1])}
    assert len(parent_by_key) == len(raw) == len(t0)
    assert {tuple(sorted(tet)) for tet in raw} == set(parent_by_key)
    # Frozen children retain their original parent's tuple positions. Reorder
    # those positions, and parent traversal, without changing any descendant.
    assert np.all(np.diff(owners) >= 0)
    stops = np.r_[0, np.cumsum(np.bincount(owners, minlength=len(t0)))]
    children = []
    for tet in raw:
        parent = parent_by_key[tuple(sorted(tet))]
        old = list(mapped[0][1][parent])
        permutation = [old.index(i) for i in tet]
        children.append(mapped[1][1][stops[parent]:stops[parent+1]][:, permutation])
    refined_tets = np.vstack(children)
    assert {tuple(sorted(tet)) for tet in refined_tets} == {tuple(sorted(tet)) for tet in mapped[1][1]}
    return ((mapped[0][0], raw), (mapped[1][0], refined_tets)), dict(
        all_original_nodes_exact=True, raw_bcc_geometric_tet_set_identical=True,
        refined_geometric_tet_set_identical=True, original_parent_tuple_order_restored=True)


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
        owners = f["original_tet"]
    faces, count, parents, _ = preflight.inventory(original[1])
    plan = preflight.selection(case, *original, faces, count)
    p, t, reproduced_owners, _, _ = preflight.subdivide(*original, faces, parents, plan)
    np.testing.assert_array_equal(p, refined[0])
    np.testing.assert_array_equal(t, refined[1])
    np.testing.assert_array_equal(reproduced_owners, owners)
    after, _, _ = preflight.metrics(*refined)
    assert after["quality_floor_met"] and after["within_caps"]
    assert after["minimum_quality"] == frozen["background"]["one_round"]["minimum_quality"]
    (baseline, one_round), ordering = reordered_background(original, refined, owners)
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
                  background=after, field_probe_count=len(probes), ordering_control=ordering,
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
    trace = record.pop("background_snap_trace", None)
    if trace is not None:
        trace_path = directory / "native-snap-trace.json"
        w.save(trace_path, trace)
        record["native_snap_trace"] = dict(path=str(trace_path), sha256=w.digest(trace_path), scans=len(trace))
        record["native_timing_scope"] = "Diagnostic trace enabled: native call includes callback allocation and trace serialization into JSON values; no speed claim"
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
    p.add_argument("--fixed-triangulation-evidence", type=Path,
                   help="Explicit ablation using a diagnosed frozen cap triangulation; keep native baseline separately")
    p.add_argument("--local-volume-guard-evidence", type=Path,
                   help="Explicit supplied-background scale correction justified by a saved determinant diagnostic")
    p.add_argument("--worker", choices=("prepare", "builtin", "baseline", "one-round"), help=argparse.SUPPRESS)
    a = p.parse_args()
    for key in ("run", "binary", "background", "fixtures", "observers"):
        setattr(a, key, getattr(a, key).resolve())
    if a.fixed_triangulation_evidence is not None:
        a.fixed_triangulation_evidence = a.fixed_triangulation_evidence.resolve()
    if a.local_volume_guard_evidence is not None:
        a.local_volume_guard_evidence = a.local_volume_guard_evidence.resolve()
        assert a.fixed_triangulation_evidence is not None
    if a.worker:
        return worker(a)
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job, result
    from job_parts import worker_part
    from microgen.shape import surface_functions

    guard_evidence = None
    guard_provenance = None
    if a.local_volume_guard_evidence is not None:
        guard_evidence = json.loads((a.local_volume_guard_evidence / "result.json").read_text())
        guard_manifest = json.loads((a.local_volume_guard_evidence / "manifest.json").read_text())
        assert guard_evidence["complete"] and guard_evidence["input_hashes_verified"] and guard_evidence["guard_unchanged"]
        guard = guard_evidence["arms"]["one-round"]
        assert guard["original_below_absolute_cutoff"] == 0
        assert guard["nonpositive_orientation"]["used_by_clipping"] == 0
        assert guard["relative_guard_failure"]["used_by_clipping"] == 0
        assert guard["positive_moderate_shrink_absolute_only"]["used_by_clipping"] > 0
        assert w.digest(a.local_volume_guard_evidence / "one-round-determinants.npz") == guard["arrays_sha256"]
        native = w.ROOT / "crates/meshers-core/src/implicit.rs"
        guard_provenance = dict(previous_implicit_sha256=guard_manifest["inputs"][str(native)],
                               current_implicit_sha256=w.digest(native),
                               retained_original_volume_fraction=.02,
                               positive_snapped_orientation_required=True,
                               ordinary_builtin_guard_unchanged=True)
        assert guard_provenance["previous_implicit_sha256"] != guard_provenance["current_implicit_sha256"]
    evidence = None
    if a.fixed_triangulation_evidence is not None:
        evidence_path = a.fixed_triangulation_evidence
        evidence = json.loads((evidence_path / "result.json").read_text())
        evidence_manifest = json.loads((evidence_path / "manifest.json").read_text())
        assert evidence["complete"] and evidence["input_hashes_verified"]
        assert evidence["all_geometric_tet_differences_are_cap_diagonal_changes"]
        assert evidence["changed_cap_pyramids"] > 0
        assert evidence["snap"]["locked_coordinates_preserved"]
        assert w.digest(evidence_path / "snap-replay.npz") == evidence["artifact_sha256"]
        for native in ("crates/meshers-core/src/implicit.rs", "crates/meshers-core/src/lib.rs",
                       "crates/meshers-core/examples/pmy_background.rs"):
            path = w.ROOT / native
            if native == "crates/meshers-core/src/implicit.rs" and guard_provenance is not None:
                assert evidence_manifest["inputs"][str(path)] == guard_provenance["previous_implicit_sha256"]
            else:
                assert w.digest(path) == evidence_manifest["inputs"][str(path)]
        diagnosed_inputs = [digest for path, digest in evidence_manifest["inputs"].items()
                            if path.endswith("/baseline/input.json")]
        assert len(diagnosed_inputs) == 1
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
    if evidence is not None:
        files += [a.fixed_triangulation_evidence / name
                  for name in ("result.json", "manifest.json", "snap-replay.npz")]
    if guard_evidence is not None:
        files += [a.local_volume_guard_evidence / name
                  for name in ("result.json", "manifest.json", "one-round-determinants.npz")]
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
        physical_specimen_unchanged=True, frozen_background_replay=True, solve_performed=False,
        baseline_control_policy="explicit-fixed-triangulation-ablation" if evidence is not None else "strict-builtin-equivalence",
        fixed_triangulation_evidence=str(a.fixed_triangulation_evidence) if evidence is not None else None,
        local_volume_guard_evidence=str(a.local_volume_guard_evidence) if guard_evidence is not None else None,
        local_volume_guard_provenance=guard_provenance,
        ablation_invariants="Frozen original vertices and geometric tets, complete physical cube, periodic caps, field and admission unchanged; refined input adds only the frozen field-selected samples"))
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
            values = dict(status=row["status"], returncode=record["returncode"], admitted=row.get("admitted", False),
                          nodes=row.get("nodes"), elements=row.get("elements"), error=row.get("error"))
            job.result(name, **values)
            w.save(a.run / "partial.json", dict(rows=rows, complete=False))
            if record["returncode"] != 0:
                break
            if name == "prepare" and evidence is not None:
                assert w.digest(a.run / "baseline/input.json") == diagnosed_inputs[0], "Fixed baseline differs from the diagnosed input"
            if name == "baseline":
                comparison = compare_baselines(a.run)
                comparison["fixed_triangulation_ablation_requested"] = evidence is not None
                comparison["diagnosis"] = str(a.fixed_triangulation_evidence) if evidence is not None else None
                # Builtin equality stays failed and visible. The separately approved
                # ablation compares this fixed background with its nested refinement.
                eligible = evidence is not None and row.get("complete", False) and row.get("field_identity", {}).get("passed", False)
                eligible = eligible and row["status"] in ("admitted", "admission-rejected")
                eligible = eligible and (a.run / "baseline/mesh.npz").exists()
                comparison["refined_ablation_eligible"] = bool(eligible)
                w.save(a.run / "baseline-equivalence.json", comparison)
                job.result("baseline-equivalence", **comparison)
                if not comparison["passed"] and not eligible:
                    break
    signal.alarm(0)
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    w.save(a.run / "result.json", dict(rows=rows, complete=True, input_hashes_verified=True,
        refined_attempted=any(r["arm"] == "one-round" for r in rows),
        admitted=any(r["arm"] == "one-round" and r["result"].get("admitted", False) for r in rows),
        fixed_triangulation_ablation=evidence is not None,
        local_volume_guard_control=guard_provenance,
        scope="Controlled native clipping experiment, not a production sampling policy or stiffness result"))
    result("experiment", refined_attempted=any(r["arm"] == "one-round" for r in rows),
           fixed_triangulation_ablation=evidence is not None)


if __name__ == "__main__":
    main()
