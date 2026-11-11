"""Replay job790 unchanged and audit borrowed linear snapshots before rejection.

Exactly one native call, with diagnostic observer allocation/serialization.
No admitted mesh, gate relaxation, parameter selection, Tet10 or speed claim.
"""

import argparse
import itertools
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import microgen_workflow as w


TIMEOUT = 240
STAGES = ("clipped", "optimized", "geometry-refined", "final-linear")


def capture(a):
    import numpy as np

    directory = a.run/"candidate"
    directory.mkdir()
    shutil.copyfile(a.source/"input.json", directory/"input.json")
    old = json.loads((a.source/"native-output.json").read_text())
    assert old["status"] == "native-rejected" and "22 elements below 0.1" in old["error"]
    subprocess.run([str(a.binary), str(directory/"input.json"), str(directory/"native-output.json")], check=True)
    row = json.loads((directory/"native-output.json").read_text())
    snapshots = row.pop("background_linear_trace")
    # Test every previous numerical field, probe, rejection and complete snap
    # trace before stripping the arrays. Only the diagnostic clock is excluded.
    equal = {k: v for k, v in row.items() if k != "native_generation_s"} == {
        k: v for k, v in old.items() if k != "native_generation_s"}
    w.save(directory/"native-reproduction.json", dict(all_previous_numerical_metadata_and_snap_arrays_equal=equal,
        ignored_fields=["native_generation_s"], added_diagnostic_field="background_linear_trace",
        prior_native_output_sha256=w.digest(a.source/"native-output.json"),
        replay_native_output_sha256=w.digest(directory/"native-output.json")))
    assert equal, "Observer changed previous native arrays/metadata; stop before causal interpretation"
    assert [s["stage"] for s in snapshots] == list(STAGES)
    assert not any(k in row for k in ("points", "tetrahedra", "surface")), "Expected preserved rejection"
    saved = []
    for state in snapshots:
        stage = state["stage"]
        path = directory/f"linear-{stage}.npz"
        arrays = {name: np.asarray(state[name], dtype=np.float64 if name == "points" else np.int64)
                  for name in ("points", "tetrahedra", "surface")}
        assert arrays["tetrahedra"].shape[1] == 4 and arrays["surface"].shape[1] == 3
        # Diagnostic snapshots are retained even if a resource/admission check
        # rejects them. The native limit is unchanged; no quadratic output exists.
        np.savez_compressed(path, **arrays)
        saved.append(dict(stage=stage, path=str(path), sha256=w.digest(path), admitted=False,
            nodes=len(arrays["points"]), elements=len(arrays["tetrahedra"])))
    trace = row.pop("background_snap_trace")
    w.save(directory/"native-snap-trace.json", trace)
    row.pop("probe_values")
    row.pop("probe_gradients")
    row.update(complete=True, native_output_sha256=w.digest(directory/"native-output.json"),
        native_snap_trace=dict(path=str(directory/"native-snap-trace.json"),
            sha256=w.digest(directory/"native-snap-trace.json"), scans=len(trace)),
        native_linear_trace=saved, previous_native_output_reproduced=True,
        native_timing_scope="Diagnostic callback allocations and JSON snapshots included; no timing comparison")
    w.save(directory/"result.json", row)


def stage_audit(a):
    import numpy as np
    from scipy.spatial import cKDTree
    from cap_quality_bound import bounds, self_check
    from facet_incidence_audit import audit
    from mmg_quality_localization import quality
    from periodic_trace import periodic_pairs

    self_check()
    directory = a.run/"candidate"
    row = json.loads((directory/"result.json").read_text())
    trace = json.loads((directory/"native-snap-trace.json").read_text())
    bits = np.asarray(json.loads((directory/"native-input-bits.json").read_text())["point_bits"], dtype=np.uint64)
    original = bits.view(np.float64)
    initial_points, initial_values = np.asarray(trace[0]["points"]), np.asarray(trace[0]["values"])
    final_points, final_values = np.asarray(trace[-1]["points"]), np.asarray(trace[-1]["values"])
    original_values = np.asarray(trace[0]["original_values"])
    changed = (initial_points.view(np.uint64) != bits).any(1) | (initial_values.view(np.uint64) != original_values.view(np.uint64))
    restored = changed & (final_points.view(np.uint64) == bits).all(1) & (final_values.view(np.uint64) == original_values.view(np.uint64))
    tree = cKDTree(original[restored]) if restored.any() else None
    rows, data = [], {}
    for frozen in row["native_linear_trace"]:
        path = Path(frozen["path"])
        assert w.digest(path) == frozen["sha256"]
        stage = frozen["stage"]
        with np.load(path) as archive:
            points, cells, surface = [archive[k].copy() for k in ("points", "tetrahedra", "surface")]
        assert max(len(points), len(cells)) <= w.CAP, "Diagnostic snapshot exceeds unchanged cap"
        det = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
        q = quality(points, cells)
        bad = np.flatnonzero(q < .1)
        topology = audit(points, cells, surface)
        try:
            _, _, periodic = periodic_pairs(points, surface, declared_bounds=([-.5]*3, [.5]*3))
            periodic["passed"] = True
        except ValueError as error:
            periodic = dict(passed=False, error=str(error))
        triangles = points[surface]
        cap = np.zeros(len(surface), dtype=bool)
        for axis in range(3):
            for side in (-.5, .5):
                cap |= np.all(abs(triangles[:, :, axis]-side) <= 1e-10, axis=1)
        ceiling, _ = bounds(triangles)
        faces = {tuple(sorted(map(int, face))): i for i, face in enumerate(surface)}
        boundary = np.zeros(len(points), dtype=bool)
        boundary[np.unique(surface)] = True
        cap_nodes = np.zeros(len(points), dtype=bool)
        cap_nodes[np.unique(surface[cap])] = True
        counts = dict(with_cap_face=0, with_wall_face=0, with_obstructing_cap_face=0,
                      with_obstructing_wall_face=0, all_vertices_on_boundary=0, with_cap_vertex=0)
        examples = []
        for cell in sorted(bad, key=lambda i: q[i]):
            ids = cells[cell]
            incident = [faces[key] for f in itertools.combinations(ids, 3)
                        if (key := tuple(sorted(map(int, f)))) in faces]
            caps = [i for i in incident if cap[i]]
            walls = [i for i in incident if not cap[i]]
            counts["with_cap_face"] += bool(caps)
            counts["with_wall_face"] += bool(walls)
            counts["with_obstructing_cap_face"] += any(ceiling[i] < .1 for i in caps)
            counts["with_obstructing_wall_face"] += any(ceiling[i] < .1 for i in walls)
            counts["all_vertices_on_boundary"] += bool(boundary[ids].all())
            counts["with_cap_vertex"] += bool(cap_nodes[ids].any())
            if len(examples) < 64:
                distances = tree.query(points[ids])[0] if tree is not None else None
                examples.append(dict(element=int(cell), nodes=ids.tolist(), quality=float(q[cell]),
                    determinant=float(det[cell]), coordinates=points[ids].tolist(),
                    cap_face_quality_ceilings=[float(ceiling[i]) for i in caps],
                    wall_face_quality_ceilings=[float(ceiling[i]) for i in walls],
                    interior_vertex_count=int((~boundary[ids]).sum()),
                    nearest_restored_vertex_distances=None if distances is None else distances.tolist(),
                    coincident_with_restored_vertex=None if distances is None else (distances <= 1e-12).tolist()))
        metrics_path = directory/f"linear-{stage}-metrics.npz"
        np.savez_compressed(metrics_path, determinant=det, quality=q, bad_elements=bad,
                            surface_quality_ceiling=ceiling, cap_faces=cap)
        record = dict(stage=stage, nodes=len(points), elements=len(cells), boundary_faces=len(surface),
            minimum_quality=float(q.min()), below_quality_gate=len(bad), minimum_determinant=float(det.min()),
            nonpositive_determinants=int((det <= 0).sum()), linear_volume=float(det.sum()/6),
            topology=topology, periodic=periodic, bad_counts=counts,
            cap_faces_below_quality_ceiling=int((cap & (ceiling < .1)).sum()),
            wall_faces_below_quality_ceiling=int((~cap & (ceiling < .1)).sum()), bad_examples=examples,
            metrics_sha256=w.digest(metrics_path), admitted=False)
        rows.append(record)
        data[stage] = (points, cells, surface, q)
    transitions = []
    for before, after in zip(STAGES[:-1], STAGES[1:], strict=True):
        pa, ca, sa, qa = data[before]
        pb, cb, sb, qb = data[after]
        same = pa.shape == pb.shape and np.array_equal(ca, cb) and np.array_equal(sa, sb)
        record = dict(before=before, after=after, complete_indexed_topology_equal=same,
                      before_bad=int((qa < .1).sum()), after_bad=int((qb < .1).sum()))
        if same:
            record.update(points_bitwise_unchanged=bool(np.array_equal(pa.view(np.uint64), pb.view(np.uint64))),
                maximum_node_move=float(np.linalg.norm(pb-pa, axis=1).max()),
                newly_below_gate=int(((qa >= .1) & (qb < .1)).sum()),
                recovered_above_gate=int(((qa < .1) & (qb >= .1)).sum()),
                persisted_below_gate=int(((qa < .1) & (qb < .1)).sum()))
        transitions.append(record)
    match = re.search(r"minimum MMG quality ([^ ]+).*\((\d+) elements below 0.1\)", row["error"])
    assert match and rows[-1]["below_quality_gate"] == int(match.group(2))
    np.testing.assert_allclose(rows[-1]["minimum_quality"], float(match.group(1)), atol=1e-12, rtol=0)
    w.save(directory/"linear-audit.json", dict(status="diagnosed", complete=True, stages=rows,
        transitions=transitions, unchanged_rejection_reproduced=True, final_quality_error_reproduced=True,
        restored_background_vertices=int(restored.sum()), native_snapshot_count=len(rows),
        spatial_relation_scope="Nearest/coincident coordinates only; mesh indices do not establish original background vertex lineage",
        fixed_face_bound_scope="Ceiling holds when that triangle stays fixed and the apex alone moves; it does not bound changes to boundary vertices or connectivity",
        scope="Rejected Tet4 intermediate states only. No quadratic construction, geometry admission, FEA or performance comparison"))


def provenance(a):
    files = [Path(__file__), a.source.parent/"manifest.json", a.prior_trace/"result.json",
             a.prior_trace/"native-input-bits.json", a.parser_binary, a.parser_receipt]
    files += [a.source/name for name in ("input.json", "result.json", "native-output.json", "native-snap-trace.json", "trace-proof.json")]
    prior = json.loads((a.source/"result.json").read_text())
    old_manifest = json.loads((a.source.parent/"manifest.json").read_text())
    assert w.digest(a.source/"native-output.json") == prior["native_output_sha256"]
    assert w.digest(a.source/"native-snap-trace.json") == prior["native_snap_trace"]["sha256"]
    prepared = a.source.parent/"prepare-result.json"
    assert w.digest(a.source/"input.json") == json.loads(prepared.read_text())["generated_input_sha256"]
    files.append(prepared)
    sources = {}
    for label, archive in (("native", a.binary.parent), ("prior-native", a.prior_binary.parent)):
        receipt = json.loads((archive/"receipt.json").read_text())
        binary = a.binary if label == "native" else a.prior_binary
        assert w.digest(binary) == receipt["binary_sha256"]
        files += [binary, archive/"receipt.json"]
        differences = []
        for name, expected in receipt["inputs"].items():
            path = archive/"sources"/name
            assert w.digest(path) == expected
            files.append(path)
            live = w.ROOT/name
            if not live.exists() or w.digest(live) != expected:
                differences.append(name)
        sources[label] = dict(binary_sha256=receipt["binary_sha256"], archived_sources_verified=True,
                              archived_source_count=len(receipt["inputs"]), live_source_differences=differences)
    receipt = json.loads(a.parser_receipt.read_text())
    assert w.digest(a.parser_binary) == receipt["binary_sha256"]
    for name, expected in receipt["inputs"].items():
        path = a.prior_binary.parent/"sources"/name
        if not path.exists():
            path = w.ROOT/name
        assert w.digest(path) == expected
        files.append(path)
    files += [w.HERE/name for name in ("pmy_severity_native_control.py", "pmy_snap_rollback_probe.py",
        "mmg_quality_localization.py", "cap_quality_bound.py", "facet_incidence_audit.py",
        "periodic_trace.py", "microgen_workflow.py", "job_parts.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs", ".toml", ".lock"):
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, sources=sources, worker_timeout_s=TIMEOUT,
        memory_bytes=w.MEMORY, numerical_threads=w.NUMERICAL_THREAD_POLICY, vtk_backend="Sequential",
        native_options=prior["options"], unchanged_gates=old_manifest["unchanged_gates"],
        native_calls=1, input_policy="Byte-identical rejected job790 candidate input",
        rejection_policy="Every old native field and snap array must match, excluding diagnostic elapsed time",
        timing_scope="Diagnostic trace allocation and serialization; no speed claim", admitted=False))
    return hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-severity-native-v1/candidate")
    parser.add_argument("--binary", type=Path, default=w.ROOT/".local-work/macbook-pmy-background-native-v4/pmy_background")
    parser.add_argument("--prior-binary", type=Path, default=w.ROOT/".local-work/macbook-pmy-background-native-v3/pmy_background")
    parser.add_argument("--prior-trace", type=Path, default=w.HERE/"runs/macbook-pmy-saved-trace-replay-v1")
    parser.add_argument("--parser-binary", type=Path, default=w.ROOT/".local-work/macbook-pmy-input-roundtrip-v1/pmy_input_roundtrip")
    parser.add_argument("--parser-receipt", type=Path, default=w.ROOT/".local-work/macbook-pmy-input-roundtrip-v1/receipt.json")
    parser.add_argument("--worker", choices=("capture", "trace", "audit"), help=argparse.SUPPRESS)
    a = parser.parse_args()
    for key, value in vars(a).items():
        if isinstance(value, Path):
            setattr(a, key, value.resolve())
    if a.worker:
        from vtkmodules.vtkCommonCore import vtkSMPTools
        assert vtkSMPTools.SetBackend("Sequential")
        vtkSMPTools.Initialize(1)
        if a.worker == "capture":
            return capture(a)
        if a.worker == "audit":
            return stage_audit(a)
        from pmy_severity_native_control import exact_trace
        return exact_trace(a)
    a.run.mkdir(parents=True, exist_ok=False)
    hashes = provenance(a)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job, result
    from job_parts import worker_part

    with Job(parts={"capture": 2, "trace": 1, "audit": 2}) as job:
        for stage in ("capture", "trace", "audit"):
            with worker_part(job, stage) as part:
                record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker", stage],
                    a.run/f"{stage}-worker.log", part, timeout=TIMEOUT)
            w.save(a.run/f"{stage}-worker.json", record)
            job.result(stage, returncode=record["returncode"], seconds=record["wall_s"],
                       peak_gib=record["peak_rss_bytes"]/1024**3)
            assert record["returncode"] == 0, f"Read {stage}-worker.log"
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    assert w.digest(a.run/"candidate/input.json") == w.digest(a.source/"input.json")
    findings = json.loads((a.run/"candidate/linear-audit.json").read_text())
    proof = json.loads((a.run/"candidate/trace-proof.json").read_text())
    findings.update(input_hashes_verified=True, native_snap_trace_verified=proof["trace_contract_verified"], admitted=False)
    w.save(a.run/"result.json", findings)
    for stage in findings["stages"]:
        result("PMY linear "+stage["stage"], nodes=stage["nodes"], elements=stage["elements"],
               minimum_q=stage["minimum_quality"], below_q=stage["below_quality_gate"],
               boundary_face_obstructions=stage["cap_faces_below_quality_ceiling"]+stage["wall_faces_below_quality_ceiling"])


if __name__ == "__main__":
    main()
