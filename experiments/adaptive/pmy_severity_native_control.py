"""One archived-native-v3 attempt on the frozen job786 severity background.

The prior776 candidate remains an immutable control. This driver generates only
one new mesh, retains all admission gates, verifies native restoration bits with
the existing standalone parser, and localizes failed forward samples afterward.
"""

import argparse
import json
from pathlib import Path
import shutil
import signal
import subprocess
import sys

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
TIMEOUT = 240


def prepare(a):
    import numpy as np
    import pmy_global_preclip_preflight as preflight

    frozen = json.loads((a.background/"result.json").read_text())
    assert frozen["complete"] and frozen["background"]["passed"] and not frozen["native_generation_performed"]
    for name, expected in frozen["frozen_decision_hashes"].items():
        assert w.digest(a.background/name) == expected
    old = json.loads((a.prior/"edge/input.json").read_text())
    with np.load(a.background/"background-edge-round.npz") as saved:
        points, cells, owners = saved["points"], saved["tetrahedra"], saved["original_tet"]
    np.testing.assert_array_equal(points[:len(old["points"])], old["points"])
    before, original_det, _ = preflight.metrics(np.asarray(old["points"]), np.asarray(old["tetrahedra"]))
    after, determinant, _ = preflight.metrics(points, cells)
    assert before == frozen["background"]["before"] and after == frozen["background"]["after"]
    error = float(np.max(abs(np.bincount(owners, weights=determinant)-original_det)/original_det))
    assert error < 1e-10 and after["quality_floor_met"] and after["within_caps"]
    directory = a.run/"candidate"
    directory.mkdir()
    for path in (a.fixtures/NAME).iterdir():
        assert w.digest(path) == frozen["inputs"][str(path)]
        shutil.copyfile(path, directory/path.name)
    data = dict(old, points=points.tolist(), tetrahedra=cells.tolist(), trace_snapping=True)
    (directory/"input.json").write_text(json.dumps(data, separators=(",", ":")))
    assert {key: value for key, value in data.items() if key not in ("points", "tetrahedra")} == {
        key: value for key, value in old.items() if key not in ("points", "tetrahedra")}
    w.save(a.run/"prepare-result.json", dict(status="prepared", complete=True,
        source_background_exact=True, old_points_exact=True, field_levels_probes_options_unchanged=True,
        before=before, after=after, maximum_parent_volume_error=error,
        generated_input_sha256=w.digest(directory/"input.json")))


def exact_trace(a):
    import numpy as np
    from pmy_snap_rollback_probe import inspect

    directory = a.run/"candidate"
    row = json.loads((directory/"result.json").read_text())
    assert w.digest(directory/"native-snap-trace.json") == row["native_snap_trace"]["sha256"]
    subprocess.run([str(a.parser_binary), str(directory/"input.json"), str(directory/"native-input-bits.json")], check=True)
    parsed = json.loads((directory/"native-input-bits.json").read_text())
    bits = np.asarray(parsed["point_bits"], dtype=np.uint64)
    points = bits.view(np.float64)
    np.testing.assert_array_equal(np.asarray(parsed["points"], dtype=np.float64).view(np.uint64), bits)
    old_proof = json.loads((a.prior_trace/"result.json").read_text())
    assert old_proof["trace_contract_verified"]
    assert w.digest(a.prior_trace/"native-input-bits.json") == old_proof["parser_output_sha256"]
    old_bits = np.asarray(json.loads((a.prior_trace/"native-input-bits.json").read_text())["point_bits"], dtype=np.uint64)
    np.testing.assert_array_equal(bits[:len(old_bits)], old_bits)
    data = json.loads((directory/"input.json").read_text())
    cells = np.asarray(data["tetrahedra"], dtype=np.int64)
    np.testing.assert_array_equal(parsed["tetrahedra"], cells)
    np.testing.assert_array_equal(np.asarray(row["options"]["levels"]).view(np.uint64), np.asarray(parsed["level_bits"], dtype=np.uint64))
    determinant = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
    trace = json.loads((directory/"native-snap-trace.json").read_text())
    assert trace == json.loads((directory/"native-output.json").read_text())["background_snap_trace"]
    original_values = np.asarray(trace[0]["original_values"], dtype=np.float64)
    initial_points, initial_values = np.asarray(trace[0]["points"]), np.asarray(trace[0]["values"])
    initial_changed = (initial_points.view(np.uint64) != bits).any(1) | (initial_values.view(np.uint64) != original_values.view(np.uint64))
    cumulative = np.zeros(len(points), dtype=bool)
    scans = []
    for state in trace:
        current, values = np.asarray(state["points"]), np.asarray(state["values"])
        np.testing.assert_array_equal(np.asarray(state["original_values"]).view(np.uint64), original_values.view(np.uint64))
        cumulative[state["restored_vertices"]] = True
        np.testing.assert_array_equal(current[cumulative].view(np.uint64), bits[cumulative])
        np.testing.assert_array_equal(values[cumulative].view(np.uint64), original_values[cumulative].view(np.uint64))
        summary, arrays = inspect(current, cells, values, state["threshold"], determinant)
        agrees = bool(np.array_equal(np.flatnonzero(arrays["bad"]), state["bad_tetrahedra"]))
        restored = initial_changed & (current.view(np.uint64) == bits).all(1) & (values.view(np.uint64) == original_values.view(np.uint64))
        scans.append(dict(topology_attempt=state["topology_attempt"], validity_round=state["validity_round"],
            native_bad_tets=len(state["bad_tetrahedra"]), explicit_restored_vertices=len(state["restored_vertices"]),
            all_observed_restored_vertices=int(restored.sum()), independent_predicate_matches=agrees, independent=summary))
    total = row.get("diagnostics", {}).get("reverted_snap_vertices")
    if total is not None:
        assert total == scans[-1]["all_observed_restored_vertices"]
    assert all(scan["independent_predicate_matches"] for scan in scans)
    assert scans[-1]["native_bad_tets"] == 0
    w.save(directory/"trace-proof.json", dict(status="verified", complete=True,
        native_input_bits_sha256=w.digest(directory/"native-input-bits.json"), trace_sha256=w.digest(directory/"native-snap-trace.json"),
        preserved_prior_input_bits=True, restored_point_and_value_bits_exact=True, numerical_tolerance_used=False,
        python_rust_original_coordinate_differences=int((np.asarray(data["points"]).view(np.uint64) != bits).sum()),
        scans=scans, native_total_reverted_vertices=total, trace_contract_verified=True,
        scope="Actual native arrays compared to authoritative parser u64 bits; geometry admission remains independent"))


def worker(a):
    from vtkmodules.vtkCommonCore import vtkSMPTools

    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    if a.worker == "prepare":
        return prepare(a)
    if a.worker == "trace":
        return exact_trace(a)
    if a.worker == "forward":
        from pmy_forward_localization import run
        return run(a)
    import pmy_native_background_control as native
    native.worker(a)


def manifest(a):
    from microgen.shape import surface_functions
    from pmy_second_edge_preflight import provenance

    native_source, files = provenance(a)
    native_archive = a.binary.parent
    receipt = json.loads((native_archive/"receipt.json").read_text())
    assert w.digest(a.binary) == receipt["binary_sha256"] == native_source["binary_sha256"]
    parser_receipt = json.loads(a.parser_receipt.read_text())
    assert w.digest(a.parser_binary) == parser_receipt["binary_sha256"]
    parser_sources = []
    for name, expected in parser_receipt["inputs"].items():
        archived = native_archive/"sources"/name
        path = archived if archived.exists() else w.ROOT/name
        assert w.digest(path) == expected
        parser_sources.append(path)
    prior_manifest = json.loads((a.prior/"manifest.json").read_text())
    prior = json.loads((a.prior/"edge/result.json").read_text())
    prepared = json.loads((a.prior/"prepare-result.json").read_text())
    assert w.digest(a.prior/"edge/input.json") == prepared["generated_input_hashes"]["edge"]
    for name, expected in (("mesh.npz", prior["mesh_sha256"]), ("native-output.json", prior["native_output_sha256"]),
                           ("native-snap-trace.json", prior["native_snap_trace"]["sha256"])):
        assert w.digest(a.prior/"edge"/name) == expected
    for path in [Path(surface_functions.__file__), w.HERE/"pmy_neck_control.py", *list((a.fixtures/NAME).iterdir())]:
        assert w.digest(path) == prior_manifest["inputs"][str(path)]
    files += [Path(__file__), a.parser_binary, a.parser_receipt, *parser_sources, Path(surface_functions.__file__),
        a.background/"result.json", a.background/"manifest.json", a.background/"background-edge-round.npz",
        a.background/"selection.npz", a.background/"frozen-schedule.json", a.background/"edge-plan.npz",
        a.prior/"manifest.json", a.prior/"prepare-result.json", a.prior/"result.json",
        a.prior_trace/"result.json", a.prior_trace/"native-input-bits.json", a.observers/"result.json",
        a.fixtures/"result.json", *list((a.fixtures/NAME).iterdir())]
    files += [a.prior/"edge"/name for name in ("input.json", "result.json", "mesh.npz", "native-output.json", "native-snap-trace.json")]
    files += [w.HERE/name for name in ("pmy_forward_localization.py", "pmy_native_background_control.py",
        "pmy_second_edge_preflight.py", "pmy_snap_rollback_probe.py", "pmy_global_preclip_preflight.py",
        "microgen_workflow.py", "job_parts.py", "pmy_neck_control.py", "tpms_family_fixtures.py",
        "comparison_mesh_audit.py", "periodic_trace.py", "scorecard.py", "quadratic_edge_removal_probe.py",
        "facet_incidence_audit.py", "nested_refinement.py", "mmg_quality_localization.py", "pmy_failure_probe.py",
        "pmy_barycenter_control.py", "box_projection.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs", ".toml", ".lock"):
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, vtk_backend="Sequential", node_cap=w.CAP, tetrahedron_cap=w.CAP,
        native_measurement_provenance=native_source, parser_receipt=parser_receipt,
        native_options=prior["options"], unchanged_gates=prior_manifest["unchanged_gates"],
        previous_control=str(a.prior/"edge"), native_mesh_attempts=1,
        input_change="Exactly job786 background arrays; prior levels/probes/flags and archived native options unchanged",
        trace_timing_scope="Diagnostic native call includes observer callback allocations and JSON-value serialization; file writing/audits additional; no speed claim",
        forward_diagnostic="Replay standard samples with identities; up to12 postmesh worst witnesses, feasible-foot upper bounds, PMY Taylor-ball exclusions and spatial-only rollback associations",
        solve_performed=False, native_default_changed=False))
    return hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--background", type=Path, default=w.HERE/"runs/macbook-pmy-severity-edge-preflight-v1")
    parser.add_argument("--prior", type=Path, default=w.HERE/"runs/macbook-pmy-rollback-native-v1")
    parser.add_argument("--prior-trace", type=Path, default=w.HERE/"runs/macbook-pmy-saved-trace-replay-v1")
    parser.add_argument("--binary", type=Path, default=w.ROOT/".local-work/macbook-pmy-background-native-v3/pmy_background")
    parser.add_argument("--parser-binary", type=Path, default=w.ROOT/".local-work/macbook-pmy-input-roundtrip-v1/pmy_input_roundtrip")
    parser.add_argument("--parser-receipt", type=Path, default=w.ROOT/".local-work/macbook-pmy-input-roundtrip-v1/receipt.json")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--observers", type=Path, default=w.HERE/"runs/macbook-pmy-neck-control-v1")
    parser.add_argument("--worker", choices=("prepare", "candidate", "trace", "forward"), help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "background", "prior", "prior_trace", "binary", "parser_binary", "parser_receipt", "fixtures", "observers"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    hashes = manifest(args)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job, result
    from job_parts import worker_part

    def expired(*_):
        raise TimeoutError("Whole PMY severity native control deadline exceeded")

    signal.signal(signal.SIGALRM, expired)
    signal.alarm(4*(TIMEOUT+5)+30)
    rows = []
    with Job(parts={"prepare": 1, "candidate": 3, "trace": 1, "forward": 2}) as job:
        for name in ("prepare", "candidate", "trace", "forward"):
            if name == "trace" and not (args.run/"candidate/native-snap-trace.json").exists():
                break
            if name == "forward":
                candidate = json.loads((args.run/"candidate/result.json").read_text())
                if candidate.get("admitted", False) or not (args.run/"candidate/mesh.npz").exists():
                    break
            with worker_part(job, name) as part:
                record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker", name],
                    args.run/f"{name}-worker.log", part, timeout=TIMEOUT)
            w.save(args.run/f"{name}-worker.json", record)
            path = args.run/"prepare-result.json" if name == "prepare" else args.run/"candidate"/(
                "trace-proof.json" if name == "trace" else "forward-localization.json" if name == "forward" else "result.json")
            row = json.loads(path.read_text()) if path.exists() else dict(status="no-result")
            rows.append(dict(arm=name, worker=record, result=row))
            job.result(name, status=row["status"], returncode=record["returncode"], admitted=row.get("admitted", False),
                       nodes=row.get("nodes"), elements=row.get("elements"), error=row.get("error"))
            w.save(args.run/"partial.json", dict(rows=rows, complete=False))
            if record["returncode"] != 0:
                break
    signal.alarm(0)
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    if (args.run/"prepare-result.json").exists():
        assert w.digest(args.run/"candidate/input.json") == json.loads((args.run/"prepare-result.json").read_text())["generated_input_sha256"]
    candidate = next((row for row in rows if row["arm"] == "candidate"), None)
    w.save(args.run/"result.json", dict(rows=rows, complete=True, input_hashes_verified=True,
        candidate_attempted=candidate is not None, admitted=bool(candidate and candidate["result"].get("admitted", False)),
        scope="One severity-background native attempt and postmesh diagnostics; no further generation, gate change, FEA or speed claim"))
    result("PMY severity native control", candidate_attempted=candidate is not None,
           admitted=bool(candidate and candidate["result"].get("admitted", False)))


if __name__ == "__main__":
    main()
