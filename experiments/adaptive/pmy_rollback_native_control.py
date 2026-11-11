"""One native supplied-background validity-rollback control with actual trace.

Replay the unchanged v4 input first, then process the exact job769 edge background.
Compare actual native first-snap/rollback arrays with the analytic job772 probe.
All admission gates and native options remain fixed; each worker is bounded.
"""

import argparse
import json
from pathlib import Path
import shutil
import signal
import sys

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
TIMEOUT = 240


def trace_comparison(a):
    import numpy as np
    from pmy_snap_rollback_probe import inspect

    data = json.loads((a.run/"edge/input.json").read_text())
    points, cells = np.asarray(data["points"]), np.asarray(data["tetrahedra"], dtype=np.int64)
    original_det = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
    trace = json.loads((a.run/"edge/native-snap-trace.json").read_text())
    assert trace and trace[0]["topology_attempt"] == trace[0]["validity_round"] == 0
    rows = []
    for state in trace:
        summary, arrays = inspect(np.asarray(state["points"]), cells, np.asarray(state["values"]), state["threshold"], original_det)
        rows.append(dict(topology_attempt=state["topology_attempt"], validity_round=state["validity_round"],
            restored_vertices=len(state["restored_vertices"]), native_bad_tets=len(state["bad_tetrahedra"]),
            independent_predicate_matches_native=bool(np.array_equal(np.flatnonzero(arrays["bad"]), state["bad_tetrahedra"])),
            independent=summary))
    initial_stable = next((state for state in trace if state["topology_attempt"] == 0 and not state["bad_tetrahedra"]), None)
    assert initial_stable is not None, "Initial validity rollback did not reach a recorded stable state"
    first = trace[0]
    with np.load(a.probe/"first-snap.npz") as analytic_first, np.load(a.probe/"rollback-state.npz") as analytic_final:
        np.testing.assert_array_equal(points, analytic_first["original_points"])
        np.testing.assert_array_equal(cells, analytic_first["tetrahedra"])

        def compare(left, right):
            left, right = np.asarray(left), np.asarray(right)
            return dict(bitwise_identical=bool(np.array_equal(left, right)),
                maximum_absolute_difference=float(abs(left-right).max(initial=0)),
                differing_scalars=int((left != right).sum()),
                within_1e_minus12=bool(np.allclose(left, right, atol=1e-12, rtol=0)))

        restored_ids = sorted({node for state in trace if state["topology_attempt"] == 0 for node in state["restored_vertices"]})
        comparison = dict(first_points=compare(first["points"], analytic_first["snapped_points"]),
            first_values=compare(first["values"], analytic_first["snapped_values"]),
            original_values=compare(first["original_values"], analytic_first["original_values"]),
            first_stable_points=compare(initial_stable["points"], analytic_final["points"]),
            first_stable_values=compare(initial_stable["values"], analytic_final["values"]),
            restored_vertex_ids_identical=bool(np.array_equal(restored_ids, np.flatnonzero(analytic_final["reverted"]))),
            first_failed_tet_ids_identical=bool(np.array_equal(first["bad_tetrahedra"], np.flatnonzero(analytic_first["bad"]))))
    for state in trace:
        restored = state["restored_vertices"]
        np.testing.assert_array_equal(np.asarray(state["points"])[restored], points[restored])
        np.testing.assert_array_equal(np.asarray(state["values"])[restored], np.asarray(state["original_values"])[restored])
    w.save(a.run/"trace-comparison.json", dict(status="observed", complete=True,
        trace_sha256=w.digest(a.run/"edge/native-snap-trace.json"), scans=rows,
        analytic_comparison=comparison, restored_coordinates_and_values_exact=True,
        initial_rollback_rounds=initial_stable["validity_round"],
        initial_restored_vertices=len(restored_ids), native_trace=True,
        scope="Actual native trace; analytic agreement is reported with explicit bitwise and tolerance checks. Background validity is not final mesh admission."))


def worker(a):
    import pmy_edge_native_control as previous

    if a.worker == "prepare":
        previous.prepare(a)
        path = a.run/"edge/input.json"
        data = json.loads(path.read_text())
        data["trace_snapping"] = True
        path.write_text(json.dumps(data, separators=(",", ":")))
        prepared = json.loads((a.run/"prepare-result.json").read_text())
        prepared["generated_input_hashes"]["edge"] = w.digest(path)
        prepared["edge_trace_opt_in"] = True
        w.save(a.run/"prepare-result.json", prepared)
        return
    if a.worker == "trace":
        return trace_comparison(a)
    previous.worker(a)


def manifest(a):
    from microgen.shape import surface_functions

    receipt = json.loads(a.receipt.read_text())
    old_receipt = json.loads((w.ROOT/".local-work/macbook-pmy-background-native-v2/receipt.json").read_text())
    prior = json.loads((a.source.parent/"manifest.json").read_text())
    probe = json.loads((a.probe/"result.json").read_text())
    assert probe["complete"] and probe["guard_unchanged"] and probe["final"]["guard_failures"] == 0
    assert w.digest(a.binary) == receipt["binary_sha256"] != old_receipt["binary_sha256"]
    assert prior["inputs"][str(w.ROOT/"target/release/examples/pmy_background")] == old_receipt["binary_sha256"]
    assert w.digest(w.ROOT/".local-work/macbook-pmy-background-native-v2/pmy_background") == old_receipt["binary_sha256"]
    files = [Path(__file__), a.binary, a.receipt, a.background/"result.json", a.background/"manifest.json",
        a.source.parent/"manifest.json", a.source.parent/"prepare-result.json", a.probe/"result.json", a.probe/"manifest.json",
        a.probe/"first-snap.npz", a.probe/"rollback-state.npz", a.observers/"result.json",
        w.ROOT/".local-work/macbook-pmy-background-native-v2/receipt.json",
        w.ROOT/".local-work/macbook-pmy-background-native-v2/pmy_background",
        Path(surface_functions.__file__), a.fixtures/"result.json", *list((a.fixtures/NAME).iterdir())]
    files += [a.source/name for name in ("input.json", "result.json", "native-output.json", "mesh.npz")]
    files += [a.background/name for name in ("background-edge-round.npz", "edge-plan.npz", "edge-decisions.json", "quality-failure-proof.json")]
    for name, expected in receipt["inputs"].items():
        path = w.ROOT/name
        assert w.digest(path) == expected, f"Native receipt source changed: {name}"
        files.append(path)
    for name, expected in probe["artifacts"].items():
        assert w.digest(a.probe/name) == expected
    files += [w.HERE/name for name in ("pmy_edge_native_control.py", "pmy_native_background_control.py",
        "pmy_snap_rollback_probe.py", "pmy_edge_sampling_control.py", "pmy_global_preclip_preflight.py",
        "microgen_workflow.py", "job_parts.py", "pmy_neck_control.py", "tpms_family_fixtures.py",
        "comparison_mesh_audit.py", "periodic_trace.py", "scorecard.py", "quadratic_edge_removal_probe.py",
        "facet_incidence_audit.py", "nested_refinement.py", "mmg_quality_localization.py", "pmy_barycenter_control.py")]
    for path in (Path(surface_functions.__file__), w.HERE/"pmy_neck_control.py", *list((a.fixtures/NAME).iterdir())):
        assert w.digest(path) == prior["inputs"][str(path)]
    source_result = json.loads((a.source/"result.json").read_text())
    assert w.digest(a.source/"native-output.json") == source_result["native_output_sha256"]
    assert w.digest(a.source/"mesh.npz") == source_result["mesh_sha256"]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs", ".toml", ".lock"):
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, vtk_backend="Sequential", node_cap=w.CAP, tetrahedron_cap=w.CAP,
        binary_receipt=receipt, previous_binary_receipt=old_receipt, native_options=source_result["options"],
        unchanged_gates=prior["unchanged_gates"], baseline_policy="Exact v4 native arrays and metadata excluding timing only; observer off",
        edge_input="Exact job769 background; diagnostic trace flag only, no native meshing option change",
        change="Experimental supplied-background periodic validity rollback before each clipping attempt; original signed/active-volume guards retained",
        native_trace=True, solve_performed=False, workflow_timing_claim=False,
        trace_timing_scope="Diagnostic native-call timing includes observer callback allocations and serialization into JSON values. File writing and admission audits are additional. No generation-speed comparison."))
    return hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--binary", type=Path, default=w.ROOT/".local-work/macbook-pmy-background-native-v3/pmy_background")
    parser.add_argument("--receipt", type=Path, default=w.ROOT/".local-work/macbook-pmy-background-native-v3/receipt.json")
    parser.add_argument("--background", type=Path, default=w.HERE/"runs/macbook-pmy-edge-sampling-v1")
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-native-background-v4/one-round")
    parser.add_argument("--probe", type=Path, default=w.HERE/"runs/macbook-pmy-snap-rollback-probe-v1")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--observers", type=Path, default=w.HERE/"runs/macbook-pmy-neck-control-v1")
    parser.add_argument("--worker", choices=("prepare", "baseline", "edge", "trace"), help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "binary", "receipt", "background", "source", "probe", "fixtures", "observers"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    hashes = manifest(args)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job, result
    from job_parts import worker_part
    from pmy_edge_native_control import compare_baseline

    def expired(*_):
        raise TimeoutError("Whole PMY native rollback control deadline exceeded")

    signal.signal(signal.SIGALRM, expired)
    signal.alarm(4*(TIMEOUT+5)+30)
    rows, comparison = [], None
    with Job(parts={"prepare": 1, "baseline": 2, "edge": 2, "trace": 1}) as job:
        for name in ("prepare", "baseline", "edge", "trace"):
            if name == "trace" and not (args.run/"edge/native-snap-trace.json").exists():
                break
            with worker_part(job, name) as part:
                record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker", name],
                    args.run/f"{name}-worker.log", part, timeout=TIMEOUT)
            w.save(args.run/f"{name}-worker.json", record)
            path = args.run/"prepare-result.json" if name == "prepare" else args.run/"trace-comparison.json" if name == "trace" else args.run/name/"result.json"
            row = json.loads(path.read_text()) if path.exists() else dict(status="no-result")
            rows.append(dict(arm=name, worker=record, result=row))
            job.result(name, status=row["status"], returncode=record["returncode"], admitted=row.get("admitted", False),
                       nodes=row.get("nodes"), elements=row.get("elements"), error=row.get("error"))
            w.save(args.run/"partial.json", dict(rows=rows, baseline_equivalence=comparison, complete=False))
            if record["returncode"] != 0:
                break
            if name == "baseline":
                comparison = compare_baseline(args)
                w.save(args.run/"baseline-equivalence.json", comparison)
                job.result("baseline replay", passed=comparison["passed"])
                if not comparison["passed"]:
                    break
    signal.alarm(0)
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    if (args.run/"prepare-result.json").exists():
        for name, expected in json.loads((args.run/"prepare-result.json").read_text())["generated_input_hashes"].items():
            assert w.digest(args.run/name/"input.json") == expected
    edge = next((row for row in rows if row["arm"] == "edge"), None)
    w.save(args.run/"result.json", dict(rows=rows, complete=True, input_hashes_verified=True,
        baseline_equivalence=comparison, edge_attempted=edge is not None,
        admitted=bool(edge and edge["result"].get("admitted", False)),
        scope="One native validity-rollback experiment; all final gates unchanged; no FEA or generation-speed claim"))
    result("PMY native rollback control", baseline_reproduced=bool(comparison and comparison["passed"]),
           edge_attempted=edge is not None, admitted=bool(edge and edge["result"].get("admitted", False)))


if __name__ == "__main__":
    main()
