"""Resolve saved PMY trace restoration against exact native input f64 bits.

Run only a standalone serde parser, then inspect the immutable job776 trace.
No meshing, snapping, production code change, or rewritten old result is involved.
Restoration requires exact u64 equality, not a floating-point tolerance.
"""

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from types import SimpleNamespace

import microgen_workflow as w


TIMEOUT = 180


def worker(a):
    import numpy as np
    from pmy_snap_rollback_probe import inspect

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    receipt = json.loads(a.receipt.read_text())
    native_receipt_path = w.ROOT/".local-work/macbook-pmy-background-native-v3/receipt.json"
    native_receipt = json.loads(native_receipt_path.read_text())
    assert w.digest(a.binary) == receipt["binary_sha256"]
    for name, expected in receipt["inputs"].items():
        assert w.digest(w.ROOT/name) == expected
        if name in native_receipt["inputs"]:
            assert expected == native_receipt["inputs"][name]
    for required in ("Cargo.lock", "Cargo.toml", "crates/meshers-core/Cargo.toml", "crates/meshers-core/examples/pmy_background.rs"):
        assert required in receipt["inputs"] and receipt["inputs"][required] == native_receipt["inputs"][required]
    rust_parser = w.ROOT/"crates/meshers-core/examples/pmy_input_roundtrip.rs"
    rust_native = w.ROOT/"crates/meshers-core/examples/pmy_background.rs"

    def declaration(path):
        source = path.read_text().split("struct Input {", 1)[1].split("}", 1)[0]
        return re.sub(r"\s+", "", source)

    assert declaration(rust_parser) == declaration(rust_native), "Native input field types differ"
    source = a.source/"edge/input.json"
    edge = json.loads((a.source/"edge/result.json").read_text())
    prepared = json.loads((a.source/"prepare-result.json").read_text())
    assert w.digest(source) == prepared["generated_input_hashes"]["edge"]
    trace_path = a.source/"edge/native-snap-trace.json"
    assert w.digest(trace_path) == edge["native_snap_trace"]["sha256"]
    assert w.digest(a.source/"edge/native-output.json") == edge["native_output_sha256"]
    analytic = json.loads((a.probe/"result.json").read_text())
    for name, expected in analytic["artifacts"].items():
        assert w.digest(a.probe/name) == expected
    files = [Path(__file__), a.binary, a.receipt, native_receipt_path, rust_parser, rust_native,
        source, trace_path, a.source/"edge/result.json", a.source/"edge/native-output.json",
        a.source/"prepare-result.json", a.source/"trace-worker.json", a.source/"trace-worker.log",
        a.probe/"result.json", a.probe/"first-snap.npz", a.probe/"rollback-state.npz",
        w.HERE/"pmy_snap_rollback_probe.py", w.HERE/"microgen_workflow.py", w.HERE/"mmg_quality_localization.py"]
    files += [w.ROOT/name for name in receipt["inputs"]]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs", ".toml", ".lock"):
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, parser_binary_receipt=receipt,
        native_binary_receipt_sha256=w.digest(native_receipt_path), identical_typed_input_declaration=True,
        reader="Same serde_json::from_reader(BufReader::new(File::open(input))) typed path",
        comparison="u64 to_bits() integers from native parser versus saved trace coordinates parsed by Python",
        no_remeshing=True, no_old_result_changes=True, floating_tolerance_not_used_for_restoration=True))
    result = dict(inputs=hashes, native_trace=True, no_remeshing=True, candidate_generated=False, admission_claim=False)
    with Job(parts={"native-input-parse": 1, "saved-trace-proof": 2, "analytic-comparison": 1}) as job:
        with job.part("native-input-parse") as part:
            part.update(0, "Recovering exact native input floating-point bits", points=edge["nodes"])
            subprocess.run([str(a.binary), str(source), str(a.run/"native-input-bits.json")], check=True)
            native = json.loads((a.run/"native-input-bits.json").read_text())
            data = json.loads(source.read_text())
            bits = np.asarray(native["point_bits"], dtype=np.uint64)
            points = bits.view(np.float64)
            python_points = np.asarray(data["points"], dtype=np.float64)
            cells = np.asarray(data["tetrahedra"], dtype=np.int64)
            np.testing.assert_array_equal(np.asarray(native["points"], dtype=np.float64).view(np.uint64), bits)
            np.testing.assert_array_equal(native["tetrahedra"], cells)
            np.testing.assert_array_equal(np.asarray(edge["options"]["levels"], dtype=np.float64).view(np.uint64), native["level_bits"])
            mismatch = python_points.view(np.uint64) != bits
            result["native_input_parse"] = dict(nodes=len(points), tetrahedra=len(cells),
                python_native_different_scalars=int(mismatch.sum()),
                python_native_different_nodes=int(mismatch.any(1).sum()),
                maximum_coordinate_difference=float(abs(points-python_points).max()),
                native_output_float_encoding_matches_authoritative_bits=True,
                native_serialize_reparse_identical=bool(np.array_equal(native["serialize_reparse_point_bits"], bits)),
                example_differences=[dict(node=int(i), axis=int(j), python=float(python_points[i, j]),
                    rust=float(points[i, j]), python_bits=int(python_points.view(np.uint64)[i, j]), rust_bits=int(bits[i, j]))
                    for i, j in np.argwhere(mismatch)[:12]])
            part.update(1, "Verified parser output against integer bit encodings", different_scalars=int(mismatch.sum()), nodes=len(points))
        with job.part("saved-trace-proof") as part:
            trace = json.loads(trace_path.read_text())
            raw = json.loads((a.source/"edge/native-output.json").read_text())["background_snap_trace"]
            assert trace == raw, "Saved trace sidecar differs from immutable raw native output"
            original_det = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
            original_values = np.asarray(trace[0]["original_values"], dtype=np.float64)
            cumulative = np.zeros(len(points), dtype=bool)
            scans, old_comparison_mismatches = [], 0
            for index, state in enumerate(trace):
                current = np.asarray(state["points"], dtype=np.float64)
                values = np.asarray(state["values"], dtype=np.float64)
                np.testing.assert_array_equal(np.asarray(state["original_values"]).view(np.uint64), original_values.view(np.uint64))
                restored = np.asarray(state["restored_vertices"], dtype=np.int64)
                cumulative[restored] = True
                # Exact recovery proof uses native input bits, with no tolerance.
                np.testing.assert_array_equal(current[cumulative].view(np.uint64), bits[cumulative])
                np.testing.assert_array_equal(values[cumulative].view(np.uint64), original_values[cumulative].view(np.uint64))
                old_comparison_mismatches += int((current[restored].view(np.uint64) != python_points[restored].view(np.uint64)).sum())
                summary, arrays = inspect(current, cells, values, state["threshold"], original_det)
                predicate_match = bool(np.array_equal(np.flatnonzero(arrays["bad"]), state["bad_tetrahedra"]))
                scans.append(dict(topology_attempt=state["topology_attempt"], validity_round=state["validity_round"],
                    restored_vertices=len(restored), cumulative_restored=int(cumulative.sum()),
                    native_bad_tetrahedra=len(state["bad_tetrahedra"]), independent_predicate_matches_native=predicate_match,
                    independent=summary, restored_coordinate_bits_exact=True, restored_field_value_bits_exact=True))
                part.update((index+1)/len(trace), "Checked native restoration with integer bit equality", scans=index+1, restored=int(cumulative.sum()), predicate_match=predicate_match)
            result.update(scans=scans, restored_vertices=int(cumulative.sum()),
                all_restored_coordinate_bits_exact=True, all_restored_field_value_bits_exact=True,
                old_python_comparison_mismatches_reproduced=old_comparison_mismatches,
                all_native_failure_predicates_reproduced=all(row["independent_predicate_matches_native"] for row in scans))
        with job.part("analytic-comparison") as part:
            first = trace[0]
            stable = next(state for state in trace if state["topology_attempt"] == 0 and not state["bad_tetrahedra"])
            restored_ids = sorted({node for state in trace if state["topology_attempt"] == 0 for node in state["restored_vertices"]})

            def compare(left, right):
                left, right = np.asarray(left, dtype=np.float64), np.asarray(right, dtype=np.float64)
                return dict(bitwise_identical=bool(np.array_equal(left.view(np.uint64), right.view(np.uint64))),
                    maximum_absolute_difference=float(abs(left-right).max(initial=0)),
                    differing_scalars=int((left.view(np.uint64) != right.view(np.uint64)).sum()),
                    within_1e_minus12=bool(np.allclose(left, right, atol=1e-12, rtol=0)))

            with np.load(a.probe/"first-snap.npz") as analytic_first, np.load(a.probe/"rollback-state.npz") as analytic_final:
                np.testing.assert_array_equal(python_points, analytic_first["original_points"])
                np.testing.assert_array_equal(cells, analytic_first["tetrahedra"])
                result["analytic_comparison"] = dict(
                    first_points=compare(first["points"], analytic_first["snapped_points"]),
                    first_values=compare(first["values"], analytic_first["snapped_values"]),
                    original_values=compare(first["original_values"], analytic_first["original_values"]),
                    first_stable_points=compare(stable["points"], analytic_final["points"]),
                    first_stable_values=compare(stable["values"], analytic_final["values"]),
                    restored_vertex_ids_identical=bool(np.array_equal(restored_ids, np.flatnonzero(analytic_final["reverted"]))),
                    first_failed_tet_ids_identical=bool(np.array_equal(first["bad_tetrahedra"], np.flatnonzero(analytic_first["bad"]))))
            result["interpretation"] = "The exact restoration contract is checked against the f64 values actually parsed by Rust. Analytic arrays remain separately compared and are not promoted to bitwise native traces. No saved failed result was changed."
            part.update(1, "Separated exact native restoration from analytic floating-point agreement", failed_ids_match=result["analytic_comparison"]["first_failed_tet_ids_identical"], restored_ids_match=result["analytic_comparison"]["restored_vertex_ids_identical"])
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    result.update(complete=True, input_hashes_verified=True, parser_output_sha256=w.digest(a.run/"native-input-bits.json"),
        trace_contract_verified=bool(result["all_native_failure_predicates_reproduced"] and scans[-1]["native_bad_tetrahedra"] == 0))
    w.save(a.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-rollback-native-v1")
    parser.add_argument("--probe", type=Path, default=w.HERE/"runs/macbook-pmy-snap-rollback-probe-v1")
    parser.add_argument("--binary", type=Path, default=w.ROOT/".local-work/macbook-pmy-input-roundtrip-v1/pmy_input_roundtrip")
    parser.add_argument("--receipt", type=Path, default=w.ROOT/".local-work/macbook-pmy-input-roundtrip-v1/receipt.json")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "source", "probe", "binary", "receipt"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-pmy_saved_trace_replay.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **kw: None), timeout=TIMEOUT)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    findings = json.loads((args.run/"result.json").read_text())
    result("PMY native input parse", different_scalars=findings["native_input_parse"]["python_native_different_scalars"],
           maximum_difference=findings["native_input_parse"]["maximum_coordinate_difference"])
    result("PMY exact saved trace proof", restored=findings["restored_vertices"],
           original_failure_reproduced=findings["old_python_comparison_mismatches_reproduced"],
           trace_contract_verified=findings["trace_contract_verified"])


if __name__ == "__main__":
    main()
