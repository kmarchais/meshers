"""Native clipping of one frozen PMY edge-star background after exact replay.

The archived local-volume-guard binary must reproduce the old one-round native
arrays and metadata, excluding elapsed time, before receiving the new background.
No field, native options, fixture, admission gate or snapping guard is changed.
Each worker has a 240-second, 8-GiB bound. The parent owns queue submission.
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


def prepare(a):
    import numpy as np
    import pmy_global_preclip_preflight as preflight

    frozen = json.loads((a.background/"result.json").read_text())
    assert frozen["complete"] and frozen["input_hashes_verified"] and frozen["background"]["passed"]
    assert not frozen["native_clipping_performed"]
    for name, digest in frozen["frozen_decision_hashes"].items():
        assert w.digest(a.background/name) == digest
    source = a.source/"input.json"
    assert w.digest(source) == frozen["inputs"][str(source)]
    old = json.loads(source.read_text())
    assert old["builtin_background"] is False
    with np.load(a.background/"background-edge-round.npz") as saved:
        points, cells, owners = saved["points"], saved["tetrahedra"], saved["original_tet"]
    np.testing.assert_array_equal(points[:len(old["points"])], old["points"])
    before, old_det, _ = preflight.metrics(np.asarray(old["points"]), np.asarray(old["tetrahedra"]))
    after, new_det, _ = preflight.metrics(points, cells)
    assert after == frozen["background"]["after"]
    assert after["quality_floor_met"] and after["within_caps"]
    parent_error = float(np.max(abs(np.bincount(owners, weights=new_det)-old_det)/old_det))
    assert parent_error < 1e-10
    fixture = a.fixtures/NAME
    for name in ("baseline", "edge"):
        directory = a.run/name
        directory.mkdir()
        for path in fixture.iterdir():
            assert w.digest(path) == frozen["inputs"][str(path)]
            shutil.copyfile(path, directory/path.name)
        if name == "baseline":
            shutil.copyfile(source, directory/"input.json")
        else:
            data = dict(old, points=points.tolist(), tetrahedra=cells.tolist())
            (directory/"input.json").write_text(json.dumps(data, separators=(",", ":")))
            assert {key: value for key, value in data.items() if key not in ("points", "tetrahedra")} == {
                key: value for key, value in old.items() if key not in ("points", "tetrahedra")}
    assert w.digest(a.run/"baseline/input.json") == w.digest(source)
    edge = json.loads((a.run/"edge/input.json").read_text())
    np.testing.assert_array_equal(edge["points"], points)
    np.testing.assert_array_equal(edge["tetrahedra"], cells)
    w.save(a.run/"prepare-result.json", dict(status="prepared", complete=True,
        baseline_input_byte_identical=True, edge_input_arrays_exact=True,
        original_points_bitwise_preserved=True, original_background=before, edge_background=after,
        maximum_parent_volume_error=parent_error, options_levels_probes_unchanged=True,
        generated_input_hashes={name: w.digest(a.run/name/"input.json") for name in ("baseline", "edge")}))


def worker(a):
    # VTK otherwise starts a second thread pool during distance observations.
    from vtkmodules.vtkCommonCore import vtkSMPTools
    vtkSMPTools.SetBackend("Sequential")
    assert vtkSMPTools.GetBackend() == "Sequential"
    if a.worker == "prepare":
        return prepare(a)
    import pmy_native_background_control as native
    native.worker(a)
    path = a.run/a.worker/"result.json"
    row = json.loads(path.read_text())
    row["observer_thread_policy"] = dict(vtk_backend="Sequential", numerical_threads=w.NUMERICAL_THREAD_POLICY)
    w.save(path, row)


def compare_baseline(a):
    import numpy as np

    previous, current = a.source, a.run/"baseline"
    required = ("native-output.json", "mesh.npz", "result.json")
    if not all((current/name).exists() for name in required):
        return dict(passed=False, reason="Replay did not archive a complete native mesh", missing=[name for name in required if not (current/name).exists()])
    old = json.loads((previous/"native-output.json").read_text())
    new = json.loads((current/"native-output.json").read_text())
    excluded = ["native_generation_s"]
    array_fields = ["points", "tetrahedra", "surface"]
    arrays = {}
    with np.load(previous/"mesh.npz") as left, np.load(current/"mesh.npz") as right:
        for key in array_fields:
            arrays[key] = bool(np.array_equal(left[key], right[key]))
            assert np.array_equal(left[key], np.asarray(old[key]))
            assert np.array_equal(right[key], np.asarray(new[key]))
    compared = sorted((set(old) | set(new))-set(excluded)-set(array_fields))
    matches = {key: key in old and key in new and old[key] == new[key] for key in compared}
    previous_result = json.loads((previous/"result.json").read_text())
    current_result = json.loads((current/"result.json").read_text())
    # Audit results are independently recomputed; record any observer differences.
    audit_keys = ("field_identity", "geometry", "periodic", "quality", "native_and_independent_pairs_agree")
    observer_matches = {key: previous_result.get(key) == current_result.get(key) for key in audit_keys}
    eligible = current_result.get("complete", False) and current_result.get("field_identity", {}).get("passed", False)
    eligible = eligible and current_result.get("status") in ("admitted", "admission-rejected")
    return dict(passed=bool(all(arrays.values()) and all(matches.values()) and eligible),
        arrays_bitwise_identical=arrays, native_metadata_exact=matches,
        failed_native_fields=[key for key, value in matches.items() if not value],
        excluded_fields=excluded, observer_results_exact=observer_matches,
        replay_status=current_result.get("status"), old_generation_s=old["native_generation_s"],
        replay_generation_s=new["native_generation_s"],
        scope="Exact old-input native replay; elapsed time alone excluded from native equality")


def manifest(a):
    from microgen.shape import surface_functions

    receipt = json.loads(a.receipt.read_text())
    source_manifest = json.loads((a.source.parent/"manifest.json").read_text())
    source_result = json.loads((a.source/"result.json").read_text())
    old_binary = str(w.ROOT/"target/release/examples/pmy_background")
    assert w.digest(a.binary) == receipt["binary_sha256"] == source_manifest["inputs"][old_binary]
    assert source_result["complete"] and source_result["status"] == "admission-rejected"
    assert w.digest(a.source/"native-output.json") == source_result["native_output_sha256"]
    assert w.digest(a.source/"mesh.npz") == source_result["mesh_sha256"]
    prepared = json.loads((a.source.parent/"prepare-result.json").read_text())
    assert w.digest(a.source/"input.json") == prepared["generated_input_hashes"]["one-round"]
    assert source_manifest["local_volume_guard_provenance"]["retained_original_volume_fraction"] == .02
    assert source_manifest["local_volume_guard_provenance"]["positive_snapped_orientation_required"]
    files = [Path(__file__), a.binary, a.receipt, a.background/"result.json", a.background/"manifest.json",
        a.source.parent/"manifest.json", a.source.parent/"prepare-result.json",
        a.fixtures/"result.json", a.fixtures/"protocol.json", a.observers/"result.json",
        Path(surface_functions.__file__), *list((a.fixtures/NAME).iterdir())]
    files += [a.source/name for name in ("input.json", "result.json", "native-output.json", "mesh.npz")]
    files += [a.background/name for name in ("background-edge-round.npz", "edge-plan.npz", "edge-decisions.json", "quality-failure-proof.json")]
    for name, expected in receipt["inputs"].items():
        path = w.ROOT/name
        assert w.digest(path) == expected, f"Build source differs from archived receipt: {name}"
        files.append(path)
    files += [w.HERE/name for name in ("pmy_native_background_control.py", "pmy_edge_sampling_control.py",
        "microgen_workflow.py", "job_parts.py", "pmy_global_preclip_preflight.py", "pmy_neck_control.py",
        "tpms_family_fixtures.py", "comparison_mesh_audit.py", "periodic_trace.py", "scorecard.py",
        "quadratic_edge_removal_probe.py", "facet_incidence_audit.py", "nested_refinement.py",
        "mmg_quality_localization.py", "pmy_barycenter_control.py")]
    # All field and frozen specimen sources must agree with the previous native run.
    for path in [Path(surface_functions.__file__), w.HERE/"pmy_neck_control.py", *list((a.fixtures/NAME).iterdir())]:
        assert w.digest(path) == source_manifest["inputs"][str(path)]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs", ".toml", ".lock"):
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, vtk_backend="Sequential",
        node_cap=w.CAP, tetrahedron_cap=w.CAP, unchanged_gates=source_manifest["unchanged_gates"],
        binary_receipt=receipt, binary_build_job=receipt["build_job"],
        native_options=source_result["options"], local_volume_guard=source_manifest["local_volume_guard_provenance"],
        source_native_run=str(a.source), source_background=str(a.background),
        baseline_policy="Replay byte-identical old one-round input; require bitwise arrays and exact native metadata except elapsed seconds",
        ablation="Only the saved edge-split background points/tets differ between the two supplied inputs",
        native_path="Archived experimental supplied-background v2 binary; production defaults unchanged",
        timing_scope="Diagnostic native-call times only; preparation, background selection and audit excluded; no workflow speed claim",
        candidate_archival="Native JSON then mesh.npz saved before count and geometry admission gates",
        solve_performed=False, no_further_refinement=True))
    return hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--binary", type=Path, default=w.ROOT/".local-work/macbook-pmy-background-native-v2/pmy_background")
    parser.add_argument("--receipt", type=Path, default=w.ROOT/".local-work/macbook-pmy-background-native-v2/receipt.json")
    parser.add_argument("--background", type=Path, default=w.HERE/"runs/macbook-pmy-edge-sampling-v1")
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-native-background-v4/one-round")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--observers", type=Path, default=w.HERE/"runs/macbook-pmy-neck-control-v1")
    parser.add_argument("--worker", choices=("prepare", "baseline", "edge"), help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "binary", "receipt", "background", "source", "fixtures", "observers"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    hashes = manifest(args)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job, result
    from job_parts import worker_part

    rows, comparison = [], None

    def expired(*_):
        raise TimeoutError("Whole PMY native edge control deadline exceeded")

    signal.signal(signal.SIGALRM, expired)
    signal.alarm(3*(TIMEOUT+5)+30)
    with Job(parts={"prepare": 1, "baseline": 1, "edge": 1}) as job:
        for name in ("prepare", "baseline", "edge"):
            with worker_part(job, name) as part:
                record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker", name],
                    args.run/f"{name}-worker.log", part, timeout=TIMEOUT)
            w.save(args.run/f"{name}-worker.json", record)
            path = args.run/"prepare-result.json" if name == "prepare" else args.run/name/"result.json"
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
                job.result("baseline replay", passed=comparison["passed"], status=row["status"])
                if not comparison["passed"]:
                    break
    signal.alarm(0)
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    prepared = args.run/"prepare-result.json"
    if prepared.exists():
        for name, expected in json.loads(prepared.read_text())["generated_input_hashes"].items():
            assert w.digest(args.run/name/"input.json") == expected
    edge = next((row for row in rows if row["arm"] == "edge"), None)
    w.save(args.run/"result.json", dict(rows=rows, complete=True, input_hashes_verified=True,
        baseline_equivalence=comparison, edge_attempted=edge is not None,
        admitted=bool(edge and edge["result"].get("admitted", False)),
        scope="One frozen-background native clipping control with unchanged geometry, quality, periodic and resource gates; no FEA or workflow timing claim"))
    result("PMY native edge control", baseline_reproduced=bool(comparison and comparison["passed"]),
           edge_attempted=edge is not None, admitted=bool(edge and edge["result"].get("admitted", False)))


if __name__ == "__main__":
    main()
