"""Capture the rejected BCC linear mesh with the original field and options.

Run through qrun --check-start --progress report, with --interrupt none:
new output directories make an interrupted run unsafe to restart in place.
Three serial workers each have a 120 s process-tree timeout and 8 GiB RSS cap.
The baseline uses the unchanged installed binary; ordinary and observed use the
isolated package from cad_observer_build.py. No meshing options are altered.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
from unittest.mock import patch

import microgen_workflow as w
from job_parts import worker_part

SOURCE = w.HERE/"runs/macbook-cad-bcc-verified-v1"
NAME = "cad-bcc-flat-struts-radius01"


def verify_inputs(manifest):
    for path, expected in manifest["inputs"].items():
        assert w.digest(Path(path)) == expected, f"Changed input: {path}"


def worker(args):
    if args.mode != "baseline":
        sys.path.insert(0, str(args.package))
    import numpy as np
    import meshers
    import cad_lattice_fixture as fixture
    import scorecard

    manifest = json.loads((args.run/"manifest.json").read_text())
    verify_inputs(manifest)
    old = json.loads((args.source/"result.json").read_text())
    expected_binary = (old["native_binary_sha256"] if args.mode == "baseline"
                       else manifest["new_native_binary_sha256"])
    assert w.digest(Path(meshers._meshers.__file__)) == expected_binary
    if args.mode != "baseline":
        assert Path(meshers.__file__).resolve().is_relative_to(args.package)
        assert getattr(meshers._meshers, "_LINEAR_OBSERVATIONS", False)
    for name in ("__init__.py", "compiler.py", "expression.py", "intersection.py"):
        source_path = w.ROOT/"crates/meshers-python/python/meshers"/name
        assert w.digest(Path(meshers.__file__).parent/name) == old["inputs"][str(source_path)]
    recipe = json.loads((args.run/"primitive-receipt.json").read_text())
    field = fixture.field(recipe["normalized_primitives"])
    with np.load(args.run/"membership-probes.npz") as frozen:
        queries = frozen["points"]
        np.testing.assert_array_equal(np.maximum(field(*queries.T), abs(queries).max(1)-.5),
                                      frozen["normalized_field"])
    original = meshers._meshers.generate
    observations = []

    def capture(*positional, **keywords):
        return original(*positional, **keywords, _observations=observations)

    result = dict(mode=args.mode, complete=False, native_admitted=False,
        native_binary_sha256=expected_binary, native_options=old["native_options"],
        comparative_timing_claim=False, finite_element_solve_performed=False,
        status="capture-control-running", snapshots=[])
    destination = args.run/args.mode
    w.save(destination/"result.json", result)
    try:
        if args.mode == "observed":
            # Scoped to this single-thread worker. The public function still owns
            # tracing, compilation, option conversion and native dispatch.
            with patch.object(meshers._meshers, "generate", capture):
                meshers.generate(field, **old["native_options"])
        else:
            meshers.generate(field, **old["native_options"])
    except Exception as error:
        result.update(error_type=type(error).__name__, error=str(error))
    else:
        result.update(error_type=None, error=None)
    for index, row in enumerate(observations):
        arrays = {key: row[key] for key in ("points", "tetrahedra", "surface")}
        path = destination/f"{index:02d}-{row['stage']}.npz"
        np.savez_compressed(path, **arrays)
        points, cells, surface = [arrays[key] for key in ("points", "tetrahedra", "surface")]
        quality, _, volumes = scorecard.tet_metrics(points, cells)
        result["snapshots"].append(dict(stage=row["stage"], path=str(path), sha256=w.digest(path),
            nodes=len(points), elements=len(cells), facets=len(surface),
            native_admitted=False, order=1, source_order_requested=old["native_options"]["order"],
            admission=row["admission"], minimum_signed_volume=float(volumes.min()),
            volume=float(volumes.sum()), minimum_quality=float(quality.min()),
            elements_below_quality_01=int((quality < .1).sum()),
            cap_satisfied=max(len(points), len(cells)) <= w.CAP))
    result["failure_reproduced"] = (result["error_type"] == old["error_type"] and result["error"] == old["error"])
    result["expected_stages"] = args.mode != "observed" or [row["stage"] for row in observations] == [
        "clipped", "adapted", "optimized", "geometry-closed", "final-linear"]
    verify_inputs(manifest)
    result.update(complete=True, input_hashes_verified=True,
        status="rejected-linear-state-captured" if args.mode == "observed" else "original-failure-reproduced")
    if not result["failure_reproduced"] or not result["expected_stages"]:
        result["status"] = "control-failed"
    w.save(destination/"result.json", result)
    assert result["failure_reproduced"], "Original exact error was not reproduced; no equivalence claim"
    assert result["expected_stages"], "Missing expected adaptive capture stages"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--mode", choices=("baseline", "ordinary", "observed"), help=argparse.SUPPRESS)
    args = parser.parse_args()
    for key in ("run", "package", "source"):
        setattr(args, key, getattr(args, key).resolve())
    if args.mode:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    old = json.loads((args.source/"result.json").read_text())
    assert old["complete"] and old["input_hashes_verified"] and old["fixture_verified"]
    assert old["native_attempted"] and not old["native_admitted"] and old["failed_stage"] == "native"
    inputs = dict(old["inputs"])
    for path, expected in inputs.items():
        assert w.digest(Path(path)) == expected, f"Changed historical input: {path}"
    for name, key in ((NAME+"-samples.npz", "samples_sha256"), (NAME+"-specimen.json", "specimen_sha256")):
        source = args.source/name
        assert w.digest(source) == old[key]
        shutil.copyfile(source, args.run/name)
        inputs[str(source)] = inputs[str(args.run/name)] = old[key]
    for name in ("primitive-receipt.json", "membership-probes.npz"):
        shutil.copyfile(args.source/name, args.run/name)
        inputs[str(args.source/name)] = inputs[str(args.run/name)] = w.digest(args.source/name)
    receipt_path = args.package.parent/"receipt.json"
    receipt = json.loads(receipt_path.read_text())
    binaries = list((args.package/"meshers").glob("_meshers*.so"))
    assert len(binaries) == 1 and w.digest(binaries[0]) == receipt["native_binary_sha256"]
    for path in (Path(__file__).resolve(), args.source/"result.json", receipt_path, *args.package.rglob("*.py"), *binaries):
        inputs[str(path)] = w.digest(path)
    for mode in ("baseline", "ordinary", "observed"):
        (args.run/mode).mkdir()
    manifest = dict(inputs=inputs, source=str(args.source), native_options=old["native_options"],
        new_native_binary_sha256=receipt["native_binary_sha256"], worker_timeout_s=120,
        memory_bytes=w.MEMORY, memory_enforcement="process-tree-rss-watchdog",
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, maximum_nodes=w.CAP, maximum_tetrahedra=w.CAP,
        samples_sha256=old["samples_sha256"], specimen_sha256=old["specimen_sha256"],
        comparative_timing_claim=False, admission="Rejected linear diagnostic states only")
    w.save(args.run/"manifest.json", manifest)
    shutil.copyfile(__file__, args.run/"driver.py")
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    records = []
    with Job(parts=["baseline", "ordinary", "observed"]) as job:
        for mode in ("baseline", "ordinary", "observed"):
            with worker_part(job, mode) as part:
                command = [sys.executable, str(Path(__file__).resolve()), "--run", str(args.run),
                    "--package", str(args.package), "--source", str(args.source), "--mode", mode]
                bounded = w.bounded(command, args.run/mode/"worker.log", part, timeout=120)
                w.save(args.run/mode/"worker.json", bounded)
                assert bounded["returncode"] == 0, f"Read {args.run/mode/'worker.log'}"
                row = json.loads((args.run/mode/"result.json").read_text())
                records.append(row)
                job.result(mode, failure_reproduced=row["failure_reproduced"], snapshots=len(row["snapshots"]),
                    peak_gib=bounded["peak_rss_bytes"]/1024**3)
    verify_inputs(manifest)
    w.save(args.run/"result.json", dict(complete=True, input_hashes_verified=True, native_admitted=False,
        status="exact-original-failure-captured", controls=records, comparative_timing_claim=False,
        samples_sha256=old["samples_sha256"], specimen_sha256=old["specimen_sha256"],
        scope="Borrowed linear snapshots before the unchanged geometry rejection; no Tet10 mesh or FE result"))


if __name__ == "__main__":
    main()
