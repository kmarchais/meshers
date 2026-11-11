"""Localize forward error on saved job776 geometry without remeshing."""

import argparse
import json
from pathlib import Path
import shutil
import sys

import microgen_workflow as w


TIMEOUT = 240
NAME = "tpms-pmy-sheet"


def worker(a):
    from pmy_forward_localization import run

    row = json.loads((a.source/"result.json").read_text())
    proof = json.loads((a.trace_proof/"result.json").read_text())
    assert row["complete"] and row["status"] == "admission-rejected" and proof["trace_contract_verified"]
    assert w.digest(a.source/"input.json") == proof["inputs"][str(a.source/"input.json")]
    assert w.digest(a.source/"mesh.npz") == row["mesh_sha256"]
    assert w.digest(a.source/"native-snap-trace.json") == row["native_snap_trace"]["sha256"]
    assert w.digest(a.trace_proof/"native-input-bits.json") == proof["parser_output_sha256"]
    directory = a.run/"candidate"
    directory.mkdir()
    for name in ("result.json", "input.json", "mesh.npz", "native-snap-trace.json"):
        shutil.copyfile(a.source/name, directory/name)
    shutil.copyfile(a.trace_proof/"native-input-bits.json", directory/"native-input-bits.json")
    run(a)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-rollback-native-v1/edge")
    parser.add_argument("--trace-proof", type=Path, default=w.HERE/"runs/macbook-pmy-saved-trace-replay-v1")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = parser.parse_args()
    for name in ("run", "source", "trace_proof", "fixtures"):
        setattr(a, name, getattr(a, name).resolve())
    if a.worker:
        return worker(a)
    from microgen.shape import surface_functions
    from pmy_second_edge_preflight import provenance

    source_manifest = json.loads((a.source.parent/"manifest.json").read_text())
    for path in [Path(surface_functions.__file__), w.HERE/"pmy_neck_control.py", *list((a.fixtures/NAME).iterdir())]:
        assert w.digest(path) == source_manifest["inputs"][str(path)]
    native_source, files = provenance(a)
    files += [Path(__file__), Path(surface_functions.__file__), a.source.parent/"manifest.json",
        a.trace_proof/"result.json", a.trace_proof/"native-input-bits.json", *list((a.fixtures/NAME).iterdir())]
    files += [a.source/name for name in ("input.json", "result.json", "mesh.npz", "native-snap-trace.json")]
    files += [w.HERE/name for name in ("pmy_forward_localization.py", "pmy_failure_probe.py", "pmy_neck_control.py",
        "pmy_second_edge_preflight.py", "tpms_family_fixtures.py", "scorecard.py", "microgen_workflow.py", "job_parts.py", "box_projection.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in files}
    a.run.mkdir(parents=True, exist_ok=False)
    for index, path in enumerate(files):
        if path.suffix in (".py", ".rs", ".toml", ".lock"):
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=TIMEOUT, memory_bytes=w.MEMORY,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, vtk_backend="Sequential", native_measurement_provenance=native_source,
        no_meshing=True, old_artifacts_copied_and_never_modified=True,
        source=str(a.source), exact_trace_proof=str(a.trace_proof)))
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job, result
    from job_parts import worker_part

    with Job(parts={"saved-forward-localization": 1}) as job:
        with worker_part(job, "saved-forward-localization") as part:
            record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
                a.run/"worker.log", part, timeout=TIMEOUT)
    w.save(a.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    assert all(w.digest(Path(path)) == expected for path, expected in hashes.items())
    findings = json.loads((a.run/"candidate/forward-localization.json").read_text())
    findings.update(inputs=hashes, input_hashes_verified=True, prior_native_result_unchanged=True)
    w.save(a.run/"result.json", findings)
    result("PMY saved776 forward", raw_forward=findings["raw_standard_forward_max"],
        reproduced=findings["standard_observer_reproduced"], **findings["summary"])


if __name__ == "__main__":
    main()
