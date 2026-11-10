"""Apply explicit hausd/hmax controls to both public MMG passes at runtime."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import sys
import time
import traceback
from unittest.mock import patch

import cad_union_comparison as base
import cad_union_declared_cell_control as declared
import microgen_workflow as w


REQUEST = dict(hausd=.0003, hmax=.03)
CONTROL = dict(declared.CONTROL,
    name="actual-cad-size003-declared-cell-mmg-hausd0003-hmax003-both-passes",
    modifications=[*declared.CONTROL["modifications"],
                   "Request hausd=.0003 and hmax=.03 in both public MMG calls, including the first nofem pass"],
    unchanged="CAD generation, Gmsh periodic constraints, BoxMesh construction and required triangles, both MMG algorithm modes, retry policy, other MMG arguments and output snapping",
    mmg_requested_both_passes=REQUEST,
    resource_rationale="Radius .15 and hausd .0003 imply a circular chord scale near .019. This is one bounded trial, not a node-count guarantee. Each intermediate retains the 250000-node/tet cap and the worker retains 180 seconds and 8 GiB.",
    geometric_scope="MMG controls distance to its input mesh. Independent analytic geometry admission remains required.")


def worker(run, source):
    import microgen

    original = microgen.Mmg.mmg3d
    original_command = microgen.Mmg._run_mmg_command
    calls = []
    journal = run / "effective-mmg-calls.json"

    def save():
        w.save(journal, dict(control=CONTROL, requested_both_passes=REQUEST, calls=calls,
            scope="Effective keyword arguments and exact binary argv are recorded before their respective calls"))

    def record_binary_command(command):
        assert calls and calls[-1]["status"] == "started"
        for name, value in REQUEST.items():
            flag = "-" + name
            assert command.count(flag) == 1 and float(command[command.index(flag) + 1]) == value
        calls[-1]["binary_argv"] = list(command)
        save()
        return original_command(command)

    def constrained_call(**parameters):
        incoming = {key: value for key, value in parameters.items() if key not in ("input", "output")}
        assert bool(parameters.get("nofem")) != bool(parameters.get("ls")), "Unexpected public MMG pass"
        effective = {**parameters, **REQUEST}
        effective_public = {key: value for key, value in effective.items() if key not in ("input", "output")}
        assert all(effective[key] == value for key, value in parameters.items() if key not in REQUEST)
        entry = dict(call=len(calls) + 1, status="started",
            pass_kind="first-nofem" if parameters.get("nofem") else "second-level-set",
            control_requested_arguments=REQUEST, public_call_arguments=incoming,
            effective_arguments=effective_public,
            input=str(Path(parameters["input"]).resolve()), output=str(Path(parameters["output"]).resolve()))
        calls.append(entry)
        save()
        started = time.perf_counter()
        try:
            value = original(**effective)
            assert "binary_argv" in entry, "Expected MMG command was not recorded"
            entry["status"] = "returned"
            return value
        except Exception as error:
            entry.update(status="failed", error_type=type(error).__name__, error=str(error),
                         traceback=traceback.format_exc())
            raise
        finally:
            entry["command_s"] = time.perf_counter() - started
            save()

    save()
    with patch.dict(declared.CONTROL, CONTROL, clear=True), patch.dict(base.MMG, REQUEST), \
         patch.object(microgen.Mmg, "mmg3d", constrained_call), \
         patch.object(microgen.Mmg, "_run_mmg_command", record_binary_command):
        declared.worker(run, source)
    result = json.loads((run / "result.json").read_text())
    recorded = result["mmg_calls"]
    assert len(recorded) == len(calls)
    for call, effective in zip(recorded, calls, strict=True):
        assert call["arguments"] == effective["public_call_arguments"]
        call["public_call_arguments"] = call["arguments"]
        call["control_requested_arguments"] = REQUEST
        call["effective_arguments"] = effective["effective_arguments"]
        call["binary_argv"] = effective.get("binary_argv")
        call["arguments"] = effective["effective_arguments"]
        call["arguments_scope"] = "Effective arguments passed to the original Mmg.mmg3d wrapper"
    result.update(control=CONTROL, effective_mmg_journal=str(journal.resolve()),
                  effective_mmg_journal_sha256=w.digest(journal))
    for key, folder in (("direct", "actual-cad-gmsh-direct"), ("final", "actual-cad-gmsh-public-mmg")):
        if key in result:
            result[key]["control"] = CONTROL
            result[key]["route"] = ("Explicit actual CAD/Gmsh size .03 and declared-cell validation" +
                (" plus public MMG with explicit hausd=.0003 and hmax=.03 in both calls"
                 if key == "final" else ", without MMG"))
            w.save(run / folder / "result.json", result[key])
    w.save(run / "result.json", result)


def compact_candidate(row):
    geometry = row["geometry"]
    return dict(qualified=row["qualified"], nodes=row["nodes"], elements=row["elements"],
        sampled_two_sided_max=geometry.get("sampled_two_sided_max"),
        sampled_reverse_max=geometry.get("sampled_reverse_max"),
        relative_volume_error=geometry["relative_volume_error"], q_min_corner=geometry["q_min_corner"],
        mesh_file=row["mesh_file"], mesh_sha256=row["mesh_sha256"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--prior-control", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    prior_manifest = json.loads((args.prior_control / "manifest.json").read_text())
    source = Path(prior_manifest["source"])
    if args.worker:
        return worker(args.run, source)
    from job_parts import worker_part
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    assert importlib.metadata.version("cadquery-ocp-novtk") == base.CAD_PIN
    prior = json.loads((args.prior_control / "summary.json").read_text())["rows"][0]
    assert prior["direct"]["qualified"] and not prior["qualified"]
    validation = json.loads((args.prior_control / "actual/validation.json").read_text())
    assert validation["required_triangles"] and all(row["equal"] for row in validation["required_triangles"])
    hashes = dict(prior_manifest["inputs"])
    for path, digest in hashes.items():
        assert w.digest(Path(path)) == digest, "Frozen input changed: " + path
    for path in (Path(__file__), args.prior_control / "manifest.json", args.prior_control / "summary.json",
                 args.prior_control / "actual/validation.json", Path(prior["direct"]["mesh_file"]),
                 Path(prior["final"]["mesh_file"])):
        hashes[str(path.resolve())] = w.digest(path)
    for call in prior["mmg_calls"]:
        for key in ("input_artifacts", "output_artifacts"):
            for path, digest in call.get(key, {}).items():
                assert w.digest(Path(path)) == digest
                hashes[path] = digest
    old_pair = json.loads((source / "summary.json").read_text())
    native = next(row for row in old_pair["rows"] if row["arm"] == "native")["final"]
    reference = dict(saved_native=compact_candidate(native), admitted_direct_cad_gmsh=compact_candidate(prior["direct"]),
        rejected_public_mmg=compact_candidate(prior["final"]),
        observation="The admitted direct CAD/Gmsh seed uses fewer tetrahedra and nodes and has lower distance and volume error than the saved native mesh. Native minimum quality is higher. No successful speed ratio is asserted.")
    args.run.mkdir(parents=True, exist_ok=False)
    for index, path in enumerate(map(Path, hashes)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run / f"input-{index}-{path.name}")
    directory = args.run / "actual"
    directory.mkdir()
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(source / "fixture" / (base.NAME + suffix), directory / (base.NAME + suffix))
    w.save(args.run / "manifest.json", dict(inputs=hashes, control=CONTROL, prior_control=str(args.prior_control.resolve()),
        source=str(source.resolve()), references=reference, cad_version_pin=base.CAD_PIN,
        packages={name: importlib.metadata.version(name) for name in prior_manifest["packages"]},
        actual_gmsh={**base.GMSH, "size": CONTROL["gmsh_size"]}, public_mmg_request={**base.MMG, **REQUEST},
        effective_mmg_request_both_passes=REQUEST, platform=platform.platform(),
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, memory_bytes=w.MEMORY,
        worker_timeout_s=180, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        timing_scope=prior_manifest["timing_scope"] + " Effective-argument journal overhead is charged.",
        mechanical_accuracy_claim=False, successful_performance_ratio=None, command=sys.argv))
    with Job(parts=["actual-mmg-control"]) as job:
        with worker_part(job, "actual-mmg-control") as part:
            record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--worker", "--run",
                str(directory.resolve()), "--prior-control", str(args.prior_control.resolve())],
                directory / "worker.log", part, timeout=180)
            w.save(directory / "worker.json", record)
            row = (json.loads((directory / "result.json").read_text()) if (directory / "result.json").exists()
                   else dict(status="no-result", qualified=False))
            row.update(worker=record, directory=str(directory.resolve()))
            if record["returncode"] != 0:
                row.update(status="worker-failed", qualified=False)
            job.result("actual-mmg-control", qualified=row["qualified"], generation_s=row.get("generation_s"),
                       wall_s=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    w.save(args.run / "summary.json", dict(rows=[row], control=CONTROL, references=reference,
        complete=True, original_results_unchanged=True, successful_performance_ratio=None, mechanical_accuracy_claim=False))
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items()), "Inputs changed during control"


if __name__ == "__main__":
    main()
