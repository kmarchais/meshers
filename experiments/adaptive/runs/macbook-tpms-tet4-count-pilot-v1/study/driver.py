"""Bounded same-order TPMS comparison on separately frozen count selections."""
import argparse
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--selection", type=Path, required=True)
    p.add_argument("--arms", default="actual,native")
    p.add_argument("--repetitions", type=int, default=1)
    p.add_argument("--warmups", type=int, default=0)
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    a.run, a.selection = a.run.resolve(), a.selection.resolve()
    selection = json.loads(a.selection.read_text())
    arms = a.arms.split(",")
    assert arms and set(arms) <= {"actual", "native"}
    assert a.repetitions >= 1 and a.warmups >= 0
    a.run.mkdir(parents=True, exist_ok=False)
    if not a.worker:
        timeout = len(selection["cases"])*len(arms)*(a.repetitions+a.warmups)*190+60
        command = [sys.executable, str(Path(__file__).resolve()), "--run", str(a.run/"study"),
            "--selection", str(a.selection), "--arms", a.arms, "--repetitions", str(a.repetitions),
            "--warmups", str(a.warmups), "--worker"]
        record = w.bounded(command, a.run/"study-worker.log",
            SimpleNamespace(update=lambda *args, **kwargs: None), timeout=timeout)
        w.save(a.run/"worker.json", record)
        assert record["returncode"] == 0, "Read study-worker.log and study/summary.json"
        return

    _, _, meshers, microgen, *_ = w.imports()
    from job_parts import worker_part
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    files = [Path(__file__).resolve(), a.selection, Path(meshers._meshers.__file__),
        Path(meshers.__file__), Path(microgen.__file__)]
    files += [w.HERE/name for name in ("tpms_family_screen.py", "tpms_family_fixtures.py",
        "comparison_benchmark.py", "comparison_fixtures.py", "comparison_mesh_audit.py",
        "microgen_workflow.py", "scorecard.py", "periodic_trace.py", "box_projection.py",
        "facet_incidence_audit.py", "quadratic_edge_removal_probe.py", "job_parts.py")]
    import microgen.remesh
    import importlib.metadata
    files += [Path(microgen.remesh.__file__), Path(microgen.shape.tpms.__file__),
        Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))]
    fixture_files = {}
    for entry in selection["cases"]:
        directory = Path(entry["fixtures"])/entry["case"]
        fixture_files[entry["case"]] = [f for f in directory.iterdir() if f.is_file()]
        files += fixture_files[entry["case"]]
        assert all(entry[arm]["order"] == 1 for arm in arms)
    hashes = {str(f): w.digest(f) for f in files}
    shutil.copyfile(__file__, a.run/"driver.py")
    shutil.copyfile(a.selection, a.run/"selection.json")
    tasks = []
    for repeat in range(-a.warmups, a.repetitions):
        order = list(reversed(arms)) if repeat % 2 else arms
        for entry in selection["cases"]:
            for arm in order:
                tasks.append((f"repeat-{repeat}-{entry['case']}-{arm}", entry, arm, repeat))
    w.save(a.run/"manifest.json", dict(inputs=hashes, selection=selection,
        arms=arms, repetitions=a.repetitions, warmups=a.warmups, numerical_threads=w.NUMERICAL_THREAD_POLICY,
        memory_bytes=w.MEMORY, worker_timeout_s=180, maximum_tets=w.CAP, maximum_nodes=w.CAP,
        command=sys.argv, production_changed=False, gates_changed=False,
        scope="Tet4 raw comparison if either route rejects; no FEA or successful admitted-speed claim from rejected meshes"))
    rows = []
    with Job(parts=[t[0] for t in tasks]) as job:
        for name, entry, arm, repeat in tasks:
            directory = a.run/name
            directory.mkdir()
            for source in fixture_files[entry["case"]]:
                shutil.copyfile(source, directory/source.name)
            config = entry[arm]
            family = entry["fixture_kind"] == "tpms"
            script = w.HERE/("tpms_family_screen.py" if family else "comparison_benchmark.py")
            command = [sys.executable, str(script), "--run", str(directory), "--fixtures", entry["fixtures"]]
            if family:
                command += ["--worker", json.dumps(dict(case=entry["case"], config=config))]
            else:
                command += ["--case", entry["case"], "--worker", json.dumps(config)]
            with worker_part(job, name) as part:
                record = w.bounded(command, directory/"worker.log", part, timeout=180)
            w.save(directory/"worker.json", record)
            path = directory/"result.json"
            row = json.loads(path.read_text()) if path.exists() else dict(qualified=False, status="no-result")
            if record["returncode"] != 0:
                row.update(qualified=False, status=record.get("limit_failure") or "worker-failed")
            mesh = directory/"mesh.npz"
            has_mesh = mesh.exists() and row.get("mesh_sha256") == w.digest(mesh)
            row.update(case=entry["case"], config=config, arm=arm, directory=str(directory),
                repeat=repeat, warmup=repeat < 0, resources=record, returned_mesh_verified=has_mesh)
            rows.append(row)
            w.save(a.run/"summary.json", dict(complete=len(rows)==len(tasks), rows=rows,
                attempted=len(rows), planned=len(tasks), selection_sha256=w.digest(a.selection)))
            job.result(name, generation_s=row.get("generation_s"), nodes=row.get("nodes"),
                tets=row.get("elements"), admitted=row["qualified"], returned_mesh=has_mesh)
    assert all(w.digest(Path(f)) == digest for f, digest in hashes.items())


if __name__ == "__main__":
    main()
