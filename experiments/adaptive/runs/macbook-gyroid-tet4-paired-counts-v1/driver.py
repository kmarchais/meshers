"""Repeat coarse native Tet4 against actual defaults, preserving failed admission."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import microgen_workflow as w


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    args.run = args.run.resolve()
    args.run.mkdir(parents=True, exist_ok=False)
    _, _, meshers, *_ = w.imports()
    from comparison_benchmark import configs
    from job_parts import worker_part
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    prior = w.HERE/"runs/macbook-gyroid-tet4-budget-corrected-v1/result.json"
    native = json.loads(prior.read_text())["candidate"]
    assert native["elements"] == 10832 and native["config"]["order"] == 1
    actual = next(c for c in configs("gyroid") if c["name"] == "actual-default-tet4")
    actual = dict(actual, periodic_trace=True, declared_unit_cell=True)
    fixtures = w.HERE/"runs/macbook-comparison-fixtures-v1"
    paths = [Path(__file__).resolve(), prior, Path(meshers._meshers.__file__)]
    paths += [w.HERE/n for n in ("comparison_benchmark.py", "comparison_fixtures.py",
        "comparison_mesh_audit.py", "microgen_workflow.py", "scorecard.py", "periodic_trace.py",
        "box_projection.py", "facet_incidence_audit.py", "quadratic_edge_removal_probe.py", "job_parts.py")]
    samples = [fixtures/"gyroid"/("gyroid"+s) for s in ("-samples.npz", "-specimen.json")]
    paths += samples
    hashes = {str(p): w.digest(p) for p in paths}
    tasks = []
    for repeat in range(-1, 5):
        arms = [("meshers", native["config"]), ("microgen", actual)]
        if repeat % 2:
            arms.reverse()
        tasks += [(f"repeat-{repeat}-{arm}", repeat, config) for arm, config in arms]
    shutil.copyfile(__file__, args.run/"driver.py")
    w.save(args.run/"manifest.json", dict(inputs=hashes, tasks=tasks,
        numerical_threads=w.NUMERICAL_THREAD_POLICY, memory_bytes=w.MEMORY, worker_timeout_s=180,
        warmup_pairs=1, measured_pairs=5, specimen_unchanged=True, admission_gates_unchanged=True,
        scope="Fresh interleaved raw generation comparison at similar Tet4 counts; failed meshes remain rejected. No admitted-speed-ratio or FE claim."))
    rows = []
    with Job(parts=[t[0] for t in tasks]) as job:
        for name, repeat, config in tasks:
            directory = args.run/name
            directory.mkdir()
            for source in samples:
                shutil.copyfile(source, directory/source.name)
            with worker_part(job, name) as part:
                record = w.bounded([sys.executable, str(w.HERE/"comparison_benchmark.py"),
                    "--run", str(directory), "--fixtures", str(fixtures), "--case", "gyroid",
                    "--worker", json.dumps(config)], directory/"worker.log", part, timeout=180)
            w.save(directory/"worker.json", record)
            assert record["returncode"] == 0, f"Read {directory/'worker.log'}"
            row = json.loads((directory/"result.json").read_text())
            row.update(repeat=repeat, excluded_warmup=repeat < 0, source=str(directory))
            rows.append(row)
            w.save(args.run/"result.json", dict(complete=len(rows)==len(tasks), rows=rows,
                production_changed=False, admitted_comparison_established=False))
            job.result(name, tets=row.get("elements"), generation_s=row.get("generation_s"), admitted=row["qualified"])
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())


if __name__ == "__main__":
    main()
