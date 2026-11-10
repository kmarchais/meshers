"""Fresh paired native fits on frozen specimens, with independent unchanged admission."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import signal
import sys

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--cases", default="graded-gyroid-density-linear,spinodoid-anisotropic")
    p.add_argument("--repetitions", type=int, default=1)
    p.add_argument("--warmups", type=int, default=0)
    a = p.parse_args()
    np, _, meshers, *_ = w.imports()
    from comparison_benchmark import configs
    from job_parts import worker_part

    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    selection_path = w.ROOT / "handoff/macbook/face-fit-selection-20261110.json"
    selection = json.loads(selection_path.read_text())
    selected = selection["cases"]
    if a.cases != "all":
        wanted = a.cases.split(",")
        selected = [next(e for e in selected if e["case"] == n) for n in wanted]
    a.run.mkdir(parents=True, exist_ok=False)
    paths = [Path(__file__), selection_path, Path(meshers._meshers.__file__)]
    paths += list((w.ROOT / "crates/meshers-core/src").rglob("*.rs"))
    paths += list((w.ROOT / "crates/meshers-python/src").rglob("*.rs"))
    paths += list((w.ROOT / "crates/meshers-python/python/meshers").glob("*.py"))
    paths += [
        w.HERE / f
        for f in (
            "comparison_benchmark.py",
            "comparison_fixtures.py",
            "comparison_mesh_audit.py",
            "microgen_workflow.py",
            "scorecard.py",
            "periodic_trace.py",
            "box_projection.py",
            "facet_incidence_audit.py",
            "quadratic_edge_removal_probe.py",
            "gradient_tie_probe.py",
            "linear_density_fixture.py",
            "job_parts.py",
            "nested_refinement.py",
        )
    ]
    for e in selected:
        paths += list((Path(e["fixtures"]) / e["case"]).iterdir())
        paths += [Path(e["source"]) / "mesh.npz"]
        paths += [f for f in Path(e["fixtures"]).glob("*.json")]
    paths.append(w.HERE / "runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json")
    hashes = {str(f): w.digest(f) for f in paths}
    for i, f in enumerate(paths):
        if f.suffix in (".py", ".rs"):
            shutil.copyfile(f, a.run / f"input-{i}-{f.name}")
    w.save(
        a.run / "manifest.json",
        dict(
            inputs=hashes,
            selection=selected,
            platform=platform.platform(),
            architecture=platform.machine(),
            native_sha256=w.digest(Path(meshers._meshers.__file__)),
            numerical_threads=w.NUMERICAL_THREAD_POLICY,
            packages={
                n: importlib.metadata.version(n) for n in ("meshers", "numpy", "scipy", "microgen")
            },
            repetitions=a.repetitions,
            warmups=a.warmups,
            worker_timeout_s=300,
            memory_bytes=w.MEMORY,
            node_cap=w.CAP,
            tetrahedron_cap=w.CAP,
            command=sys.argv,
            timing_scope="Complete field preparation and native generation including fitting; excludes imports, fixture loading, archive and independent admission",
            scope="Opt-in native fitting; no actual comparator or transferred FEA score in this experiment",
        ),
    )
    tasks = []
    for e in selected:
        config = next(c for c in configs(e["case"]) if c["name"] == "meshers-curved-tet10")
        for repeat in range(-a.warmups, a.repetitions):
            for fit in [False, True] if repeat % 2 == 0 else [True, False]:
                name = f"{e['case']}-repeat-{repeat}-{'fit' if fit else 'baseline'}"
                tasks.append((name, e, repeat, dict(config, curved_fit=fit, periodic_trace=True)))
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("Study deadline")))
    signal.alarm(len(tasks) * 310 + 60)
    rows = []
    with Job(parts=[t[0] for t in tasks]) as job:
        for name, e, repeat, config in tasks:
            directory = a.run / name
            directory.mkdir()
            for f in (Path(e["fixtures"]) / e["case"]).iterdir():
                shutil.copyfile(f, directory / f.name)
            with worker_part(job, name) as part:
                record = w.bounded(
                    [
                        sys.executable,
                        str(w.HERE / "comparison_benchmark.py"),
                        "--worker",
                        json.dumps(config),
                        "--run",
                        str(directory),
                        "--fixtures",
                        e["fixtures"],
                        "--case",
                        e["case"],
                    ],
                    directory / "worker.log",
                    part,
                    timeout=300,
                )
                w.save(directory / "worker.json", record)
                row = (
                    json.loads((directory / "result.json").read_text())
                    if (directory / "result.json").exists()
                    else dict(status="no-result", qualified=False)
                )
                row.update(
                    directory=str(directory), worker=record, repeat=repeat, warmup=repeat < 0
                )
                if not config["curved_fit"] and (directory / "mesh.npz").exists():
                    with (
                        np.load(directory / "mesh.npz") as new,
                        np.load(Path(e["source"]) / "mesh.npz") as old,
                    ):
                        row["frozen_baseline_arrays_equal"] = all(
                            np.array_equal(new[k], old[k])
                            for k in ("points", "tetrahedra", "surface")
                        )
                rows.append(row)
                job.result(
                    name,
                    generation_s=row.get("generation_s"),
                    qualified=row["qualified"],
                    elements=row.get("elements"),
                    nodes=row.get("nodes"),
                )
                w.save(a.run / "summary.json", dict(rows=rows, complete=len(rows) == len(tasks)))
    signal.alarm(0)
    assert all(w.digest(Path(f)) == h for f, h in hashes.items()), (
        "Inputs changed during measurement"
    )


if __name__ == "__main__":
    main()
