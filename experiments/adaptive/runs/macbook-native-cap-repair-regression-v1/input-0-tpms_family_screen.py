"""One declared coverage screen using actual Microgen seeds and shared mesh admission."""

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
    p.add_argument("--fixtures", type=Path, required=True)
    p.add_argument("--fixture-kind", choices=("tpms", "primitive"), default="tpms")
    p.add_argument("--parts", default="sheet,lower skeletal,upper skeletal")
    p.add_argument("--names")
    p.add_argument("--cases", help="Exact specimen names selected after a completed screen")
    p.add_argument("--arms", default="native-fit,actual-default")
    p.add_argument("--repetitions", type=int, default=1)
    p.add_argument("--warmups", type=int, default=0)
    p.add_argument("--worker", type=json.loads)
    a = p.parse_args()
    assert a.repetitions >= 1 and a.warmups >= 0
    _, _, meshers, *_ = w.imports()
    import comparison_benchmark as benchmark

    if a.fixture_kind == "tpms":
        import tpms_family_fixtures as fixtures
    else:
        import primitive_family_fixtures as fixtures

    if a.worker:
        benchmark.fixtures.case = fixtures.case
        benchmark.fixtures.actual_shape = fixtures.actual_shape
        benchmark.fixtures.periodic = lambda _: a.fixture_kind == "tpms"
        return benchmark.worker(a.run, a.worker["case"], a.worker["config"], a.fixtures)

    from job_parts import worker_part

    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    frozen = json.loads((a.fixtures / "result.json").read_text())
    assert frozen["passed"]
    assert all(w.digest(Path(f)) == h for f, h in frozen["inputs"].items())
    selected = [e for e in frozen["cases"] if e["part"] in a.parts.split(",")]
    if a.names:
        selected = [e for e in selected if e["surface_function"] in a.names.split(",")]
    if a.cases:
        selected = [e for e in selected if e["name"] in a.cases.split(",")]
    assert selected
    arms = {
        "native": dict(
            name="native-uniform16-tet10",
            arm="meshers",
            cells=16,
            order=2,
            adapt=False,
            tolerance=1.0,
            curved_fit=False,
            periodic_trace=True,
        ),
        "native-fit": dict(
            name="native-uniform16-fit-tet10",
            arm="meshers",
            cells=16,
            order=2,
            adapt=False,
            tolerance=1.0,
            curved_fit=True,
            periodic_trace=True,
        ),
        "actual-default": dict(
            name="actual-default-tet4",
            arm="microgen",
            resolution=None,
            order=1,
            periodic_trace=True,
        ),
    }
    arms["native-adaptive"] = dict(
        arms["native-fit"],
        name="native-adaptive002-fit-tet10",
        adapt=True,
        adapt_tolerance=0.002,
        min_length=0.002,
    )
    arms["actual-modified"] = dict(
        arms["actual-default"],
        name="actual-modified40-first-pass-repair-tet4",
        resolution=40,
        first_only=True,
        first_hausd=0.0005,
        first_hmax=0.1,
        mmg=dict(hausd=0.0005, hmax=0.1, hgrad=1.3),
        repair=True,
    )
    arms["native-adaptive8"] = dict(
        arms["native-adaptive"],
        name="native-adaptive002-eight-passes-fit-tet10",
        adapt_iterations=8,
    )
    arms["native-adaptive16"] = dict(
        arms["native-adaptive"],
        name="native-adaptive002-sixteen-passes-fit-tet10",
        adapt_iterations=16,
    )
    arms["native-features"] = dict(
        arms["native-adaptive"],
        name="native-retained-adaptive002-tet4",
        order=1,
        curved_fit=False,
        primitive_retained=True,
    )
    arms["actual-modified24"] = dict(
        arms["actual-modified"],
        name="actual-modified24-first-pass-repair-tet4",
        resolution=24,
    )
    arms["actual-modified24-tet10"] = dict(
        arms["actual-modified24"],
        name="actual-modified24-first-pass-repair-straight-tet10",
        order=2,
    )
    arms = {name: dict(config, declared_unit_cell=True) for name, config in arms.items()}
    if a.fixture_kind == "primitive":
        arms = {
            name: dict(
                config,
                nonperiodic_reason="Complete bounded basic solid; no periodic identification",
                actual_generate_options={}
                if config.get("resolution") is None
                else dict(resolution=config["resolution"]),
            )
            for name, config in arms.items()
        }
    tasks = []
    repeated = a.repetitions > 1 or a.warmups > 0
    for repeat in range(-a.warmups, a.repetitions):
        names = a.arms.split(",")
        if repeated and repeat % 2:
            names.reverse()
        tasks.extend(
            (
                (f"repeat-{repeat}-" if repeated else "") + e["name"] + "-" + name,
                e,
                arms[name],
                repeat,
            )
            for e in selected
            for name in names
        )
    a.run.mkdir(parents=True, exist_ok=False)
    paths = [
        Path(__file__),
        Path(fixtures.__file__),
        Path(benchmark.__file__),
        Path(meshers._meshers.__file__),
        a.fixtures / "result.json",
    ]
    paths += [
        w.HERE / n
        for n in (
            "microgen_workflow.py",
            "scorecard.py",
            "comparison_mesh_audit.py",
            "periodic_trace.py",
            "box_projection.py",
            "facet_incidence_audit.py",
            "quadratic_edge_removal_probe.py",
        )
    ]
    paths += [Path(meshers.__file__), Path(meshers.expression.__file__)]
    import microgen.remesh

    paths += [
        Path(microgen.remesh.__file__),
        Path(microgen.shape.tpms.__file__),
        Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3")),
    ]
    hashes = {str(f): w.digest(f) for f in paths}
    w.save(
        a.run / "manifest.json",
        dict(
            inputs=hashes,
            fixture_hashes=frozen["inputs"],
            selected=selected,
            repetitions=a.repetitions,
            warmups=a.warmups,
            arms={n: arms[n] for n in a.arms.split(",")},
            worker_timeout_s=180,
            memory_bytes=w.MEMORY,
            node_cap=w.CAP,
            tetrahedron_cap=w.CAP,
            numerical_threads=w.NUMERICAL_THREAD_POLICY,
            platform=platform.platform(),
            command=sys.argv,
            coverage_screen_only=True,
            claim_limits="Single generation observations, differing element orders, no FEA or matched-accuracy superiority claim. All admission gates unchanged; defaults and failures retained.",
        ),
    )
    for i, f in enumerate(paths):
        if f.suffix == ".py":
            shutil.copyfile(f, a.run / f"input-{i}-{f.name}")
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("Screen deadline")))
    signal.alarm(len(tasks) * 190 + 60)
    rows = []
    with Job(parts=[t[0] for t in tasks]) as job:
        for name, e, config, repeat in tasks:
            directory = a.run / name
            directory.mkdir()
            for f in (a.fixtures / e["name"]).iterdir():
                shutil.copyfile(f, directory / f.name)
            with worker_part(job, name) as part:
                record = w.bounded(
                    [
                        sys.executable,
                        __file__,
                        "--run",
                        str(directory),
                        "--fixtures",
                        str(a.fixtures),
                        "--fixture-kind",
                        a.fixture_kind,
                        "--worker",
                        json.dumps(dict(case=e["name"], config=config)),
                    ],
                    directory / "worker.log",
                    part,
                    timeout=180,
                )
                w.save(directory / "worker.json", record)
                row = (
                    json.loads((directory / "result.json").read_text())
                    if (directory / "result.json").exists()
                    else dict(status="no-result", qualified=False)
                )
                row.update(
                    worker=record, directory=str(directory), repeat=repeat, warmup=repeat < 0
                )
                if record["returncode"] != 0:
                    row["status"] = record.get("limit_failure") or "worker-failed"
                    row["qualified"] = False
                rows.append(row)
                job.result(
                    name,
                    qualified=row["qualified"],
                    seconds=record["wall_s"],
                    generation_s=row.get("generation_s"),
                    elements=row.get("elements"),
                )
                w.save(
                    a.run / "summary.json",
                    dict(
                        rows=rows,
                        complete=len(rows) == len(tasks),
                        attempted=len(rows),
                        planned=len(tasks),
                    ),
                )
    signal.alarm(0)
    assert all(w.digest(Path(f)) == h for f, h in {**hashes, **frozen["inputs"]}.items())


if __name__ == "__main__":
    main()
