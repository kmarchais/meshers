"""Serial, fresh-process measurements with retained outputs and failures."""

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

from .cases import CASES, WORKFLOWS, Case


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def failure_catalog(document):
    """Export reproductions without credentials, environments or mesh payloads."""
    failures = []
    config = document["provenance"]["config"]
    for row in document["rows"]:
        if row["repeat"] < 0 or row["status"] not in {"failed", "timeout", "geometry-warning"}:
            continue
        failures.append(
            {
                **{key: value for key, value in row.items() if key != "preview"},
                "reproduce": (
                    "uv run --project benchmarks/microgen --locked microgen-benchmark"
                    f" --cases {row['case']} --workflows {row['workflow']} --cells {row['cells']}"
                    f" --repeats 1 --warmups 0 --threads {config['threads']}"
                    f" --timeout {config['timeout_s']} --geometry-tolerance {config['geometry_tolerance']}"
                    f" --volume-tolerance {config['volume_tolerance']} --minimum-quality {config['minimum_quality']}"
                    f" --max-tetrahedra {config['max_tetrahedra']} --pair-cap {config['pair_cap']}"
                    f" --mmg-boundary {config['mmg_boundary']}"
                ),
            }
        )
    return {
        "manifest_sha256": document["manifest_sha256"],
        "provenance": document["provenance"],
        "cases": [{k: v for k, v in c.items() if k != "target_preview"} for c in document["cases"]],
        "failures": failures,
        "scope": "Pipeline requiring investigation, not a confirmed root cause. Timeouts are retained limits, not proof of library bugs.",
    }


def verify_sources(expected):
    root = Path(__file__).parent
    for name, sha in expected.items():
        if digest(root / name) != sha:
            raise RuntimeError(f"Benchmark source changed during the run: {name}")


def worker(request_path):
    request = json.loads(request_path.read_text(encoding="utf-8"))
    out = request_path.parent
    result = {"status": "failed", "stage": "imports", "failure_owner": "Benchmark dependencies"}
    phase_start = time.perf_counter()

    def checkpoint(**values):
        nonlocal phase_start
        if "stage" in values:
            phase_start = time.perf_counter()
        result.update(values)
        write_json(out / "progress.json", result)

    try:
        import numpy as np

        from .backends import generate
        from .elasticity import solve
        from .geometry import canonical_mesh, validate
        from .preview import surface_preview

        case = Case(**request["case"])
        result["stage"] = "generation"
        print(f"Generating {case.id}: {request['workflow']} / {request['cells']} cells", flush=True)
        points, tets, timing = generate(
            case, request["workflow"], request["cells"], request["config"], checkpoint
        )
        result.update(timing)
        pipeline = (
            "Meshers"
            if request["workflow"].startswith("meshers")
            else (
                "Microgen / VTK"
                if request["workflow"] == "vtk"
                else "Microgen / MMG pipeline"
                if request["workflow"] == "mmg-microgen"
                else "MMG binding pipeline"
            )
        )
        checkpoint(stage="canonicalization", failure_owner=pipeline)
        points, tets, flips, unused = canonical_mesh(points, tets)
        np.savez_compressed(out / "mesh.npz", points=points, tetrahedra=tets)
        result.update(
            mesh_sha256=digest(out / "mesh.npz"),
            orientation_flips=flips,
            unused_nodes=unused,
            tetrahedra=len(tets),
            nodes=len(points),
        )
        if request["warmup"]:
            result.update(status="warmup", stage="complete")
        else:
            checkpoint(stage="geometry validation", failure_owner=pipeline)
            print("Checking topology, tetrahedral overlap and sampled geometry", flush=True)
            metrics, surface = validate(points, tets, case, request["config"]["pair_cap"])
            result.update(metrics)
            target_fraction = request["target_volume_fraction"]
            result["relative_volume_error"] = (
                abs(metrics["solid_fraction"] - target_fraction) / target_fraction
            )
            checkpoint(stage="elasticity", failure_owner="Benchmark FEA solver")
            print("Solving 18 quadratic fields and six affine patch tests", flush=True)
            fea, fields = solve(points, tets, surface, case.scale)
            result.update(fea)
            np.savez_compressed(out / "fields.npz", **fields)
            result["fields_sha256"] = digest(out / "fields.npz")
            # Display simplification leaves mesh arrays and numerical checks intact.
            try:
                result["preview"] = surface_preview(points, surface, case.scale)
            except Exception as exc:
                result["preview_error"] = f"{type(exc).__name__}: {exc}"
            distance = result.get("sampled_geometry_distance")
            result["comparison_eligible"] = (
                distance is not None
                and distance <= request["config"]["geometry_tolerance"]
                and result["relative_volume_error"] <= request["config"]["volume_tolerance"]
            )
            result["status"] = "ok" if result["comparison_eligible"] else "geometry-warning"
            result["failure_owner"] = pipeline
            result["admission_reasons"] = []
            if distance is None:
                result["admission_reasons"].append(
                    "Boundary distance unavailable: near-zero field gradient"
                )
            elif distance > request["config"]["geometry_tolerance"]:
                result["admission_reasons"].append(
                    f"Sampled boundary distance {distance:.6g}L exceeds {request['config']['geometry_tolerance']:.6g}L"
                )
            if result["relative_volume_error"] > request["config"]["volume_tolerance"]:
                result["admission_reasons"].append(
                    f"Relative volume error {result['relative_volume_error']:.2%} exceeds {request['config']['volume_tolerance']:.2%}"
                )
            result["stage"] = "complete"
            if result["status"] == "ok":
                result.pop("failure_owner", None)
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["failed_stage_elapsed_s"] = time.perf_counter() - phase_start
        result["failure_kind"] = (
            "native rejection" if type(exc).__name__ == "MeshingFailure" else "pipeline failure"
        )
        if (
            result["failure_kind"] == "native rejection"
            and request["cells"] == 8
            and case.id in {"gyroid", "schwarz-p", "neovius", "graded-gyroid", "thin-gyroid"}
        ):
            result["issue_url"] = "https://github.com/kmarchais/meshers/issues/3"
        traceback.print_exc()
    write_json(out / "result.json", result)
    print(f"Result: {result['status']} / {result.get('error', '')}", flush=True)


def stop_tree(process):
    """Kill only descendants of this still-owned Popen process, children first."""
    import psutil

    if process.poll() is not None:
        return
    try:
        parent = psutil.Process(process.pid)
        children = parent.children(recursive=True)
        for child in reversed(children):
            try:
                child.kill()
            except psutil.NoSuchProcess:
                pass
        parent.kill()
        psutil.wait_procs(children + [parent], timeout=5)
    except psutil.NoSuchProcess:
        pass
    process.wait(timeout=10)


def progress(done, total, message, start, failures):
    elapsed = time.monotonic() - start
    eta = elapsed / done * (total - done) if done else None
    try:
        from pueue_progress import report
    except ImportError:
        helper = Path.home() / "bin" / "pueue_progress.py"
        if helper.is_file():
            import importlib.util

            spec = importlib.util.spec_from_file_location("pueue_progress", helper)
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            report = module.report
        else:
            report = None
    if report:
        report(done / total, message, eta_s=eta, completed=done, total=total, failures=failures)
    print(
        f"\r\x1b]9;4;1;{int(100 * done / total)}\x07[{done}/{total}] {message} | elapsed {elapsed:.0f}s",
        end="\n",
        flush=True,
    )


def run_child(request, directory, timeout, heartbeat):
    directory.mkdir(parents=True)
    request_path = directory / "request.json"
    write_json(request_path, request)
    env = os.environ.copy()
    for name in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        env[name] = str(request["config"]["threads"])
    env["PYTHONUNBUFFERED"] = "1"
    # Works both from an installed project and a source checkout.
    source_root = str(Path(__file__).resolve().parent.parent)
    env["PYTHONPATH"] = source_root + os.pathsep + env.get("PYTHONPATH", "")
    with (directory / "worker.log").open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            [sys.executable, "-u", "-m", "microgen_benchmark", "--worker", str(request_path)],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )

        def drain():
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end="", flush=True)

        reader = threading.Thread(target=drain, daemon=True)
        reader.start()
        start = time.monotonic()
        timed_out = False
        try:
            while process.poll() is None:
                heartbeat()
                if time.monotonic() - start > timeout:
                    timed_out = True
                    stop_tree(process)
                    break
                time.sleep(2)
        except BaseException:
            stop_tree(process)
            raise
        finally:
            reader.join(timeout=10)
    path = directory / "result.json"
    if timed_out:
        partial = directory / "progress.json"
        return {
            "stage": "worker",
            "failure_owner": "Owned worker pipeline",
            **(json.loads(partial.read_text(encoding="utf-8")) if partial.is_file() else {}),
            "status": "timeout",
            "error": f"Owned worker tree killed after {timeout}s",
            "failure_kind": "time limit; not a proven library bug",
            "worker_exit_code": process.returncode,
        }
    if process.returncode != 0 or not path.is_file():
        return {
            "status": "failed",
            "stage": "worker",
            "failure_owner": "Owned worker pipeline",
            "error": "Nonzero worker exit or missing completed result",
            "worker_exit_code": process.returncode,
        }
    result = json.loads(path.read_text(encoding="utf-8"))
    result["worker_exit_code"] = process.returncode
    return result


def provenance(config):
    import psutil

    packages = {}
    for name in ["meshers", "microgen", "mmgpy", "vtk", "pyvista", "numpy", "scipy", "psutil"]:
        try:
            distribution = importlib.metadata.distribution(name)
            binary_hashes = {
                str(f): digest(Path(distribution.locate_file(f)))
                for f in distribution.files or []
                if (str(f).endswith((".pyd", ".so", ".dll", ".exe", ".dylib")) or ".so." in str(f))
                and name in {"meshers", "mmgpy"}
            }
            packages[name] = {"version": distribution.version, "native_sha256": binary_hashes}
        except importlib.metadata.PackageNotFoundError:
            packages[name] = {"version": "unavailable"}
    root = Path(__file__).parent
    try:
        git = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        git = None
    return {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "machine": platform.node(),
        "platform": platform.platform(),
        "cpu": platform.processor(),
        "logical_cpus": psutil.cpu_count(),
        "ram_bytes": psutil.virtual_memory().total,
        "python": sys.version,
        "packages": packages,
        "git_head": git,
        "source_sha256": {
            p.name: digest(p) for p in root.iterdir() if p.suffix in {".py", ".html"}
        },
        "config": config,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Run portable Microgen meshing and exact-field elasticity comparisons."
    )
    parser.add_argument("--suite", choices=["smoke", "standard"], default="standard")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cases", default=",".join(c.id for c in CASES))
    parser.add_argument("--workflows", default=",".join(WORKFLOWS))
    parser.add_argument(
        "--cells", help="Comma-separated interval counts, default smoke 8, standard 12,20,28"
    )
    parser.add_argument("--repeats", type=int)
    parser.add_argument("--warmups", type=int)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=float, help="Seconds per worker: smoke 60, standard 300")
    parser.add_argument("--geometry-tolerance", type=float, default=0.01)
    parser.add_argument("--volume-tolerance", type=float, default=0.1)
    parser.add_argument("--mmg-boundary", choices=["caps", "fixed"], default="caps")
    parser.add_argument("--minimum-quality", type=float, default=0.1)
    parser.add_argument("--max-tetrahedra", type=int, default=250_000)
    parser.add_argument("--pair-cap", type=int, default=2_000_000)
    parser.add_argument(
        "--render", type=Path, help="Regenerate HTML from a retained results.json without compute"
    )
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker.resolve())
        return
    from .report import render

    if args.render:
        target = args.output or args.render.with_name("report.html")
        target.write_text(
            render(json.loads(args.render.read_text(encoding="utf-8"))), encoding="utf-8"
        )
        print(target.resolve())
        return
    case_ids, workflows = args.cases.split(","), args.workflows.split(",")
    if len(set(case_ids)) != len(case_ids) or not set(case_ids) <= {c.id for c in CASES}:
        parser.error("Unknown or repeated case IDs")
    if len(set(workflows)) != len(workflows) or not set(workflows) <= set(WORKFLOWS):
        parser.error("Unknown or repeated workflows")
    cells = [
        int(n) for n in (args.cells or ("8" if args.suite == "smoke" else "12,20,28")).split(",")
    ]
    repeats = args.repeats if args.repeats is not None else (1 if args.suite == "smoke" else 3)
    warmups = args.warmups if args.warmups is not None else (0 if args.suite == "smoke" else 1)
    timeout = args.timeout if args.timeout is not None else (60 if args.suite == "smoke" else 300)
    if min(cells) < 2 or repeats < 1 or warmups < 0 or args.threads < 1 or timeout <= 0:
        parser.error("Positive cells >=2, repeats, threads and timeout required; warmups >=0")
    if (
        len(set(cells)) != len(cells)
        or args.geometry_tolerance <= 0
        or args.volume_tolerance <= 0
        or not 0 <= args.minimum_quality <= 1
        or args.pair_cap < 1
    ):
        parser.error("Distinct cells and valid geometry, quality and validation limits required")
    output = (
        args.output or Path("benchmarks/microgen/output") / datetime.now().strftime("%Y%m%d-%H%M%S")
    ).resolve()
    if output.exists() and any(output.iterdir()):
        parser.error("Output directory is not empty; use a new directory")
    output.mkdir(parents=True, exist_ok=True)
    config = {
        "case_ids": case_ids,
        "cells": cells,
        "repeats": repeats,
        "warmups": warmups,
        "threads": args.threads,
        "geometry_tolerance": args.geometry_tolerance,
        "volume_tolerance": args.volume_tolerance,
        "mmg_boundary": args.mmg_boundary,
        "minimum_quality": args.minimum_quality,
        "max_tetrahedra": args.max_tetrahedra,
        "pair_cap": args.pair_cap,
        "timeout_s": timeout,
        "workflows": workflows,
        "suite": args.suite,
        "matched_count_tolerance": 0.05,
        "error_targets": [0.05, 0.1, 0.2],
    }
    cases = [c for c in CASES if c.id in case_ids]
    document = {
        "schema_version": 1,
        "state": "running",
        "provenance": provenance(config),
        "cases": [c.record() for c in cases],
        "rows": [],
    }
    write_json(output / "manifest.json", document["provenance"])
    document["manifest_sha256"] = digest(output / "manifest.json")
    tasks = []
    for c in cases:
        for n in cells:
            for r in range(-warmups, repeats):
                # Rotate backend order, including warmups, to distribute order effects.
                shift = (r + warmups) % len(workflows)
                for w in workflows[shift:] + workflows[:shift]:
                    tasks.append((c, n, r, w))
    start, failures = time.monotonic(), 0
    progress(0, len(tasks), "Preparing benchmark", start, failures)
    from .reference import target

    for c in cases:
        progress(0, len(tasks), f"Independent target diagnostics: {c.id}", start, failures)
        record = next(record for record in document["cases"] if record["id"] == c.id)
        record.update(target(c))
    write_json(output / "results.json", document)
    (output / "report.html").write_text(render(document), encoding="utf-8")
    try:
        for index, (c, n, r, w) in enumerate(tasks):
            verify_sources(document["provenance"]["source_sha256"])
            identity = f"{c.id}/{w}/n{n}/{'warmup' if r < 0 else 'repeat'}{abs(r)}"
            request = {
                "case": c.record(),
                "workflow": w,
                "cells": n,
                "repeat": r,
                "warmup": r < 0,
                "config": config,
                "target_volume_fraction": next(
                    record["target_volume_fraction"]
                    for record in document["cases"]
                    if record["id"] == c.id
                ),
            }
            result = run_child(
                request,
                output / identity,
                timeout,
                lambda i=index, msg=identity, f=failures: progress(i, len(tasks), msg, start, f),
            )
            row = {
                "case": c.id,
                "workflow": w,
                "cells": n,
                "repeat": r,
                "artifact": identity,
                **result,
            }
            document["rows"].append(row)
            failures += result["status"] not in {"ok", "warmup", "geometry-warning"}
            write_json(output / "results.json", document)
            (output / "report.html").write_text(render(document), encoding="utf-8")
    except BaseException:
        document["state"] = "interrupted"
        write_json(output / "results.json", document)
        (output / "report.html").write_text(render(document), encoding="utf-8")
        raise
    document["state"] = "complete"
    document["elapsed_s"] = time.monotonic() - start
    write_json(output / "results.json", document)
    write_json(output / "failures.json", failure_catalog(document))
    (output / "report.html").write_text(render(document), encoding="utf-8")
    progress(len(tasks), len(tasks), "Benchmark complete; failures retained", start, failures)
    print(f"Report: {output / 'report.html'}", flush=True)
