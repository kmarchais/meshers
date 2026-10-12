"""Serial, fresh-process measurements with retained outputs and failures."""

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import math
import os
import platform
import shlex
import shutil
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

from .cases import CASES, DEFAULT_WORKFLOWS, SMOKE_WORKFLOWS, SURFACE_WORKFLOWS, WORKFLOWS, Case
from .qualification import admission, summarize


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")
    deadline = time.perf_counter() + 5
    while True:
        try:
            temporary.replace(path)
            return
        except PermissionError as error:
            # Windows readers without delete sharing temporarily block rename.
            # Keep the old complete JSON visible; persistent errors still fail.
            if (
                os.name != "nt"
                or error.winerror not in {5, 32, 33}
                or time.perf_counter() >= deadline
            ):
                raise
            time.sleep(0.05)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def meshers_native_identity():
    """Hash the imported extension even when an editable wheel has no binary RECORD."""
    native = importlib.import_module("meshers._meshers")
    path = Path(native.__file__).resolve()
    return {"path": str(path), "sha256": digest(path)}


def failure_catalog(document):
    """Export reproductions without credentials, environments or mesh payloads."""
    failures = []
    config = document["provenance"]["config"]
    mmg_args = (
        ["--mmg-bin-directory", config["mmg_bin_directory"]]
        if config.get("mmg_bin_directory")
        else []
    )
    mmg_command = subprocess.list2cmdline(mmg_args) if os.name == "nt" else shlex.join(mmg_args)
    for row in document["rows"]:
        if row["repeat"] < 0 or row["status"] not in {
            "failed",
            "timeout",
            "memory-limit",
            "geometry-warning",
        }:
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
                    f" --mmg-hausd {config.get('mmg_hausd', config['geometry_tolerance'])}"
                    f" --meshers-solid-region {config.get('meshers_solid_region', 'band')}"
                    f" --purpose {config.get('purpose', 'fea')} --fea-target {config.get('fea_target', 0.05)}"
                    f" --error-targets {','.join(map(str, config.get('error_targets', [0.01, 0.02, 0.05])))}"
                    f" --memory-limit-mib {config.get('memory_limit_mib', 8192)}"
                    f" --linear-solver {config.get('linear_solver', 'direct')}"
                    + (f" {mmg_command}" if mmg_command else "")
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


def snapshot_sources(output, expected):
    """Retain the exact runner and lockfile, including uncommitted experiments."""
    root = Path(__file__).parent
    destination = output / "source" / root.name
    destination.mkdir(parents=True)
    for name, sha in expected.items():
        shutil.copy2(root / name, destination / name)
        if digest(destination / name) != sha:
            raise RuntimeError(f"Source snapshot hash disagrees: {name}")
    dependencies = {}
    for name in ["pyproject.toml", "uv.lock", ".python-version"]:
        path = root.parent / name
        shutil.copy2(path, destination.parent / name)
        dependencies[name] = digest(path)
    return dependencies


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

        from .backends import generate, generate_surface
        from .elasticity import solve
        from .geometry import canonical_mesh, canonical_surface, validate
        from .preview import surface_preview

        native_identity = meshers_native_identity()
        result["meshers_native_sha256"] = native_identity["sha256"]
        expected_native = request.get("meshers_native_sha256")
        if expected_native is not None and native_identity["sha256"] != expected_native:
            raise RuntimeError("Loaded Meshers native binary differs from the run manifest")
        case = Case(**request["case"])
        surface_mode = request["config"].get("purpose") == "surface"
        result["stage"] = "generation"
        print(f"Generating {case.id}: {request['workflow']} / {request['cells']} cells", flush=True)
        points, tets, timing = (generate_surface if surface_mode else generate)(
            case, request["workflow"], request["cells"], request["config"], checkpoint
        )
        result.update(timing)
        pipeline = (
            "Meshers"
            if request["workflow"].startswith("meshers")
            else (
                "VTK Flying Edges control"
                if request["workflow"] == "vtk-flying-edges"
                else "Microgen / VTK"
                if request["workflow"].startswith("vtk")
                else "Microgen / MMG pipeline"
                if request["workflow"] == "mmg-microgen"
                else "MMG binding pipeline"
            )
        )
        checkpoint(stage="canonicalization", failure_owner=pipeline)
        if surface_mode:
            points, tets, flips, unused = canonical_surface(points, tets)
            np.savez_compressed(out / "mesh.npz", points=points, triangles=tets)
        else:
            points, tets, flips, unused = canonical_mesh(points, tets)
            np.savez_compressed(out / "mesh.npz", points=points, tetrahedra=tets)
        result.update(
            mesh_sha256=digest(out / "mesh.npz"),
            orientation_flips=flips,
            unused_nodes=unused,
            nodes=len(points),
        )
        result["triangles" if surface_mode else "tetrahedra"] = len(tets)
        if request["warmup"]:
            result.update(status="warmup", stage="complete")
        else:
            checkpoint(stage="geometry validation", failure_owner=pipeline)
            print("Checking topology, tetrahedral overlap and sampled geometry", flush=True)
            target_points = None
            if request.get("target_boundary_path"):
                reference_path = Path(request["target_boundary_path"])
                if digest(reference_path) != request["target_boundary_sha256"]:
                    raise RuntimeError("Independent boundary reference hash changed")
                target_points = np.load(reference_path)["boundary_points"]
            if surface_mode:
                from .surface import validate_surface

                metrics = validate_surface(
                    points, tets, case, request["config"]["pair_cap"], target_points
                )
                surface = tets
            else:
                metrics, surface = validate(
                    points, tets, case, request["config"]["pair_cap"], target_points
                )
            result.update(metrics)
            target_fraction = request["target_volume_fraction"]
            result["relative_volume_error"] = (
                abs(metrics["solid_fraction"] - target_fraction) / target_fraction
            )
            if request["config"].get("purpose", "fea") == "fea":
                checkpoint(stage="elasticity", failure_owner="Benchmark FEA solver")
                print("Solving 18 quadratic fields and six affine patch tests", flush=True)
                fea, fields = solve(
                    points,
                    tets,
                    surface,
                    case.scale,
                    linear_solver=request["config"].get("linear_solver", "direct"),
                    checkpoint=checkpoint,
                )
                result.update(fea)
                np.savez_compressed(out / "fields.npz", **fields)
                result["fields_sha256"] = digest(out / "fields.npz")
            else:
                result["fea_skipped_reason"] = "Geometry-only run; no elasticity solve requested"
            # Display simplification leaves mesh arrays and numerical checks intact.
            try:
                result["preview"] = surface_preview(points, surface, case.scale)
            except Exception as exc:
                result["preview_error"] = f"{type(exc).__name__}: {exc}"
            result.update(
                admission(
                    result, request["config"], request.get("relative_reference_uncertainty", 0)
                )
            )
            result["status"] = "ok" if result["comparison_eligible"] else "geometry-warning"
            result["failure_owner"] = pipeline
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
            and case.id
            in {
                "gyroid",
                "schwarz-p",
                "neovius",
                "graded-gyroid",
                "thin-gyroid",
                "implicit-lattice",
                "spinodoid",
            }
        ):
            number = 4 if case.implicit else 3
            result["issue_url"] = f"https://github.com/kmarchais/meshers/issues/{number}"
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


def progress(done, total, message, start, failures, **values):
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
        report(done / total, message, eta_s=eta, failures=failures, **values)
    print(
        f"\r\x1b]9;4;1;{int(100 * done / total)}\x07[{done}/{total}] {message} | elapsed {elapsed:.0f}s",
        end="\n",
        flush=True,
    )


def run_child(request, directory, timeout, heartbeat):
    import psutil

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
        memory_limited = False
        peak_rss = 0
        memory_limit = request["config"].get("memory_limit_mib", 8192) * 2**20
        monitored = psutil.Process(process.pid)
        try:
            while process.poll() is None:
                live = {}
                try:
                    live = json.loads((directory / "progress.json").read_text(encoding="utf-8"))
                except (FileNotFoundError, json.JSONDecodeError):
                    pass
                try:
                    rss = monitored.memory_info().rss
                    for child in monitored.children(recursive=True):
                        try:
                            rss += child.memory_info().rss
                        except psutil.NoSuchProcess:
                            pass
                    peak_rss = max(peak_rss, rss)
                    live["worker_rss_mib"] = rss / 2**20
                    if rss > memory_limit:
                        memory_limited = True
                        stop_tree(process)
                        break
                except psutil.NoSuchProcess:
                    pass
                heartbeat(live)
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
    if timed_out or memory_limited:
        partial = directory / "progress.json"
        return {
            "stage": "worker",
            "failure_owner": "Owned worker pipeline",
            **(json.loads(partial.read_text(encoding="utf-8")) if partial.is_file() else {}),
            "status": "memory-limit" if memory_limited else "timeout",
            "error": f"Owned worker tree exceeded {memory_limit / 2**20:g} MiB RSS"
            if memory_limited
            else f"Owned worker tree killed after {timeout}s",
            "failure_kind": "resource limit; not a proven library bug",
            "peak_worker_tree_rss_mib": peak_rss / 2**20,
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
    result["peak_worker_tree_rss_mib"] = peak_rss / 2**20
    return result


def provenance(config):
    import psutil

    packages = {}
    for name in [
        "meshers",
        "microgen",
        "mmgpy",
        "vtk",
        "pyvista",
        "numpy",
        "scipy",
        "psutil",
        "pyamg",
    ]:
        try:
            distribution = importlib.metadata.distribution(name)
            binary_hashes = {
                str(f): digest(Path(distribution.locate_file(f)))
                for f in distribution.files or []
                if (str(f).endswith((".pyd", ".so", ".dll", ".exe", ".dylib")) or ".so." in str(f))
                and name in {"meshers", "mmgpy"}
            }
            packages[name] = {"version": distribution.version, "native_sha256": binary_hashes}
            if name == "meshers":
                packages[name]["loaded_native"] = meshers_native_identity()
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
    parser.add_argument("--purpose", choices=["fea", "geometry", "surface"], default="fea")
    parser.add_argument(
        "--linear-solver",
        choices=["direct", "amg"],
        help="Smoke: direct; standard: AMG-preconditioned CG",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cases", default=",".join(c.id for c in CASES))
    parser.add_argument(
        "--workflows", help="Comma-separated workflows; standard adds cleaned controls and one pass"
    )
    parser.add_argument(
        "--mmg-bin-directory",
        type=Path,
        help="Directory containing native mmg3d_O3 for Microgen's CLI path; default: bundled MMG binary",
    )
    parser.add_argument(
        "--meshers-solid-region",
        choices=["band", "negative"],
        default="band",
        help="Native implicit solids: legacy inactive-lower band or equivalent negative region; sheets keep both walls",
    )
    parser.add_argument(
        "--cells", help="Comma-separated interval counts, default smoke 8, standard 24,36,54,80"
    )
    parser.add_argument("--repeats", type=int)
    parser.add_argument("--warmups", type=int)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=float, help="Seconds per worker: smoke 60, standard 1800")
    parser.add_argument(
        "--geometry-tolerance", type=float, help="Boundary distance / L: smoke .01, standard .002"
    )
    parser.add_argument(
        "--volume-tolerance", type=float, help="Relative volume error: smoke .1, standard .01"
    )
    parser.add_argument(
        "--fea-target",
        type=float,
        help="Required worst relative energy error: smoke .05, standard .02",
    )
    parser.add_argument(
        "--error-targets", default="0.01,0.02,0.05", help="Measured FEA error targets"
    )
    parser.add_argument("--mmg-boundary", choices=["caps", "fixed"], default="caps")
    parser.add_argument(
        "--mmg-hausd",
        type=float,
        help="MMG approximation distance / L; default: geometry tolerance, independent output gates stay fixed",
    )
    parser.add_argument(
        "--minimum-quality", type=float, help="Native minimum q: FEA .1, geometry/surface 0"
    )
    parser.add_argument(
        "--max-tetrahedra", type=int, help="Native budget: smoke 250000, standard 2000000"
    )
    parser.add_argument(
        "--pair-cap", type=int, help="Validation cap: smoke 2000000, standard 60000000"
    )
    parser.add_argument(
        "--render", type=Path, help="Regenerate HTML from a retained results.json without compute"
    )
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--memory-limit-mib",
        type=float,
        help="Worker tree RSS ceiling, checked every 2s: smoke 8192 MiB, standard 16384 MiB",
    )
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
    defaults = DEFAULT_WORKFLOWS if args.suite == "standard" else SMOKE_WORKFLOWS
    if args.purpose == "surface":
        defaults = [*defaults, "vtk-flying-edges"]
    case_ids = args.cases.split(",")
    workflows = args.workflows.split(",") if args.workflows else list(defaults)
    if len(set(case_ids)) != len(case_ids) or not set(case_ids) <= {c.id for c in CASES}:
        parser.error("Unknown or repeated case IDs")
    available = SURFACE_WORKFLOWS if args.purpose == "surface" else WORKFLOWS
    if len(set(workflows)) != len(workflows) or not set(workflows) <= set(available):
        parser.error("Unknown or repeated workflows")
    cells = [
        int(n) for n in (args.cells or ("8" if args.suite == "smoke" else "24,36,54,80")).split(",")
    ]
    repeats = args.repeats if args.repeats is not None else (1 if args.suite == "smoke" else 3)
    warmups = args.warmups if args.warmups is not None else (0 if args.suite == "smoke" else 1)
    timeout = args.timeout if args.timeout is not None else (60 if args.suite == "smoke" else 1800)
    args.geometry_tolerance = (
        args.geometry_tolerance
        if args.geometry_tolerance is not None
        else (0.01 if args.suite == "smoke" else 0.002)
    )
    args.volume_tolerance = (
        args.volume_tolerance
        if args.volume_tolerance is not None
        else (0.1 if args.suite == "smoke" else 0.01)
    )
    args.mmg_hausd = args.mmg_hausd if args.mmg_hausd is not None else args.geometry_tolerance
    args.fea_target = (
        args.fea_target
        if args.fea_target is not None
        else (0.05 if args.suite == "smoke" else 0.02)
    )
    args.memory_limit_mib = (
        args.memory_limit_mib
        if args.memory_limit_mib is not None
        else (8192 if args.suite == "smoke" else 16384)
    )
    args.minimum_quality = (
        args.minimum_quality
        if args.minimum_quality is not None
        else (0.1 if args.purpose == "fea" else 0.0)
    )
    args.max_tetrahedra = (
        args.max_tetrahedra
        if args.max_tetrahedra is not None
        else (250_000 if args.suite == "smoke" else 2_000_000)
    )
    args.pair_cap = (
        args.pair_cap
        if args.pair_cap is not None
        else (2_000_000 if args.suite == "smoke" else 60_000_000)
    )
    try:
        error_targets = [float(t) for t in args.error_targets.split(",")]
    except ValueError:
        parser.error("Error targets must be comma-separated numbers")
    if min(cells) < 2 or repeats < 1 or warmups < 0 or args.threads < 1 or timeout <= 0:
        parser.error("Positive cells >=2, repeats, threads and timeout required; warmups >=0")
    if (
        len(set(cells)) != len(cells)
        or args.geometry_tolerance <= 0
        or args.volume_tolerance <= 0
        or args.mmg_hausd <= 0
        or not 0 <= args.minimum_quality <= 1
        or args.pair_cap < 1
        or args.max_tetrahedra < 1
        or not 0 < args.fea_target < 1
        or any(not 0 < t < 1 for t in error_targets)
        or len(set(error_targets)) != len(error_targets)
        or not 0 < args.memory_limit_mib < float("inf")
        or not all(
            math.isfinite(v)
            for v in [args.geometry_tolerance, args.volume_tolerance, args.mmg_hausd, timeout]
        )
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
        "mmg_hausd": args.mmg_hausd,
        "meshers_solid_region": args.meshers_solid_region,
        "mmg_bin_directory": str(args.mmg_bin_directory.resolve())
        if args.mmg_bin_directory
        else None,
        "minimum_quality": args.minimum_quality,
        "max_tetrahedra": args.max_tetrahedra,
        "pair_cap": args.pair_cap,
        "timeout_s": timeout,
        "workflows": workflows,
        "suite": args.suite,
        "matched_count_tolerance": 0.05,
        "error_targets": sorted(error_targets),
        "fea_target": args.fea_target,
        "purpose": args.purpose,
        "linear_solver": args.linear_solver or ("direct" if args.suite == "smoke" else "amg"),
        "memory_limit_mib": args.memory_limit_mib,
    }
    cases = [c for c in CASES if c.id in case_ids]
    document = {
        "schema_version": 2,
        "state": "running",
        "provenance": provenance(config),
        "cases": [c.record() for c in cases],
        "rows": [],
    }
    document["provenance"]["dependency_sha256"] = snapshot_sources(
        output, document["provenance"]["source_sha256"]
    )
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
    total_units = len(tasks) + len(cases)
    progress(0, total_units, "Preparing benchmark", start, failures)
    from .reference import target

    for reference_index, c in enumerate(cases):
        message = f"Independent target diagnostics: {c.id}"
        progress(reference_index, total_units, message, start, failures)
        record = next(record for record in document["cases"] if record["id"] == c.id)
        stop_heartbeat = threading.Event()

        def reference_heartbeat(index=reference_index, msg=message, signal=stop_heartbeat):
            while not signal.wait(2):
                progress(index, total_units, msg, start, failures)

        reporter = threading.Thread(target=reference_heartbeat, daemon=True)
        reporter.start()
        try:
            record.update(target(c, output))
        finally:
            stop_heartbeat.set()
            reporter.join()
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
                "meshers_native_sha256": document["provenance"]["packages"]["meshers"][
                    "loaded_native"
                ]["sha256"],
                "target_volume_fraction": next(
                    record["target_volume_fraction"]
                    for record in document["cases"]
                    if record["id"] == c.id
                ),
                "relative_reference_uncertainty": next(
                    record.get("target_volume_uncertainty", 0) / record["target_volume_fraction"]
                    for record in document["cases"]
                    if record["id"] == c.id
                ),
            }
            reference_record = next(record for record in document["cases"] if record["id"] == c.id)
            request.update(
                target_boundary_path=str(output / reference_record["target_boundary_artifact"]),
                target_boundary_sha256=reference_record["target_boundary_sha256"],
            )
            result = run_child(
                request,
                output / identity,
                timeout,
                lambda live, i=index, msg=identity, f=failures: progress(
                    i + len(cases),
                    total_units,
                    f"{msg}: {live.get('stage', 'imports')}",
                    start,
                    f,
                    **{
                        key: live[key]
                        for key in ["worker_rss_mib", "fea_cg_relative_residual"]
                        if key in live
                    },
                ),
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
            document["qualification"] = summarize(document)
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
    write_json(output / "qualification.json", document["qualification"])
    (output / "report.html").write_text(render(document), encoding="utf-8")
    progress(total_units, total_units, "Benchmark complete; failures retained", start, failures)
    print(f"Report: {output / 'report.html'}", flush=True)
