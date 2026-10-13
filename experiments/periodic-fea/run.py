"""Serial periodic comparison with frozen sources and bounded owned workers."""

import argparse
import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    deadline = time.monotonic() + 5
    while True:
        try:
            temporary.replace(path)
            return
        except PermissionError:
            if time.monotonic() >= deadline:
                raise
            time.sleep(0.05)


def prepare(output, native_package, task_specs=None, mmg_hausd=0.002, worker_timeout=1200,
            max_tetrahedra=2_000_000, mmg_seed_cells=None, mmg_size_map=None, surface_seed_cells=None,
            native_build_manifest=None, specimens_path=None, geometry_tolerance=0.002,
            native_geometry_tolerance=None):
    """Freeze the executable experiment before it enters the shared queue."""
    output.mkdir(parents=True, exist_ok=False)
    source = output / "source"
    source.mkdir()
    for path in HERE.glob("*.py"):
        shutil.copy2(path, source / path.name)
    shutil.copytree(
        ROOT / "benchmarks/microgen/microgen_benchmark",
        source / "microgen_benchmark",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    shutil.copytree(
        native_package, source / "meshers", ignore=shutil.ignore_patterns("__pycache__")
    )
    for name in ["pyproject.toml", "uv.lock"]:
        shutil.copy2(ROOT / "benchmarks/microgen" / name, source / name)
    if mmg_size_map is not None:
        shutil.copy2(mmg_size_map, source / "size-map.npz")
    if specimens_path is not None:
        specimens_document = json.loads(specimens_path.read_text())
        if "spinodoid_spectrum" in specimens_document:
            spectrum = specimens_path.parent / specimens_document["spinodoid_spectrum"]
            if digest(spectrum) != specimens_document["spinodoid_spectrum_sha256"]:
                raise ValueError("Changed spinodoid spectrum")
            shutil.copy2(spectrum, source / "spinodoid-spectrum.npz")
        shutil.copy2(specimens_path, source / "specimens.json")
    native = next(
        path
        for path in (source / "meshers").iterdir()
        if path.name.startswith("_meshers") and path.suffix in {".pyd", ".so"}
    )
    native_hash = digest(native)
    native_commit = "66b800cf80e791788b980c01a8438f0e7f8ec8c6"
    if native_build_manifest is not None:
        build = json.loads(native_build_manifest.read_text())
        if build["state"] != "complete" or native_hash != build["native_sha256"] or any(
            command["exit_code"] != 0 for command in build["commands"]
        ):
            raise RuntimeError("Experimental native build is incomplete or does not match package")
        native_commit = build["source_commit"]
        shutil.copy2(native_build_manifest, source / "native-build.json")
    elif native_hash != "4ee6b2e90ae9e82a5dcbb40ba0bed7ac5bc1fc0cbe3567a678ec92797d28438a":
        raise RuntimeError("Expected the previously verified packaged PR #5 native build")
    mmg_executable = Path(
        importlib.metadata.distribution("mmgpy").locate_file(
            "mmgpy/bin/mmg3d_O3.exe" if os.name == "nt" else "mmgpy/bin/mmg3d_O3"
        )
    )
    dependency_sources = {str(mmg_executable.resolve()): digest(mmg_executable)}
    for package in ["microgen", "fedoo"]:
        distribution = importlib.metadata.distribution(package)
        for name in distribution.files:
            if str(name).endswith(".py") and not str(name).startswith(".."):
                path = Path(distribution.locate_file(name)).resolve()
                dependency_sources[str(path)] = digest(path)
    manifest = {
        "state": "prepared",
        "source_sha256": {
            str(path.relative_to(source)): digest(path)
            for path in source.rglob("*")
            if path.is_file()
        },
        "native_path": str(native.resolve()),
        "native_sha256": native_hash,
        "dependency_sources_sha256": dependency_sources,
        "native_source_commit": native_commit,
        "packages": {
            name: importlib.metadata.version(name)
            for name in [
                "meshers",
                "microgen",
                "fedoo",
                "numpy",
                "scipy",
                "pyamg",
                "pyvista",
                "vtk",
                "mmgpy",
                "meshio",
                "simcoon",
            ]
        },
        "interpreter": sys.executable,
        "config": {
            "periodic": True,
            "geometry_tolerance": geometry_tolerance,
            "native_geometry_tolerance": native_geometry_tolerance or geometry_tolerance,
            "volume_tolerance": 0.01,
            "minimum_quality": 0.1,
            "threads": 1,
            "max_tetrahedra": max_tetrahedra,
            "max_validation_pairs": 60_000_000,
            "memory_limit_mib": 16384,
            "mmg_boundary": "caps",
            "mmg_hausd": mmg_hausd,
            "periodic_pair_tolerance": 1e-10,
            "worker_timeout_s": worker_timeout,
        },
        "specimens": {
            "solid-cell": {"function": "solid", "offset": 1},
            "periodic-gyroid": {"function": "gyroid", "offset": 1.1},
            "periodic-schwarz-p": {"function": "schwarz_p", "offset": 1.2},
        },
        "measurements": "Single calibration samples, no warmup; no performance claim",
        "reference_status": "No converged continuum stiffness reference yet",
    }
    tasks = [
        {"case": "solid-cell", "workflow": workflow, "cells": 4}
        for workflow in ["meshers-fast", "mmg-microgen"]
    ]
    if specimens_path is not None:
        manifest["specimens"] = specimens_document["specimens"]
    tasks += [
        {"case": "periodic-gyroid", "workflow": workflow, "cells": 16}
        for workflow in ["meshers-fast", "meshers-pass1", "meshers-balanced", "mmg-microgen"]
    ]
    tasks += [
        {"case": "periodic-gyroid", "workflow": workflow, "cells": 24}
        for workflow in ["meshers-pass1", "mmg-microgen"]
    ]
    if task_specs:
        known = {
            "meshers-fast",
            "meshers-pass1",
            "meshers-balanced",
            "mmg-microgen",
            "mmg-microgen-first-distance",
            "mmg-native-seed",
            "mmg-native-seed-metric",
            "gmsh-hxt-native-surface",
        }
        tasks = []
        for spec in task_specs:
            descriptor, resolution = spec.split(":")
            case_id, workflow = (
                descriptor.split("/") if "/" in descriptor else ("periodic-gyroid", descriptor)
            )
            cells = int(resolution)
            if workflow not in known or cells < 2 or case_id not in manifest["specimens"]:
                raise ValueError(f"Unknown workflow or invalid resolution: {spec}")
            tasks.append({"case": case_id, "workflow": workflow, "cells": cells})
    manifest["tasks"] = tasks
    if any(task["workflow"] == "gmsh-hxt-native-surface" for task in tasks):
        import gmsh
        manifest["packages"]["gmsh"] = importlib.metadata.version("gmsh")
        for path in [Path(gmsh.__file__), Path(gmsh.libpath)]:
            manifest["dependency_sources_sha256"][str(path.resolve())] = digest(path)
    if mmg_seed_cells is not None:
        manifest["config"]["mmg_seed_cells"] = mmg_seed_cells
    if surface_seed_cells is not None:
        manifest["config"]["surface_seed_cells"] = surface_seed_cells
    if any(task["workflow"] == "mmg-native-seed-metric" for task in tasks):
        if mmg_size_map is None:
            raise ValueError("The metric workflow requires --mmg-size-map")
        manifest["config"].update(mmg_size_map_path=str((source / "size-map.npz").resolve()),
                                  mmg_size_map_sha256=digest(source / "size-map.npz"))
    write_json(output / "manifest.json", manifest)
    print(
        json.dumps(
            {
                "prepared": str(output),
                "tasks": len(tasks),
                "native_sha256": manifest["native_sha256"],
            }
        ),
        flush=True,
    )


def native_seed_control(case, cells, config, directory, checkpoint):
    """Standard MMG adaptation on a valid native seed, with required box caps."""
    from contextlib import nullcontext
    from unittest.mock import patch

    import numpy as np
    import pyvista as pv
    from microgen import BoxMesh, Mmg
    from microgen_benchmark.backends import (
        configure_mmg_cli,
        generate,
        reopenable_temporary_file,
        tetrahedra,
    )

    remesh = importlib.import_module("microgen.remesh")
    seed_cells = config.get("mmg_seed_cells", cells)
    points, tets, seed = generate(case, "meshers-pass1", seed_cells, config, checkpoint)
    np.savez_compressed(directory / "native-seed.npz", points=points, tetrahedra=tets)
    cli = configure_mmg_cli(config.get("mmg_bin_directory"))
    started = time.perf_counter()
    grid = pv.UnstructuredGrid(
        np.c_[np.full(len(tets), 4), tets].ravel(), np.full(len(tets), pv.CellType.TETRA), points
    )
    input_path, raw_path, output_path = [
        directory / name
        for name in [
            "native-seed-required.mesh",
            "standard-mmg-raw.mesh",
            "standard-mmg-volume.mesh",
        ]
    ]
    context = (
        patch("microgen.remesh.NamedTemporaryFile", reopenable_temporary_file)
        if os.name == "nt"
        else nullcontext()
    )
    with context:
        remesh._generate_mesh_with_required_triangles(BoxMesh.from_pyvista(grid), str(input_path))
    # A standard adaptation uses the valid volume directly. It needs neither
    # nofem nor a level-set pass. This is an experimental algorithm control,
    # separate from the unchanged public two-pass Microgen baseline.
    sizing, metric_metadata = {"hsiz": case.scale / cells}, {}
    if config.get("mmg_size_map_path"):
        from make_size_map import write_mmg_metric
        metric_path, metric_metadata = write_mmg_metric(
            Path(config["mmg_size_map_path"]), input_path, points, tets,
            case.scale / cells, config["mmg_size_map_sha256"],
        )
        sizing = {"metric": str(metric_path)}
    Mmg.mmg3d(
        input=str(input_path),
        output=str(raw_path),
        hausd=config["mmg_hausd"] * case.scale,
        **sizing,
        hgrad=1.3,
        nr=True,
    )
    remesh._remove_unnecessary_fields_from_mesh_file(str(raw_path), str(output_path), 2, 3)
    output = pv.UnstructuredGrid(str(output_path))
    remesh._snap_to_expected_bounds(output, points, 1e-8)
    adapted_points, adapted_tets = tetrahedra(output)
    adaptation_s = time.perf_counter() - started
    return (
        adapted_points,
        adapted_tets,
        {
            **seed,
            **cli,
            **metric_metadata,
            "mesh_s": seed["mesh_s"] + adaptation_s,
            "generation_s": seed["generation_s"] + adaptation_s,
            "native_seed_mesh_s": seed["mesh_s"],
            "native_seed_cells": seed_cells,
            "mmg_sizing_cells": cells,
            "native_seed_tetrahedra": len(tets),
            "standard_adaptation_s_including_capture": adaptation_s,
            "mmg_standard_adaptation_hausd_over_L": config["mmg_hausd"],
            "mmg_algorithm_control": "single standard adaptation of a valid Meshers seed",
            "mmg_periodic_policy": "Microgen required cap triangles, independently validated",
        },
    )


def worker(request_path):
    """Generate, independently admit, then audit six periodic elasticity loads."""
    import importlib
    from types import SimpleNamespace
    from unittest.mock import patch

    import homogenization
    import numpy as np
    from microgen_benchmark.backends import generate
    from fixtures import make_case
    from microgen_benchmark.geometry import canonical_mesh, validate
    from periodic import periodic_pairs

    request = json.loads(request_path.read_text(encoding="utf-8"))
    directory = request_path.parent
    manifest = request["manifest"]
    task, config = request["task"], manifest["config"]
    result = {**task, "status": "failed", "stage": "startup"}
    stage_start = time.monotonic()
    live = {}

    def checkpoint(**values):
        live.update(values)
        write_json(directory / "progress.json", live)

    try:
        native = Path(importlib.import_module("meshers._meshers").__file__)
        if digest(native) != manifest["native_sha256"]:
            raise RuntimeError("Imported native binary differs from prepared manifest")
        result["native_sha256"] = digest(native)
        for name, version in manifest["packages"].items():
            if importlib.metadata.version(name) != version:
                raise RuntimeError(f"Dependency changed since preparation: {name}")
        for path, expected in manifest["dependency_sources_sha256"].items():
            if digest(path) != expected:
                raise RuntimeError(f"Dependency source changed since preparation: {path}")
        specimen = manifest["specimens"][task["case"]]
        case = make_case(task["case"], specimen)
        result["stage"] = "generation"
        checkpoint(stage="generation")
        calls, capture_seconds = [], 0.0
        from microgen.external import Mmg

        original = Mmg.mmg3d

        def recorded_mmg(**parameters):
            nonlocal capture_seconds
            index = len(calls)
            record = {"parameters": parameters.copy()}
            calls.append(record)
            start = time.perf_counter()
            source = Path(parameters["input"])
            shutil.copy2(source, directory / f"mmg-{index}-input.mesh")
            record["input_sha256"] = digest(source)
            checkpoint(stage=f"MMG pass {index + 1}", mmg_calls=len(calls))
            capture_seconds += time.perf_counter() - start
            start = time.perf_counter()
            try:
                return original(**parameters)
            finally:
                record["command_seconds"] = time.perf_counter() - start
                start = time.perf_counter()
                output = Path(parameters["output"])
                if output.is_file():
                    shutil.copy2(output, directory / f"mmg-{index}-output.mesh")
                    record["output_sha256"] = digest(output)
                write_json(directory / "mmg-calls.json", calls)
                capture_seconds += time.perf_counter() - start

        with patch.object(Mmg, "mmg3d", recorded_mmg):
            generation_config = dict(config)
            generation_config["geometry_tolerance"] = config.get("native_geometry_tolerance", config["geometry_tolerance"])
            result["requested_native_geometry_tolerance"] = generation_config["geometry_tolerance"]
            if task["workflow"] in {"mmg-native-seed", "mmg-native-seed-metric"}:
                control_config = dict(generation_config)
                if task["workflow"] == "mmg-native-seed":
                    control_config.pop("mmg_size_map_path", None)
                points, tets, metrics = native_seed_control(
                    case, task["cells"], control_config, directory, checkpoint
                )
            elif task["workflow"] == "gmsh-hxt-native-surface":
                from gmsh_surface import generate as generate_hxt
                points, tets, metrics = generate_hxt(
                    case, task["cells"], generation_config, directory, checkpoint
                )
            else:
                points, tets, metrics = generate(
                    case, task["workflow"], task["cells"], generation_config, checkpoint
                )
        result.update(metrics)
        result["mmg_diagnostic_capture_seconds"] = capture_seconds
        for key in ["generation_s", "mesh_s"]:
            result[key] -= capture_seconds
        if "standard_adaptation_s_including_capture" in result:
            result["standard_adaptation_s"] = (
                result["standard_adaptation_s_including_capture"] - capture_seconds
            )
        np.savez_compressed(directory / "raw-mesh.npz", points=points, tetrahedra=tets)
        result["stage"] = "independent geometry validation"
        checkpoint(stage=result["stage"])
        points, tets, flipped, unused = canonical_mesh(points, tets)
        result.update(reoriented_tets=flipped, unused_nodes_removed=unused)
        reference = np.load(Path(request["reference_path"]))
        metrics, surface = validate(
            points, tets, case, config["max_validation_pairs"], reference["boundary_points"]
        )
        result.update(metrics)
        np.savez_compressed(directory / "mesh.npz", points=points, tetrahedra=tets, surface=surface)
        result["mesh_sha256"] = digest(directory / "mesh.npz")
        volume = request["volume_reference"]
        result["relative_volume_error"] = (
            abs(result["solid_fraction"] - volume["fraction"]) / volume["fraction"]
        )
        result["relative_volume_error_with_uncertainty"] = (
            result["relative_volume_error"] + volume["uncertainty"] / volume["fraction"]
        )
        result["quality_eligible"] = result["q_min"] >= config["minimum_quality"]
        result["geometry_eligible"] = (
            result["sampled_geometry_distance"] is not None
            and result["sampled_geometry_distance"] <= config["geometry_tolerance"]
            and result["sampled_reverse_distance"] <= config["geometry_tolerance"]
            and result["sampled_cap_overshoot"] is not None
            and result["sampled_cap_overshoot"] <= config["geometry_tolerance"]
            and result["sampled_outward_orientation"]
            and result["relative_volume_error_with_uncertainty"] <= config["volume_tolerance"]
        )
        result["stage"] = "periodic matching"
        checkpoint(stage=result["stage"])
        pairs, classes, diagnostics = periodic_pairs(
            points, surface, config["periodic_pair_tolerance"]
        )
        result.update(diagnostics, periodic_eligible=True)
        np.savez_compressed(
            directory / "periodic.npz",
            classes=classes,
            **{f"pairs_{i}": pair for i, pair in enumerate(pairs)},
        )
        mesh = SimpleNamespace(
            points=points + 0.5,
            tetrahedra=tets,
            surface=surface,
            periodic_pairs=pairs,
            diagnostics={"minimum_mmg_quality": result["q_min"]},
        )
        result["stage"] = "periodic homogenization"
        checkpoint(stage=result["stage"], load=0)
        original_solver = homogenization.AuditedCG

        solver_count = 0

        class ProgressCG(original_solver):
            def __init__(self):
                nonlocal solver_count
                super().__init__()
                self.load_offset = 6 * solver_count
                solver_count += 1

            def __call__(self, matrix, rhs, **kwargs):
                checkpoint(
                    stage="periodic elasticity solve",
                    load=self.load_offset + len(self.records) + 1,
                    free_dofs=len(rhs),
                )
                solution = super().__call__(matrix, rhs, **kwargs)
                checkpoint(
                    load=self.load_offset + len(self.records),
                    residual=self.records[-1]["relative_residual"],
                )
                return solution

        with patch.object(homogenization, "AuditedCG", ProgressCG):
            result["homogenization"] = homogenization.homogenize(mesh, directory)
        result["homogenization_audited"] = True
        if task["case"] == "solid-cell":
            tensor = np.array(result["homogenization"]["stiffness_engineering_voigt"])
            error = float(
                np.linalg.norm(tensor - homogenization.MATERIAL)
                / np.linalg.norm(homogenization.MATERIAL)
            )
            result["analytic_stiffness_relative_error"] = error
            if error > 1e-8:
                raise ValueError(f"Solid-cell stiffness recovery failed: {error}")
        result["status"] = "ok"
        result["stage"] = "complete"
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        traceback.print_exc()
    result["worker_seconds"] = time.monotonic() - stage_start
    write_json(directory / "result.json", result)
    print(
        json.dumps(
            {k: result[k] for k in ["case", "workflow", "cells", "status", "stage"]}
            | {"error": result.get("error")}
        ),
        flush=True,
    )


def run(output):
    """Run frozen tests and one coarse batch in the serial heavy-compute slot."""
    import psutil

    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import report

    started = time.monotonic()
    report(0, "Verify frozen periodic experiment", eta_s=1800)
    source = output / "source"
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    for name, sha in manifest["source_sha256"].items():
        if digest(source / name) != sha:
            raise RuntimeError(f"Frozen source changed: {name}")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(source)
    env["PYTHONUNBUFFERED"] = "1"
    for name in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        env[name] = "1"
    tests = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", *[str(path) for path in sorted(source.glob("test_*.py"))]],
        env=env,
        capture_output=True,
        text=True,
    )
    (output / "tests.log").write_text(tests.stdout + tests.stderr, encoding="utf-8")
    print(tests.stdout + tests.stderr, flush=True)
    if tests.returncode:
        raise RuntimeError("Periodic matching/tensor regression checks failed")
    # Fixed analytic probes, independent of each generated candidate.
    sys.path.insert(0, str(source))
    import numpy as np
    from fixtures import make_case
    from microgen_benchmark.reference import boundary_samples
    from scipy.stats import qmc

    references = {}
    selected_cases = {task["case"] for task in manifest["tasks"]}
    for name, definition in manifest["specimens"].items():
        if name not in selected_cases:
            continue
        report(
            0.03, f"Analytic reference samples: {name}", eta_s=1800, cases=len(manifest["tasks"])
        )
        case = make_case(name, definition)
        fractions, nested = [], 0.0
        for seed in [20261010, 20261011, 20261012, 20261013]:
            sample = qmc.Sobol(3, scramble=True, seed=seed).random_base2(18) - 0.5
            material = case.material(*sample.T)
            fractions.append(float(material.mean()))
            nested = max(nested, abs(fractions[-1] - float(material[: 2**17].mean())))
        average = float(np.mean(fractions))
        uncertainty = max(nested, max(abs(fraction - average) for fraction in fractions))
        reference_path = output / f"reference-{name}.npz"
        np.savez_compressed(reference_path, boundary_points=boundary_samples(case))
        references[name] = {
            "fraction": average,
            "uncertainty": uncertainty,
            "samples": 4 * 2**18,
            "scramble_fractions": fractions,
            "scope": "Empirical spread/nested-sample change; not a confidence bound",
            "path": str(reference_path),
            "sha256": digest(reference_path),
        }
    document = {
        "state": "running",
        "manifest_sha256": digest(output / "manifest.json"),
        "references": references,
        "rows": [],
    }
    write_json(output / "results.json", document)
    tasks = manifest["tasks"]
    for index, task in enumerate(tasks):
        directory = output / f"{index:02d}-{task['case']}-{task['workflow']}-{task['cells']}"
        directory.mkdir()
        request = {
            "manifest": manifest,
            "task": task,
            "reference_path": references[task["case"]]["path"],
            "volume_reference": references[task["case"]],
        }
        write_json(directory / "request.json", request)
        worker_started = time.monotonic()
        live, peak_rss, stopped = {}, 0, None
        with (directory / "worker.log").open("w", encoding="utf-8") as log:
            process = subprocess.Popen(
                [
                    sys.executable,
                    "-u",
                    str(source / "run.py"),
                    "--worker",
                    str(directory / "request.json"),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
            )
            monitored = psutil.Process(process.pid)
            try:
                while process.poll() is None:
                    # Do not hold a Windows handle against atomically replaced JSON.
                    try:
                        with (directory / "progress.json").open(encoding="utf-8") as stream:
                            live = json.load(stream)
                    except (FileNotFoundError, PermissionError, json.JSONDecodeError):
                        pass
                    try:
                        rss = monitored.memory_info().rss
                        for child in monitored.children(recursive=True):
                            try:
                                rss += child.memory_info().rss
                            except psutil.NoSuchProcess:
                                pass
                        peak_rss = max(peak_rss, rss)
                    except psutil.NoSuchProcess:
                        pass
                    elapsed = time.monotonic() - worker_started
                    phase = live.get("stage", "starting worker")
                    substep = min(live.get("load", 0), 12) / 16
                    progress = 0.05 + 0.95 * (index + substep) / len(tasks)
                    eta = ((time.monotonic() - started) / max(index, 1)) * (len(tasks) - index)
                    report(
                        progress,
                        f"{task['workflow']} {task['cells']}: {phase}",
                        eta_s=eta,
                        rss_mib=round(peak_rss / 2**20),
                        seconds=round(elapsed),
                        load=live.get("load", 0),
                    )
                    if elapsed > manifest["config"]["worker_timeout_s"]:
                        stopped = "timeout"
                    if peak_rss > manifest["config"]["memory_limit_mib"] * 2**20:
                        stopped = "memory-limit"
                    if stopped:
                        break
                    time.sleep(2)
            finally:
                if process.poll() is None:
                    try:
                        children = monitored.children(recursive=True)
                    except psutil.NoSuchProcess:
                        children = []
                    for child in children:
                        try:
                            child.kill()
                        except psutil.NoSuchProcess:
                            pass
                    process.kill()
                    process.wait()
        path = directory / "result.json"
        if stopped or process.returncode or not path.is_file():
            result = {
                **task,
                "status": stopped or "worker-failed",
                "stage": live.get("stage"),
                "worker_exit_code": process.returncode,
            }
        else:
            result = json.loads(path.read_text(encoding="utf-8"))
        result["peak_rss_mib"] = peak_rss / 2**20
        result["artifact"] = directory.name
        document["rows"].append(result)
        write_json(output / "results.json", document)
        print(
            json.dumps(
                {
                    key: result.get(key)
                    for key in [
                        "case",
                        "workflow",
                        "cells",
                        "status",
                        "stage",
                        "error",
                        "mesh_s",
                        "q_min",
                        "geometry_eligible",
                        "quality_eligible",
                    ]
                }
            ),
            flush=True,
        )
        if task["case"] == "solid-cell" and result["status"] != "ok":
            document["state"] = "failed-foundation"
            write_json(output / "results.json", document)
            raise RuntimeError("Analytic solid-cell control failed; further measurements skipped")
    # Record coarse/fine disagreement only; no converged-reference claim.
    from periodic import tensor_difference

    for workflow in {task["workflow"] for task in tasks}:
        rows = [
            row
            for row in document["rows"]
            if row["case"] == "periodic-gyroid"
            and row["workflow"] == workflow
            and row.get("homogenization_audited")
        ]
        rows.sort(key=lambda row: row["cells"])
        for previous, current in zip(rows, rows[1:], strict=False):
            previous["difference_to_next_measured_level"] = {
                "reference_cells": current["cells"],
                **tensor_difference(
                    previous["homogenization"]["stiffness_engineering_voigt"],
                    current["homogenization"]["stiffness_engineering_voigt"],
                ),
            }
    document.update(
        state="verifying-artifacts",
        elapsed_s=time.monotonic() - started,
        accuracy_scope="Coarse/fine stability only; no converged reference or 2% qualification",
    )
    write_json(output / "results.json", document)
    from verify import verify_run

    report(0.99, "Reconstruct stiffness and equilibrium from saved fields", eta_s=60)
    document["artifact_verification"] = verify_run(output)
    document["state"] = "complete"
    write_json(output / "results.json", document)
    report(
        1,
        "Periodic coarse comparison complete",
        eta_s=0,
        audited=sum(bool(row.get("homogenization_audited")) for row in document["rows"]),
        failures=sum(row["status"] != "ok" for row in document["rows"]),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", type=Path)
    mode.add_argument("--run", type=Path)
    mode.add_argument("--worker", type=Path)
    parser.add_argument(
        "--native-package",
        type=Path,
        help="Previously verified packaged Meshers Python/native directory",
    )
    parser.add_argument("--tasks", nargs="+", help="Periodic Gyroid workflow:cells pairs")
    parser.add_argument("--mmg-hausd", type=float, default=0.002)
    parser.add_argument("--worker-timeout", type=float, default=1200)
    parser.add_argument("--max-tetrahedra", type=int, default=2_000_000)
    parser.add_argument("--mmg-seed-cells", type=int, help="Native seed grid, independent of MMG target sizing")
    parser.add_argument("--mmg-size-map", type=Path, help="Saved-point scalar sizing map for the metric workflow")
    parser.add_argument("--surface-seed-cells", type=int, help="Native surface seed grid, independent of HXT volume sizing")
    parser.add_argument("--native-build-manifest", type=Path, help="Verified isolated experimental build provenance")
    parser.add_argument("--specimens", type=Path, help="Frozen Microgen TPMS and spectral-field definitions")
    parser.add_argument("--geometry-tolerance", type=float, default=0.002)
    parser.add_argument("--native-geometry-tolerance", type=float, help="Native refinement request, independent of the shared geometry admission gate")
    args = parser.parse_args()
    if args.prepare:
        if args.native_package is None:
            parser.error("--prepare requires --native-package")
        if args.mmg_hausd <= 0 or args.worker_timeout <= 0 or args.max_tetrahedra <= 0 or args.geometry_tolerance <= 0:
            parser.error("--mmg-hausd, --worker-timeout, --max-tetrahedra and --geometry-tolerance must be positive")
        if args.mmg_seed_cells is not None and not 4 <= args.mmg_seed_cells <= 128:
            parser.error("--mmg-seed-cells must be between 4 and 128")
        if args.native_geometry_tolerance is not None and args.native_geometry_tolerance <= 0:
            parser.error("--native-geometry-tolerance must be positive")
        if args.surface_seed_cells is not None and not 4 <= args.surface_seed_cells <= 128:
            parser.error("--surface-seed-cells must be between 4 and 128")
        prepare(
            args.prepare.resolve(),
            args.native_package.resolve(),
            args.tasks,
            args.mmg_hausd,
            args.worker_timeout,
            args.max_tetrahedra,
            args.mmg_seed_cells,
            args.mmg_size_map,
            args.surface_seed_cells,
            args.native_build_manifest,
            args.specimens,
            args.geometry_tolerance,
            args.native_geometry_tolerance,
        )
    elif args.run:
        run(args.run.resolve())
    else:
        worker(args.worker.resolve())
