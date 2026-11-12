"""Same-host fixture comparison through Microgen's actual backend selectors.

A returned mesh must still pass the unchanged independent admission gates.
Constructor and backend-generation costs are charged; imports and audits are not.
"""

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import traceback
from unittest.mock import patch

import microgen_workflow as w

AUDIT_SOURCES = (
    "microgen_workflow.py",
    "comparison_mesh_audit.py",
    "periodic_trace.py",
    "facet_incidence_audit.py",
    "scorecard.py",
    "comparison_fixtures.py",
    "tpms_family_fixtures.py",
    "primitive_family_fixtures.py",
    "linear_density_fixture.py",
)


def factory(entry):
    module = importlib.import_module(
        {
            "study": "comparison_fixtures",
            "tpms": "tpms_family_fixtures",
            "primitive": "primitive_family_fixtures",
        }[entry["fixture_kind"]]
    )
    return module, Path(entry["fixtures"])


def worker(run, entry, arm, package, reuse=None):
    microgen_root = w.ROOT.parent / "microgen"
    sys.path[:0] = [str(package), str(microgen_root)]
    from vtkmodules.vtkCommonCore import vtkSMPTools

    vtkSMPTools.SetBackend("Sequential")
    np, _, meshers, microgen, pv, scorecard = w.imports()
    from comparison_mesh_audit import observe
    from periodic_trace import periodic_pairs
    from functools import partial

    assert Path(meshers.__file__).resolve().is_relative_to(package)
    assert Path(microgen.__file__).resolve().is_relative_to(microgen_root / "microgen")
    assert vtkSMPTools.GetBackend() == "Sequential"
    module, fixtures = factory(entry)
    name = entry["case"]
    target = module.case(name, fixtures)
    periodic = entry["periodic"]
    config = entry[arm]
    assert config.get("order", 1) == 1, "This benchmark compares Tet4 only"
    protocol = {name: w.digest(w.HERE / name) for name in AUDIT_SOURCES}
    protocol["driver"] = w.digest(Path(__file__))
    protocol.update({
        str(path.relative_to(microgen_root)): w.digest(path)
        for path in (microgen_root / "microgen").rglob("*.py")
    })
    row = dict(
        case=name,
        arm=arm,
        config=config,
        status="started",
        qualified=False,
        periodic_required=periodic,
        audit_protocol_sha256=hashlib.sha256(
            json.dumps(protocol, sort_keys=True).encode()
        ).hexdigest(),
        backend_selector="meshers" if arm == "native" else "mmg",
        mmg_calls=[],
        constructor_s=None,
        generation_s=None,
        workflow_s=None,
        instrumentation_s=0.0,
        imports_charged=False,
        audits_charged=False,
        fea_run=False,
        versions={
            n: importlib.metadata.version(n)
            for n in ["meshers", "microgen", "numpy", "scipy", "pyvista", "vtk", "mmgpy"]
        },
        native_binary_sha256=w.digest(Path(meshers._meshers.__file__)),
        microgen_module=str(Path(microgen.__file__).resolve()),
        vtk_backend=vtkSMPTools.GetBackend(),
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY,
    )
    run.mkdir(parents=True, exist_ok=False)
    for source in (fixtures / name).iterdir():
        if source.is_file():
            shutil.copyfile(source, run / source.name)
    additional = []
    if name == "spinodoid-anisotropic":
        additional.append(fixtures / "spinodoid-anisotropic-coefficients.json")
    if name == "graded-gyroid-density-linear":
        additional.append(fixtures / "linear-density-calibration.json")
    if name == "spinodoid":
        additional.append(w.HERE / "runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json")
    for source in additional:
        if source.is_file():
            shutil.copyfile(source, run / source.name)
    row["fixture_hashes"] = {p.name: w.digest(p) for p in run.iterdir()}
    w.save(run / "result.json", row)
    overall = time.perf_counter()
    try:
        start = time.perf_counter()
        shape, part = module.actual_shape(name, fixtures, config.get("resolution"))
        row["constructor_s"] = time.perf_counter() - start
        # Confirm the unchanged physical specimen outside the charged interval.
        probes = np.random.default_rng(20261112).uniform(-0.5, 0.5, (4096, 3))
        expected = target.function(*probes.T)
        inside = (expected >= target.levels[0]) & (expected <= target.levels[1])
        if isinstance(shape, microgen.Tpms):
            values = shape.raw_field(*probes.T)
            half = 0.5 * (
                shape._offset_func(*probes.T) if shape._offset_func is not None else shape.offset
            )
            actual_inside = (
                np.abs(values) <= half
                if part == "sheet"
                else values <= -half
                if part == "lower skeletal"
                else values >= half
            )
        else:
            local = probes + 0.5 if name.startswith("spinodoid") else probes
            values = shape.require_func()(*local.T)
            actual_inside = values <= 0
        np.testing.assert_array_equal(inside, actual_inside)
        row["field_identity_points"] = len(probes)
        original_mmg = microgen.Mmg.mmg3d
        remesh_module = importlib.import_module("microgen.remesh")
        original_remesh = remesh_module.remesh_keeping_boundaries_for_fem
        binary = Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))
        os.environ["PATH"] = str(binary.parent) + os.pathsep + os.environ["PATH"]
        row["mmg_binary_sha256"] = w.digest(binary)

        def checkpoint():
            begin = time.perf_counter()
            w.save(run / "result.json", row)
            row["instrumentation_s"] += time.perf_counter() - begin

        def mmg(**parameters):
            call = {k: v for k, v in parameters.items() if k not in ("input", "output")}
            row["mmg_calls"].append(call)
            checkpoint()
            result = original_mmg(**parameters)
            begin = time.perf_counter()
            stage = pv.UnstructuredGrid(str(parameters["output"]))
            call.update(
                nodes=stage.n_points,
                tetrahedra=int(np.count_nonzero(stage.celltypes == pv.CellType.TETRA)),
            )
            row["instrumentation_s"] += time.perf_counter() - begin
            checkpoint()
            if max(call["nodes"], call["tetrahedra"]) > w.CAP:
                raise ValueError("Intermediate MMG node/tet cap exceeded")
            return result

        def remesh(seed, **options):
            row.update(seed_nodes=seed.n_points, seed_cells=seed.n_cells, mmg_options=options)
            checkpoint()
            if max(seed.n_points, seed.n_cells) > w.CAP:
                raise ValueError("Actual VTK seed node/tet cap exceeded")
            return original_remesh(seed, **options)

        options = dict(type_part=part) if isinstance(shape, microgen.Tpms) else {}
        if arm == "native":
            options.update(
                cells=config["cells"],
                order=1,
                adapt=config.get("adapt", False),
                geometry_tolerance=config.get("tolerance", 1.0),
                minimum_quality=0.1,
                max_tetrahedra=w.CAP,
                threads=1,
                optimize_passes=config.get("optimization_passes", 1),
                periodic=(periodic,) * 3,
            )
            for key in [
                "min_length",
                "max_length",
                "gradation",
                "adapt_tolerance",
                "adapt_iterations",
            ]:
                if key in config:
                    options[key] = config[key]
            if not isinstance(shape, (microgen.Tpms, microgen.Spinodoid)):
                options["bounds"] = w.BOUNDS
        elif not periodic:
            options["periodic"] = False
        start = time.perf_counter()
        with (
            patch.object(microgen.Mmg, "mmg3d", mmg),
            patch.object(remesh_module, "remesh_keeping_boundaries_for_fem", remesh),
        ):
            grid = shape.generate_volume_mesh(backend=row["backend_selector"], **options)
        row["generation_s"] = time.perf_counter() - start - row["instrumentation_s"]
        row["workflow_s"] = row["constructor_s"] + row["generation_s"]
        assert row["generation_s"] >= 0
        assert (len(row["mmg_calls"]) == 0) if arm == "native" else (len(row["mmg_calls"]) > 0)
        if name.startswith("spinodoid"):
            grid.translate((-0.5,) * 3, inplace=True)
        points = np.asarray(grid.points).copy()
        cells = np.asarray(grid.cells_dict[pv.CellType.TETRA]).copy()
        row.update(nodes=len(points), elements=len(cells))
        if max(len(points), len(cells)) > w.CAP:
            raise ValueError("Final node/tet cap exceeded")
        if arm == "native":
            surface = np.asarray(grid.field_data["meshers_boundary_faces"]).copy()
        else:
            poly = grid.extract_surface(
                pass_pointid=True, algorithm="dataset_surface"
            ).triangulate()
            surface = np.asarray(poly.point_data["vtkOriginalPointIds"])[
                poly.faces.reshape(-1, 4)[:, 1:]
            ]
        cells, row["audit_reoriented_tets"] = w.reorient_linear(points, cells)
        np.savez_compressed(run / "mesh.npz", points=points, tetrahedra=cells, surface=surface)
        digest = hashlib.sha256()
        for array in [points, cells, surface]:
            digest.update(str((array.shape, array.dtype.str)).encode())
            digest.update(array.tobytes())
        row["arrays_sha256"] = digest.hexdigest()
        row["mesh_sha256"] = w.digest(run / "mesh.npz")
        start = time.perf_counter()
        previous = json.loads(reuse.read_text()) if reuse else None
        auditor = (
            partial(periodic_pairs, declared_bounds=([-0.5] * 3, [0.5] * 3)) if periodic else None
        )
        if (
            previous
            and previous.get("status") in ("admitted", "admission-rejected")
            and previous.get("arrays_sha256") == row["arrays_sha256"]
            and previous.get("case") == name
            and previous.get("arm") == arm
            and previous.get("fixture_hashes") == row["fixture_hashes"]
            and previous.get("audit_protocol_sha256") == row["audit_protocol_sha256"]
            and previous.get("periodic_required") == periodic
        ):
            for key in ["geometry", "periodic", "quality_stats", "qualified"]:
                row[key] = previous[key]
            row["audit_reused_from"] = str(reuse)
            row["audit_reused_sha256"] = w.digest(reuse)
        else:
            row["geometry"], row["periodic"], quantiles = observe(
                run, name, target, points, cells, surface, periodic, periodic_auditor=auditor
            )
            q, _, _ = scorecard.tet_metrics(points, cells)
            row["quality_stats"] = {**quantiles, "mean": float(q.mean())}
            row["qualified"] = row["geometry"]["geometry_admitted"]
        if not periodic:
            row["periodic"]["reason"] = (
                "Declared bounded nonperiodic specimen; no periodic matching required"
            )
        if arm == "native" and periodic and "error" not in row["periodic"]:
            independently, _, _ = auditor(points, surface)
            for axis, pairs in zip("xyz", independently, strict=True):
                key = f"periodic_pairs_{axis}"
                original = (
                    np.asarray(grid.field_data[key])
                    if key in grid.field_data
                    else np.empty((0, 2), dtype=np.int64)
                )
                assert {tuple(x) for x in pairs} == {tuple(x) for x in original}
            row["native_pairs_match_independent_audit"] = True
        row["audit_s"] = time.perf_counter() - start
        row["status"] = "admitted" if row["qualified"] else "admission-rejected"
    except Exception as error:
        row.update(
            status="failed",
            qualified=False,
            error_type=type(error).__name__,
            error=str(error),
            traceback=traceback.format_exc(),
        )
    row["elapsed_s"] = time.perf_counter() - overall
    w.save(run / "result.json", row)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--selection", type=Path, required=True)
    p.add_argument("--package", type=Path, required=True)
    p.add_argument("--worker", type=int)
    p.add_argument("--arm", choices=["native", "actual"])
    p.add_argument("--reuse", type=Path)
    p.add_argument("--repetitions", type=int, default=1)
    p.add_argument("--warmups", type=int, default=0)
    a = p.parse_args()
    assert a.repetitions > 0 and a.warmups >= 0
    entries = json.loads(a.selection.read_text())["cases"]
    package = a.package.resolve()
    if a.worker is not None:
        return worker(a.run, entries[a.worker], a.arm, package, a.reuse)
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job
    from job_parts import worker_part

    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run / "driver.py")
    shutil.copyfile(a.selection, a.run / "selection.json")
    source_dir = a.run / "audit-source"
    source_dir.mkdir()
    for name in AUDIT_SOURCES:
        shutil.copyfile(w.HERE / name, source_dir / name)
    manifest = dict(
        platform=platform.platform(),
        machine=platform.machine(),
        package=str(package),
        microgen_commit=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=w.ROOT.parent / "microgen", text=True
        ).strip(),
        microgen_source_sha256={
            str(path.relative_to(w.ROOT.parent / "microgen")): w.digest(path)
            for path in (w.ROOT.parent / "microgen/microgen").rglob("*.py")
        },
        selection_sha256=w.digest(a.selection),
        script_sha256=w.digest(Path(__file__)),
        audit_source_sha256={name: w.digest(source_dir / name) for name in AUDIT_SOURCES},
        repetitions=a.repetitions,
        warmups=a.warmups,
        case_count=len(entries),
        caps=dict(nodes=w.CAP, tetrahedra=w.CAP),
        admission=dict(distance=0.002, relative_volume=0.01, quality=0.1, periodic_tolerance=1e-10),
        all_four_original_cases_retained=all(n in {x["case"] for x in entries} for n in w.CASES),
        all_24_lattice_models_retained=True,
        lattice_models_evaluated_here=0,
        qualified_scope="Mesh admission only; no FEA accuracy claim",
        timing_rule="Successful generation and worker; absolute realtime-minus-monotonic gap <= 0.1 seconds",
        new_fea=False,
        reference_promoted=False,
    )
    w.save(a.run / "manifest.json", manifest)
    rows = []
    audit_cache = {}
    tasks = [
        (rep, index, arm)
        for rep in range(-a.warmups, a.repetitions)
        for index in range(len(entries))
        for arm in (["native", "actual"] if rep % 2 == 0 else ["actual", "native"])
    ]
    parts = [f"{rep}-{index}-{arm}" for rep, index, arm in tasks]
    with Job(parts=parts) as job:
        for task, (rep, index, arm) in zip(parts, tasks, strict=True):
            entry = entries[index]
            directory = a.run / f"repeat-{rep}-{entry['case']}-{arm}"
            base = audit_cache.get((index, arm))
            if base is None and a.reuse:
                base = a.reuse / f"repeat-0-{entry['case']}-{arm}/result.json"
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--run",
                str(directory),
                "--selection",
                str(a.selection.resolve()),
                "--package",
                str(package),
                "--worker",
                str(index),
                "--arm",
                arm,
            ]
            if base and base.exists():
                command += ["--reuse", str(base)]
            with worker_part(job, task, memory_bytes=12 * 1024**3) as part:
                result = w.bounded(
                    command, a.run / (task + ".log"), part, timeout=300, memory_bytes=12 * 1024**3
                )
                row = (
                    json.loads((directory / "result.json").read_text())
                    if (directory / "result.json").exists()
                    else dict(case=entry["case"], arm=arm, status="worker-failed", qualified=False)
                )
                row.update(
                    repetition=rep,
                    warmup=rep < 0,
                    resource=result,
                    directory=str(directory.resolve()),
                )
                if result["returncode"]:
                    row.update(status="worker-failed", qualified=False)
                row["timing_valid"] = (
                    result["returncode"] == 0
                    and row.get("workflow_s") is not None
                    and abs(result["suspension_or_clock_adjustment_s"]) <= 0.1
                )
                rows.append(row)
                if row["status"] in ("admitted", "admission-rejected"):
                    audit_cache[index, arm] = directory / "result.json"
                w.save(
                    a.run / "result.json",
                    dict(complete=False, rows=rows, new_fea=False, reference_promoted=False),
                )
                job.result(
                    task,
                    status=row["status"],
                    qualified=row["qualified"],
                    seconds=result["wall_s"],
                    peak_gib=result["peak_rss_bytes"] / 1024**3,
                )
        w.save(
            a.run / "result.json",
            dict(
                complete=True,
                rows=rows,
                returned=sum("elements" in r for r in rows),
                admitted=sum(r["qualified"] for r in rows),
                new_fea=False,
                reference_promoted=False,
            ),
        )


if __name__ == "__main__":
    main()
