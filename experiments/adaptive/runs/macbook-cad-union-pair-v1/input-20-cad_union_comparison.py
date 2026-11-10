"""One bounded Tet4 pair on the exact crossed-cylinder CustomLattice.

The actual arm records CAD/Gmsh as a separate candidate before the unmodified
public two-pass MMG workflow. Admission uses the existing shared observers.
"""

import argparse
from contextlib import contextmanager
from functools import partial
import importlib
import importlib.metadata
import inspect
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import tempfile
import time
import traceback
from unittest.mock import patch

import cad_coverage_preflight as fixture
import microgen_workflow as w


NAME = "cad-crossed-flat-cylinders"
CAD_PIN = "7.9.3.1.1"
NATIVE = dict(bounds=fixture.BOUNDS, cells=16, periodic=(True, True, True),
              threads=1, order=1, adapt=True, adapt_iterations=4,
              optimize_passes=1, geometry_tolerance=.002, adapt_tolerance=.002,
              min_length=.002, max_length=.1, gradation=.3, quality_guard=.3,
              minimum_quality=.1, max_tetrahedra=w.CAP)
GMSH = dict(size=.06, order=1, msh_file_version=4, tol=1e-8)
MMG = dict(periodic=True, mesh_version=2, dimension=3, tol=1e-8,
           hausd=None, hgrad=None, hmax=None, hmin=None, hsiz=None)


def target():
    import numpy as np
    import scorecard

    volume = 2 * np.pi * fixture.RADIUS**2 - 16 * fixture.RADIUS**3 / 3
    return scorecard.Case(NAME, fixture.periodic_field, (-np.inf, 0), volume)


def prepare(run, preflight):
    """Freeze analytic probes, independent of either generated mesh."""
    import numpy as np
    from scipy.stats import qmc

    reference = json.loads((preflight / "result.json").read_text())
    assert reference["cad"]["status"] == "verified"
    assert reference["cad"]["valid"] and reference["cad"]["membership_points"] == 512
    r = fixture.RADIUS
    uv = qmc.Sobol(2, scramble=True, seed=20261110).random_base2(13)
    angle = 2 * np.pi * uv[:, 0]
    transverse, z = r * np.cos(angle), r * np.sin(angle)
    start = np.abs(transverse)
    axial = start + uv[:, 1] * (.5 - start)
    walls = []
    for axis in (0, 1):
        for sign in (-1, 1):
            points = np.column_stack((sign * axial, transverse, z))
            if axis == 1:
                points = points[:, [1, 0, 2]]
            walls.append(points)
    angle = np.linspace(0, 2 * np.pi, 2048, endpoint=False)
    ridges = np.vstack([np.column_stack((r * np.cos(angle), sign * r * np.cos(angle),
                                        r * np.sin(angle))) for sign in (-1, 1)])
    uv = qmc.Sobol(2, scramble=True, seed=20261111).random_base2(11)
    radial, cap_angle = r * np.sqrt(uv[:, 0]), 2 * np.pi * uv[:, 1]
    caps, rims = [], []
    for axis in (0, 1):
        for sign in (-1, 1):
            cap = np.column_stack((np.full(len(uv), sign * .5),
                                   radial * np.cos(cap_angle), radial * np.sin(cap_angle)))
            rim = np.column_stack((np.full(len(angle), sign * .5),
                                   r * np.cos(angle), r * np.sin(angle)))
            if axis == 1:
                cap, rim = cap[:, [1, 0, 2]], rim[:, [1, 0, 2]]
            caps.append(cap)
            rims.append(rim)
    walls = np.vstack([*walls, ridges, *rims])
    caps = np.vstack(caps)
    assert np.max(np.abs(fixture.periodic_field(*walls.T))) < 1e-14
    assert np.all(fixture.periodic_field(*caps.T) < 0)
    assert np.max(np.abs(np.vstack((walls, caps)))) <= .5
    np.testing.assert_allclose(target().volume, reference["fixture"]["analytic_volume"], rtol=0, atol=0)
    sample_file = run / (NAME + "-samples.npz")
    np.savez_compressed(sample_file, walls=walls, caps=caps, ridges=ridges, rims=np.vstack(rims))
    w.save(run / (NAME + "-specimen.json"), dict(volume=target().volume,
        volume_reference="analytic: 2*pi*r^2*L - 16*r^3/3", volume_scramble_spread=0,
        wall_samples=len(walls), cap_samples=len(caps), periodic_axes=[True] * 3,
        declared_bounds=[[-.5] * 3, [.5] * 3], definition=reference["fixture"]))
    step = Path(reference["cad"]["step_file"])
    assert w.digest(step) == reference["cad"]["step_sha256"]
    shutil.copyfile(step, run / "verified-cad.step")
    w.save(run / "result.json", dict(prepared=True, samples_sha256=w.digest(sample_file),
        verified_cad_sha256=w.digest(run / "verified-cad.step"), wall_samples=len(walls),
        cap_samples=len(caps), ridge_samples=len(ridges), rim_samples=sum(map(len, rims)),
        analytic_volume=target().volume, original_samples_unchanged=True))


@contextmanager
def phase(run, row, name):
    row["active_phase"] = name
    w.save(run / "result.json", row)
    started = time.perf_counter()
    try:
        yield
    finally:
        row["phase_seconds"][name] = time.perf_counter() - started
        w.save(run / "result.json", row)


def extract(grid):
    import numpy as np
    import pyvista as pv

    if set(grid.cells_dict) != {int(pv.CellType.TETRA)}:
        raise ValueError("Actual output is not exclusively Tet4")
    points = np.asarray(grid.points).copy()
    cells = np.asarray(grid.cells_dict[pv.CellType.TETRA]).copy()
    poly = grid.extract_surface(pass_pointid=True, algorithm="dataset_surface").triangulate()
    surface = np.asarray(poly.point_data["vtkOriginalPointIds"])[poly.faces.reshape(-1, 4)[:, 1:]]
    return points, cells, surface


def archive_and_audit(run, label, arrays, *, reorient=False, native_pairs=None):
    """Save arrays before any admission checks, then use the shared observer."""
    import numpy as np
    from comparison_mesh_audit import observe
    from periodic_trace import periodic_pairs

    points, cells, surface = arrays
    directory = run / label
    directory.mkdir()
    for suffix in ("-specimen.json", "-samples.npz"):
        shutil.copyfile(run / (NAME + suffix), directory / (NAME + suffix))
    mesh_file = directory / "mesh.npz"
    np.savez_compressed(mesh_file, points=points, tetrahedra=cells, surface=surface)
    result = dict(nodes=len(points), elements=len(cells), order=1, qualified=False,
                  mesh_file=str(mesh_file.resolve()), raw_mesh_sha256=w.digest(mesh_file))
    w.save(directory / "result.json", result)
    if max(len(points), len(cells)) > w.CAP:
        result.update(status="cap-rejected", error="Node or tetrahedron cap exceeded")
        w.save(directory / "result.json", result)
        return result
    if reorient:
        shutil.copyfile(mesh_file, directory / "raw-mesh.npz")
        started = time.perf_counter()
        cells, result["reoriented_tetrahedra"] = w.reorient_linear(points, cells)
        result["reorientation_s"] = time.perf_counter() - started
        np.savez_compressed(mesh_file, points=points, tetrahedra=cells, surface=surface)
    result["mesh_sha256"] = w.digest(mesh_file)
    started = time.perf_counter()
    try:
        check = partial(periodic_pairs, declared_bounds=([-.5] * 3, [.5] * 3))
        result["geometry"], result["periodic"], result["quality_quantiles"] = observe(
            directory, NAME, target(), points, cells, surface, True, periodic_auditor=check)
        result["qualified"] = result["geometry"]["geometry_admitted"]
        if native_pairs is not None:
            try:
                independent, _, _ = check(points, surface)
                assert all({tuple(map(int, p)) for p in supplied} ==
                           {tuple(map(int, p)) for p in checked}
                           for supplied, checked in zip(native_pairs, independent, strict=True))
                result["native_and_independent_pairs_agree"] = True
            except (ValueError, AssertionError) as error:
                result.update(qualified=False, native_pair_error=str(error))
        result["status"] = "admitted" if result["qualified"] else "admission-rejected"
    except Exception as error:
        result.update(status="observer-failed", error=str(error), traceback=traceback.format_exc())
    result["audit_s"] = time.perf_counter() - started
    w.save(directory / "result.json", result)
    return result


def actual_worker(run, row):
    import numpy as np
    import gmsh
    import microgen
    import pyvista as pv
    from microgen.remesh import remesh_keeping_boundaries_for_fem

    with phase(run, row, "shape_setup"):
        shape = fixture.make_shape()
    with phase(run, row, "cad_generate"):
        cad = shape.generate_cad()
    with phase(run, row, "cad_volume_check"):
        row["cad_volume"] = cad.volume()
        np.testing.assert_allclose(row["cad_volume"], target().volume, rtol=1e-8, atol=1e-10)
    step = run / "actual-cad.step"
    with phase(run, row, "step_export"):
        cad.export_step(str(step))
    row["cad_step_sha256"] = w.digest(step)
    seed = run / "gmsh-seed.vtk"
    original_initialize = gmsh.initialize

    def initialize_one_thread(*args, **kwargs):
        result = original_initialize(*args, **kwargs)
        for option in ("General.NumThreads", "Mesh.MaxNumThreads1D", "Mesh.MaxNumThreads2D", "Mesh.MaxNumThreads3D"):
            gmsh.option.setNumber(option, 1)
        return result

    public_error = None
    with phase(run, row, "public_cad_gmsh"):
        try:
            with patch.object(gmsh, "initialize", initialize_one_thread):
                microgen.mesh_periodic(mesh_file=str(step), rve=shape.rve,
                    list_phases=[microgen.Phase(cad)], output_file=str(seed), **GMSH)
        except Exception as error:
            public_error = error
            row["gmsh_public_error"] = dict(type=type(error).__name__, message=str(error),
                                             traceback=traceback.format_exc())
        finally:
            if gmsh.isInitialized():
                gmsh.finalize()
    if not seed.exists() or not seed.stat().st_size:
        raise public_error or ValueError("Public CAD/Gmsh produced no mesh file")
    row["gmsh_seed_sha256"] = w.digest(seed)
    with phase(run, row, "gmsh_read"):
        grid = pv.read(seed)
        row.update(seed_nodes=grid.n_points, seed_elements=grid.n_cells)
        if max(grid.n_points, grid.n_cells) > w.CAP:
            raise ValueError("Gmsh seed cap exceeded; original VTK file retained")
    row["cad_gmsh_read_s"] = sum(row["phase_seconds"].values())
    with phase(run, row, "gmsh_seed_extract"):
        arrays = extract(grid)
    row["direct_generation_s"] = sum(row["phase_seconds"].values())
    with phase(run, row, "direct_archive_audit"):
        row["direct"] = archive_and_audit(run, "actual-cad-gmsh-direct", arrays, reorient=True)
        row["direct_generation_s"] += row["direct"].get("reorientation_s", 0.)
        row["direct"].update(route="Actual Microgen CAD/STEP/Gmsh, without MMG",
                             generation_s=row["direct_generation_s"])
        w.save(run / "actual-cad-gmsh-direct/result.json", row["direct"])
    if public_error is not None:
        row["direct"]["qualified"] = False
        row["direct"]["status"] = "public-route-failed"
        row["direct"]["public_route_error"] = row["gmsh_public_error"]
        w.save(run / "actual-cad-gmsh-direct/result.json", row["direct"])
        raise public_error

    binary = Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))
    os.environ["PATH"] = str(binary.parent) + os.pathsep + os.environ["PATH"]
    original = microgen.Mmg.mmg3d

    def record_mmg(**parameters):
        call = dict(arguments={k: v for k, v in parameters.items() if k not in ("input", "output")},
                    status="started")
        row["mmg_calls"].append(call)
        saved = run / ("mmg-call-" + str(len(row["mmg_calls"])))
        saved.mkdir()

        def retain(which):
            source = Path(parameters[which])
            files = {}
            for path in (source, source.with_suffix(".sol")):
                if path.exists() and path.stat().st_size:
                    dest = saved / (which + path.suffix)
                    shutil.copyfile(path, dest)
                    files[str(dest.resolve())] = w.digest(dest)
            call[which + "_artifacts"] = files

        started = time.perf_counter()
        retain("input")
        w.save(run / "result.json", row)
        row["mmg_instrumentation_s"] += time.perf_counter() - started
        started = time.perf_counter()
        failed = None
        try:
            result = original(**parameters)
            call["status"] = "returned"
        except Exception as error:
            failed = error
            call.update(status="failed", error=str(error), traceback=traceback.format_exc())
            output = getattr(error.__cause__, "output", b"")
            if output:
                (saved / "command.log").write_bytes(output if isinstance(output, bytes) else output.encode())
        call["command_s"] = time.perf_counter() - started
        started = time.perf_counter()
        retain("output")
        output = Path(parameters["output"])
        if output.exists() and output.stat().st_size:
            try:
                stage = pv.UnstructuredGrid(str(output))
                call.update(output_nodes=stage.n_points,
                    output_tetrahedra=int(np.count_nonzero(stage.celltypes == pv.CellType.TETRA)))
            except Exception as error:
                call["output_parse_error"] = str(error)
        w.save(run / "result.json", row)
        row["mmg_instrumentation_s"] += time.perf_counter() - started
        # Stop on an intermediate cap violation even if public MMG would retry.
        call["within_node_tetrahedron_cap"] = (
            max(call.get("output_nodes", 0), call.get("output_tetrahedra", 0)) <= w.CAP
            if "output_nodes" in call else None)
        if call["within_node_tetrahedron_cap"] is False:
            raise ValueError("Intermediate MMG cap exceeded; files retained")
        if failed is not None:
            raise failed
        if "output_nodes" not in call:
            raise ValueError("MMG output could not be inspected; files retained")
        return result

    with phase(run, row, "public_mmg"):
        with patch.object(microgen.Mmg, "mmg3d", record_mmg):
            grid = remesh_keeping_boundaries_for_fem(grid, **MMG)
    row["mmg_generation_s"] = row["phase_seconds"]["public_mmg"] - row["mmg_instrumentation_s"]
    with phase(run, row, "final_grid_archive"):
        grid.save(run / "mmg-final.vtu")
        row["final_grid_sha256"] = w.digest(run / "mmg-final.vtu")
    with phase(run, row, "final_extract"):
        if max(grid.n_points, grid.n_cells) > w.CAP:
            raise ValueError("Final MMG cap exceeded; VTU file retained")
        arrays = extract(grid)
    row["generation_s"] = row["cad_gmsh_read_s"] + row["mmg_generation_s"] + row["phase_seconds"]["final_extract"]
    with phase(run, row, "final_archive_audit"):
        row["final"] = archive_and_audit(run, "actual-cad-gmsh-public-mmg", arrays, reorient=True)
        row["generation_s"] += row["final"].get("reorientation_s", 0.)
        row["final"].update(route="Actual Microgen CAD/STEP/Gmsh plus public MMG, both passes and retries unchanged",
                            generation_s=row["generation_s"])
        w.save(run / "actual-cad-gmsh-public-mmg/result.json", row["final"])
    row["qualified"] = row["final"]["qualified"]


def worker(run, arm):
    np, _, meshers, _, _, _ = w.imports()
    from vtkmodules.vtkCommonCore import vtkSMPTools

    assert vtkSMPTools.GetBackend() == "Sequential" and vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    (run / "temporary").mkdir()
    tempfile.tempdir = str((run / "temporary").resolve())
    row = dict(arm=arm, status="started", qualified=False, phase_seconds={}, mmg_calls=[],
               mmg_instrumentation_s=0., vtk_backend="Sequential", numerical_threads=1,
               solve_performed=False)
    w.save(run / "result.json", row)
    started = time.perf_counter()
    try:
        if arm == "native":
            with phase(run, row, "native_generate"):
                mesh = meshers.generate(fixture.periodic_field, **NATIVE)
            row["diagnostics"] = mesh.diagnostics
            row["generation_s"] = row["phase_seconds"]["native_generate"]
            with phase(run, row, "native_archive_audit"):
                row["final"] = archive_and_audit(run, "native-tet4", (
                    np.asarray(mesh.points), np.asarray(mesh.tetrahedra), np.asarray(mesh.surface)),
                    native_pairs=mesh.periodic_pairs)
                row["final"]["generation_s"] = row["generation_s"]
                w.save(run / "native-tet4/result.json", row["final"])
            row["qualified"] = row["final"]["qualified"]
        else:
            actual_worker(run, row)
        row["status"] = "admitted" if row["qualified"] else "admission-rejected"
    except Exception as error:
        row.update(status="failed", error_type=type(error).__name__, error=str(error),
                   traceback=traceback.format_exc(), qualified=False)
        if arm == "native" and "final" not in row:
            row["mesh_retention"] = "Native API raised without returning a mesh; no partial mesh is exposed"
    row["wall_s"] = time.perf_counter() - started
    w.save(run / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--worker", choices=("fixture", "native", "actual"))
    args = parser.parse_args()
    if args.worker == "fixture":
        return prepare(args.run, args.preflight)
    if args.worker:
        return worker(args.run, args.worker)
    from job_parts import worker_part
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    assert importlib.metadata.version("cadquery-ocp-novtk") == CAD_PIN
    preflight = json.loads((args.preflight / "result.json").read_text())
    assert preflight["cad"]["status"] == "verified", "Exact actual CAD preflight must pass first"
    for path, digest in preflight["inputs"].items():
        assert w.digest(Path(path)) == digest, "Preflight source changed: " + path
    args.run.mkdir(parents=True, exist_ok=False)
    paths = {Path(__file__).resolve(), Path(fixture.__file__).resolve(), Path(w.__file__).resolve(),
             args.preflight / "result.json", Path(preflight["cad"]["step_file"])}
    paths.update(Path(p) for p in preflight["inputs"])
    paths.update(w.HERE / name for name in ("comparison_mesh_audit.py", "scorecard.py", "periodic_trace.py",
        "facet_incidence_audit.py", "quadratic_edge_removal_probe.py", "nested_refinement.py",
        "box_projection.py", "job_parts.py"))
    for module in ("microgen.mesh", "microgen.remesh", "microgen.external", "microgen.operations",
                   "microgen.phase", "microgen.rve", "microgen.box_mesh", "microgen.single_mesh", "meshers._meshers"):
        paths.add(Path(inspect.getfile(importlib.import_module(module))))
    binary = Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))
    assert os.access(binary, os.X_OK), "Installed MMG executable is not executable"
    paths.add(binary)
    hashes = {str(p.resolve()): w.digest(p) for p in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run / f"input-{index}-{path.name}")
    w.save(args.run / "manifest.json", dict(inputs=hashes, preflight=str(args.preflight.resolve()),
        cad_version_pin=CAD_PIN, packages={n: importlib.metadata.version(n) for n in
            ("microgen", "cadquery-ocp-novtk", "gmsh", "mmgpy", "numpy", "scipy", "pyvista", "vtk", "meshers")},
        platform=platform.platform(), architecture=platform.machine(), command=sys.argv,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, native=NATIVE, actual_gmsh=GMSH, actual_mmg=MMG,
        worker_timeout_s=180, memory_bytes=w.MEMORY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        timing_scope="One fresh pilot per arm. Generation includes actual CAD construction, volume check, STEP export, public Gmsh, public MMG, array extraction and reorientation. Imports, fixture preparation, archive, observer and measured MMG instrumentation are excluded.",
        actual_route="Public Microgen CAD/STEP/Gmsh seed followed by unmodified public two-pass MMG. Direct Gmsh candidate is reported separately. No repair, remesh parameter tuning or periodic-validator override.",
        admission="Existing comparison_mesh_audit.observe with declared unit-cell periodic_trace observer; .002 distance, 1% absolute relative volume, .1 quality, positive Jacobians, coherent facets, matching full caps and node classes.",
        mechanical_accuracy_claim=False, repetitions=1, warmups=0))
    rows = []
    fixture_dir = args.run / "fixture"
    fixture_hashes = {}
    with Job(parts={"fixture": 1, "native": 3, "actual": 3}) as job:
        for arm, timeout in (("fixture", 60), ("native", 180), ("actual", 180)):
            directory = args.run / arm
            directory.mkdir()
            if arm != "fixture":
                assert all(w.digest(Path(p)) == digest for p, digest in fixture_hashes.items())
                for suffix in ("-samples.npz", "-specimen.json"):
                    shutil.copyfile(fixture_dir / (NAME + suffix), directory / (NAME + suffix))
            with worker_part(job, arm) as part:
                record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--run", str(directory.resolve()),
                    "--preflight", str(args.preflight.resolve()), "--worker", arm],
                    directory / "worker.log", part, timeout=timeout)
                w.save(directory / "worker.json", record)
                row = (json.loads((directory / "result.json").read_text()) if (directory / "result.json").exists()
                       else dict(status="no-result", qualified=False))
                row.update(directory=str(directory.resolve()), worker=record)
                if record["returncode"] != 0:
                    row.update(status="worker-failed", qualified=False)
                job.result(arm, wall_s=record["wall_s"], generation_s=row.get("generation_s"),
                           qualified=row.get("qualified", row.get("prepared", False)),
                           peak_gib=record["peak_rss_bytes"] / 1024**3)
                if arm == "fixture":
                    assert record["returncode"] == 0 and row.get("prepared"), "Read fixture worker log"
                    fixture_hashes = {str(p.resolve()): w.digest(p) for p in (
                        fixture_dir / (NAME + "-samples.npz"), fixture_dir / (NAME + "-specimen.json"),
                        fixture_dir / "verified-cad.step")}
                    manifest = json.loads((args.run / "manifest.json").read_text())
                    manifest["fixture_hashes"] = fixture_hashes
                    w.save(args.run / "manifest.json", manifest)
                else:
                    rows.append(row)
                    w.save(args.run / "summary.json", dict(rows=rows, complete=len(rows) == 2,
                        mechanical_accuracy_claim=False, successful_performance_ratio=None,
                        scope="Single same-domain Tet4 pilot. Failed seeds and MMG stages retained. Direct CAD/Gmsh is a separate actual route."))
    assert all(w.digest(Path(p)) == digest for p, digest in hashes.items()), "Inputs changed during the run"
    assert all(w.digest(Path(p)) == digest for p, digest in fixture_hashes.items()), "Frozen fixture changed"


if __name__ == "__main__":
    main()
