"""One explicit actual control: Gmsh size .03 and declared-cell validation.

The installed CAD, Gmsh meshing and MMG algorithms are unchanged. Reject the
control if BoxMesh's required triangles differ from the physical cell caps.
"""

import argparse
from contextlib import contextmanager
import importlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import sys
import tempfile
import time
import traceback
from unittest.mock import patch

import cad_union_comparison as base
import microgen_workflow as w


CONTROL = dict(name="actual-cad-size003-declared-cell-validation", gmsh_size=.03,
    declared_bounds=[[-.5] * 3, [.5] * 3], tolerance=1e-10,
    modifications=["Gmsh size .06 to .03", "Replace periodic predicates with full declared-cell node/facet/class validation"],
    unchanged="CAD generation, Gmsh periodic constraints, BoxMesh construction and required triangles, both public MMG passes, MMG arguments/retries and output snapping",
    required_triangle_condition="BoxMesh-selected required triangles must equal declared-cell material-cap triangles exactly",
    size_rationale="Halving chord length predicts about one quarter of the coarse .004369 reverse distance and 1.678% volume error; admission must be measured")


def diagnose(source, run):
    """Replay the saved public rejection and record its per-axis reason."""
    import numpy as np
    import gmsh
    from periodic_trace import periodic_pairs

    public = importlib.import_module("microgen.mesh")
    path = source / "actual/actual-cad-gmsh-direct/mesh.npz"
    with np.load(path) as data:
        points, surface = data["points"], data["surface"]
    _, _, physical = periodic_pairs(points, surface, declared_bounds=CONTROL["declared_bounds"])
    lo, hi = public._get_bounding_box(points)
    faces = public._sort_adjacent_faces_3d(public._extract_face_nodes(points, lo, hi, 1e-8, 3), points, 1e-8)
    axes = []
    for axis, label in enumerate("xyz"):
        low, high = faces[label + "-"], faces[label + "+"]
        row = dict(axis=axis, inferred_planes=[float(lo[axis]), float(hi[axis])],
                   low_nodes=len(low), high_nodes=len(high))
        if len(low) == len(high):
            other = [i for i in range(3) if i != axis]
            difference = points[high][:, other] - points[low][:, other]
            row.update(max_signed_transverse_difference=float(np.max(difference, initial=0)),
                       max_absolute_transverse_difference=float(np.max(abs(difference), initial=0)),
                       public_rejects=bool(np.any(difference > 1e-8)))
        else:
            row["public_rejects"] = True
        axes.append(row)
    original_rejected = False
    try:
        public._check_output_mesh_periodicity(str(source / "actual/gmsh-seed.vtk"), 1e-8)
    except public.OutputMeshNotPeriodicError:
        original_rejected = True
    finally:
        if gmsh.isInitialized():
            gmsh.finalize()
    result = dict(original_public_call_rejected=original_rejected,
        public_coordinate_check=bool(public.is_periodic(points)), axes=axes,
        declared_cell=physical, source_mesh_sha256=w.digest(path),
        hypotheses=["Mesh extrema substitute for physical empty Z faces", "Physical X/Y caps mismatch", "Public sorting rejects otherwise paired X/Y nodes"],
        recorded_failure="OutputMeshNotPeriodicError from the unmodified public post-Gmsh check")
    w.save(run / "diagnosis.json", result)
    assert original_rejected and not result["public_coordinate_check"]
    assert [row["public_rejects"] for row in axes] == [False, False, True], "Diagnosis differs from the declared-cell hypothesis"
    assert physical["periodic_pair_counts"][2] == 0 and physical["periodic_triangles_match"]


@contextmanager
def declared_validation(run):
    import numpy as np
    import pyvista as pv
    from periodic_trace import periodic_pairs

    mesh_module = importlib.import_module("microgen.mesh")
    remesh_module = importlib.import_module("microgen.remesh")
    original_remesh = remesh_module.remesh_keeping_boundaries_for_fem
    original_pipeline = remesh_module._run_mmg_pipeline
    original_required = remesh_module._generate_mesh_with_required_triangles
    records = dict(control=CONTROL, validators=[], required_triangles=[])

    def save():
        w.save(run / "validation.json", records)

    def check(grid, label):
        entry = dict(phase=label, passed=False)
        records["validators"].append(entry)
        file = run / f"validation-{len(records['validators'])}-{label}.vtu"
        grid.save(file)
        entry.update(mesh_file=str(file.resolve()), mesh_sha256=w.digest(file))
        try:
            if max(grid.n_points, grid.n_cells) > w.CAP:
                raise ValueError("Validation-stage mesh exceeds node/tetrahedron cap")
            points, _, surface = base.extract(grid)
            _, _, entry["periodic"] = periodic_pairs(points, surface,
                tolerance=CONTROL["tolerance"], declared_bounds=CONTROL["declared_bounds"])
            entry["passed"] = True
        except ValueError as error:
            entry["error"] = str(error)
        save()
        return entry["passed"]

    def gmsh_check(output_mesh_file, tol=1e-8):
        if not check(pv.read(output_mesh_file), "gmsh-output"):
            raise mesh_module.OutputMeshNotPeriodicError("Declared-cell full-cap validation failed; see validation.json")

    def required_guard(input_mesh, output_mesh="merged_reqtri.mesh"):
        # Inspection only. The public writer runs with the original BoxMesh,
        # inferred RVE and selected triangles if this equality holds.
        selected_mesh, selected_ids = input_mesh.boundary_elements(input_mesh.rve)
        surface = input_mesh.surface
        faces = surface.faces.reshape(-1, 4)
        assert np.all(faces[:, 0] == 3)
        corners = surface.points[faces[:, 1:]]
        physical = np.any(np.all(abs(corners + .5) <= CONTROL["tolerance"], axis=1) |
                          np.all(abs(corners - .5) <= CONTROL["tolerance"], axis=1), axis=1)
        expected = np.flatnonzero(physical)
        selected = np.asarray(selected_ids, dtype=int)
        entry = dict(boxmesh_bounds=[input_mesh.rve.min_point.tolist(), input_mesh.rve.max_point.tolist()],
            declared_bounds=CONTROL["declared_bounds"], selected_triangles=len(selected),
            declared_cap_triangles=len(expected), duplicate_selections=len(selected) - len(np.unique(selected)),
            extra_selected_ids=np.setdiff1d(selected, expected).tolist(),
            missing_cap_ids=np.setdiff1d(expected, selected).tolist(),
            equal=bool(np.array_equal(np.sort(selected), expected)))
        file = run / "boxmesh-required-triangles.vtp"
        selected_mesh.save(file)
        entry.update(mesh_file=str(file.resolve()), mesh_sha256=w.digest(file))
        records["required_triangles"].append(entry)
        save()
        if not entry["equal"]:
            raise ValueError("Validator-only override is insufficient: BoxMesh required triangles differ from physical caps")
        return original_required(input_mesh, output_mesh)

    def remesh_with_validation(input_mesh, **kwargs):
        # The public predicate receives only coordinates. Associate them with
        # its actual input/output grid so the replacement checks full facets.
        context = dict(grid=input_mesh if isinstance(input_mesh, pv.UnstructuredGrid) else input_mesh.to_pyvista(),
                       phase="mmg-input")

        def predicate(coords, tol=1e-8):
            assert np.array_equal(coords, context["grid"].points), "Validator lost its matching volume mesh"
            return check(context["grid"], context["phase"])

        def pipeline(*args, **options):
            grid = original_pipeline(*args, **options)
            context.update(grid=grid, phase="mmg-output")
            # The original caller still applies _snap_to_expected_bounds to
            # this same object before invoking the replacement predicate.
            return grid

        with patch.object(remesh_module, "is_periodic", predicate), \
             patch.object(remesh_module, "_run_mmg_pipeline", pipeline), \
             patch.object(remesh_module, "_generate_mesh_with_required_triangles", required_guard):
            return original_remesh(input_mesh, **kwargs)

    save()
    with patch.object(mesh_module, "_check_output_mesh_periodicity", gmsh_check), \
         patch.object(remesh_module, "remesh_keeping_boundaries_for_fem", remesh_with_validation):
        yield


def worker(run, source):
    w.imports()
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.GetBackend() == "Sequential" and vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    (run / "temporary").mkdir()
    tempfile.tempdir = str((run / "temporary").resolve())
    row = dict(arm="actual", control=CONTROL, status="started", qualified=False,
        phase_seconds={}, mmg_calls=[], mmg_instrumentation_s=0., vtk_backend="Sequential",
        numerical_threads=1, solve_performed=False)
    w.save(run / "result.json", row)
    started = time.perf_counter()
    try:
        diagnosis = time.perf_counter()
        diagnose(source, run)
        row["diagnosis_s"] = time.perf_counter() - diagnosis
        with patch.dict(base.GMSH, size=CONTROL["gmsh_size"]), declared_validation(run):
            base.actual_worker(run, row)
        row["status"] = "admitted" if row["qualified"] else "admission-rejected"
    except Exception as error:
        row.update(status="failed", qualified=False, error_type=type(error).__name__,
                   error=str(error), traceback=traceback.format_exc())
    for key, folder in (("direct", "actual-cad-gmsh-direct"), ("final", "actual-cad-gmsh-public-mmg")):
        if key in row:
            row[key]["control"] = CONTROL
            row[key]["route"] = ("Explicit actual CAD/Gmsh size .03 with declared-cell validation" +
                (" plus unchanged public MMG passes, arguments and retries" if key == "final" else ", without MMG"))
            w.save(run / folder / "result.json", row[key])
    row["wall_s"] = time.perf_counter() - started
    w.save(run / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True, help="Completed original paired run")
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        return worker(args.run, args.source)
    from job_parts import worker_part
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    assert importlib.metadata.version("cadquery-ocp-novtk") == base.CAD_PIN
    original = json.loads((args.source / "manifest.json").read_text())
    # Reuse the original frozen samples. The actual control does not execute
    # or claim a fresh native result from whichever binary is installed now.
    hashes = {p: digest for p, digest in original["inputs"].items() if not p.endswith("_meshers.abi3.so")}
    hashes.update(original["fixture_hashes"])
    for path, digest in hashes.items():
        assert w.digest(Path(path)) == digest, "Original source or fixture changed: " + path
    for path in (Path(__file__), args.source / "manifest.json", args.source / "summary.json",
                 args.source / "actual/gmsh-seed.vtk", args.source / "actual/actual-cad-gmsh-direct/mesh.npz"):
        hashes[str(path.resolve())] = w.digest(path)
    args.run.mkdir(parents=True, exist_ok=False)
    for index, path in enumerate(map(Path, hashes)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run / f"input-{index}-{path.name}")
    directory = args.run / "actual"
    directory.mkdir()
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(args.source / "fixture" / (base.NAME + suffix), directory / (base.NAME + suffix))
    w.save(args.run / "manifest.json", dict(inputs=hashes, control=CONTROL,
        source=str(args.source.resolve()), original_native_binary_sha256=original["inputs"].get(
            str(w.ROOT / "crates/meshers-python/python/meshers/_meshers.abi3.so")),
        cad_version_pin=base.CAD_PIN, packages={n: importlib.metadata.version(n) for n in original["packages"]},
        actual_gmsh={**base.GMSH, "size": CONTROL["gmsh_size"]}, actual_mmg=base.MMG,
        platform=platform.platform(), numerical_thread_policy=w.NUMERICAL_THREAD_POLICY,
        memory_bytes=w.MEMORY, worker_timeout_s=180, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        timing_scope="One fresh actual control. CAD construction, volume check, STEP export, Gmsh, MMG, output extraction, reorientation, declared-cell validation and required-triangle equality checks are charged. Imports, saved-seed diagnosis, archive, independent admission and measured MMG instrumentation are excluded.",
        mechanical_accuracy_claim=False, successful_performance_ratio=None, command=sys.argv))
    with Job(parts=["actual-control"]) as job:
        with worker_part(job, "actual-control") as part:
            record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--worker", "--run",
                str(directory.resolve()), "--source", str(args.source.resolve())], directory / "worker.log", part, timeout=180)
            w.save(directory / "worker.json", record)
            row = (json.loads((directory / "result.json").read_text()) if (directory / "result.json").exists()
                   else dict(status="no-result", qualified=False))
            row.update(worker=record, directory=str(directory.resolve()))
            if record["returncode"] != 0:
                row.update(status="worker-failed", qualified=False)
            job.result("actual-control", qualified=row["qualified"], generation_s=row.get("generation_s"),
                       wall_s=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    w.save(args.run / "summary.json", dict(rows=[row], control=CONTROL, complete=True,
        original_results_unchanged=True, successful_performance_ratio=None, mechanical_accuracy_claim=False))
    assert all(w.digest(Path(p)) == digest for p, digest in hashes.items()), "Inputs changed during control"


if __name__ == "__main__":
    main()
