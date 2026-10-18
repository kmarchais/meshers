"""One bounded tuned TPMS Microgen/MMG candidate per invocation."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import time

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


class IntermediateCapError(ValueError):
    pass


def worker(args):
    np, fd, _, microgen, pv, _ = workflow.imports()
    from microgen.shape import surface_functions
    from microgen.remesh import remesh_keeping_boundaries_for_fem
    from periodic import periodic_pairs, tensor_difference
    from homogenization import AuditedCG
    from cap_quality_bound import bounds, caps
    from unittest.mock import patch

    output = args.run / "result.json"
    row = {"case": args.case, "route": "microgen-vtk-mmg-tuned",
           "resolution": args.resolution,
           "mmg_arguments": {"hausd": args.hausd, "hmax": 0.06, "hgrad": 1.3},
           "status": "started", "mmg_calls": []}
    row["stage_audit_enabled"] = args.audit_stages
    row["timing_comparable"] = not args.audit_stages
    row["first_pass_hausd_override"] = args.first_pass_hausd
    row["first_pass_hmax_override"] = args.first_pass_hmax if args.first_pass_hausd is not None else None
    row["first_pass_only"] = args.first_pass_only
    row["optimize_second_pass"] = args.optimize_second_pass
    row["post_optimize"] = args.post_optimize
    row["post_optimization_method"] = args.post_optimization_method if args.post_optimize else None
    row["solver_formulation"] = "prescribed-strain" if args.strain_controlled else "stress-controlled"
    if args.first_pass_hausd is not None:
        row["route"] = "microgen-vtk-mmg-first-pass-override"
    if args.first_pass_only:
        row["route"] = "microgen-vtk-single-mmg-diagnostic-override"
    if args.optimize_second_pass:
        row["route"] = "microgen-vtk-mmg-interior-optimization-override"
    if args.post_optimize:
        row["route"] = "microgen-vtk-mmg-third-pass-optimization-override"
    def audit_stage(label, stage_grid):
        tetrahedra_count = int(np.count_nonzero(stage_grid.celltypes == int(pv.CellType.TETRA)))
        if tetrahedra_count > workflow.CAP:
            raise IntermediateCapError(f"{label}: {tetrahedra_count} tetrahedra exceed the intermediate element cap")
        stage_grid = stage_grid.extract_cells(stage_grid.celltypes == int(pv.CellType.TETRA))
        stage_points = np.asarray(stage_grid.points).copy()
        stage_tets, _ = workflow.reorient_linear(
            stage_points, stage_grid.cells_dict[pv.CellType.TETRA].copy())
        stage_poly = stage_grid.extract_surface(pass_pointid=True, algorithm="dataset_surface").triangulate()
        stage_surface = np.asarray(stage_poly.point_data["vtkOriginalPointIds"])[
            stage_poly.faces.reshape(-1, 4)[:, 1:]]
        np.savez_compressed(args.run / (label + ".npz"), points=stage_points,
                            tetrahedra=stage_tets, surface=stage_surface)
        metrics = workflow.geometry(args.run, args.case, stage_points, stage_tets, stage_surface)
        workflow.save(args.run / (label + ".json"), {
            "nodes": len(stage_points), "elements": len(stage_tets), "geometry": metrics})
    started = time.perf_counter()
    function, offset, part = ((surface_functions.gyroid, 1.1, "sheet") if args.case == "gyroid"
        else (surface_functions.schwarz_p, 0, "lower skeletal"))
    shape = microgen.Tpms(function, offset=offset, resolution=args.resolution)
    row["specimen_setup_s"] = time.perf_counter() - started
    grid = shape.generate_volume_mesh(type_part=part).triangulate()
    row.update(seed_nodes=grid.n_points, seed_elements=grid.n_cells,
               vtk_extraction_s=time.perf_counter() - started)
    workflow.save(output, row)
    if grid.n_cells > workflow.CAP:
        row.update(status="seed-element-cap", total_workflow_s=time.perf_counter() - started)
        workflow.save(output, row)
        return
    if args.audit_stages:
        audit_stage("vtk-seed", grid)
    binary = Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))
    os.environ["PATH"] = str(binary.parent) + os.pathsep + os.environ["PATH"]
    original = microgen.Mmg.mmg3d
    def record_call(**parameters):
        if args.optimize_second_pass and parameters.get("ls"):
            row["replaced_second_pass"] = {key: value for key, value in parameters.items()
                if key not in {"input", "output"}}
            metric = Path(parameters["input"]).with_suffix(".sol")
            if metric.exists():
                row["optimization_input_metric_sha256"] = workflow.digest(metric)
                metric.unlink()
                row["optimization_input_metric_removed"] = True
            parameters = {"input": parameters["input"], "output": parameters["output"],
                "optim": True, "nosurf": True, "noinsert": True}
        if args.first_pass_only and parameters.get("ls"):
            row["omitted_second_pass"] = {key: value for key, value in parameters.items()
                if key not in {"input", "output"}}
            workflow.save(output, row)
            shutil.copyfile(parameters["input"], parameters["output"])
            return
        if args.first_pass_hausd is not None and parameters.get("nofem"):
            parameters.update(hausd=args.first_pass_hausd, hmax=args.first_pass_hmax, hgrad=1.3)
        row["mmg_calls"].append({key: value for key, value in parameters.items()
                                 if key not in {"input", "output"}})
        workflow.save(output, row)
        try:
            result = original(**parameters)
            if args.audit_stages:
                label = f"mmg-stage-{len(row['mmg_calls'])}"
                captured = args.run / (label + ".mesh")
                shutil.copyfile(parameters["output"], captured)
                from microgen.remesh import _remove_unnecessary_fields_from_mesh_file
                cleaned = args.run / (label + "-tetra.mesh")
                _remove_unnecessary_fields_from_mesh_file(str(captured), str(cleaned), 2, 3)
                audit_stage(label, pv.UnstructuredGrid(str(cleaned)))
            return result
        except Exception as error:
            captured = getattr(error.__cause__, "output", b"")
            (args.run / f"mmg-call-{len(row['mmg_calls'])}.log").write_bytes(captured or b"")
            raise
    try:
        with patch.object(microgen.Mmg, "mmg3d", record_call):
            grid = remesh_keeping_boundaries_for_fem(grid, **row["mmg_arguments"])
            if args.post_optimize:
                from microgen.remesh import (_generate_mesh_with_required_triangles,
                    _remove_unnecessary_fields_from_mesh_file)
                optimization_input = args.run / "postopt-input.mesh"
                optimization_output = args.run / "postopt-output.mesh"
                optimization_tetrahedra = args.run / "postopt-tetra.mesh"
                _generate_mesh_with_required_triangles(microgen.BoxMesh.from_pyvista(grid), str(optimization_input))
                microgen.Mmg.mmg3d(input=str(optimization_input), output=str(optimization_output),
                    **{args.post_optimization_method: True}, nosurf=True, noinsert=True)
                _remove_unnecessary_fields_from_mesh_file(str(optimization_output), str(optimization_tetrahedra), 2, 3)
                grid = pv.UnstructuredGrid(str(optimization_tetrahedra))
    except Exception as error:
        status = "intermediate-element-cap" if isinstance(error, IntermediateCapError) else "microgen-mmg-failed"
        row.update(status=status, error=str(error),
                   total_workflow_s=time.perf_counter() - started)
        workflow.save(output, row)
        return
    row["generation_s"] = time.perf_counter() - started
    row.update(nodes=grid.n_points, elements=grid.n_cells)
    if grid.n_cells > workflow.CAP or grid.n_points > 100_000:
        row.update(status="output-size-cap", total_workflow_s=time.perf_counter() - started)
        workflow.save(output, row)
        return
    assert set(grid.cells_dict) == {int(pv.CellType.TETRA)}
    points = np.asarray(grid.points).copy()
    tetrahedra, changed = workflow.reorient_linear(points, grid.cells_dict[pv.CellType.TETRA].copy())
    poly = grid.extract_surface(pass_pointid=True, algorithm="dataset_surface").triangulate()
    surface = np.asarray(poly.point_data["vtkOriginalPointIds"])[poly.faces.reshape(-1, 4)[:, 1:]]
    np.savez_compressed(args.run / "mesh.npz", points=points, tetrahedra=tetrahedra, surface=surface)
    row["reoriented_tetrahedra"] = changed
    _, triangles = caps({"points": points, "surface": surface})
    ceiling, _ = bounds(triangles)
    row["cap_minimum_quality_ceiling"] = float(ceiling.min())
    row["caps_below_quality_gate"] = int((ceiling < 0.1).sum())
    audit_started = time.perf_counter()
    row["geometry"] = workflow.geometry(args.run, args.case, points, tetrahedra, surface)
    try:
        _, _, row["periodic"] = periodic_pairs(points, surface)
    except ValueError as error:
        row["periodic_error"] = str(error)
    row["audit_s"] = time.perf_counter() - audit_started
    row["qualified"] = row["geometry"]["geometry_admitted"] and "periodic" in row
    row["status"] = "geometry-rejected"
    if row["qualified"]:
        solve_started = time.perf_counter()
        if args.strain_controlled:
            from stress_recovery import solve
            fields, stiffness, records = solve(points, tetrahedra)
        else:
            fd.ModelingSpace("3D")
            fem = fd.Mesh(points.copy(), tetrahedra.copy(), "tet4")
            assembly = fd.Assembly.create(fd.weakform.StressEquilibrium(
                fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)), fem)
            solver = AuditedCG()
            stiffness = np.asarray(fd.homogen.get_homogenized_stiffness(
                assembly, meshperio=True, solver=solver, rigid_body_constraint="pin"))
            records = solver.records
        solve_seconds = time.perf_counter() - solve_started
        assert len(records) == 6
        if args.strain_controlled:
            np.savez_compressed(args.run / "fields.npz", displacement=fields)
        reference = json.loads((args.source / "summary.json").read_text())["references"][args.case]
        assert reference["qualified_inputs"]
        row.update(stiffness=stiffness.tolist(), solver_records=records,
                   solve_s=solve_seconds, status="solved",
                   reference_disagreement=tensor_difference(stiffness, np.array(reference["stiffness"])))
    row["total_workflow_s"] = time.perf_counter() - started
    workflow.save(output, row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--case", choices=("gyroid", "schwarz-solid"), default="schwarz-solid")
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--source", type=Path, default=workflow.HERE / "runs/g3-microgen-defaults-v4")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--audit-stages", action="store_true")
    parser.add_argument("--hausd", type=float, default=0.0005)
    parser.add_argument("--first-pass-hausd", type=float)
    parser.add_argument("--first-pass-hmax", type=float, default=0.06)
    parser.add_argument("--first-pass-only", action="store_true")
    parser.add_argument("--optimize-second-pass", action="store_true")
    parser.add_argument("--post-optimize", action="store_true")
    parser.add_argument("--post-optimization-method", choices=("optim", "optimLES"), default="optim")
    parser.add_argument("--strain-controlled", action="store_true")
    args = parser.parse_args()
    assert args.hausd > 0
    assert args.first_pass_hmax > 0
    assert args.first_pass_hmax == 0.06 or args.first_pass_hausd is not None
    assert not (args.first_pass_only and args.optimize_second_pass)
    assert not (args.post_optimize and (args.first_pass_only or args.optimize_second_pass))
    assert args.post_optimize or args.post_optimization_method == "optim"
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    inputs = [args.source / name for name in
              (args.case + "-specimen.json", args.case + "-samples.npz", "summary.json")]
    hashes = {str(path): workflow.digest(path) for path in inputs}
    for path in inputs[:2]:
        shutil.copyfile(path, args.run / path.name)
    workflow.save(args.run / "manifest.json", {"script_sha256": workflow.digest(Path(__file__)),
        "inputs": hashes, "limits": {"rss_bytes": workflow.MEMORY, "timeout_s": 600,
        "tetrahedra": workflow.CAP, "nodes": 100_000, "threads": 1},
        "packages": {name: importlib.metadata.version(name) for name in
                     ("microgen", "mmgpy", "fedoo", "numpy", "scipy", "vtk", "pyvista")},
        "baseline_manifest_sha256": workflow.digest(args.source / "manifest.json")})
    with Job(parts=["tuned-candidate"]) as job:
        with job.part("tuned-candidate") as part:
            command = [sys.executable, __file__, "--run", str(args.run),
                "--source", str(args.source), "--case", args.case,
                "--resolution", str(args.resolution), "--hausd", str(args.hausd), "--worker"]
            if args.audit_stages:
                command.append("--audit-stages")
            if args.first_pass_hausd is not None:
                command.extend(["--first-pass-hausd", str(args.first_pass_hausd)])
                command.extend(["--first-pass-hmax", str(args.first_pass_hmax)])
            if args.first_pass_only:
                command.append("--first-pass-only")
            if args.optimize_second_pass:
                command.append("--optimize-second-pass")
            if args.post_optimize:
                command.extend(["--post-optimize", "--post-optimization-method", args.post_optimization_method])
            if args.strain_controlled:
                command.append("--strain-controlled")
            result = workflow.bounded(command,
                args.run / "worker.log", part, timeout=600)
            workflow.save(args.run / "worker.json", result)
            assert result["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("tuned-candidate", qualified=int(row.get("qualified", False)),
                       nodes=row.get("nodes", 0), seconds=row["total_workflow_s"],
                       quality=row.get("geometry", {}).get("q_min_corner"))
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
