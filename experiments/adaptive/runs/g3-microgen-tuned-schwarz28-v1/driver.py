"""One bounded tuned Schwarz-solid Microgen/MMG candidate per invocation."""

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


def worker(args):
    np, fd, _, microgen, pv, _ = workflow.imports()
    from microgen.shape import surface_functions
    from microgen.remesh import remesh_keeping_boundaries_for_fem
    from periodic import periodic_pairs, tensor_difference
    from homogenization import AuditedCG
    from cap_quality_bound import bounds, caps
    from unittest.mock import patch

    output = args.run / "result.json"
    row = {"case": "schwarz-solid", "route": "microgen-vtk-mmg-tuned",
           "resolution": args.resolution,
           "mmg_arguments": {"hausd": 0.0005, "hmax": 0.06, "hgrad": 1.3},
           "status": "started", "mmg_calls": []}
    started = time.perf_counter()
    shape = microgen.Tpms(surface_functions.schwarz_p, offset=0, resolution=args.resolution)
    row["specimen_setup_s"] = time.perf_counter() - started
    grid = shape.generate_volume_mesh(type_part="lower skeletal").triangulate()
    row.update(seed_nodes=grid.n_points, seed_elements=grid.n_cells,
               vtk_extraction_s=time.perf_counter() - started)
    workflow.save(output, row)
    if grid.n_cells > workflow.CAP:
        row.update(status="seed-element-cap", total_workflow_s=time.perf_counter() - started)
        workflow.save(output, row)
        return
    binary = Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))
    os.environ["PATH"] = str(binary.parent) + os.pathsep + os.environ["PATH"]
    original = microgen.Mmg.mmg3d
    def record_call(**parameters):
        row["mmg_calls"].append({key: value for key, value in parameters.items()
                                 if key not in {"input", "output"}})
        workflow.save(output, row)
        return original(**parameters)
    with patch.object(microgen.Mmg, "mmg3d", record_call):
        grid = remesh_keeping_boundaries_for_fem(grid, **row["mmg_arguments"])
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
    row["geometry"] = workflow.geometry(args.run, "schwarz-solid", points, tetrahedra, surface)
    try:
        _, _, row["periodic"] = periodic_pairs(points, surface)
    except ValueError as error:
        row["periodic_error"] = str(error)
    row["audit_s"] = time.perf_counter() - audit_started
    row["qualified"] = row["geometry"]["geometry_admitted"] and "periodic" in row
    row["status"] = "geometry-rejected"
    if row["qualified"]:
        solve_started = time.perf_counter()
        fd.ModelingSpace("3D")
        fem = fd.Mesh(points.copy(), tetrahedra.copy(), "tet4")
        assembly = fd.Assembly.create(fd.weakform.StressEquilibrium(
            fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)), fem)
        solver = AuditedCG()
        stiffness = np.asarray(fd.homogen.get_homogenized_stiffness(
            assembly, meshperio=True, solver=solver, rigid_body_constraint="pin"))
        assert len(solver.records) == 6
        reference = json.loads((args.source / "summary.json").read_text())["references"]["schwarz-solid"]
        assert reference["qualified_inputs"]
        row.update(stiffness=stiffness.tolist(), solver_records=solver.records,
                   solve_s=time.perf_counter() - solve_started, status="solved",
                   reference_disagreement=tensor_difference(stiffness, np.array(reference["stiffness"])))
    row["total_workflow_s"] = time.perf_counter() - started
    workflow.save(output, row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--source", type=Path, default=workflow.HERE / "runs/g3-microgen-defaults-v4")
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    inputs = [args.source / name for name in
              ("schwarz-solid-specimen.json", "schwarz-solid-samples.npz", "summary.json")]
    hashes = {str(path): workflow.digest(path) for path in inputs}
    for path in inputs[:2]:
        shutil.copyfile(path, args.run / path.name)
    workflow.save(args.run / "manifest.json", {"script_sha256": workflow.digest(Path(__file__)),
        "inputs": hashes, "limits": {"rss_bytes": workflow.MEMORY, "timeout_s": 600,
        "tetrahedra": workflow.CAP, "nodes": 100_000, "threads": 1},
        "baseline_manifest_sha256": workflow.digest(args.source / "manifest.json")})
    with Job(parts=["tuned-candidate"]) as job:
        with job.part("tuned-candidate") as part:
            result = workflow.bounded([sys.executable, __file__, "--run", str(args.run),
                "--source", str(args.source), "--resolution", str(args.resolution), "--worker"],
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
