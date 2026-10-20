"""Paired fixed-resolution Microgen extraction-order diagnostic, no MMG."""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import time

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

import microgen_workflow as workflow
from cap_quality_bound import bounds, caps
from facet_incidence_audit import audit
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def sampled_inputs(grid):
    arrays = {"points": np.asarray(grid.points), **{name: np.asarray(grid.point_data[name]) for name in grid.point_data}}
    return {name: {"shape": list(values.shape), "dtype": str(values.dtype),
        "sha256": hashlib.sha256(values.tobytes(order="C")).hexdigest()} for name, values in arrays.items()}


def triangulate_background(grid, cap=workflow.CAP, tetra_per_cell=6):
    if tetra_per_cell not in (5, 6):
        raise ValueError("Only five- or six-tet backgrounds are supported")
    if tetra_per_cell * grid.n_cells > cap:
        raise ValueError("Background tetrahedron upper bound exceeds staged cap")
    points = np.asarray(grid.points).copy()
    values = {name: np.asarray(grid.point_data[name]).copy() for name in grid.point_data}
    if tetra_per_cell == 6:
        converted = grid.cast_to_unstructured_grid()
        converted.points = points
        converted = converted.triangulate()
    else:
        import pyvista as pv

        rectilinear = pv.RectilinearGrid(*(np.unique(points[:, axis]) for axis in range(3)))
        np.testing.assert_array_equal(rectilinear.points, points)
        converted = rectilinear.to_tetrahedra(tetra_per_cell=5, pass_data=False, pass_cell_ids=False)
        np.testing.assert_array_equal(converted.points, points.astype(converted.points.dtype))
        converted.points = points
        for name, expected in values.items():
            converted.point_data[name] = expected
    assert converted.n_cells <= cap and set(map(int, converted.celltypes)) == {10}
    np.testing.assert_array_equal(converted.points, points)
    for name, expected in values.items():
        np.testing.assert_array_equal(converted.point_data[name], expected)
    return converted


def worker(run, mode, tetra_per_cell):
    _, _, _, microgen, pv, _ = workflow.imports()
    from microgen.shape import surface_functions
    from periodic import periodic_pairs

    started = time.monotonic()
    shape = microgen.Tpms(surface_functions.gyroid, offset=1.1, resolution=31)
    sampled = sampled_inputs(shape.grid)
    background_cells = shape.grid.n_cells
    assert background_cells <= workflow.CAP and shape._grid_sheet is None
    if mode == "pretriangulated":
        shape.grid = triangulate_background(shape.grid, tetra_per_cell=tetra_per_cell)
        assert sampled_inputs(shape.grid) == sampled
    background_tets = shape.grid.n_cells if mode == "pretriangulated" else None
    grid = shape.generate_volume_mesh(type_part="sheet").triangulate()
    row = {"case": "gyroid", "resolution": 31, "offset": 1.1, "mode": mode,
        "microgen_generate_volume_mesh_called": True, "mmg_calls": [], "qualified": False,
        "background_cells_original": background_cells, "preclipped_background_tetrahedra": background_tets,
        "background_tetrahedra_per_cell": tetra_per_cell if mode == "pretriangulated" else None,
        "sampled_inputs": sampled, "sample_values_and_coordinates_unchanged": True,
        "nodes": grid.n_points, "elements": grid.n_cells, "generation_s": time.monotonic() - started,
        "route": "microgen-vtk-unmodified-extraction-at31" if mode == "unmodified" else "microgen-vtk-explicit-preclipping-background-triangulation-at31",
        "default_microgen_pipeline": mode == "unmodified", "default_resolution": False,
        "status": "seed-element-cap", "stiffness_solve_performed": False}
    output = run / mode
    if grid.n_cells <= workflow.CAP:
        assert set(grid.cells_dict) == {int(pv.CellType.TETRA)}
        points = np.asarray(grid.points).copy()
        original = grid.cells_dict[pv.CellType.TETRA].copy()
        determinants = np.linalg.det(points[original[:, 1:]] - points[original[:, :1]])
        cells, reordered = workflow.reorient_linear(points, original)
        poly = grid.extract_surface(pass_pointid=True, algorithm="dataset_surface").triangulate()
        surface = np.asarray(poly.point_data["vtkOriginalPointIds"])[poly.faces.reshape(-1, 4)[:, 1:]]
        np.savez_compressed(output / "mesh.npz", points=points, tetrahedra=cells, surface=surface)
        _, triangles = caps({"points": points, "surface": surface})
        ceiling, _ = bounds(triangles)
        finite = np.isfinite(ceiling)
        row.update(reoriented_tets=reordered, zero_volume_tets=int((determinants == 0).sum()),
            cap_minimum_quality_ceiling=float(ceiling[finite].min()) if finite.any() else None,
            cap_faces_below_quality_gate=int((ceiling < .1).sum()), nonfinite_cap_ceilings=int((~finite).sum()),
            facet_consistency=audit(points, cells, surface))
        try:
            _, _, row["periodic"] = periodic_pairs(points, surface)
        except ValueError as error:
            row["periodic_error"] = str(error)
        if row["zero_volume_tets"]:
            row["status"] = "degenerate-seed"
        else:
            row["geometry"] = workflow.geometry(output, "gyroid", points, cells, surface)
            row["status"] = "seed-audited-no-mmg"
        row["necessary_fixed_cap_feasibility"] = bool(row["facet_consistency"]["oriented_facet_consistency_passed"]
            and not row["cap_faces_below_quality_gate"] and not row["nonfinite_cap_ceilings"] and "periodic" in row
            and not row["zero_volume_tets"])
    row.update(total_diagnostic_s=time.monotonic() - started,
        scope="Fixed analytic Gyroid and identical sampled field; one extraction-order variant, not default resolution or a Meshers seed; seed checks do not establish an admitted MMG or matched-accuracy workflow")
    workflow.save(output / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--worker", choices=("unmodified", "pretriangulated"))
    parser.add_argument("--background-tets", type=int, choices=(5, 6), default=6)
    args = parser.parse_args()
    if args.worker:
        worker(args.run, args.worker, args.background_tets)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    source = workflow.HERE / "runs/g3-microgen-defaults-v4"
    paths = [Path(__file__)]
    paths.extend(workflow.HERE / name for name in ("test_preclip_probe.py", "microgen_workflow.py",
        "facet_incidence_audit.py", "cap_quality_bound.py", "job_parts.py", "test_job_parts.py"))
    import microgen.shape.tpms
    paths.append(Path(microgen.shape.tpms.__file__))
    for mode in ("unmodified", "pretriangulated"):
        output = args.run / mode
        output.mkdir()
        for suffix in ("-samples.npz", "-specimen.json"):
            path = source / ("gyroid" + suffix)
            paths.append(path)
            shutil.copyfile(path, output / path.name)
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "resolution": 31, "offset": 1.1,
        "background_tet_cap": workflow.CAP, "material_tet_cap": workflow.CAP,
        "background_tetrahedra_per_cell": args.background_tets,
        "no_mmg_or_fea": True, "one_extraction_order_variant": True,
        "packages": {name: importlib.metadata.version(name) for name in ("microgen", "vtk", "pyvista", "numpy")}})
    rows = {}
    with Job(parts=["controls", "unmodified", "pretriangulated"]) as job:
        for stage in ("controls", "unmodified", "pretriangulated"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_preclip_probe.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker", stage,
                        "--background-tets", str(args.background_tets)])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                values = {"seconds": record["wall_s"], "peak_gib": record["peak_rss_bytes"] / 1024**3}
                if stage != "controls":
                    rows[stage] = json.loads((args.run / stage / "result.json").read_text())
                    values.update(elements=rows[stage]["elements"], blocking_caps=rows[stage].get("cap_faces_below_quality_gate"))
                job.result(stage, **values)
    assert rows["unmodified"]["sampled_inputs"] == rows["pretriangulated"]["sampled_inputs"]
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())
    workflow.save(args.run / "summary.json", {"inputs_unchanged": True, "sampled_field_identical": True, "rows": rows})


if __name__ == "__main__":
    main()
