"""Screen actual VTK seeds for retained-cap feasibility before MMG tuning."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    np, _, _, microgen, pv, _ = workflow.imports()
    from cap_quality_bound import bounds, caps
    from microgen.shape import surface_functions

    field_verification = None
    if args.case == "spinodoid":
        case, shape, _ = workflow.specimen(args.case)
        frep = shape._frep
        coefficients = [(*map(float, mode / frep.cell_size), float(amplitude), float(phase))
            for mode, amplitude, phase in zip(frep.modes, frep.amplitudes, frep.phases, strict=True)]
        frozen = json.loads((args.run / "fourier-coefficients.json").read_text())
        np.testing.assert_array_equal(coefficients, frozen["coefficients"])
        assert frep.threshold == frozen["threshold"] and shape.resolution == 16
        queries = np.random.default_rng(20261018).uniform(-0.5, 0.5, (4096, 3))
        before = case.function(*queries.T)
        shape.resolution = args.resolution
        grid = shape.generate_volume_mesh().triangulate()
        grid.translate((-0.5,) * 3, inplace=True)
        np.testing.assert_array_equal(before, case.function(*queries.T))
        after = [(*map(float, mode / frep.cell_size), float(amplitude), float(phase))
            for mode, amplitude, phase in zip(frep.modes, frep.amplitudes, frep.phases, strict=True)]
        np.testing.assert_array_equal(after, frozen["coefficients"])
        field_verification = {"fft_resolution": 16, "sampling_resolution_cells": args.resolution,
            "fourier_modes": len(coefficients), "unchanged_evaluation_points": len(queries),
            "field_values_bitwise_unchanged": True, "coefficients_unchanged": True}
    elif args.case == "smooth-lattice":
        case, shape, _ = workflow.specimen(args.case)
        queries = np.random.default_rng(20261018).uniform(-0.5, 0.5, (4096, 3))
        before = shape.evaluate(*queries.T)
        expected = case.function(*queries.T)
        np.testing.assert_allclose(before, expected, rtol=1e-12, atol=1e-12)
        grid = shape.generate_volume_mesh(resolution=args.resolution).triangulate()
        np.testing.assert_array_equal(before, shape.evaluate(*queries.T))
        field_verification = {"sampling_resolution_points": args.resolution,
            "unchanged_evaluation_points": len(queries), "field_values_bitwise_unchanged": True,
            "shape_case_maximum_value_difference": float(abs(before - expected).max())}
    else:
        function, offset, part = ((surface_functions.gyroid, 1.1, "sheet") if args.case == "gyroid"
            else (surface_functions.schwarz_p, 0, "lower skeletal"))
        grid = microgen.Tpms(function, offset=offset,
                            resolution=args.resolution).generate_volume_mesh(type_part=part).triangulate()
    row = {"case": args.case, "resolution": args.resolution, "nodes": grid.n_points, "elements": grid.n_cells}
    if field_verification:
        row["field_verification"] = field_verification
    if grid.n_cells > workflow.CAP:
        row["status"] = "seed-element-cap"
    else:
        points = np.asarray(grid.points).copy()
        tetrahedra, _ = workflow.reorient_linear(points, grid.cells_dict[pv.CellType.TETRA].copy())
        poly = grid.extract_surface(pass_pointid=True, algorithm="dataset_surface").triangulate()
        surface = np.asarray(poly.point_data["vtkOriginalPointIds"])[poly.faces.reshape(-1, 4)[:, 1:]]
        default_resolution = {"spinodoid": 16, "smooth-lattice": 50}.get(args.case)
        if args.resolution == default_resolution:
            previous = np.load(args.run / (args.case + "-vtk-default-0.npz"))
            for name, actual in (("points", points), ("tetrahedra", tetrahedra), ("surface", surface)):
                np.testing.assert_array_equal(actual, previous[name])
            row["default_seed_arrays_identical"] = True
        _, triangles = caps({"points": points, "surface": surface})
        ceiling, _ = bounds(triangles)
        valid_ceiling = np.isfinite(ceiling)
        row.update(cap_minimum_quality_ceiling=float(ceiling[valid_ceiling].min()) if valid_ceiling.any() else None,
                   caps_below_quality_gate=int((ceiling < 0.1).sum()),
                   degenerate_caps=int((~valid_ceiling).sum()))
        determinants = np.linalg.det(points[tetrahedra[:, 1:]] - points[tetrahedra[:, :1]])
        row["zero_volume_tetrahedra"] = int((determinants == 0).sum())
        if row["zero_volume_tetrahedra"]:
            row["status"] = "degenerate-seed"
        elif args.cap_only:
            row["status"] = "cap-screened-no-geometry"
        else:
            row["geometry"] = workflow.geometry(args.run, args.case, points, tetrahedra, surface)
            row["status"] = "screened"
    workflow.save(args.run / f"seed-{args.resolution}.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--case", choices=("gyroid", "schwarz-solid", "spinodoid", "smooth-lattice"),
                        default="schwarz-solid")
    parser.add_argument("--cap-only", action="store_true")
    parser.add_argument("--resolution", type=int)
    parser.add_argument("--resolutions", type=int, nargs="+", default=(21, 22, 23, 25, 26, 27))
    args = parser.parse_args()
    if args.resolution is not None:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    source = workflow.HERE / "runs/g3-microgen-defaults-v4"
    inputs = [source / (args.case + suffix) for suffix in ("-samples.npz", "-specimen.json")]
    if args.case == "spinodoid":
        coefficients = workflow.HERE / "runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json"
        inputs.append(coefficients)
    default_resolution = {"spinodoid": 16, "smooth-lattice": 50}.get(args.case)
    if default_resolution in args.resolutions:
        inputs.append(source / (args.case + "-vtk-default-0.npz"))
    hashes = {str(path): workflow.digest(path) for path in inputs}
    for path in inputs:
        shutil.copyfile(path, args.run / path.name)
    resolutions = args.resolutions
    assert all(resolution >= 2 for resolution in resolutions) and len(set(resolutions)) == len(resolutions)
    workflow.save(args.run / "manifest.json", {"script_sha256": workflow.digest(Path(__file__)),
        "inputs": hashes, "case": args.case, "cap_only": args.cap_only,
        "resolutions": resolutions, "memory_bytes": workflow.MEMORY,
        "workflow_source_sha256": workflow.digest(Path(workflow.__file__)),
        "baseline_manifest_sha256": workflow.digest(source / "manifest.json"),
        "threads": 1, "worker_timeout_s": 120, "tetrahedra_cap": workflow.CAP})
    rows = []
    with Job(parts=[str(resolution) for resolution in resolutions]) as job:
        for resolution in resolutions:
            with job.part(str(resolution)) as part:
                command = [sys.executable, __file__, "--run", str(args.run),
                    "--case", args.case, "--resolution", str(resolution)]
                if args.cap_only:
                    command.append("--cap-only")
                result = workflow.bounded(command,
                    args.run / f"seed-{resolution}.log", part, timeout=120)
                workflow.save(args.run / f"worker-{resolution}.json", result)
                assert result["returncode"] == 0, f"Read seed-{resolution}.log"
                row = json.loads((args.run / f"seed-{resolution}.json").read_text())
                rows.append(row)
                job.result(str(resolution), nodes=row["nodes"],
                    blocking_caps=row.get("caps_below_quality_gate"), seconds=result["wall_s"])
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())
    workflow.save(args.run / "summary.json", {"inputs_unchanged": True, "rows": rows})


if __name__ == "__main__":
    main()
