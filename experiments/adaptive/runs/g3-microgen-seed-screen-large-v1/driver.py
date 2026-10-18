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

    grid = microgen.Tpms(surface_functions.schwarz_p, offset=0,
                        resolution=args.resolution).generate_volume_mesh(type_part="lower skeletal").triangulate()
    row = {"resolution": args.resolution, "nodes": grid.n_points, "elements": grid.n_cells}
    if grid.n_cells > workflow.CAP:
        row["status"] = "seed-element-cap"
    else:
        points = np.asarray(grid.points).copy()
        tetrahedra, _ = workflow.reorient_linear(points, grid.cells_dict[pv.CellType.TETRA].copy())
        poly = grid.extract_surface(pass_pointid=True, algorithm="dataset_surface").triangulate()
        surface = np.asarray(poly.point_data["vtkOriginalPointIds"])[poly.faces.reshape(-1, 4)[:, 1:]]
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
        else:
            row["geometry"] = workflow.geometry(args.run, "schwarz-solid", points, tetrahedra, surface)
            row["status"] = "screened"
    workflow.save(args.run / f"seed-{args.resolution}.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--resolution", type=int)
    parser.add_argument("--resolutions", type=int, nargs="+", default=(21, 22, 23, 25, 26, 27))
    args = parser.parse_args()
    if args.resolution is not None:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    source = workflow.HERE / "runs/g3-microgen-defaults-v4"
    inputs = [source / ("schwarz-solid" + suffix) for suffix in ("-samples.npz", "-specimen.json")]
    hashes = {str(path): workflow.digest(path) for path in inputs}
    for path in inputs:
        shutil.copyfile(path, args.run / path.name)
    resolutions = args.resolutions
    assert all(resolution >= 2 for resolution in resolutions) and len(set(resolutions)) == len(resolutions)
    workflow.save(args.run / "manifest.json", {"script_sha256": workflow.digest(Path(__file__)),
        "inputs": hashes, "resolutions": resolutions, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "tetrahedra_cap": workflow.CAP})
    rows = []
    with Job(parts=[str(resolution) for resolution in resolutions]) as job:
        for resolution in resolutions:
            with job.part(str(resolution)) as part:
                result = workflow.bounded([sys.executable, __file__, "--run", str(args.run),
                    "--resolution", str(resolution)], args.run / f"seed-{resolution}.log", part, timeout=120)
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
