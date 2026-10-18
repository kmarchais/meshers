"""Recheck saved candidate geometry with strictly absolute box-plane tolerance."""

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


def source_files(source):
    return (source.with_suffix(".json"), source) if source.is_file() else (source / "result.json", source / "mesh.npz")


def worker(args):
    np, *_ = workflow.imports()
    record_path, mesh_path = source_files(args.source)
    row = json.loads(record_path.read_text())
    mesh = np.load(mesh_path)
    points, surface = mesh["points"], mesh["surface"]
    probes = points[surface[:, :3]].reshape(-1, 3)
    loose = np.isclose(abs(probes), 0.5, atol=1e-12)
    strict = workflow.on_box_plane(abs(probes), 0.5)
    assert np.isclose(0.5 - 1e-6, 0.5, atol=1e-12)
    assert not workflow.on_box_plane(0.5 - 1e-6, 0.5)
    assert workflow.on_box_plane(0.5 + 1e-13, 0.5)
    metrics = workflow.geometry(mesh_path.parent, row["case"], points, mesh["tetrahedra"], surface)
    workflow.save(args.run / (args.source.stem + ".json"), {"source": str(args.source),
        "case": row["case"], "previous_geometry": row["geometry"], "geometry": metrics,
        "near_plane_vertex_coordinate_locks_removed": int((loose & ~strict).sum())})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--legacy-references", action="store_true")
    args = parser.parse_args()
    if args.source:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    sources = [workflow.HERE / "runs" / name for name in (
        "g3-spinodoid-linear32-solve-v1", "g3-spinodoid-curved24-v1",
        "g3-lattice-void-clamp32-v1", "g3-lattice-curved16-v1")]
    if args.legacy_references:
        sources = [workflow.HERE / "runs/g3-microgen-defaults-v4" / f"{case}-meshers-curved-{cells}.npz"
                   for case in ("gyroid", "schwarz-solid") for cells in (16, 20)]
    inputs = {}
    for source in sources:
        record_path, mesh_path = source_files(source)
        case = json.loads(record_path.read_text())["case"]
        for path in (record_path, mesh_path, mesh_path.parent / (case + "-samples.npz"),
                     mesh_path.parent / (case + "-specimen.json")):
            inputs[str(path)] = workflow.digest(path)
    workflow.save(args.run / "manifest.json", {"script_sha256": workflow.digest(Path(__file__)),
        "geometry_source_sha256": workflow.digest(Path(workflow.__file__)), "inputs": inputs,
        "plane_absolute_tolerance": 1e-12, "plane_relative_tolerance": 0})
    with Job(parts=[source.stem for source in sources]) as job:
        for source in sources:
            with job.part(source.stem) as part:
                result = workflow.bounded([sys.executable, __file__, "--run", str(args.run),
                    "--source", str(source)], args.run / (source.name + ".log"), part, timeout=180)
                workflow.save(args.run / (source.name + "-worker.json"), result)
                assert result["returncode"] == 0, f"Read {source.name}.log"
                row = json.loads((args.run / (source.stem + ".json")).read_text())
                job.result(source.stem, admitted=int(row["geometry"]["geometry_admitted"]),
                    distance=row["geometry"]["sampled_two_sided_max"], seconds=result["wall_s"])
    assert all(workflow.digest(Path(path)) == value for path, value in inputs.items())


if __name__ == "__main__":
    main()
