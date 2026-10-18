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


def worker(args):
    np, *_ = workflow.imports()
    row = json.loads((args.source / "result.json").read_text())
    mesh = np.load(args.source / "mesh.npz")
    points, surface = mesh["points"], mesh["surface"]
    probes = points[surface[:, :3]].reshape(-1, 3)
    loose = np.isclose(abs(probes), 0.5, atol=1e-12)
    strict = np.isclose(abs(probes), 0.5, atol=1e-12, rtol=0)
    assert np.isclose(0.5 - 1e-6, 0.5, atol=1e-12)
    assert not np.isclose(0.5 - 1e-6, 0.5, atol=1e-12, rtol=0)
    metrics = workflow.geometry(args.source, row["case"], points, mesh["tetrahedra"], surface)
    workflow.save(args.run / (args.source.name + ".json"), {"source": str(args.source),
        "case": row["case"], "previous_geometry": row["geometry"], "geometry": metrics,
        "near_plane_vertex_coordinate_locks_removed": int((loose & ~strict).sum())})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    args = parser.parse_args()
    if args.source:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    sources = [workflow.HERE / "runs" / name for name in (
        "g3-spinodoid-linear32-solve-v1", "g3-spinodoid-curved24-v1",
        "g3-lattice-void-clamp32-v1", "g3-lattice-curved16-v1")]
    inputs = {str(source / name): workflow.digest(source / name) for source in sources
              for name in ("result.json", "mesh.npz")}
    workflow.save(args.run / "manifest.json", {"script_sha256": workflow.digest(Path(__file__)),
        "geometry_source_sha256": workflow.digest(Path(workflow.__file__)), "inputs": inputs,
        "plane_absolute_tolerance": 1e-12, "plane_relative_tolerance": 0})
    with Job(parts=[source.name for source in sources]) as job:
        for source in sources:
            with job.part(source.name) as part:
                result = workflow.bounded([sys.executable, __file__, "--run", str(args.run),
                    "--source", str(source)], args.run / (source.name + ".log"), part, timeout=180)
                workflow.save(args.run / (source.name + "-worker.json"), result)
                assert result["returncode"] == 0, f"Read {source.name}.log"
                row = json.loads((args.run / (source.name + ".json")).read_text())
                job.result(source.name, admitted=int(row["geometry"]["geometry_admitted"]),
                    distance=row["geometry"]["sampled_two_sided_max"], seconds=result["wall_s"])
    assert all(workflow.digest(Path(path)) == value for path, value in inputs.items())


if __name__ == "__main__":
    main()
