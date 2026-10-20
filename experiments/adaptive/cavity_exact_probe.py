"""Enumerate only the two unresolved six-vertex cap-flip cavities."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

import microgen_workflow as workflow
from cap_flip_probe import SOURCE
from cap_quality_bound import caps
from cavity_triangulation import enumerate_cavity, replacement_boundary
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

CAVITIES = workflow.HERE / "runs/g3-gyroid-cap-cavity-preflight-v2/result.json"
CAPS = workflow.HERE / "runs/g3-gyroid-cap-flip-preflight-v2/result.json"


def worker(run):
    mesh = np.load(SOURCE / "mesh.npz")
    points, cells = mesh["points"], mesh["tetrahedra"]
    cap_faces, _ = caps(mesh)
    unresolved = [row for row in json.loads(CAVITIES.read_text())["candidates"] if not row["cavity_passed"]]
    assert len(unresolved) == 2
    candidates = json.loads(CAPS.read_text())["candidates"]
    rows = []
    for row in unresolved:
        candidate = next(candidate for candidate in candidates if candidate["shared_edge"] == row["edge"])
        original = cells[row["incident_elements"]]
        target = replacement_boundary(points, original, cap_faces[row["cap_triangles"]], np.asarray(candidate["alternative_triangles"]))
        alternative, result = enumerate_cavity(points, original, target)
        result.update(edge=row["edge"], incident_elements=row["incident_elements"],
            alternative_cells=None if alternative is None else alternative.tolist())
        rows.append(result)
    workflow.save(run / "result.json", {"source": str(SOURCE), "cavities": rows,
        "solved_cavities": sum(row["passed"] for row in rows), "source_unchanged": True,
        "qualified": False, "mmg_calls": [], "stiffness_solve_performed": False,
        "scope": "Exhaustive floating-point signed-face matching for these six existing vertices only; not an obstruction to larger cavities, Steiner vertices, or general remeshing"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), CAVITIES, CAPS, SOURCE / "mesh.npz"]
    paths.extend(workflow.HERE / name for name in ("cavity_triangulation.py", "test_cavity_triangulation.py",
        "bistellar_probe.py", "cap_flip_probe.py", "cap_quality_bound.py", "microgen_workflow.py", "job_parts.py", "test_job_parts.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "maximum_vertices": 6, "maximum_subsets": 32767})
    with Job(parts=["controls", "enumeration"]) as job:
        for stage in ("controls", "enumeration"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_cavity_triangulation.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker"])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
