"""One neighboring tet at a time around two measured unresolved stars."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

import microgen_workflow as workflow
from bistellar_probe import oriented_boundary
from cap_flip_probe import SOURCE
from cap_quality_bound import caps
from cavity_milp import solve_cavity
from cavity_triangulation import replacement_boundary
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

CAVITIES = workflow.HERE / "runs/g3-gyroid-cap-cavity-preflight-v2/result.json"
CAPS = workflow.HERE / "runs/g3-gyroid-cap-flip-preflight-v2/result.json"


def worker(run, full_ring):
    mesh = np.load(SOURCE / "mesh.npz")
    points, cells = mesh["points"], mesh["tetrahedra"]
    assert len(cells) <= workflow.CAP
    cap_faces, _ = caps(mesh)
    previous = json.loads(CAVITIES.read_text())["candidates"]
    unresolved = [row for row in previous if not row["cavity_passed"]]
    assert len(unresolved) == 2
    used = {index for row in previous if row["cavity_passed"] for index in row["incident_elements"]}
    candidates = json.loads(CAPS.read_text())["candidates"]
    rows = []
    for row in unresolved:
        candidate = next(candidate for candidate in candidates if candidate["shared_edge"] == row["edge"])
        original_indices = row["incident_elements"]
        boundary = oriented_boundary(cells[original_indices])
        assert boundary is not None
        neighbors = set()
        for face in boundary:
            selected = np.ones(len(cells), dtype=bool)
            for vertex in face:
                selected &= (cells == vertex).any(axis=1)
            neighbors.update(map(int, np.flatnonzero(selected)))
        neighbors -= set(original_indices)
        assert len(neighbors) <= 8
        report = {"edge": row["edge"], "neighbor_candidates": sorted(neighbors), "attempts": [], "passed": False}
        expansions = [sorted(neighbors)] if full_ring else [[neighbor] for neighbor in sorted(neighbors)]
        for expansion in expansions:
            indices = sorted(set(original_indices) | set(expansion))
            if set(indices) & used:
                report["attempts"].append({"neighbors": expansion, "status": "conflicts-with-prior-cavity"})
                continue
            original = cells[indices]
            target = replacement_boundary(points, original, cap_faces[row["cap_triangles"]], np.asarray(candidate["alternative_triangles"]))
            alternative, outcome = solve_cavity(points, original, target, maximum_vertices=9 if full_ring else 7)
            outcome["neighbors"] = expansion
            report["attempts"].append(outcome)
            if alternative is not None:
                report.update(passed=True, incident_elements=indices, alternative_cells=alternative.tolist(),
                    new_cap_faces=candidate["alternative_triangles"], old_cap_faces=cap_faces[row["cap_triangles"]].tolist())
                used.update(indices)
                break
        rows.append(report)
    workflow.save(run / "result.json", {"source": str(SOURCE), "cavities": rows,
        "solved_cavities": sum(row["passed"] for row in rows), "source_unchanged": True,
        "qualified": False, "replacement_applied": False, "mmg_calls": [], "stiffness_solve_performed": False,
        "full_ring": full_ring,
        "scope": "Measured face-neighbor expansions, at most nine existing vertices, independent integer replay; no mesh replacement or workflow accuracy claim"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--full-ring", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.run, args.full_ring)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), CAVITIES, CAPS, SOURCE / "mesh.npz"]
    paths.extend(workflow.HERE / name for name in ("cavity_milp.py", "test_cavity_milp.py", "cavity_triangulation.py",
        "test_cavity_triangulation.py", "bistellar_probe.py", "cap_flip_probe.py", "cap_quality_bound.py",
        "microgen_workflow.py", "job_parts.py", "test_job_parts.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "maximum_vertices": 9 if args.full_ring else 7,
        "full_ring": args.full_ring,
        "maximum_neighbor_expansions_per_star": 8, "binary_solve_time_limit_s": 2,
        "binary_solve_node_limit": 1000, "unchanged_source": True})
    with Job(parts=["controls", "expansions"]) as job:
        for stage in ("controls", "expansions"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_cavity_milp.py", "experiments/adaptive/test_cavity_triangulation.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker"])
                if stage != "controls" and args.full_ring:
                    arguments.append("--full-ring")
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
