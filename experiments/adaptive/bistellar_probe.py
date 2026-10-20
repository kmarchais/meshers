"""Bounded 2-to-3 cavity alternatives next to a checked Gyroid obstruction."""

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import shutil
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

import microgen_workflow as workflow
from job_parts import worker_part
from mmg_quality_localization import quality

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

SOURCE = workflow.HERE / "runs/g3-gyroid-firstpass-fem43-interior-relocation-v1"
WITNESS = workflow.HERE / "runs/g3-gyroid-convex-star-witness-audit-v1.json"
FACES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))


def oriented_boundary(cells):
    counts = defaultdict(int)
    for cell in cells:
        for face in FACES:
            vertices = list(map(int, cell[list(face)]))
            inversions = sum(vertices[first] > vertices[second] for first in range(3) for second in range(first + 1, 3))
            counts[tuple(sorted(vertices))] += -1 if inversions % 2 else 1
    if any(abs(count) > 1 for count in counts.values()):
        return None
    return {face: count for face, count in counts.items() if count}


def flip_candidate(points, cells, indices):
    original = cells[list(indices)]
    common = sorted(set(original[0]) & set(original[1]))
    assert len(common) == 3 and len(np.unique(original)) == 5
    opposite = [int(next(iter(set(cell) - set(common)))) for cell in original]
    row = {"old_elements": list(map(int, indices)), "shared_face": list(map(int, common)),
        "proposed_edge": opposite, "old_minimum_quality": float(quality(points, original).min()),
        "local_target_passed": False, "status": "new-edge-already-present"}
    if ((cells == opposite[0]).any(axis=1) & (cells == opposite[1]).any(axis=1)).any():
        return None, row
    alternative = np.array([[*opposite, common[first], common[second]] for first, second in ((0, 1), (1, 2), (2, 0))])
    determinants = np.linalg.det(points[alternative[:, 1:]] - points[alternative[:, :1]])
    if np.any(determinants == 0):
        row["status"] = "degenerate-alternative"
        return None, row
    negative = determinants < 0
    alternative[negative, 1], alternative[negative, 2] = alternative[negative, 2].copy(), alternative[negative, 1].copy()
    original_boundary = oriented_boundary(original)
    assert original_boundary is not None, "Source cavity has incoherent oriented incidence; do not apply flips"
    if original_boundary != oriented_boundary(alternative):
        row["status"] = "oriented-cavity-boundary-mismatch"
        return None, row
    old_volume = float(np.linalg.det(points[original[:, 1:]] - points[original[:, :1]]).sum() / 6)
    new_volume = float(abs(determinants).sum() / 6)
    assert old_volume > 0
    relative = abs(new_volume - old_volume) / old_volume
    assert relative < 1e-10
    measured = quality(points, alternative)
    row.update(status="local-target-passed" if measured.min() >= .1 else "local-quality-rejected",
        local_target_passed=bool(measured.min() >= .1), alternative_cells=alternative.tolist(),
        new_minimum_quality=float(measured.min()), minimum_signed_determinant=float(abs(determinants).min()),
        relative_cavity_volume_difference=relative, oriented_cavity_boundary_identical=True)
    return alternative, row


def worker(run):
    witness = json.loads(WITNESS.read_text())
    assert witness["passed"] and witness["checked_fixed_neighbor_obstruction"]
    assert not witness["other_interior_vertices_in_active_constraints"]
    recorded = json.loads(Path(witness["source"]).joinpath("result.json").read_text())
    assert workflow.digest(Path(witness["source"]) / "result.json") == witness["source_result_sha256"]
    assert Path(recorded["source"]).resolve() == SOURCE.resolve()
    mesh = np.load(SOURCE / "mesh.npz")
    points, cells = mesh["points"], mesh["tetrahedra"]
    assert len(cells) + 1 <= workflow.CAP
    active = [record["element"] for record in witness["active_constraints"]]
    targets = {tuple(sorted(map(int, cells[index, list(face)]))) for index in active for face in FACES}
    adjacency = {face: [] for face in targets}
    for index, cell in enumerate(cells):
        for face in FACES:
            key = tuple(sorted(map(int, cell[list(face)])))
            if key in adjacency:
                adjacency[key].append(index)
    candidates, visited = [], set()
    for indices in adjacency.values():
        assert len(indices) in (1, 2)
        if len(indices) != 2 or tuple(indices) in visited:
            continue
        visited.add(tuple(indices))
        _, report = flip_candidate(points, cells, indices)
        candidates.append(report)
    assert len(candidates) <= 12
    workflow.save(run / "result.json", {"case": "gyroid", "source": str(SOURCE), "witness": str(WITNESS),
        "active_elements": active, "candidates": candidates, "candidates_evaluated": len(candidates),
        "local_target_candidates": sum(record["local_target_passed"] for record in candidates),
        "source_mesh_unchanged": True, "source_boundary_unchanged": True,
        "qualified": False, "stiffness_solve_performed": False,
        "scope": "Single five-vertex 2-to-3 alternatives next to checked obstruction; no flips applied, no whole-mesh admission or continuum claim"})


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
    paths = [Path(__file__), SOURCE / "result.json", SOURCE / "mesh.npz", WITNESS]
    paths.extend(workflow.HERE / name for name in ("test_bistellar_probe.py", "job_parts.py", "test_job_parts.py",
        "microgen_workflow.py", "mmg_quality_localization.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "quality_gate": .1,
        "maximum_cavities": 12, "cavity_boundary_orientation_preserved": True, "source_unmodified": True})
    with Job(parts=["controls", "two-to-three"]) as job:
        for stage in ("controls", "two-to-three"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_bistellar_probe.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker"])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
