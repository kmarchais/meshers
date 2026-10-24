"""Planar cap-edge alternatives only; no surface or volume mesh is modified."""

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
from cap_quality_bound import bounds, caps
from facet_incidence_audit import audit
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

SOURCE = workflow.HERE / "runs/g3-gyroid-preclip31-five-v2/pretriangulated"


def triangle_boundary(triangles):
    edges = defaultdict(int)
    for triangle in triangles:
        for first, second in ((0, 1), (1, 2), (2, 0)):
            start, end = map(int, triangle[[first, second]])
            edges[tuple(sorted((start, end)))] += 1 if start < end else -1
    return {edge: count for edge, count in edges.items() if count}


def flip_pair(points, triangles, indices):
    original = triangles[list(indices)].copy()
    shared = sorted(map(int, set(original[0]) & set(original[1])))
    assert len(shared) == 2 and len(np.unique(original)) == 4
    opposite = [int(next(iter(set(triangle) - set(shared)))) for triangle in original]
    row = {"triangles": list(map(int, indices)), "shared_edge": shared,
        "alternative_edge": opposite, "local_target_passed": False}
    axis = next((axis for axis in range(3) if any(
        np.all(abs(points[original, axis] - value) < 1e-10) for value in (-.5, .5))), None)
    if axis is None:
        return None, dict(row, status="different-cap-planes")
    if ((triangles == opposite[0]).any(axis=1) & (triangles == opposite[1]).any(axis=1)).any():
        return None, dict(row, status="alternative-edge-already-present")
    transverse = [coordinate for coordinate in range(3) if coordinate != axis]

    def signed_area(cells):
        coordinates = points[cells][:, :, transverse]
        return np.linalg.det(coordinates[:, 1:] - coordinates[:, :1]) / 2

    areas = signed_area(original)
    if np.any(areas == 0):
        return None, dict(row, status="degenerate-source")
    negative = areas < 0
    original[negative, 1], original[negative, 2] = original[negative, 2].copy(), original[negative, 1].copy()
    alternative = np.array([[*opposite, vertex] for vertex in shared])
    new_areas = signed_area(alternative)
    if np.any(new_areas == 0):
        return None, dict(row, status="degenerate-alternative")
    negative = new_areas < 0
    alternative[negative, 1], alternative[negative, 2] = alternative[negative, 2].copy(), alternative[negative, 1].copy()
    if triangle_boundary(original) != triangle_boundary(alternative):
        return None, dict(row, status="nonconvex-or-boundary-mismatch")
    difference = abs(abs(areas).sum() - abs(new_areas).sum()) / abs(areas).sum()
    assert difference < 1e-10
    ceiling, _ = bounds(points[alternative])
    passed = bool(ceiling.min() >= .1)
    row.update(status="cap-ceiling-passed" if passed else "cap-ceiling-rejected",
        axis=axis, side=float(points[original[0, 0], axis]), alternative_triangles=alternative.tolist(),
        old_minimum_ceiling=float(bounds(points[original])[0].min()),
        new_minimum_ceiling=float(ceiling.min()), local_target_passed=passed,
        relative_area_difference=float(difference), oriented_patch_boundary_identical=True)
    return alternative, row


def worker(run, source=SOURCE, maximum_blocked_faces=100):
    workflow.imports()
    from periodic import periodic_pairs

    source_row = json.loads((source / "result.json").read_text())
    mesh = np.load(source / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert len(cells) <= workflow.CAP
    assert audit(points, cells, surface)["oriented_facet_consistency_passed"]
    _, _, periodic = periodic_pairs(points, surface)
    cap_faces, triangles = caps(mesh)
    ceiling, _ = bounds(triangles)
    blocked = np.flatnonzero(ceiling < .1)
    assert len(blocked) <= maximum_blocked_faces, "Bound the local candidate count before enumeration"
    adjacency = defaultdict(list)
    for index, face in enumerate(cap_faces):
        for first, second in ((0, 1), (1, 2), (2, 0)):
            adjacency[tuple(sorted(map(int, face[[first, second]])))].append(index)
    candidates, visited, outer_edges = [], set(), []
    for index in blocked:
        face = cap_faces[index]
        for first, second in ((0, 1), (1, 2), (2, 0)):
            edge = tuple(sorted(map(int, face[[first, second]])))
            neighbors = adjacency[edge]
            assert len(neighbors) <= 2
            if len(neighbors) == 1:
                outer_edges.append(list(edge))
                continue
            key = tuple(sorted(neighbors))
            if key in visited:
                continue
            visited.add(key)
            _, row = flip_pair(points, cap_faces, key)
            candidates.append(row)
    feasible = [row for row in candidates if row["local_target_passed"]]
    covered = {index for row in feasible for index in row["triangles"]}
    workflow.save(run / "result.json", {"source": str(source), "case": source_row["case"],
        "blocked_cap_faces": list(map(int, blocked)), "blocked_count": len(blocked),
        "candidate_count": len(candidates), "locally_feasible_count": len(feasible),
        "blocked_faces_with_local_alternative": len(covered & set(blocked)),
        "blocked_faces_without_local_alternative": sorted(map(int, set(blocked) - covered)),
        "outer_edges": outer_edges, "candidates": candidates, "periodic_source": periodic,
        "source_mesh_unchanged": True, "qualified": False, "mmg_calls": [],
        "stiffness_solve_performed": False,
        "scope": "Individual convex planar edge-flip cap-ceiling feasibility only; alternatives are not applied, jointly compatible, periodic replacements, or valid volume cavities"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--maximum-blocked-faces", type=int, choices=(100, 256), default=100)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    args.run, args.source = args.run.resolve(), args.source.resolve()
    if args.worker:
        worker(args.run, args.source, args.maximum_blocked_faces)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    shutil.copyfile(workflow.HERE / "test_cap_flip_probe.py", args.run / "test_cap_flip_probe.py")
    paths = [Path(__file__), args.source / "mesh.npz", args.source / "result.json"]
    paths.extend(workflow.HERE / name for name in ("test_cap_flip_probe.py", "cap_quality_bound.py",
        "facet_incidence_audit.py", "microgen_workflow.py", "job_parts.py", "test_job_parts.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY,
        "platform": sys.platform,
        "threads": 1, "worker_timeout_s": 120, "tet_cap": workflow.CAP,
        "maximum_blocked_faces": args.maximum_blocked_faces, "no_mmg_or_fea": True})
    with Job(parts=["controls", "cap-flips"]) as job:
        for stage in ("controls", "cap-flips"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_cap_flip_probe.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker",
                        "--source", str(args.source), "--maximum-blocked-faces", str(args.maximum_blocked_faces)])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
