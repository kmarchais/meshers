"""Audit oriented tet-face incidence without changing frozen qualification flags."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

import microgen_workflow as workflow
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

FACES = np.array(((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)))
RUNS = workflow.HERE / "runs"
MESHES = {
    "gyroid-vtk43": RUNS / "g3-gyroid-first-pass-fem43-v1/vtk-seed.npz",
    "gyroid-mmg-stage1": RUNS / "g3-gyroid-first-pass-fem43-v1/mmg-stage-1.npz",
    "gyroid-fem-final": RUNS / "g3-gyroid-first-pass-fem43-v1/mesh.npz",
    "gyroid-partial-repair": RUNS / "g3-gyroid-firstpass-fem43-interior-relocation-v1/mesh.npz",
    "gyroid-two-pass": RUNS / "g3-gyroid-tuned43-v3/mesh.npz",
    "schwarz-microgen-linear": RUNS / "g3-schwarz-matched-v1/repeat-1-microgen/mesh.npz",
    "schwarz-microgen-tet10-parent": RUNS / "g3-schwarz-microgen-hmax01-seed-v1/mesh.npz",
    "schwarz-native": RUNS / "g3-schwarz-tet10-live-amg-repeated-v1/repeat-1/meshers-tet10/mesh.npz",
    "gyroid-native": RUNS / "g3-gyroid-curved9-v1/mesh.npz",
    "spinodoid-native": RUNS / "g3-spinodoid-partial-settle-curved24-tol003-solve-v1/mesh.npz",
    "lattice-native": RUNS / "g3-lattice-adapt24-tol0028-v1/mesh.npz",
    "lattice-microgen": RUNS / "g3-lattice-tuned51-tet10-seed-hausd001-v1/mesh.npz",
}


def audit(points, cells, surface):
    cells = cells[:, :4]
    surface = surface[:, :3]
    assert cells.min() >= 0 and cells.max() < len(points)
    assert not (np.diff(np.sort(cells, axis=1), axis=1) == 0).any()
    oriented = cells[:, FACES].reshape(-1, 3)
    inversions = (oriented[:, 0] > oriented[:, 1]).astype(int)
    inversions += oriented[:, 0] > oriented[:, 2]
    inversions += oriented[:, 1] > oriented[:, 2]
    signs = np.where(inversions % 2, -1, 1)
    sorted_faces = np.sort(oriented, axis=1)
    ordering = np.lexsort(sorted_faces[:, ::-1].T)
    ordered = sorted_faces[ordering]
    starts = np.r_[0, np.flatnonzero((ordered[1:] != ordered[:-1]).any(axis=1)) + 1]
    counts = np.diff(np.r_[starts, len(ordered)])
    sums = np.add.reduceat(signs[ordering], starts)
    unique = ordered[starts]
    conflicts = np.flatnonzero((counts == 2) & (sums != 0))
    extracted = {tuple(map(int, face)) for face in unique[counts == 1]}
    provided = {tuple(sorted(map(int, face))) for face in surface}
    examples = []
    for group in conflicts[:10]:
        rows = ordering[starts[group]:starts[group] + counts[group]]
        indices = rows // 4
        face = unique[group]
        triangle = points[face]
        normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
        opposite = [int(next(iter(set(cells[index]) - set(face)))) for index in indices]
        distances = [float(np.dot(points[node] - triangle[0], normal)) for node in opposite]
        examples.append({"face": face.tolist(), "incident_elements": indices.tolist(),
            "orientation_sum": int(sums[group]), "opposite_vertices": opposite,
            "signed_plane_products": distances, "same_physical_side": bool(distances[0] * distances[1] > 0)})
    determinants = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])
    duplicate_cells = len(cells) - len(np.unique(np.sort(cells, axis=1), axis=0))
    duplicate_surface = len(surface) - len(provided)
    matched = extracted == provided
    coherent = not len(conflicts) and not (counts > 2).any() and not duplicate_cells and not duplicate_surface and matched
    return {"nodes": len(points), "elements": len(cells), "unique_faces": len(unique),
        "interior_faces": int((counts == 2).sum()), "extracted_boundary_faces": len(extracted),
        "same_orientation_shared_faces": len(conflicts), "nonmanifold_faces": int((counts > 2).sum()),
        "duplicate_cells": duplicate_cells, "duplicate_surface_faces": duplicate_surface,
        "provided_boundary_matches_extraction": matched, "missing_provided_faces": len(extracted - provided),
        "extra_provided_faces": len(provided - extracted), "nonpositive_corner_determinants": int((determinants <= 0).sum()),
        "minimum_corner_determinant": float(determinants.min()), "oriented_facet_consistency_passed": bool(coherent),
        "conflict_examples": examples,
        "scope": "Oriented corner-face incidence and same-side neighbor diagnostics; necessary consistency, not a certificate of curved geometry or absence of nonlocal intersections"}


def worker(run, name):
    path = MESHES[name]
    archive = np.load(path)
    assert len(archive["tetrahedra"]) <= workflow.CAP
    result = audit(archive["points"], archive["tetrahedra"], archive["surface"])
    result.update(source=str(path), source_sha256=workflow.digest(path))
    result_path = path.parent / "result.json"
    if path.name == "mesh.npz" and result_path.exists():
        previous = json.loads(result_path.read_text())
        result.update(previous_qualified_flag=previous.get("qualified"), previous_case=previous.get("case"))
    workflow.save(run / (name + ".json"), result)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", choices=tuple(MESHES))
    args = parser.parse_args()
    if args.worker:
        worker(args.run, args.worker)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), *MESHES.values()]
    paths.extend(workflow.HERE / name for name in ("test_facet_incidence_audit.py", "job_parts.py", "test_job_parts.py", "microgen_workflow.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "meshes": {name: str(path) for name, path in MESHES.items()},
        "memory_bytes": workflow.MEMORY, "threads": 1, "worker_timeout_s": 120,
        "frozen_source_flags_unchanged": True, "no_meshing_or_solves": True})
    rows = {}
    with Job(parts=["controls", *MESHES]) as job:
        for stage in ("controls", *MESHES):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_facet_incidence_audit.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker", stage])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                values = {"seconds": record["wall_s"], "peak_gib": record["peak_rss_bytes"] / 1024**3}
                if stage != "controls":
                    rows[stage] = json.loads((args.run / (stage + ".json")).read_text())
                    values.update(conflicts=rows[stage]["same_orientation_shared_faces"], coherent=int(rows[stage]["oriented_facet_consistency_passed"]))
                job.result(stage, **values)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())
    workflow.save(args.run / "summary.json", {"inputs_unchanged": True, "meshes": rows,
        "scope": "Supplemental validity audit; no frozen flags or numerical outputs reclassified in place"})


if __name__ == "__main__":
    main()
