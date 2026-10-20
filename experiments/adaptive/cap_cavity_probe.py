"""Check complete bounded boundary stars before applying cap-edge alternatives."""

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
from cap_quality_bound import bounds, caps
from facet_incidence_audit import audit
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

PREFLIGHT = workflow.HERE / "runs/g3-gyroid-cap-flip-preflight-v2/result.json"


def cavity(points, cells, old_faces, new_faces, edge, *, minimum_relative_determinant=0):
    indices = np.flatnonzero((cells == edge[0]).any(axis=1) & (cells == edge[1]).any(axis=1))
    row = {"edge": list(map(int, edge)), "incident_elements": indices.tolist(), "cavity_passed": False}
    if len(indices) not in (2, 3):
        return None, dict(row, status="not-two-or-three-tet-edge-star")
    original = cells[indices]
    if set(new_faces.ravel()) != set(old_faces.ravel()):
        return None, dict(row, status="cap-patch-vertex-mismatch")
    apex = set(original.ravel()) - set(old_faces.ravel())
    if len(apex) != len(indices) - 1:
        return None, dict(row, status="unexpected-edge-star-vertices")
    before = oriented_boundary(original)
    assert before is not None, "Source cavity has inconsistent orientation"
    old_keys = {tuple(sorted(map(int, face))) for face in old_faces}
    new_keys = {tuple(sorted(map(int, face))) for face in new_faces}
    if not old_keys <= set(before):
        return None, dict(row, status="cap-faces-not-cavity-boundary")
    old_volume = np.linalg.det(points[original[:, 1:]] - points[original[:, :1]]).sum()
    assert old_volume > 0
    outer = {key: value for key, value in before.items() if key not in old_keys}
    attempts = []
    for vertex in sorted(map(int, apex)):
        faces = [face for face in [*outer, *map(tuple, new_faces.tolist())] if vertex not in face]
        alternative = np.array([[*face, vertex] for face in faces])
        determinants = np.linalg.det(points[alternative[:, 1:]] - points[alternative[:, :1]])
        if np.any(abs(determinants) <= minimum_relative_determinant * old_volume):
            attempts.append({"apex": vertex, "status": "degenerate-alternative"})
            continue
        negative = determinants < 0
        alternative[negative, 1], alternative[negative, 2] = alternative[negative, 2].copy(), alternative[negative, 1].copy()
        after = oriented_boundary(alternative)
        if after is None or not new_keys <= set(after) or outer != {key: value for key, value in after.items() if key not in new_keys}:
            attempts.append({"apex": vertex, "status": "noncap-oriented-boundary-mismatch"})
            continue
        difference = float(abs(abs(determinants).sum() - old_volume) / old_volume)
        if difference >= 1e-10:
            attempts.append({"apex": vertex, "status": "cavity-volume-mismatch"})
            continue
        row.update(status="cavity-passed", cavity_passed=True, alternative_cells=alternative.tolist(),
            relative_volume_difference=difference, noncap_oriented_boundary_identical=True,
            minimum_corner_determinant=float(abs(determinants).min()), apex=vertex, prior_attempts=attempts)
        return alternative, row
    return None, dict(row, status="no-valid-boundary-cone", attempts=attempts)


def worker(run):
    workflow.imports()
    from periodic import periodic_pairs

    mesh = np.load(SOURCE / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert len(cells) <= workflow.CAP
    assert audit(points, cells, surface)["oriented_facet_consistency_passed"]
    cap_faces, _ = caps(mesh)
    candidates = [row for row in json.loads(PREFLIGHT.read_text())["candidates"] if row["local_target_passed"]]
    assert len(candidates) <= 100
    rows, plans, used_cells, used_faces = [], [], set(), set()
    for candidate in candidates:
        old_faces = cap_faces[candidate["triangles"]]
        new_faces = np.asarray(candidate["alternative_triangles"])
        alternative, row = cavity(points, cells, old_faces, new_faces, candidate["shared_edge"])
        row["cap_triangles"] = candidate["triangles"]
        if alternative is not None:
            cell_indices = set(row["incident_elements"])
            face_keys = {tuple(sorted(map(int, face))) for face in old_faces}
            row["disjoint_from_prior_cavities"] = not bool(cell_indices & used_cells or face_keys & used_faces)
            if row["disjoint_from_prior_cavities"]:
                plans.append((row["incident_elements"], alternative, face_keys, new_faces))
                used_cells |= cell_indices
                used_faces |= face_keys
        rows.append(row)
    result = {"source": str(SOURCE), "preflight": str(PREFLIGHT), "candidates": rows,
        "candidate_count": len(rows), "cavities_passed": sum(row["cavity_passed"] for row in rows),
        "compatible_cavities": len(plans), "qualified": False, "mmg_calls": [],
        "stiffness_solve_performed": False, "source_unchanged": True, "replacement_applied": False}
    if len(plans) == len(candidates):
        replaced = [cell for index, cell in enumerate(cells) if index not in used_cells]
        kept = [face for face in surface if tuple(sorted(map(int, face))) not in used_faces]
        for indices, alternative, _, new_faces in plans:
            replaced.extend(alternative)
            kept.extend(new_faces)
        replacement_surface = np.asarray(kept)
        replaced = np.asarray(replaced)
        assert len(replaced) <= workflow.CAP
        consistency = audit(points, replaced, replacement_surface)
        result["facet_consistency"] = consistency
        assert consistency["oriented_facet_consistency_passed"] and not consistency["nonpositive_corner_determinants"]
        _, _, result["periodic"] = periodic_pairs(points, replacement_surface)
        _, triangles = caps({"points": points, "surface": replacement_surface})
        ceiling, _ = bounds(triangles)
        result.update(cap_faces_below_quality_gate=int((ceiling < .1).sum()),
            cap_minimum_quality_ceiling=float(ceiling.min()), replacement_applied=True)
        for suffix in ("-samples.npz", "-specimen.json"):
            shutil.copyfile(SOURCE / ("gyroid" + suffix), run / ("gyroid" + suffix))
        result["geometry"] = workflow.geometry(run, "gyroid", points, replaced, replacement_surface)
        np.savez_compressed(run / "mesh.npz", points=points, tetrahedra=replaced, surface=replacement_surface)
    result["scope"] = "Only complete two- or three-tet boundary-star replacements; no MMG or FE and no admitted volume quality claim"
    workflow.save(run / "result.json", result)


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
    shutil.copyfile(workflow.HERE / "test_cap_cavity_probe.py", args.run / "test_cap_cavity_probe.py")
    paths = [Path(__file__), PREFLIGHT, SOURCE / "mesh.npz", SOURCE / "result.json"]
    paths.extend(SOURCE / ("gyroid" + suffix) for suffix in ("-samples.npz", "-specimen.json"))
    paths.extend(workflow.HERE / name for name in ("test_cap_cavity_probe.py", "cap_flip_probe.py", "bistellar_probe.py",
        "cap_quality_bound.py", "facet_incidence_audit.py", "microgen_workflow.py", "job_parts.py", "test_job_parts.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "tet_cap": workflow.CAP, "no_mmg_or_fea": True})
    with Job(parts=["controls", "cavities"]) as job:
        for stage in ("controls", "cavities"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_cap_cavity_probe.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker"])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
