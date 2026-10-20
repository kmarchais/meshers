"""Apply the complete measured cap plan to a new copy, then audit it."""

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
from cap_quality_bound import bounds, caps
from facet_incidence_audit import audit
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

RUNS = workflow.HERE / "runs"
CAVITIES = RUNS / "g3-gyroid-cap-cavity-preflight-v2/result.json"
CAPS = RUNS / "g3-gyroid-cap-flip-preflight-v2/result.json"
EXPANSIONS = RUNS / "g3-gyroid-cap-cavity-fullring-v1/result.json"


def replace_plans(cells, surface, plans):
    used_cells, old_faces, new_cells, new_faces = set(), set(), [], []
    original_faces = {tuple(sorted(map(int, face))) for face in surface}
    for plan in plans:
        indices = set(map(int, plan["incident_elements"]))
        patch = {tuple(sorted(map(int, face))) for face in plan["old_cap_faces"]}
        if indices & used_cells or patch & old_faces:
            raise ValueError("Replacement plans overlap")
        if not indices or min(indices) < 0 or max(indices) >= len(cells) or not patch <= original_faces:
            raise ValueError("Replacement plan does not match source")
        used_cells |= indices
        old_faces |= patch
        new_cells.extend(plan["alternative_cells"])
        new_faces.extend(plan["new_cap_faces"])
    replaced = [cell for index, cell in enumerate(cells) if index not in used_cells]
    replaced.extend(new_cells)
    boundary = [face for face in surface if tuple(sorted(map(int, face))) not in old_faces]
    boundary.extend(new_faces)
    return np.asarray(replaced), np.asarray(boundary)


def worker(run, plans_path):
    workflow.imports()
    from periodic import periodic_pairs

    mesh = np.load(SOURCE / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert audit(points, cells, surface)["oriented_facet_consistency_passed"]
    cap_faces, _ = caps(mesh)
    candidates = json.loads(CAPS.read_text())["candidates"]
    plans = []
    if plans_path is not None:
        prepared = json.loads(plans_path.read_text())
        assert prepared["source"] == str(SOURCE) and prepared["minimum_relative_determinant"] == 1e-12
        plans = prepared["complete_plans"]
        assert plans is not None
    else:
        for row in json.loads(CAVITIES.read_text())["candidates"]:
            if row["cavity_passed"]:
                candidate = next(item for item in candidates if item["shared_edge"] == row["edge"])
                plans.append({"incident_elements": row["incident_elements"], "alternative_cells": row["alternative_cells"],
                    "old_cap_faces": cap_faces[row["cap_triangles"]].tolist(), "new_cap_faces": candidate["alternative_triangles"]})
        expanded = json.loads(EXPANSIONS.read_text())
        assert expanded["solved_cavities"] == 2
        plans.extend(expanded["cavities"])
    assert len(plans) == 24
    replaced, boundary = replace_plans(cells, surface, plans)
    minimum_relative = min(float(np.linalg.det(points[np.asarray(plan["alternative_cells"])[:, 1:]]
        - points[np.asarray(plan["alternative_cells"])[:, :1]]).min() /
        np.linalg.det(points[cells[plan["incident_elements"]][:, 1:]]
        - points[cells[plan["incident_elements"]][:, :1]]).sum()) for plan in plans)
    if plans_path is not None:
        assert minimum_relative > 1e-12
    assert len(replaced) <= workflow.CAP
    consistency = audit(points, replaced, boundary)
    assert consistency["oriented_facet_consistency_passed"] and not consistency["nonpositive_corner_determinants"]
    _, _, periodic = periodic_pairs(points, boundary)
    _, triangles = caps({"points": points, "surface": boundary})
    ceiling, _ = bounds(triangles)
    assert (ceiling >= .1).all()
    old_caps, _ = caps(mesh)
    new_caps, _ = caps({"points": points, "surface": boundary})
    old_keys = {tuple(sorted(map(int, face))) for face in old_caps}
    new_keys = {tuple(sorted(map(int, face))) for face in new_caps}
    old_walls = {tuple(sorted(map(int, face))) for face in surface} - old_keys
    new_walls = {tuple(sorted(map(int, face))) for face in boundary} - new_keys
    assert old_walls == new_walls
    volume = float(np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]]).sum())
    new_volume = float(np.linalg.det(points[replaced[:, 1:]] - points[replaced[:, :1]]).sum())
    relative_volume = abs(volume - new_volume) / volume
    assert relative_volume < 1e-12
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(SOURCE / ("gyroid" + suffix), run / ("gyroid" + suffix))
    geometry = workflow.geometry(run, "gyroid", points, replaced, boundary)
    np.savez_compressed(run / "mesh.npz", points=points, tetrahedra=replaced, surface=boundary)
    workflow.save(run / "result.json", {"source": str(SOURCE), "case": "gyroid", "nodes": len(points),
        "elements": len(replaced), "plans_applied": len(plans), "plans": plans, "geometry": geometry,
        "periodic": periodic, "facet_consistency": consistency, "cap_faces_below_quality_gate": 0,
        "cap_minimum_quality_ceiling": float(ceiling.min()), "wall_triangles_identical": True,
        "original_points_identical": True, "relative_volume_change": relative_volume,
        "minimum_replacement_relative_determinant": minimum_relative,
        "source_unchanged": True, "qualified": False, "mmg_calls": [], "stiffness_solve_performed": False,
        "scope": "Explicit cap retriangulation of a non-default pretriangulated Microgen seed; saved-seed preparation, not an admitted MMG mesh or live workflow"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--plans", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.run, args.plans)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), CAVITIES, CAPS, EXPANSIONS, SOURCE / "mesh.npz"]
    if args.plans is not None:
        paths.append(args.plans)
    paths.extend(SOURCE / ("gyroid" + suffix) for suffix in ("-samples.npz", "-specimen.json"))
    paths.extend(workflow.HERE / name for name in ("test_cap_replacement_trial.py", "cap_flip_probe.py",
        "cap_quality_bound.py", "facet_incidence_audit.py", "microgen_workflow.py", "job_parts.py", "test_job_parts.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "tet_cap": workflow.CAP,
        "no_mmg_or_fea": True, "explicit_nondefault_seed_preparation": True})
    with Job(parts=["controls", "replacement"]) as job:
        for stage in ("controls", "replacement"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_cap_replacement_trial.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker"])
                if stage != "controls" and args.plans is not None:
                    arguments.extend(["--plans", str(args.plans)])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
