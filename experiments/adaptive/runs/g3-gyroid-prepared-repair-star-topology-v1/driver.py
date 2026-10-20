"""Quality-constrained fixed-boundary retriangulation of one complete star."""

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
from cap_quality_bound import bounds
from cavity_milp import solve_cavity
from facet_incidence_audit import audit
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

SOURCE = workflow.HERE / "runs/g3-gyroid-prepared-mmg-interior-repair-v1"
WITNESS = workflow.HERE / "runs/g3-gyroid-prepared-repair-boundary-witness-audit-v1.json"


def worker(run):
    witness = json.loads(WITNESS.read_text())
    assert witness["passed"] and witness["checked_fixed_neighbor_obstruction"]
    assert not witness["other_interior_vertices_in_active_constraints"]
    diagnosis_path = Path(witness["source"]) / "result.json"
    assert workflow.digest(diagnosis_path) == witness["source_result_sha256"]
    diagnosis = json.loads(diagnosis_path.read_text())
    assert Path(diagnosis["source"]).resolve() == SOURCE.resolve()
    node = diagnosis["diagnostic"]["node"]
    mesh = np.load(SOURCE / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert len(cells) <= workflow.CAP and audit(points, cells, surface)["oriented_facet_consistency_passed"]
    selected = np.flatnonzero((cells == node).any(axis=1))
    original = cells[selected]
    target = oriented_boundary(original)
    assert target is not None and len(selected) <= 32 and node not in surface
    vertices = np.unique(original)
    row = {"source": str(SOURCE), "node": node, "incident_elements": selected.tolist(),
        "original_tetrahedra": len(original), "vertices": vertices.tolist(), "passed": False,
        "qualified": False, "source_unchanged": True, "stiffness_solve_performed": False,
        "maximum_vertices": 16, "quality_target": .100001}
    ceiling, _ = bounds(points[np.asarray(list(target))])
    row["minimum_fixed_star_face_ceiling"] = float(ceiling.min())
    if len(vertices) > 16:
        row["status"] = "vertex-cap"
    elif ceiling.min() < .1:
        row["status"] = "fixed-star-boundary-face-obstruction"
    else:
        alternative, outcome = solve_cavity(points, original, target, maximum_vertices=16,
            minimum_relative_determinant=1e-12, quality_target=.100001)
        row.update(status="local-topology-passed" if outcome["passed"] else "local-topology-not-found",
            passed=outcome["passed"], binary_solve=outcome,
            alternative_cells=None if alternative is None else alternative.tolist())
        if alternative is not None:
            assert len(cells) - len(original) + len(alternative) <= workflow.CAP
            row["unused_old_star_vertices"] = sorted(map(int, set(vertices) - set(alternative.ravel())))
    row["scope"] = "One complete fixed-boundary star, finite existing-vertex candidates and bounded binary solve; no applied changes or whole-mesh admission"
    workflow.save(run / "result.json", row)


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
    paths.extend(workflow.HERE / name for name in ("cavity_milp.py", "test_cavity_milp.py", "bistellar_probe.py",
        "cap_quality_bound.py", "facet_incidence_audit.py", "mmg_quality_localization.py", "microgen_workflow.py",
        "job_parts.py", "test_job_parts.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 120, "maximum_vertices": 16, "maximum_candidates": 1820,
        "maximum_original_tets": 32, "binary_time_limit_s": 2, "binary_node_limit": 1000,
        "quality_target": .100001, "source_unchanged": True})
    with Job(parts=["controls", "star-topology"]) as job:
        for stage in ("controls", "star-topology"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_cavity_milp.py", "experiments/adaptive/test_job_parts.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--worker"])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
