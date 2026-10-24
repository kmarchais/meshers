"""One bounded fixed-boundary quality repair of saved genuine Microgen output."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import numpy as np

import microgen_workflow as workflow
from facet_incidence_audit import audit
from interior_star_probe import relocate_star
from job_parts import worker_part
from mmg_quality_localization import quality

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

SOURCE = workflow.HERE / "runs/g3-gyroid-first-pass-fem43-v1"


def repair(points, cells, surface, *, maximum_attempts=400):
    if not audit(points, cells, surface)["oriented_facet_consistency_passed"]:
        raise ValueError("Source fails oriented facet consistency; relocation is unsupported")
    moved = points.copy()
    boundary = np.zeros(len(points), dtype=bool)
    boundary[np.unique(surface)] = True
    attempts, accepted = [], []
    initial = quality(moved, cells)
    measured = initial.copy()
    for sweep in range(2):
        accepted_before = len(accepted)
        for element in np.argsort(measured):
            if measured[element] >= .1 or len(attempts) >= maximum_attempts:
                continue
            free = cells[element][~boundary[cells[element]]]
            if not len(free):
                attempts.append({"element": int(element), "sweep": sweep, "status": "all-boundary", "target_attained": False})
            for node in free:
                if len(attempts) >= maximum_attempts:
                    break
                candidate, record = relocate_star(moved, cells, int(node), target=.100001)
                record.update(element=int(element), sweep=sweep)
                attempts.append(record)
                if record["target_attained"]:
                    bad_before = int((measured < .1).sum())
                    moved[node] = candidate
                    measured = quality(moved, cells)
                    assert int((measured < .1).sum()) < bad_before
                    accepted.append(int(node))
                    break
        if not (measured < .1).any() or len(accepted) == accepted_before or len(attempts) >= maximum_attempts:
            break
    np.testing.assert_array_equal(moved[boundary], points[boundary])
    old_determinants = np.linalg.det(points[cells[:, 1:]] - points[cells[:, :1]])
    new_determinants = np.linalg.det(moved[cells[:, 1:]] - moved[cells[:, :1]])
    assert new_determinants.min() > 0
    volume_difference = float((new_determinants.sum() - old_determinants.sum()) / 6)
    assert abs(volume_difference) < 1e-12
    return moved, {"initial_bad_elements": int((initial < .1).sum()), "remaining_bad_elements": int((measured < .1).sum()),
        "initial_minimum_quality": float(initial.min()), "final_minimum_quality": float(measured.min()),
        "accepted_moves": len(accepted), "moved_vertices": len(set(accepted)), "attempts": attempts,
        "maximum_attempts": maximum_attempts, "maximum_sweeps": 2, "boundary_coordinates_unchanged": True,
        "signed_volume_difference": volume_difference, "minimum_signed_determinant": float(new_determinants.min())}


def worker(run, source_path):
    worker_started = time.monotonic()
    workflow.imports()
    from periodic import periodic_pairs

    source = json.loads((source_path / "result.json").read_text())
    case_name = source["case"]
    assert case_name in ("gyroid", "schwarz-solid")
    old_route = source["route"].startswith("microgen-vtk") and len(source["mmg_calls"]) == 1
    if old_route:
        assert source["first_pass_only"]
        if case_name == "gyroid":
            assert source["first_pass_fem_override"]
        else:
            assert source["resolution"] == 30 and source["first_pass_hmax_override"] == .1
    else:
        assert case_name == "gyroid"
        assert source["route"] == "microgen-explicit-preclipping-and-cap-retriangulation-mmg-tuned"
        assert source["saved_seed_replay"] and len(source["mmg_calls"]) == 2
        assert source["geometry"]["sampled_two_sided_max"] <= .002 and abs(source["geometry"]["relative_volume_error"]) <= .01
    archive = np.load(source_path / "mesh.npz")
    points, cells, surface = (archive[name] for name in ("points", "tetrahedra", "surface"))
    assert len(cells) <= workflow.CAP and len(points) <= 100000
    topology = audit(points, cells, surface)
    if not topology["oriented_facet_consistency_passed"]:
        workflow.save(run / "result.json", {"case": case_name, "source": str(source_path),
            "status": "source-incidence-rejected", "qualified": False, "facet_consistency": topology,
            "repair_performed": False, "stiffness_solve_performed": False,
            "scope": "Retired invalid one-pass source; no relocation, new mesh or FEA is attempted"})
        return
    started = time.monotonic()
    moved, diagnostics = repair(points, cells, surface)
    repair_s = time.monotonic() - started
    np.savez_compressed(run / "mesh.npz", points=moved, tetrahedra=cells, surface=surface)
    geometry = workflow.geometry(run, case_name, moved, cells, surface)
    _, _, periodic = periodic_pairs(moved, surface)
    row = {"case": case_name, "route": source["route"] + "-explicit-slsqp-interior-relocation",
        "parent_route": source["route"], "source": str(source_path), "parent_mmg_calls": source["mmg_calls"],
        "nodes": len(moved), "elements": len(cells), "order": 1, "qualified": geometry["geometry_admitted"],
        "status": "admitted-needs-reference" if geometry["geometry_admitted"] else "repair-incomplete",
        "geometry": geometry, "periodic": periodic, "repair_s": repair_s, "repair": diagnostics,
        "total_saved_mesh_study_s": time.monotonic() - worker_started,
        "live_workflow_timing": False, "stiffness_solve_performed": False,
        "scope": "Saved genuine Microgen VTK/MMG output plus explicitly added SciPy interior relocation; fixed boundary and topology, not default/public Microgen pipeline or a Meshers seed; no matched-accuracy claim"}
    workflow.save(run / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--source", type=Path, default=SOURCE)
    args = parser.parse_args()
    if args.worker:
        worker(args.run, args.source)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    case_name = json.loads((args.source / "result.json").read_text())["case"]
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), args.source / "result.json", args.source / "mesh.npz"]
    paths.extend(workflow.HERE / name for name in ("job_parts.py", "test_job_parts.py", "interior_star_probe.py",
        "test_interior_star_probe.py", "test_interior_relocation_trial.py", "microgen_workflow.py",
        "cap_quality_bound.py", "mmg_quality_localization.py", "facet_incidence_audit.py", "test_facet_incidence_audit.py"))
    for suffix in ("-specimen.json", "-samples.npz"):
        path = args.source / (case_name + suffix)
        paths.append(path)
        shutil.copyfile(path, args.run / path.name)
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY,
        "platform": sys.platform,
        "threads": 1, "worker_timeout_s": 600, "quality_gate": .1, "local_constraint_target": .100001,
        "maximum_attempts": 400, "maximum_sweeps": 2, "live_workflow_timing": False})
    with Job(parts=["controls", "repair-and-admit"]) as job:
        for stage in ("controls", "repair-and-admit"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_job_parts.py",
                    "experiments/adaptive/test_interior_star_probe.py", "experiments/adaptive/test_interior_relocation_trial.py"]
                    if stage == "controls" else [__file__, "--run", str(args.run), "--source", str(args.source), "--worker"])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"), part,
                    timeout=120 if stage == "controls" else 600)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
