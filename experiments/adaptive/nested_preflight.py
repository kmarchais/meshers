"""Bounded fixed-domain refinement preflight, without an equilibrium solve."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import microgen_workflow as workflow
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    np, *_ = workflow.imports()
    from nested_refinement import interpolation_check, refine

    row = json.loads((args.source / "result.json").read_text())
    assert row["case"] == "smooth-lattice" and row["qualified"] and row["order"] == 2
    assert len(row["solver_records"]) == 6 and row["geometry"]["geometry_admitted"]
    mesh = np.load(args.source / "mesh.npz")
    fields = np.load(args.source / "fields.npz")["displacement"]
    result = refine(mesh["points"], mesh["tetrahedra"], mesh["surface"], fields)
    checks = interpolation_check(mesh["points"], mesh["tetrahedra"], fields, result)
    assert checks["elements_checked"] == len(result["tetrahedra"])
    np.savez_compressed(args.run / "mesh.npz", **{key: result[key] for key in ("points", "tetrahedra", "surface")})
    np.savez_compressed(args.run / "prolongation.npz", **{key: result[key]
        for key in ("displacement", "parents", "parent_barycentric_maps")})
    geometry = workflow.geometry(args.run, "smooth-lattice", result["points"], result["tetrahedra"], result["surface"])
    volume_change = float(geometry["volume"] / row["geometry"]["volume"] - 1)
    assert abs(volume_change) < 1e-10
    workflow.save(args.run / "result.json", {"case": row["case"], "source": str(args.source),
        "nodes": len(result["points"]), "elements": len(result["tetrahedra"]), "order": 2,
        "qualified": geometry["geometry_admitted"], "geometry": geometry, "relative_volume_change": volume_change,
        "refinement": result["diagnostics"], "interpolation": checks,
        "equilibrium_solve_performed": False, "reference_promoted": False,
        "status": "nested-preflight-passed" if geometry["geometry_admitted"] else "geometry-rejected",
        "scope": "One conforming edge-star bisection stage with exact parent Tet10 map restriction and field prolongation; no convergence or continuum claim"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(args.source / ("smooth-lattice" + suffix), args.run / ("smooth-lattice" + suffix))
    paths = [args.source / name for name in ("result.json", "mesh.npz", "fields.npz",
        "smooth-lattice-samples.npz", "smooth-lattice-specimen.json")]
    paths.extend(workflow.HERE / name for name in ("nested_preflight.py", "nested_refinement.py",
        "test_nested_refinement.py", "job_parts.py", "test_job_parts.py", "microgen_workflow.py",
        "facet_incidence_audit.py", "scorecard.py"))
    paths.append(workflow.ROOT.parent / "periodic-microgen/experiments/periodic-fea/periodic.py")
    inputs = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": inputs, "memory_bytes": workflow.MEMORY,
        "threads": 1, "controls_timeout_s": 120, "worker_timeout_s": 300,
        "maximum_tetrahedra": workflow.CAP, "maximum_nodes": 250000,
        "minimum_corner_quality": .1, "maximum_bisections_per_parent": 1,
        "selection": "Descending physical edge length, complete periodic orbits and disjoint incident stars"})
    with Job(parts=["controls", "nested-preflight"]) as job:
        for stage in ("controls", "nested-preflight"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_nested_refinement.py",
                    "experiments/adaptive/test_job_parts.py"] if stage == "controls"
                    else [__file__, *sys.argv[1:], "--worker"])
                record = workflow.bounded([sys.executable, *arguments], args.run / (stage + ".log"),
                    part, timeout=120 if stage == "controls" else 300)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                values = {"seconds": record["wall_s"], "peak_gib": record["peak_rss_bytes"] / 1024**3}
                if stage == "nested-preflight":
                    row = json.loads((args.run / "result.json").read_text())
                    values.update(nodes=row["nodes"], qualified=int(row["qualified"]))
                job.result(stage, **values)
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
