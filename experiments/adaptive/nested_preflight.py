"""Bounded fixed-domain refinement preflight, without an equilibrium solve."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

import microgen_workflow as workflow
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def fixture_inputs(directory):
    """Use saved samples only when their hashes match their frozen provenance."""
    manifest = json.loads((directory / "manifest.json").read_text())
    paths = [directory / ("smooth-lattice" + suffix)
        for suffix in ("-samples.npz", "-specimen.json")]
    for path in paths:
        hashes = {value for name, value in manifest["inputs"].items() if Path(name).name == path.name}
        assert hashes == {workflow.digest(path)}, f"Fixture provenance mismatch: {path}"
    return paths


def worker(args):
    np, *_ = workflow.imports()
    from nested_refinement import interpolation_check, refine

    row = json.loads((args.source / "result.json").read_text())
    assert row["case"] == "smooth-lattice" and row["qualified"] and row["order"] == 2
    assert len(row["solver_records"]) == 6 and row["geometry"]["geometry_admitted"]
    mesh = np.load(args.source / "mesh.npz")
    fields = np.load(args.source / "fields.npz")["displacement"]
    preferred = None
    if args.prefer_box_curves:
        from periodic_crease_audit import features
        preferred = features(mesh)["curve_edges"]
    result = refine(mesh["points"], mesh["tetrahedra"], mesh["surface"], fields, preferred_edges=preferred)
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
    parser.add_argument("--fixture-source", type=Path,
        help="Original parent/prepared directory with hash-identical fixed samples, when absent from solve output")
    parser.add_argument("--prefer-box-curves", action="store_true",
        help="Prioritize actual wall/cap intersection edge orbits before other edges, with every existing guard retained")
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    args.run, args.source = args.run.resolve(), args.source.resolve()
    if args.fixture_source:
        args.fixture_source = args.fixture_source.resolve()
    if args.worker:
        worker(args)
        return
    fixtures = fixture_inputs(args.fixture_source or args.source)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    for path in fixtures:
        shutil.copyfile(path, args.run / path.name)
    paths = [args.source / name for name in ("result.json", "mesh.npz", "fields.npz")]
    paths.extend([*fixtures, (args.fixture_source or args.source) / "manifest.json"])
    paths.extend(workflow.HERE / name for name in ("nested_preflight.py", "nested_refinement.py",
        "test_nested_refinement.py", "job_parts.py", "test_job_parts.py", "microgen_workflow.py",
        "facet_incidence_audit.py", "scorecard.py"))
    if args.prefer_box_curves:
        paths.append(workflow.HERE / "periodic_crease_audit.py")
    paths.append(workflow.ROOT.parent / "periodic-microgen/experiments/periodic-fea/periodic.py")
    inputs = {str(path): workflow.digest(path) for path in paths}
    import meshers._meshers as native
    (args.run / "requirements.txt").write_text(subprocess.check_output(
        ["uv", "pip", "freeze", "--python", sys.executable], text=True))
    workflow.save(args.run / "manifest.json", {"inputs": inputs, "memory_bytes": workflow.MEMORY,
        "platform": platform.platform(), "architecture": platform.machine(), "python": sys.version,
        "memory_enforcement": "rlimit-as-and-process-tree-rss" if sys.platform == "linux"
            else "process-tree-rss-watchdog", "rss_poll_interval_s": .25,
        "native_path": native.__file__, "native_sha256": workflow.digest(Path(native.__file__)),
        "versions": {name: importlib.metadata.version(name) for name in
            ("numpy", "scipy", "fedoo", "microgen", "mmgpy", "pyamg", "pyvista", "vtk")},
        "fixture_source": str(args.fixture_source or args.source),
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY,
        "threads": 1, "controls_timeout_s": 120, "worker_timeout_s": 300,
        "maximum_tetrahedra": workflow.CAP, "maximum_nodes": 250000,
        "minimum_corner_quality": .1, "maximum_bisections_per_parent": 1,
        "prefer_box_curves": args.prefer_box_curves,
        "selection": ("Wall/cap intersection orbits first, then descending physical length; complete periodic orbits and disjoint incident stars"
            if args.prefer_box_curves else "Descending physical edge length, complete periodic orbits and disjoint incident stars")})
    with Job(parts=["controls", "nested-preflight"]) as job:
        for stage in ("controls", "nested-preflight"):
            with worker_part(job, stage) as part:
                arguments = (["-m", "pytest", "-q", "experiments/adaptive/test_nested_refinement.py",
                    "experiments/adaptive/test_job_parts.py", "experiments/adaptive/test_bounded.py"] if stage == "controls"
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
