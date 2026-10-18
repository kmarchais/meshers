"""Unranked stiffness diagnostics for returned, quality-rejected Microgen meshes."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(source, run, reference_path):
    np, *_ = workflow.imports()
    from homogenization import strain_tensor
    from periodic import periodic_pairs, tensor_difference
    from stress_recovery import solve
    from tet10_integration import integrated_tensors

    original = json.loads((source / "result.json").read_text())
    geometry = original["geometry"]
    assert original["qualified"] is False and geometry["geometry_admitted"] is False
    assert 0 < geometry["q_min_corner"] < 0.1
    assert geometry["sampled_two_sided_max"] <= 0.002 and abs(geometry["relative_volume_error"]) <= 0.01
    assert geometry["minimum_sampled_jacobian"] > 0 and geometry["projection_failed"] == 0
    assert original["periodic"]["periodic_triangles_match"] and len(original["mmg_calls"]) == 2
    assert original["route"] == "microgen-vtk-mmg-first-pass-override"
    mesh = np.load(source / "mesh.npz")
    points, tetrahedra = mesh["points"], mesh["tetrahedra"]
    assert len(points) <= 100_000 and len(tetrahedra) <= workflow.CAP and tetrahedra.shape[1] == 4
    pairs, _, _ = periodic_pairs(points, mesh["surface"])
    reference = json.loads(reference_path.read_text())["references"][original["case"]]
    assert reference["qualified_inputs"]
    started = time.perf_counter()
    fields, tensor, records = solve(points, tetrahedra)
    solve_seconds = time.perf_counter() - started
    assert len(records) == 6 and max(record["relative_residual"] for record in records) < 2e-9
    assert np.isfinite(tensor).all() and np.linalg.eigvalsh(tensor).min() > 0
    stress, energy, volume = integrated_tensors(points, tetrahedra, fields)
    reaction_error = float(np.linalg.norm(stress - tensor) / np.linalg.norm(tensor))
    energy_error = float(np.linalg.norm(energy - tensor) / np.linalg.norm(tensor))
    assert max(reaction_error, energy_error) < 1e-7
    assert abs(volume / geometry["volume"] - 1) < 1e-10
    jump_max = 0.0
    for load, engineering in enumerate(np.eye(6)):
        for axis, matching in enumerate(pairs):
            jump = fields[load, matching[:, 1]] - fields[load, matching[:, 0]]
            jump_max = max(jump_max, float(abs(jump - strain_tensor(engineering)[:, axis]).max(initial=0)))
    assert jump_max < 1e-9
    np.savez_compressed(run / (source.name + "-fields.npz"), displacement=fields)
    workflow.save(run / (source.name + ".json"), {"case": original["case"], "source": str(source),
        "qualified": False, "status": "unqualified-stiffness-diagnostic", "source_geometry": geometry,
        "nodes": len(points), "elements": len(tetrahedra), "stiffness": tensor.tolist(),
        "reference_disagreement": tensor_difference(tensor, np.asarray(reference["stiffness"])),
        "solver_records": records, "solve_s": solve_seconds, "reaction_relative_error": reaction_error,
        "hill_mandel_relative_error": energy_error, "periodic_jump_max": jump_max,
        "scope": "Saved quality-rejected mesh diagnostic only; excluded from admitted rankings and matched-workflow timing; no continuum or local-stress certification"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--sources", nargs="+", type=Path)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--reference", type=Path, default=workflow.HERE / "runs/g3-microgen-defaults-v4/summary.json")
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.run, args.reference)
        return
    assert args.sources and len({source.name for source in args.sources}) == len(args.sources)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    inputs = {str(source / name): workflow.digest(source / name) for source in args.sources
        for name in ("result.json", "mesh.npz")}
    inputs[str(args.reference)] = workflow.digest(args.reference)
    sources = {str(path): workflow.digest(path) for path in [Path(__file__), Path(workflow.__file__),
        workflow.HERE / "stress_recovery.py", workflow.HERE / "tet10_integration.py"]}
    workflow.save(args.run / "manifest.json", {"inputs": inputs, "source_sha256": sources,
        "memory_bytes": workflow.MEMORY, "worker_timeout_s": 600, "threads": 1,
        "admitted_ranking_eligible": False, "live_workflow_timing": False})
    with Job(parts=[source.name for source in args.sources]) as job:
        for source in args.sources:
            with job.part(source.name) as part:
                record = workflow.bounded([sys.executable, __file__, "--run", str(args.run),
                    "--worker", str(source), "--reference", str(args.reference)],
                    args.run / (source.name + ".log"), part, timeout=600)
                workflow.save(args.run / (source.name + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {source.name}.log"
                row = json.loads((args.run / (source.name + ".json")).read_text())
                assert row["qualified"] is False
                job.result(source.name, admitted=0, solve_seconds=row["solve_s"],
                    stiffness_disagreement=row["reference_disagreement"]["worst_relative_strain_energy_difference"],
                    peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == value for path, value in {**inputs, **sources}.items())


if __name__ == "__main__":
    main()
