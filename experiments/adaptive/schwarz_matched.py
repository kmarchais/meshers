"""Repeat two live Schwarz workflows at the fixed one-percent stiffness target."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import statistics
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job
import microgen_workflow as workflow


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    source = workflow.HERE / "runs/g3-microgen-defaults-v4"
    source_files = [workflow.HERE / name for name in (
        "schwarz_matched.py", "microgen_workflow.py", "microgen_tuned.py",
        "admission_pilot.py", "stress_recovery.py", "gradient_tie_probe.py")]
    inputs = [source / name for name in (
        "summary.json", "schwarz-solid-specimen.json", "schwarz-solid-samples.npz")]
    hashes = {str(path): workflow.digest(path) for path in source_files + inputs}
    np, *_ = workflow.imports()
    import meshers._meshers
    from periodic import tensor_difference

    native = Path(meshers._meshers.__file__)
    hashes[str(native)] = workflow.digest(native)
    reference = json.loads((source / "summary.json").read_text())["references"]["schwarz-solid"]
    assert reference["qualified_inputs"]
    workflow.save(args.run / "reference.json", reference)
    workflow.save(args.run / "manifest.json", {"inputs_and_sources": hashes,
        "packages": {name: importlib.metadata.version(name) for name in (
            "meshers", "microgen", "mmgpy", "fedoo", "numpy", "scipy", "vtk", "pyvista")},
        "memory_bytes": workflow.MEMORY, "worker_timeout_s": 600, "threads": 1,
        "tetrahedron_cap": workflow.CAP, "solve_node_cap": 100000,
        "stiffness_target": 0.01, "solver_formulation": "six prescribed engineering strains",
        "repetitions": 5, "warmups_per_arm": 1,
        "microgen_recipe": "30-point VTK; actual first MMG call at hausd .0005 hmax .05 hgrad 1.3; explicit second-call omission",
        "meshers_recipe": "8 cells, Tet10, no adaptation, one optimization pass, native linear tolerance 1; independent admission .002"})
    schedule = [("warmup", arm) for arm in ("microgen", "meshers")]
    for repetition in range(1, 6):
        arms = ("microgen", "meshers") if repetition % 2 else ("meshers", "microgen")
        schedule.extend((f"repeat-{repetition}", arm) for arm in arms)
    labels = {f"{phase}-{arm}": 5 if arm == "microgen" else 1 for phase, arm in schedule}
    rows, workers = [], {}
    with Job(parts=labels) as job:
        for phase, arm in schedule:
            label = f"{phase}-{arm}"
            output = args.run / label
            output.mkdir()
            for path in inputs[1:]:
                shutil.copyfile(path, output / path.name)
            script = workflow.HERE / ("microgen_tuned.py" if arm == "microgen" else "admission_pilot.py")
            shutil.copyfile(script, output / "driver.py")
            command = [sys.executable, str(script), "--run", str(output), "--worker"]
            if arm == "microgen":
                command.extend(["--source", str(source), "--resolution", "30",
                    "--first-pass-hausd", "0.0005", "--first-pass-hmax", "0.05",
                    "--first-pass-only", "--strain-controlled"])
            else:
                command.extend(["--case", "schwarz-solid", "--cells", "8", "--order", "2",
                    "--native-tolerance", "1", "--no-adapt", "--solve"])
            with job.part(label) as part:
                workers[label] = workflow.bounded(command, output / "worker.log", part, timeout=600)
                workflow.save(args.run / "workers.json", workers)
                assert workers[label]["returncode"] == 0, f"Read {output / 'worker.log'}"
                row = json.loads((output / "result.json").read_text())
                assert row["qualified"] and row["geometry"]["geometry_admitted"]
                assert len(row["solver_records"]) == 6
                assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
                tensor = np.asarray(row["stiffness"])
                assert tensor.shape == (6, 6) and np.isfinite(tensor).all()
                assert np.linalg.eigvalsh(tensor).min() > 0
                disagreement = tensor_difference(tensor, np.asarray(reference["stiffness"]))
                error = disagreement["worst_relative_strain_energy_difference"]
                assert error <= 0.01, f"{label} misses stiffness target: {error}"
                rows.append({"phase": phase, "arm": arm, "source": str(output),
                    "nodes": row["nodes"], "elements": row["elements"], "geometry": row["geometry"],
                    "stiffness_disagreement": disagreement, "generation_s": row["generation_s"],
                    "solve_s": row["solve_s"], "total_workflow_s": row["total_workflow_s"]})
                workflow.save(args.run / "checkpoint.json", {"completed": rows})
                job.result(label, seconds=row["total_workflow_s"], nodes=row["nodes"],
                    stiffness_error=error, peak_gib=workers[label]["peak_rss_bytes"] / 1024**3)
    measured = [row for row in rows if row["phase"] != "warmup"]
    statistics_by_arm = {}
    for arm in ("microgen", "meshers"):
        chosen = [row for row in measured if row["arm"] == arm]
        statistics_by_arm[arm] = {}
        for metric in ("generation_s", "solve_s", "total_workflow_s"):
            values = [row[metric] for row in chosen]
            statistics_by_arm[arm][metric] = {"median": statistics.median(values),
                "minimum": min(values), "maximum": max(values)}
    paired = []
    for repetition in range(1, 6):
        pair = {row["arm"]: row for row in measured if row["phase"] == f"repeat-{repetition}"}
        paired.append(pair["microgen"]["total_workflow_s"] / pair["meshers"]["total_workflow_s"])
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())
    workflow.save(args.run / "summary.json", {"case": "schwarz-solid", "rows": rows,
        "statistics": statistics_by_arm, "paired_workflow_speedups": paired,
        "all_pairs_meshers_faster": min(paired) > 1,
        "scope": "One-specimen matched 1 percent target against explicit one-pass Microgen variant; not default pipeline or four-fixture victory"})


if __name__ == "__main__":
    main()
