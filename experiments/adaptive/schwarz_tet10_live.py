"""One fresh Schwarz Tet10 workflow pair with the same explicit AMG preset."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def comparison_summary(parent, elevated, native, comparisons):
    for row in (parent, elevated, native):
        assert row["qualified"] and row["geometry"]["geometry_admitted"]
        assert row["case"] == "schwarz-solid" and row["elements"] <= workflow.CAP
        assert "periodic" in row
    assert parent["first_pass_only"] and len(parent["mmg_calls"]) == 1
    assert parent["resolution"] == 30 and parent["first_pass_hmax_override"] == .1
    for row in (elevated, native):
        assert row["order"] == 2 and row["nodes"] <= 250000
        assert row["solver"] == "cached-elastic-block-amg-cg"
        assert len(row["solver_records"]) == 6
        assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
    assert elevated["quadrature_sensitivity_passed"]
    assert elevated["quadrature_sensitivity_limit"] <= 1e-4
    assert elevated["solver_quadrature_points"] == 4 and native["solver_quadrature_points"] == 15
    assert all(report["stiffness_target"] == .001 for report in comparisons.values())
    admitted = all(report["passes_target_under_all_declared_models"] for report in comparisons.values())
    microgen_total = parent["total_workflow_s"] + elevated["total_saved_mesh_study_s"]
    native_total = native["total_workflow_s"]
    microgen_phases = parent["generation_s"] + elevated["lift_and_conversion_s"] + elevated["solve_s"]
    native_phases = native["generation_s"] + native["solve_s"]
    assert 0 < microgen_phases <= microgen_total and 0 < native_phases <= native_total
    return {"case": "schwarz-solid", "stiffness_target": .001,
        "both_pass_declared_reference_models": admitted,
        "live_workflow_timing": True, "repetitions": 1, "warmups": 0,
        "science_audited_workflow_s": {"microgen": microgen_total, "meshers": native_total},
        "generation_lift_solve_s": {"microgen": microgen_phases, "meshers": native_phases},
        "single_pair_workflow_ratio_microgen_over_meshers": microgen_total / native_total if admitted else None,
        "single_pair_phase_ratio_microgen_over_meshers": microgen_phases / native_phases if admitted else None,
        "comparison_reports": comparisons,
        "timing_scope": "Fresh VTK/MMG parent plus immediately elevated fresh parent, versus fresh native mesh; science-audited total charges parent and final audits; subtotal uses recorded generation/lift/solve timers and retains checkpoint I/O within those phases",
        "scope": "Single G3 pair at supplementary 0.1 percent target, explicit one-pass Microgen override plus manual straight Tet10 and direct fedoo, not default/public Microgen Tet10 or a repeated/four-fixture advantage"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    source = workflow.HERE / "runs/g3-microgen-defaults-v4"
    directories = {name: args.run / name for name in ("microgen-parent", "microgen-tet10", "meshers-tet10", "audits")}
    for directory in directories.values():
        directory.mkdir()
        for suffix in ("-specimen.json", "-samples.npz"):
            shutil.copyfile(source / ("schwarz-solid" + suffix), directory / ("schwarz-solid" + suffix))
    import meshers._meshers
    import microgen.remesh
    paths = [Path(__file__), Path(meshers._meshers.__file__), Path(microgen.remesh.__file__)]
    paths.extend(workflow.HERE / name for name in ("microgen_workflow.py", "microgen_tuned.py",
        "microgen_tet10_control.py", "admission_pilot.py", "stress_recovery.py", "tet10_integration.py",
        "amg_reference.py", "stiffness_pair.py", "test_cubic_reference.py", "test_reference_admission.py",
        "test_schwarz_tet10_live.py", "cap_quality_bound.py", "gradient_tie_probe.py"))
    paths.extend(workflow.ROOT / "crates/meshers-python/python/meshers" / name for name in ("__init__.py", "expression.py"))
    paths.extend(workflow.ROOT.parent / "periodic-microgen/experiments/periodic-fea" / name
                 for name in ("homogenization.py", "periodic.py"))
    paths.extend(source / name for name in ("summary.json", "schwarz-solid-specimen.json", "schwarz-solid-samples.npz"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "solve_node_cap": 250000, "tetrahedron_cap": workflow.CAP,
        "stiffness_target": .001, "repetitions": 1, "live_workflow_timing": True,
        "solver": "cached-elastic-block-amg-cg", "cg_rtol": 1e-12, "cg_atol": 1e-14,
        "packages": {name: importlib.metadata.version(name) for name in (
            "meshers", "microgen", "mmgpy", "fedoo", "numpy", "scipy", "vtk", "pyvista", "pyamg")},
        "microgen_recipe": "Fresh VTK30 and actual MMG nofem hausd.0005 hmax.1 hgrad1.3; explicit second-call omission; manual straight Tet10; four-point exact affine rule",
        "meshers_recipe": "Fresh curved16, no adaptation, one optimization pass, minimumq.1, native linear tolerance1; fifteen-point curved rule",
        "independent_admission": "H.002 volume1percent q.1 periodic caps and positive sampled Jacobians",
        "reference_models": "Previously retained Schwarz16/20 sensitivity, two proxies and assumed orders1,2,3,4"})
    parent, elevated, native, audits = (directories[name] for name in ("microgen-parent", "microgen-tet10", "meshers-tet10", "audits"))
    commands = [
        ("controls", 120, ["-m", "pytest", "-q", "experiments/adaptive/test_schwarz_tet10_live.py",
            "experiments/adaptive/test_reference_admission.py", "experiments/adaptive/test_cubic_reference.py"]),
        ("microgen-parent", 600, [str(workflow.HERE / "microgen_tuned.py"), "--worker", "--run", str(parent),
            "--source", str(source), "--resolution", "30", "--first-pass-hausd", ".0005",
            "--first-pass-hmax", ".1", "--first-pass-only", "--geometry-only"]),
        ("microgen-tet10", 1800, [str(workflow.HERE / "microgen_tet10_control.py"), "--worker", "--run", str(elevated),
            "--source", str(parent), "--amg"]),
        ("meshers-tet10", 600, [str(workflow.HERE / "admission_pilot.py"), "--worker", "--run", str(native),
            "--case", "schwarz-solid", "--cells", "16", "--order", "2", "--native-tolerance", "1",
            "--no-adapt", "--solve", "--amg"]),
    ]
    for arm, directory in (("microgen", elevated), ("meshers", native)):
        commands.append((arm + "-audit", 300, [str(workflow.HERE / "tet10_integration.py"), "--run", str(audits),
            "--saved-worker", str(directory)]))
        commands.append((arm + "-reference", 60, [str(workflow.HERE / "stiffness_pair.py"),
            "--first", str(directory / "result.json"), "--reference-summary", str(source / "summary.json"),
            "--orders", "1", "2", "3", "4", "--target", ".001", "--output", str(args.run / (arm + "-reference.json"))]))
    records = {}
    with Job(parts={name: timeout for name, timeout, _ in commands}) as job:
        for name, timeout, arguments in commands:
            with job.part(name) as part:
                records[name] = workflow.bounded([sys.executable, *arguments], args.run / (name + ".log"), part, timeout=timeout)
                workflow.save(args.run / "workers.json", records)
                assert records[name]["returncode"] == 0, f"Read {name}.log"
                if name in ("microgen-parent", "microgen-tet10", "meshers-tet10"):
                    row = json.loads((directories[name] / "result.json").read_text())
                    assert row["qualified"] and row["geometry"]["geometry_admitted"]
                    if name != "microgen-parent":
                        assert row["status"] == "solved-needs-reference"
                job.result(name, seconds=records[name]["wall_s"], peak_gib=records[name]["peak_rss_bytes"] / 1024**3)
    comparisons = {arm: json.loads((args.run / (arm + "-reference.json")).read_text()) for arm in ("microgen", "meshers")}
    report = comparison_summary(*(json.loads((directory / "result.json").read_text()) for directory in (parent, elevated, native)), comparisons)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())
    report.update(inputs_unchanged=True, workers=records)
    report["artifact_sha256"] = {str(directory / name): workflow.digest(directory / name)
        for directory in (parent, elevated, native) for name in ("result.json", "mesh.npz", "fields.npz")
        if (directory / name).exists()}
    workflow.save(args.run / "summary.json", report)


if __name__ == "__main__":
    main()
