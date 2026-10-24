"""One fresh Schwarz Tet10 workflow pair with the same explicit AMG preset."""

import argparse
from contextlib import nullcontext
import importlib.metadata
import json
import os
import platform
from pathlib import Path
import shutil
import sys
import time

import microgen_workflow as workflow
from job_parts import worker_part

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

STAGE_TIMEOUTS = {"controls": 120, "microgen-parent": 600, "microgen-tet10": 1800,
    "meshers-tet10": 600, "microgen-audit": 300, "microgen-reference": 60,
    "meshers-audit": 300, "meshers-reference": 60}


def stage_names(native_first=False, repair_interior=False):
    names = list(STAGE_TIMEOUTS)
    if repair_interior:
        names.insert(names.index("microgen-tet10"), "microgen-repair")
    if native_first:
        names.remove("meshers-tet10")
        names.insert(1, "meshers-tet10")
    return names


def comparison_summary(parent, elevated, native, comparisons, repair=None):
    for row in (repair or parent, elevated, native):
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
    if repair is not None:
        assert repair["repair"]["boundary_coordinates_unchanged"] and repair["repair"]["remaining_bad_elements"] == 0
        assert abs(repair["repair"]["signed_volume_difference"]) < 1e-12
        microgen_total += repair["total_saved_mesh_study_s"]
        microgen_phases += repair["repair_s"]
    native_phases = native["generation_s"] + native["solve_s"]
    assert 0 < microgen_phases <= microgen_total and 0 < native_phases <= native_total
    return {"case": "schwarz-solid", "stiffness_target": .001,
        "both_pass_declared_reference_models": admitted,
        "live_workflow_timing": True, "repetitions": 1, "warmups": 0,
        "science_audited_workflow_s": {"microgen": microgen_total, "meshers": native_total},
        "generation_lift_solve_s": {"microgen": microgen_phases, "meshers": native_phases},
        "single_pair_workflow_ratio_microgen_over_meshers": microgen_total / native_total if admitted else None,
        "single_pair_phase_ratio_microgen_over_meshers": microgen_phases / native_phases if admitted else None,
        "comparison_reports": comparisons, "explicit_interior_repair": repair is not None,
        "timing_scope": "Fresh VTK/MMG parent plus immediately elevated fresh parent, versus fresh native mesh; science-audited total charges parent and final audits; subtotal uses recorded generation/lift/solve timers and retains checkpoint I/O within those phases",
        "scope": "Single same-host pair at supplementary 0.1 percent target, explicit one-pass Microgen override plus manual straight Tet10 and direct fedoo, not default/public Microgen Tet10 or a repeated/four-fixture advantage"}


def run_pair(run, *, native_first=False, job=None, prefix="", deadline=None, repair_interior=False):
    wall_start, monotonic_start = time.time(), time.monotonic()
    run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, run / "driver.py")
    source = workflow.HERE / "runs/g3-microgen-defaults-v4"
    directories = {name: run / name for name in ("microgen-parent", "microgen-tet10", "meshers-tet10", "audits")}
    if repair_interior:
        directories["microgen-repair"] = run / "microgen-repair"
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
        "test_schwarz_tet10_live.py", "schwarz_tet10_repeated.py", "test_schwarz_tet10_repeated.py",
        "cap_quality_bound.py", "gradient_tie_probe.py", "job_parts.py", "test_job_parts.py", "facet_incidence_audit.py",
        "interior_relocation_trial.py", "interior_star_probe.py"))
    paths.extend(workflow.ROOT / "crates/meshers-python/python/meshers" / name for name in ("__init__.py", "expression.py"))
    paths.extend(workflow.ROOT.parent / "periodic-microgen/experiments/periodic-fea" / name
                 for name in ("homogenization.py", "periodic.py"))
    paths.extend(source / name for name in ("summary.json", "schwarz-solid-specimen.json", "schwarz-solid-samples.npz"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    timeouts = dict(STAGE_TIMEOUTS)
    if repair_interior:
        timeouts["microgen-repair"] = 600
    workflow.save(run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "platform": platform.platform(), "architecture": platform.machine(),
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY,
        "memory_enforcement": "process-tree-rss-watchdog" if sys.platform == "darwin" else "rlimit-as-and-process-tree-rss",
        "threads": 1, "solve_node_cap": 250000, "tetrahedron_cap": workflow.CAP,
        "stiffness_target": .001, "repetitions": 1, "live_workflow_timing": True,
        "native_first": native_first, "stage_timeouts_s": timeouts,
        "explicit_interior_repair": repair_interior,
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
            "experiments/adaptive/test_reference_admission.py", "experiments/adaptive/test_cubic_reference.py",
            "experiments/adaptive/test_schwarz_tet10_repeated.py", "experiments/adaptive/test_job_parts.py"]),
        ("microgen-parent", 600, [str(workflow.HERE / "microgen_tuned.py"), "--worker", "--run", str(parent),
            "--source", str(source), "--resolution", "30", "--first-pass-hausd", ".0005",
            "--first-pass-hmax", ".1", "--first-pass-only", "--geometry-only"]),
        ("microgen-tet10", 1800, [str(workflow.HERE / "microgen_tet10_control.py"), "--worker", "--run", str(elevated),
            "--source", str(parent), "--amg"]),
        ("meshers-tet10", 600, [str(workflow.HERE / "admission_pilot.py"), "--worker", "--run", str(native),
            "--case", "schwarz-solid", "--cells", "16", "--order", "2", "--native-tolerance", "1",
            "--no-adapt", "--solve", "--amg"]),
    ]
    if repair_interior:
        commands.insert(2, ("microgen-repair", 600, [str(workflow.HERE / "interior_relocation_trial.py"),
            "--worker", "--run", str(directories["microgen-repair"]), "--source", str(parent)]))
        for name, _, arguments in commands:
            if name == "microgen-tet10":
                arguments[arguments.index("--source") + 1] = str(directories["microgen-repair"])
    for arm, directory in (("microgen", elevated), ("meshers", native)):
        commands.append((arm + "-audit", 300, [str(workflow.HERE / "tet10_integration.py"), "--run", str(audits),
            "--saved-worker", str(directory)]))
        commands.append((arm + "-reference", 60, [str(workflow.HERE / "stiffness_pair.py"),
            "--first", str(directory / "result.json"), "--reference-summary", str(source / "summary.json"),
            "--orders", "1", "2", "3", "4", "--target", ".001", "--output", str(run / (arm + "-reference.json"))]))
    assert {name: timeout for name, timeout, _ in commands} == timeouts
    commands = {name: arguments for name, timeout, arguments in commands}
    records = {}
    context = Job(parts=stage_names(native_first, repair_interior)) if job is None else nullcontext(job)
    with context as active_job:
        for name in stage_names(native_first, repair_interior):
            timeout = timeouts[name]
            if deadline is not None:
                timeout = min(timeout, deadline - time.monotonic())
                assert timeout > 0, "Whole-experiment deadline reached before the next worker"
            with worker_part(active_job, prefix + name) as part:
                records[name] = workflow.bounded([sys.executable, *commands[name]], run / (name + ".log"), part, timeout=timeout)
                workflow.save(run / "workers.json", records)
                assert records[name]["returncode"] == 0, f"Read {name}.log"
                if name in ("microgen-parent", "microgen-repair", "microgen-tet10", "meshers-tet10"):
                    row = json.loads((directories[name] / "result.json").read_text())
                    if name == "microgen-parent" and repair_interior:
                        prospective = dict(row["geometry"], q_min_corner=.1)
                        assert workflow.geometry_is_admitted(prospective) and "periodic" in row
                    else:
                        assert row["qualified"] and row["geometry"]["geometry_admitted"]
                    if name in ("microgen-tet10", "meshers-tet10"):
                        assert row["status"] == "solved-needs-reference"
                active_job.result(prefix + name, seconds=records[name]["wall_s"], peak_gib=records[name]["peak_rss_bytes"] / 1024**3)
    comparisons = {arm: json.loads((run / (arm + "-reference.json")).read_text()) for arm in ("microgen", "meshers")}
    repair = json.loads((directories["microgen-repair"] / "result.json").read_text()) if repair_interior else None
    report = comparison_summary(*(json.loads((directory / "result.json").read_text()) for directory in (parent, elevated, native)), comparisons, repair)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())
    report.update(inputs_unchanged=True, workers=records)
    report["artifact_sha256"] = {str(directory / name): workflow.digest(directory / name)
        for directory in directories.values() for name in ("result.json", "mesh.npz", "fields.npz")
        if (directory / name).exists()}
    wall_elapsed, monotonic_elapsed = time.time() - wall_start, time.monotonic() - monotonic_start
    report.update(native_first=native_first, elapsed_wall_s=wall_elapsed, elapsed_monotonic_s=monotonic_elapsed,
        wall_minus_monotonic_s=wall_elapsed - monotonic_elapsed,
        timing_clock_consistent=abs(wall_elapsed - monotonic_elapsed) <= 5)
    workflow.save(run / "summary.json", report)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--native-first", action="store_true")
    parser.add_argument("--repair-interior", action="store_true")
    args = parser.parse_args()
    run_pair(args.run, native_first=args.native_first, repair_interior=args.repair_interior)


if __name__ == "__main__":
    main()
