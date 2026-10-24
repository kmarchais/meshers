"""Five interleaved fresh Tet10 workflow pairs after one warmup per arm."""

import argparse
import json
from pathlib import Path
import platform
import shutil
import statistics
import sys
import time

import microgen_workflow as workflow
from schwarz_tet10_live import Job, STAGE_TIMEOUTS, run_pair, stage_names

PILOT = workflow.HERE / "runs/g3-schwarz-tet10-live-amg-pilot-v1"
QUADRATURE = workflow.HERE / "runs/g3-schwarz-tet10-live-amg-pilot-quadrature-v1"


def schedule():
    return [("warmup", False)] + [(f"repeat-{number}", bool(number % 2)) for number in range(1, 6)]


def summarize(pairs):
    assert [pair["phase"] for pair in pairs] == [phase for phase, _ in schedule()]
    assert all(pair["report"]["both_pass_declared_reference_models"] for pair in pairs)
    assert all(pair["report"]["timing_clock_consistent"] for pair in pairs)
    assert all(pair["report"]["stiffness_target"] == .001 for pair in pairs)
    assert all(pair["native_first"] == native_first for pair, (_, native_first) in zip(pairs, schedule(), strict=True))
    repaired = pairs[0]["report"].get("explicit_interior_repair", False)
    assert all(pair["report"].get("explicit_interior_repair", False) == repaired for pair in pairs)
    measured = pairs[1:]
    statistics_by_metric, ratios = {}, {}
    for metric in ("science_audited_workflow_s", "generation_lift_solve_s"):
        statistics_by_metric[metric] = {}
        for arm in ("microgen", "meshers"):
            values = [pair["report"][metric][arm] for pair in measured]
            assert min(values) > 0
            statistics_by_metric[metric][arm] = {"median": statistics.median(values),
                "minimum": min(values), "maximum": max(values)}
        ratios[metric] = [pair["report"][metric]["microgen"] / pair["report"][metric]["meshers"] for pair in measured]
    return {"case": "schwarz-solid", "stiffness_target": .001, "repetitions": 5, "warmups_per_arm": 1,
        "all_pairs_qualified": True, "statistics": statistics_by_metric, "explicit_interior_repair": repaired,
        "paired_ratios_microgen_over_meshers": ratios,
        "meshers_faster_in_every_measured_pair": all(min(values) > 1 for values in ratios.values()),
        "pairs": pairs,
        "scope": "Repeated same-host common-target Schwarz comparison against explicit one-pass actual Microgen plus manual straight Tet10/direct fedoo, with any declared interior repair charged; common AMG; not default pipeline, identical measured errors, or four-fixture victory"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--pilot", type=Path, default=PILOT)
    parser.add_argument("--quadrature", type=Path, default=QUADRATURE)
    parser.add_argument("--repair-interior", action="store_true")
    args = parser.parse_args()
    started = time.monotonic()
    deadline = started + 3600
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    pilot = json.loads((args.pilot / "summary.json").read_text())
    assert pilot["both_pass_declared_reference_models"] and pilot["stiffness_target"] == .001
    assert pilot.get("explicit_interior_repair", False) == args.repair_interior
    checked = json.loads((args.quadrature / "meshers-tet10.json").read_text())
    checked_source = Path(checked["source"].replace("/home/kmarchais/dev/", str(workflow.ROOT.parent) + "/", 1))
    assert checked_source.resolve() == (args.pilot / "meshers-tet10").resolve()
    assert checked["rules"][-1]["energy_vs_solved_tensor"]["worst_relative_strain_energy_difference"] < 1e-4
    assert checked["rules"][-1]["energy_vs_previous_rule"]["worst_relative_strain_energy_difference"] < 1e-4
    verified = json.loads((args.quadrature / "manifest.json").read_text())["inputs"]
    artifacts = pilot.get("artifact_sha256")
    if artifacts is None:
        artifacts = json.loads((workflow.HERE / "runs/g3-schwarz-tet10-live-amg-pilot-artifacts-v1.json").read_text())["artifact_sha256"]
    evidence = {}
    for path, digest in {**verified, **artifacts}.items():
        resolved = Path(path.replace("/home/kmarchais/dev/", str(workflow.ROOT.parent) + "/", 1))
        assert workflow.digest(resolved) == digest
        evidence[str(resolved)] = digest
    paths = [Path(__file__), workflow.HERE / "schwarz_tet10_live.py", workflow.HERE / "test_schwarz_tet10_repeated.py",
        args.pilot / "summary.json", args.quadrature / "meshers-tet10.json", args.quadrature / "manifest.json"]
    hashes = {**evidence, **{str(path): workflow.digest(path) for path in paths}}
    timeouts = dict(STAGE_TIMEOUTS)
    if args.repair_interior:
        timeouts["microgen-repair"] = 600
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "whole_experiment_timeout_s": 3600, "stage_timeouts_s": timeouts,
        "platform": platform.platform(), "architecture": platform.machine(),
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY,
        "memory_enforcement": "process-tree-rss-watchdog" if sys.platform == "darwin" else "rlimit-as-and-process-tree-rss",
        "pilot": str(args.pilot), "quadrature": str(args.quadrature), "explicit_interior_repair": args.repair_interior,
        "repetitions": 5, "warmups_per_arm": 1, "schedule": schedule(),
        "stiffness_target": .001, "solve_node_cap": 250000, "tetrahedron_cap": workflow.CAP,
        "maximum_wall_minus_monotonic_s": 5,
        "timing_scope": "Both recorded costs retained; warmup excluded; alternate native-first and Microgen-first; no previous-run mesh, load or preconditioner replay"})
    parts = {phase + "/" + name: timeouts[name] for phase, native_first in schedule()
        for name in stage_names(native_first, args.repair_interior)}
    pairs = []
    with Job(parts=parts) as job:
        for phase, native_first in schedule():
            report = run_pair(args.run / phase, native_first=native_first, job=job, prefix=phase + "/",
                deadline=deadline, repair_interior=args.repair_interior)
            pairs.append({"phase": phase, "native_first": native_first, "source": str(args.run / phase), "report": report})
            workflow.save(args.run / "checkpoint.json", {"pairs": pairs, "completed_pairs": len(pairs),
                "elapsed_s": time.monotonic() - started})
            assert report["both_pass_declared_reference_models"], "Accuracy failure: read the completed pair before any retry"
            assert report["timing_clock_consistent"], "Clock interruption: preserve and exclude this timing experiment"
    report = summarize(pairs)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())
    report.update(inputs_unchanged=True, elapsed_s=time.monotonic() - started)
    workflow.save(args.run / "summary.json", report)


if __name__ == "__main__":
    main()
