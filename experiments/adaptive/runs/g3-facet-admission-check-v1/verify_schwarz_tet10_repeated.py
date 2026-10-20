"""Independently recompute the completed repeated study from frozen records."""

import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess

import numpy as np
from scipy.linalg import eigvalsh

import microgen_workflow as workflow


def model_errors(candidate, coarse, fine, ratios):
    errors = []
    for ratio in ratios:
        assert 0 < ratio < 1
        for order in range(1, 5):
            model = fine + (fine - coarse) / (ratio**(-order) - 1)
            assert eigvalsh(model).min() > 0
            errors.append(float(abs(eigvalsh(candidate, model) - 1).max()))
    return errors


def verify_hash(path, expected, run, source_ref=None):
    if workflow.digest(path) == expected:
        return False
    archives = {"schwarz_tet10_live.py": run / "warmup/driver.py",
        "schwarz_tet10_repeated.py": run / "driver.py"}
    archived = archives.get(path.name)
    if archived is not None and workflow.digest(archived) == expected:
        return True
    if source_ref and path.suffix == ".py" and path.resolve().is_relative_to(workflow.ROOT):
        relative = path.resolve().relative_to(workflow.ROOT).as_posix()
        content = subprocess.check_output(["git", "show", f"{source_ref}:{relative}"], cwd=workflow.ROOT)
        assert hashlib.sha256(content).hexdigest() == expected, f"Measured Git source differs: {path}"
        return True
    raise AssertionError(f"Hash mismatch: {path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-ref")
    args = parser.parse_args()
    assert not args.output.exists()
    report = json.loads((args.run / "summary.json").read_text())
    assert report["repetitions"] == 5 and report["warmups_per_arm"] == 1
    assert report["stiffness_target"] == .001 and report["all_pairs_qualified"]
    source = workflow.HERE / "runs/g3-microgen-defaults-v4/summary.json"
    baseline = json.loads(source.read_text())
    levels = {row["cells"]: row for row in baseline["rows"]
        if row["case"] == "schwarz-solid" and row["route"] == "meshers-curved"}
    coarse, fine = (np.asarray(levels[count]["stiffness"]) for count in (16, 20))
    ratios = [16 / 20, (levels[16]["elements"] / levels[20]["elements"])**(1 / 3)]
    manifest = json.loads((args.run / "manifest.json").read_text())
    verified_hashes, archived_sources = 0, 0
    for path, digest in manifest["inputs"].items():
        archived_sources += verify_hash(Path(path), digest, args.run, args.source_ref)
        verified_hashes += 1
    times = {metric: {arm: [] for arm in ("microgen", "meshers")}
        for metric in ("science_audited_workflow_s", "generation_lift_solve_s")}
    checked = []
    maximum_reaction, maximum_energy, maximum_jump = 0., 0., 0.
    for index, pair in enumerate(report["pairs"]):
        phase = "warmup" if index == 0 else f"repeat-{index}"
        assert pair["phase"] == phase and pair["native_first"] == bool(index % 2)
        directory = args.run / phase
        saved = json.loads((directory / "summary.json").read_text())
        assert saved == pair["report"] and saved["timing_clock_consistent"]
        assert abs(saved["wall_minus_monotonic_s"]) <= 5
        inputs = json.loads((directory / "manifest.json").read_text())["inputs"]
        assert inputs[str(source)] == workflow.digest(source)
        for path, digest in {**inputs, **saved["artifact_sha256"]}.items():
            archived_sources += verify_hash(Path(path), digest, args.run, args.source_ref)
            verified_hashes += 1
        parent = json.loads((directory / "microgen-parent/result.json").read_text())
        rows = {arm: json.loads((directory / (arm + "-tet10/result.json")).read_text())
            for arm in ("microgen", "meshers")}
        ranges = {}
        for arm, row in rows.items():
            assert row["qualified"] and row["geometry"]["geometry_admitted"]
            geometry = row["geometry"]
            assert geometry["sampled_two_sided_max"] <= .002 and abs(geometry["relative_volume_error"]) <= .01
            assert geometry["q_min_corner"] >= .1 and geometry["minimum_sampled_jacobian"] > 0
            assert geometry["projection_failed"] == 0 and "periodic" in row
            assert row["nodes"] <= 250000 and row["elements"] <= workflow.CAP
            assert len(row["solver_records"]) == 6
            assert max(load["relative_residual"] for load in row["solver_records"]) < 2e-9
            candidate = np.asarray(row["stiffness"])
            np.testing.assert_allclose(candidate, candidate.T, atol=1e-10, rtol=0)
            errors = model_errors((candidate + candidate.T) / 2, coarse, fine, ratios)
            assert len(errors) == 8 and max(errors) <= .001
            comparison = json.loads((directory / (arm + "-reference.json")).read_text())
            np.testing.assert_allclose([min(errors), max(errors)], comparison["declared_order_difference_range"], atol=1e-12, rtol=0)
            ranges[arm] = [min(errors), max(errors)]
            audit = json.loads((directory / "audits" / (arm + "-tet10.json")).read_text())
            assert max(audit["reaction_relative_error"], audit["hill_mandel_relative_error"]) < 1e-7
            assert audit["periodic_jump_max"] < 1e-9
            maximum_reaction = max(maximum_reaction, audit["reaction_relative_error"])
            maximum_energy = max(maximum_energy, audit["hill_mandel_relative_error"])
            maximum_jump = max(maximum_jump, audit["periodic_jump_max"])
        totals = {"microgen": parent["total_workflow_s"] + rows["microgen"]["total_saved_mesh_study_s"],
            "meshers": rows["meshers"]["total_workflow_s"]}
        phases = {"microgen": parent["generation_s"] + rows["microgen"]["lift_and_conversion_s"] + rows["microgen"]["solve_s"],
            "meshers": rows["meshers"]["generation_s"] + rows["meshers"]["solve_s"]}
        for metric, values in (("science_audited_workflow_s", totals), ("generation_lift_solve_s", phases)):
            assert values == saved[metric]
            if index:
                for arm, value in values.items():
                    times[metric][arm].append(value)
        checked.append({"phase": phase, "model_difference_ranges": ranges})
    assert len(checked) == 6
    verified_faster = True
    for metric, arms in times.items():
        for arm, values in arms.items():
            expected = {"median": statistics.median(values), "minimum": min(values), "maximum": max(values)}
            assert expected == report["statistics"][metric][arm]
        paired = np.asarray(arms["microgen"]) / np.asarray(arms["meshers"])
        np.testing.assert_allclose(paired, report["paired_ratios_microgen_over_meshers"][metric], atol=1e-12, rtol=0)
        verified_faster = verified_faster and bool((paired > 1).all())
    assert verified_faster == report["meshers_faster_in_every_measured_pair"]
    workflow.save(args.output, {"source": str(args.run), "source_summary_sha256": workflow.digest(args.run / "summary.json"),
        "verifier_sha256": workflow.digest(Path(__file__)), "passed": True, "verified_hash_records": verified_hashes,
        "archived_source_hash_records": archived_sources,
        "measurement_source_ref": args.source_ref,
        "meshers_faster_in_every_measured_pair": verified_faster,
        "independent_generalized_eigenvalue_comparisons": 6 * 2 * 8, "pairs": checked,
        "maximum_reaction_relative_error": maximum_reaction, "maximum_hill_mandel_relative_error": maximum_energy,
        "maximum_periodic_jump": maximum_jump,
        "scope": "Independent frozen-record tensor, timing and provenance checks; not re-solving fields or continuum certification"})


if __name__ == "__main__":
    main()
