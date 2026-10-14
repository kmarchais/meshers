"""Read-only saved-field audits and reference comparisons from frozen sources."""

import json
import os
import subprocess
import sys
from pathlib import Path

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[name] = "1"
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report
from run import digest, write_json
from verify import verify_run
import compare_accuracy_references as comparison


def main():
    output = Path(__file__).resolve().parents[1]
    source = output / "source"
    manifest = json.loads((output / "manifest.json").read_text())
    report(0, "Verify frozen accuracy-analysis inputs", eta_s=60)
    for name, expected in manifest["source_sha256"].items():
        assert digest(source / name) == expected, name
    for name, expected in manifest["input_sha256"].items():
        assert digest(name) == expected, name
    env = os.environ.copy()
    env["PYTHONPATH"] = str(source)
    for label, interpreter, tests in [
        ("regression", sys.executable, sorted(source.glob("test_*.py"))),
        ("pardiso-smoke", manifest["pardiso_interpreter"], [source / "test_pardiso_solver.py"]),
    ]:
        report(.05, f"Run {label} checks", eta_s=60)
        checked = subprocess.run([str(interpreter), "-m", "pytest", "-q", *map(str, tests)],
                                 env=env, capture_output=True, text=True)
        (output / f"{label}.log").write_text(checked.stdout + checked.stderr, encoding="utf-8")
        print(checked.stdout + checked.stderr, flush=True)
        if checked.returncode:
            raise RuntimeError(f"{label} checks failed")
    revalidation = {"manifest_sha256": digest(output / "manifest.json")}
    if manifest.get("body_run"):
        body_run = Path(manifest["body_run"])
        audit_path = body_run / "artifact-verification.json"
        audit_hash = digest(audit_path)
        report(.10, "Revalidate six body fields without rewriting their audit", eta_s=45)
        audit = verify_run(body_run, write_report=False)
        assert digest(audit_path) == audit_hash
        body_row = json.loads((body_run / "results.json").read_text())["rows"][0]
        assert body_row["status"] == "ok"
        assert abs(body_row["solid_fraction"] - body_row["seed_volume"]) <= 1e-10
        assert max(body_row["nested_refinement_energy_change"]["relative_energy_difference_eigenvalues"]) <= 1e-8
        revalidation.update({
            "audit": audit, "unchanged_audit_sha256": audit_hash,
            "manifest_sha256": digest(output / "manifest.json"),
            "volume_change": body_row["solid_fraction"] - body_row["seed_volume"],
            "energy_change": body_row["nested_refinement_energy_change"],
        })
        if manifest.get("prototype_result"):
            # Cell order differs between implementations; test the same coordinates and
            # unordered tetrahedra before comparing their retained nodal solutions.
            import numpy as np
            from periodic import tensor_difference
            prototype_path = Path(manifest["prototype_result"])
            prototype_row = json.loads(prototype_path.read_text())
            native_path = body_run / body_row["artifact"] / "result.json"
            native_row = json.loads(native_path.read_text())
            with np.load(prototype_path.parent / "mesh.npz") as old, np.load(native_path.parent / "mesh.npz") as new:
                def cells(t):
                    t = np.sort(t, axis=1)
                    return t[np.lexsort(t.T[::-1])]
                identical_points = np.array_equal(old["points"], new["points"])
                identical_tetrahedra = identical_points and np.array_equal(cells(old["tetrahedra"]), cells(new["tetrahedra"]))
            agreement = {"exact_coordinates": identical_points, "same_unordered_tetrahedra": identical_tetrahedra,
                         "stiffness": tensor_difference(native_row["homogenization"]["stiffness_engineering_voigt"],
                                                         prototype_row["homogenization"]["stiffness_engineering_voigt"])}
            if identical_tetrahedra:
                differences = []
                for a, b in zip(native_row["homogenization"]["displacement_fields"], prototype_row["homogenization"]["displacement_fields"]):
                    with np.load(native_path.parent / a["artifact"]) as n, np.load(prototype_path.parent / b["artifact"]) as p:
                        assert np.array_equal(n["engineering"], p["engineering"])
                        delta = n["displacement"] - p["displacement"]
                        differences.append({"maximum_absolute": float(np.abs(delta).max()),
                                            "relative_l2": float(np.linalg.norm(delta) / np.linalg.norm(p["displacement"]))})
                agreement["six_displacement_field_differences"] = differences
                assert max(x["relative_l2"] for x in differences) < 1e-7
            write_json(output / "native-prototype-equivalence.json", agreement)
    if manifest.get("audit_runs"):
        audits = []
        for index, run_path in enumerate(manifest["audit_runs"]):
            run_path = Path(run_path)
            audit_path = run_path / "artifact-verification.json"
            audit_hash = digest(audit_path)
            report(.10, f"Read-only saved-field audit {index + 1}/{len(manifest['audit_runs'])}", eta_s=120)
            audit = verify_run(run_path, write_report=False)
            assert digest(audit_path) == audit_hash
            assert audit["verified_rows"] > 0
            audits.append({"run": str(run_path), "audit": audit,
                           "unchanged_audit_sha256": audit_hash})
        revalidation["audits"] = audits
    write_json(output / "readonly-revalidation.json", revalidation)
    comparison.report = lambda fraction, message, **kwargs: report(.20 + .80 * fraction, message, **kwargs)
    comparison.compare(list(map(Path, manifest["references"])), list(map(Path, manifest["candidates"])),
                       output / "comparisons", manifest["geometry_gate"], manifest["power"])
    for name, expected in manifest["input_sha256"].items():
        assert digest(name) == expected, name
    write_json(output / "completion.json", {"state": "complete", "manifest_sha256": digest(output / "manifest.json"),
                                           "input_hashes_unchanged": True, "pardiso_smoke_passed": True})
    report(1, "Body accuracy comparisons and optional direct-solver checks complete", eta_s=0)


if __name__ == "__main__":
    main()
