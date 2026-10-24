"""Independent constrained-basis replay and unchanged 24-model assessment."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as workflow
import numpy as np
from scipy.linalg import eigvalsh


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from cubic_refinement import MODEL_REPORTS

    report = json.loads((args.source / "result.json").read_text())
    inherited = json.loads((args.source / "analysis-inputs.json").read_text())
    assert all(workflow.digest(Path(path)) == digest for path, digest in inherited.items())
    paths = [Path(__file__), args.source / "result.json", args.source / "analysis-inputs.json"]
    paths.extend(workflow.HERE / "runs" / name for name in MODEL_REPORTS)
    inputs = {**inherited, **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    constraint = np.hstack((np.eye(6), np.eye(6)))
    gains, checks = [], []
    for row in report["rules"]:
        energy = np.asarray(row["energy_matrix"])
        np.testing.assert_allclose(energy, energy.T, atol=1e-12)
        assert eigvalsh(energy).min() > 0
        saddle = np.block([[energy, constraint.T], [constraint, np.zeros((6, 6))]])
        rhs = np.vstack((np.zeros((12, 6)), np.eye(6)))
        solved = np.linalg.solve(saddle, rhs)
        basis = solved[:12]
        np.testing.assert_allclose(constraint @ basis, np.eye(6), atol=1e-12)
        residual = float(np.linalg.norm(saddle @ solved - rhs))
        assert residual < 1e-11
        trial = basis.T @ energy @ basis
        np.testing.assert_allclose(trial, row["trial_energy_tensor"], atol=1e-12, rtol=1e-10)
        gain = energy[:6, :6] - trial
        np.testing.assert_allclose(gain, row["gain_tensor"], atol=1e-12)
        assert eigvalsh(gain).min() >= -1e-12
        assert eigvalsh(energy[6:, 6:] - trial).min() >= -1e-12
        maximum = float(eigvalsh(gain, energy[:6, :6]).max())
        np.testing.assert_allclose(maximum, row["maximum_relative_gain"], atol=1e-12)
        gains.append(gain)
        checks.append({"positive_rule_points": row["positive_rule_points"], "saddle_residual_norm": residual,
            "maximum_relative_gain": maximum})
    precision = float(np.linalg.norm(gains[-1] - gains[0]) / np.linalg.norm(gains[-1]))
    np.testing.assert_allclose(precision, report["gain_rule_relative_change"], atol=1e-10)
    trial = np.asarray(report["rules"][-1]["trial_energy_tensor"])
    models = []
    for name in MODEL_REPORTS:
        path = workflow.HERE / "runs" / name
        reference = json.loads(path.read_text())
        assert reference["case"] == "smooth-lattice" and reference["stiffness_target"] == .01
        for index, model in enumerate(reference["models"]):
            matrix = np.asarray(model["stiffness"])
            values = eigvalsh(trial - matrix, matrix)
            models.append({"source": str(path), "index": index, "signed_relative_energy_eigenvalues": values.tolist(),
                "worst_relative_difference": float(abs(values).max())})
    assert len(models) == 24
    worst = max(row["worst_relative_difference"] for row in models)
    workflow.save(args.output, {"passed": True, "source": str(args.source), "rule_checks": checks,
        "gain_rule_relative_change": precision, "precision_passed": precision <= .01,
        "signal_passed": checks[-1]["maximum_relative_gain"] >= .0005,
        "retained_model_count": len(models), "models": models,
        "trial_worst_retained_model_disagreement": worst, "trial_passes_target_all_models": worst <= .01,
        "reference_promoted": False, "native_candidate_promoted": False, "new_equilibrium_solve": False,
        "scope": "Independent 18-by-18 constrained-basis replay of fixed trial fields and all original model tensors; no new physical mesh, equilibrium reference or workflow timing claim."})
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())
    print(f"Constrained-basis replay passes; trial all-model disagreement {worst:.9g}")


if __name__ == "__main__":
    main()
