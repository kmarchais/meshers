"""Test a declared leading-error model on three independently audited cubic levels."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

LEVELS = (
    (23, "g3-lattice-cubic23-reference-v1", "g3-lattice-cubic23-assessment-v1"),
    (24, "g3-lattice-cubic-reference-cache-release-v1", "g3-lattice-cubic-assessment-v1"),
    (25, "g3-lattice-cubic25-reference-blocked-symmetry-v2", "g3-lattice-cubic25-assessment-v1"),
)
MODEL_REPORTS = (
    "g3-lattice-curve-grading-reference-comparison-v1.json",
    "g3-lattice-reference-sensitivity-v2.json",
    "g3-lattice-microgen-reference-comparison-v1.json",
)


def fit_refinement(tensors, lengths, quadrature_shift):
    import numpy as np
    from scipy.optimize import brentq

    tensors = np.asarray(tensors)
    lengths = np.asarray(lengths)
    assert tensors.shape == (3, 6, 6) and np.isfinite(tensors).all()
    assert lengths.shape == (3,) and np.isfinite(lengths).all() and (lengths > 0).all()
    assert (np.diff(lengths) < 0).all() and np.isfinite(quadrature_shift) and quadrature_shift >= 0
    for tensor in tensors:
        np.testing.assert_allclose(tensor, tensor.T, atol=1e-12, rtol=0)
        assert np.linalg.eigvalsh(tensor).min() > 0
    factor = np.linalg.cholesky(tensors[-1])
    differences = [tensors[0] - tensors[1], tensors[1] - tensors[2]]
    whitened = [np.linalg.solve(factor, np.linalg.solve(factor, difference).T).T
                for difference in differences]
    norms = [float(np.linalg.norm(difference)) for difference in whitened]
    spectra = [np.linalg.eigvalsh(difference) for difference in whitened]
    signal = min(float(abs(spectrum).max()) for spectrum in spectra)
    result = {"lengths": lengths.tolist(), "quadrature_shift": quadrature_shift,
        "correction_frobenius_norms": norms, "correction_spectra": [spectrum.tolist() for spectrum in spectra],
        "minimum_mixed_load_correction": signal,
        "signal_above_ten_times_quadrature": signal >= 10 * quadrature_shift,
        "declared_order_interval": [.25, 6.], "matrix_model_relative_residual_limit": .1,
        "compatible": False, "observed_order": None, "extrapolated_stiffness": None}
    if min(norms) <= 1e-30:
        result["reason"] = "Insufficient correction signal"
        return result
    coefficient = float(np.sum(whitened[0] * whitened[1]) / norms[1]**2)
    residual = float(np.linalg.norm(whitened[0] - coefficient * whitened[1]) / norms[0])
    cosine = float(np.sum(whitened[0] * whitened[1]) / (norms[0] * norms[1]))
    scaled = lengths / lengths[-1]
    def correction_ratio(order):
        return float((scaled[0]**order - scaled[1]**order) / (scaled[1]**order - 1))
    bounds = [correction_ratio(order) for order in (.25, 6.)]
    result.update(best_correction_ratio=coefficient, correction_direction_cosine=cosine,
                  matrix_model_relative_residual=residual, positive_order_ratio_interval=bounds)
    reasons = []
    if coefficient <= 0:
        reasons.append("Correction directions reverse")
    if residual > .1:
        reasons.append("Corrections are not proportional under one leading matrix model")
    if not min(bounds) <= coefficient <= max(bounds):
        reasons.append("Correction ratio is outside the declared positive-order interval")
    else:
        result["observed_order"] = float(brentq(lambda order: correction_ratio(order) - coefficient, .25, 6.))
    if not result["signal_above_ten_times_quadrature"]:
        reasons.append("Correction is not resolved above the empirical quadrature-sensitivity screen")
    if not reasons:
        order = result["observed_order"]
        extrapolated = tensors[-1] - differences[-1] / (scaled[1]**order - 1)
        if np.linalg.eigvalsh(extrapolated).min() <= 0:
            reasons.append("Extrapolated tensor is not positive definite")
        else:
            result.update(compatible=True, extrapolated_stiffness=extrapolated.tolist())
    result["reasons"] = reasons
    return result


def analyze(run):
    import numpy as np

    canonical = []
    for name in MODEL_REPORTS:
        report = json.loads((workflow.HERE / "runs" / name).read_text())
        assert report["case"] == "smooth-lattice" and len(report["models"]) == 8
        canonical.extend(np.asarray(model["stiffness"], dtype=np.float64).tobytes() for model in report["models"])
    tensors, cells, elements, sensitivity, records = [], [], [], [], []
    for level, source_name, audit_name in LEVELS:
        source, audit = workflow.HERE / "runs" / source_name, workflow.HERE / "runs" / audit_name
        row = json.loads((source / "result.json").read_text())
        checked = json.loads((audit / "result.json").read_text())
        assert row["qualified"] and row["order"] == 3 and row["completed_loads"] == 6
        assert row["quadrature_sensitivity_passed"] and row["quadrature_sensitivity_limit"] <= 1e-4
        assert max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
        assert checked["fresh_process_audit_passed"] and checked["retained_model_count"] == 24
        assert Path(checked["source"]).resolve() == source.resolve()
        frozen = {Path(path).resolve(): value for path, value in
                  json.loads((audit / "manifest.json").read_text())["inputs"].items()}
        for path in (source / "result.json", source / "fields.npz"):
            assert workflow.digest(path) == frozen[path.resolve()]
        actual = []
        for model in checked["models"]:
            path = Path(model["source"])
            path = path if path.is_absolute() else workflow.ROOT / path
            assert workflow.digest(path) == frozen[path.resolve()]
            matrix = json.loads(path.read_text())["models"][model["index"]]["stiffness"]
            actual.append(np.asarray(matrix, dtype=np.float64).tobytes())
        assert sorted(actual) == sorted(canonical), "Audited model tensors differ from the retained 24-model set"
        parent = Path(row["geometry_source"])
        parent = parent if parent.is_absolute() else workflow.ROOT / parent
        parent_row = json.loads((parent / "result.json").read_text())
        provenance = {Path(path).resolve(): value for path, value in
                      json.loads((source / "manifest.json").read_text())["sources"].items()}
        assert workflow.digest(parent / "result.json") == provenance[(parent / "result.json").resolve()]
        assert row["case"] == parent_row["case"] == "smooth-lattice"
        assert parent_row["cells"] == level and parent_row["adaptive"] is False
        assert parent_row["geometry"]["geometry_admitted"]
        assert parent_row["void_clamp"]["positive_ceiling"] == .02 and parent_row["symmetric_union"]
        tensors.append(row["stiffness"])
        cells.append(level)
        elements.append(row["elements"])
        sensitivity.append(row["higher_rule_fixed_field_difference"]["worst_relative_strain_energy_difference"])
        records.append({"cells": level, "source": str(source), "elements": row["elements"],
                        "nodes": row["nodes"], "quadrature_shift": sensitivity[-1]})
    shift = max(sensitivity)
    screens = {"seed_spacing": fit_refinement(tensors, 1 / np.asarray(cells), shift),
               "element_density": fit_refinement(tensors, np.asarray(elements, dtype=float)**(-1 / 3), shift)}
    workflow.save(run / "result.json", {"records": records, "screens": screens,
        "retained_model_count": 24, "retained_models_unchanged": True,
        "scope": "Compatibility with one assumed common leading-error matrix, not a continuum certificate or a replacement for existing benchmark admission. Nonnested domains and spaces may violate this model."})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        analyze(args.run)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), workflow.HERE / "test_cubic_refinement.py", Path(workflow.__file__)]
    paths.extend(workflow.HERE / "runs" / name for name in MODEL_REPORTS)
    for _, source_name, audit_name in LEVELS:
        source, audit = workflow.HERE / "runs" / source_name, workflow.HERE / "runs" / audit_name
        paths.extend((source / "result.json", source / "fields.npz", source / "manifest.json",
                      audit / "result.json", audit / "manifest.json"))
        for model in json.loads((audit / "result.json").read_text())["models"]:
            path = Path(model["source"])
            paths.append(path if path.is_absolute() else workflow.ROOT / path)
        parent = Path(json.loads((source / "result.json").read_text())["geometry_source"])
        parent = parent if parent.is_absolute() else workflow.ROOT / parent
        paths.append(parent / "result.json")
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "threads": 1, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": 120, "leading_matrix_model": "C(h) = C_limit + A*h**p",
        "order_interval": [.25, 6.], "maximum_relative_matrix_residual": .1,
        "minimum_correction_to_quadrature_ratio": 10, "no_reference_admission_change": True})
    with Job(parts=["synthetic-controls", "three-level-screen"]) as job:
        for stage in ("synthetic-controls", "three-level-screen"):
            with job.part(stage) as part:
                command = ([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_cubic_refinement.py"]
                    if stage == "synthetic-controls" else [sys.executable, __file__, "--run", str(args.run), "--worker"])
                record = workflow.bounded(command, args.run / (stage + ".log"), part, timeout=120)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
