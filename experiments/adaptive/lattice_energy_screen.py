"""Screen a non-oracle all-load energy indicator against full cubic enrichment."""

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


def worker(args):
    np, fd, *_ = workflow.imports()
    from cubic_reference import register_cubic
    from quadrature_audit import duffy_rule
    from scipy.linalg import eigvalsh
    from stress_recovery import stresses, gram
    from tet10_integration import audit_saved

    audit_saved(args.parent, args.run / "parent-field-consistency.json")
    parent = json.loads((args.parent / "result.json").read_text())
    cubic = json.loads((args.cubic / "result.json").read_text())
    assert parent["qualified"] and parent["geometry"]["geometry_admitted"]
    assert cubic["qualified"] and cubic["quadrature_sensitivity_passed"]
    assert Path(cubic["geometry_source"]).resolve() == args.parent.resolve()
    original = np.load(args.parent / "mesh.npz")
    enriched = np.load(args.prepared / "mesh.npz")
    parent_points, parent_cells = original["points"], original["tetrahedra"]
    points, cells = enriched["points"], enriched["tetrahedra"]
    parent_fields = np.load(args.parent / "fields.npz")["displacement"]
    cubic_fields = np.load(args.cubic / "fields.npz")["displacement"]
    assert len(cells) == len(parent_cells) <= workflow.CAP and len(points) <= 250000
    assert parent_fields.shape == (6, len(parent_points), 3)
    assert cubic_fields.shape == (6, len(points), 3)
    element = register_cubic()()
    quadratic = fd.lib_elements.element_list.get_element("tet10")(15)
    reconstructed = np.einsum("qn,tnd->tqd", quadratic.shape_function(element.xi_nd),
                              parent_points[parent_cells])
    np.testing.assert_allclose(points[cells], reconstructed, atol=1e-12, rtol=0)
    parent_tensor = np.asarray(parent["stiffness"])
    expected = parent_tensor - np.asarray(cubic["stiffness"])
    inverse = np.linalg.inv(parent_tensor)
    assert np.linalg.eigvalsh(expected).min() > 0
    rotation, _ = np.linalg.qr(np.random.default_rng(20261019).normal(size=(6, 6)))
    rotated_inverse = np.linalg.inv(rotation.T @ parent_tensor @ rotation)
    count = int(np.ceil(.2 * len(cells)))
    results, previous_mass, previous_error = [], None, None
    for label, quadrature in (("own15", (element.xi_pg, element.w_pg)),
                             ("positive64", duffy_rule(4)), ("positive125", duffy_rule(5))):
        locations, weights = quadrature
        errors = np.zeros((len(cells), 6, 6))
        magnitudes = np.zeros_like(errors)
        volumes = np.zeros(len(cells))
        for parent_derivative, derivative, weight in zip(
            quadratic.shape_function_derivative(locations),
            element.shape_function_derivative(locations), weights, strict=True):
            parent_jacobian = np.einsum("an,tnb->tab", parent_derivative, parent_points[parent_cells])
            jacobian = np.einsum("an,tnb->tab", derivative, points[cells])
            np.testing.assert_allclose(jacobian, parent_jacobian, atol=1e-12, rtol=1e-10)
            determinant = np.linalg.det(jacobian)
            assert determinant.min() > 0 and weight > 0
            parent_stress = stresses(parent_points, parent_cells, parent_fields, parent_derivative)
            cubic_stress = stresses(points, cells, cubic_fields, derivative)
            errors += (weight * determinant)[:, None, None] * gram(parent_stress - cubic_stress)
            magnitudes += (weight * determinant)[:, None, None] * gram(parent_stress)
            volumes += weight * determinant
        total_error = errors.sum(axis=0)
        identity_error = float(np.linalg.norm(total_error - expected) / np.linalg.norm(expected))
        if label == "own15":
            assert identity_error < 1e-6
            np.testing.assert_allclose(magnitudes.sum(axis=0), parent_tensor, atol=1e-9, rtol=1e-8)
        observed = np.einsum("ij,tji->t", inverse, errors)
        magnitude = np.einsum("ij,tji->t", inverse, magnitudes)
        assert observed.min() >= -1e-12 and magnitude.min() > 0 and volumes.min() > 0
        rotated = np.einsum("ij,tji->t", rotated_inverse,
                           np.einsum("ai,tab,bj->tij", rotation, magnitudes, rotation))
        invariance_error = float(np.linalg.norm(rotated - magnitude) / np.linalg.norm(magnitude))
        assert invariance_error < 1e-12
        captures = {}
        coverage = {}
        for name, indicator in (("all_load_stress_mass", magnitude),
                                ("all_load_stress_density", magnitude / volumes),
                                ("volume", volumes), ("oracle", observed)):
            selected = np.argsort(indicator)[-count:]
            captures[name] = float(observed[selected].sum() / observed.sum())
            selected_error = errors[selected].sum(axis=0)
            coverage[name] = float(eigvalsh(selected_error, total_error).min())
        results.append({"quadrature": label, "galerkin_matrix_relative_error": identity_error,
            "basis_invariance_relative_error": invariance_error,
            "total_whitened_enrichment": float(observed.sum()),
            "top_twenty_percent_trace_capture": captures,
            "minimum_mixed_load_capture": coverage,
            "stress_mass_relative_change_from_previous_rule": None if previous_mass is None else
                float(np.linalg.norm(magnitude - previous_mass, ord=1) / magnitude.sum()),
            "enrichment_relative_change_from_previous_rule": None if previous_error is None else
                abs(float(observed.sum() / previous_error - 1))})
        previous_mass, previous_error = magnitude, float(observed.sum())
        np.savez_compressed(args.run / (label + ".npz"), enrichment_gram=errors,
                            stress_gram=magnitudes, volume=volumes, stress_mass=magnitude,
                            whitened_enrichment=observed)
    measured = results[-1]
    captures = measured["top_twenty_percent_trace_capture"]
    passed = (captures["all_load_stress_mass"] >= .8
        and captures["all_load_stress_mass"] - captures["volume"] >= .2
        and measured["minimum_mixed_load_capture"]["all_load_stress_mass"] >= .5
        and measured["stress_mass_relative_change_from_previous_rule"] <= .01
        and measured["enrichment_relative_change_from_previous_rule"] <= .01)
    workflow.save(args.run / "result.json", {"parent": str(args.parent), "cubic": str(args.cubic),
        "elements": len(cells), "nodes": len(parent_points), "rules": results,
        "predeclared_energy_baseline_screen_passed": bool(passed),
        "scope": "All-six-load compliance-whitened stress-energy ranking versus a same-domain full cubic correction; no oracle strain in the indicator; no recovery promotion, stiffness gain, continuum bound or local-stress certificate"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--cubic", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), args.prepared / "mesh.npz"]
    paths.extend(args.parent / name for name in ("result.json", "mesh.npz", "fields.npz"))
    paths.extend(args.cubic / name for name in ("result.json", "fields.npz"))
    paths.extend(workflow.HERE / name for name in ("cubic_reference.py", "quadrature_audit.py",
                 "stress_recovery.py", "tet10_integration.py", "microgen_workflow.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 600, "minimum_trace_capture": .8,
        "minimum_gain_over_volume": .2, "minimum_capture_every_mixed_load": .5,
        "maximum_final_quadrature_change": .01, "maximum_own_rule_galerkin_matrix_error": 1e-6,
        "primary_indicator": "trace(C_parent^-1 * per_element_stress_energy_gram)",
        "oracle_role": "Validation only; not an input to the primary indicator"})
    with Job(parts=["all-load-energy-screen"]) as job:
        with job.part("all-load-energy-screen") as part:
            record = workflow.bounded([sys.executable, __file__, *sys.argv[1:], "--worker"],
                                      args.run / "worker.log", part, timeout=600)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("all-load-energy-screen", seconds=record["wall_s"],
                peak_gib=record["peak_rss_bytes"] / 1024**3,
                screen_passed=int(row["predeclared_energy_baseline_screen_passed"]))
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
