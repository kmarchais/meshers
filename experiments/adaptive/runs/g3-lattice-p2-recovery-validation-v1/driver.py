"""Validate periodic Tet10 stress recovery against saved cubic enrichment energy."""

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


def recover_vertices(samples, volumes, corners, classes):
    import numpy as np

    assert samples.shape == (len(corners), 4, 3, 3)
    assert np.isfinite(samples).all() and np.isfinite(volumes).all() and (volumes > 0).all()
    keys = classes[corners]
    totals = np.zeros((len(classes), 3, 3))
    weights = np.zeros(len(classes))
    np.add.at(totals, keys.ravel(), (samples * volumes[:, None, None, None]).reshape(-1, 3, 3))
    np.add.at(weights, keys.ravel(), np.repeat(volumes, 4))
    assert (weights[keys] > 0).all()
    return totals[keys] / weights[keys, None, None]


def worker(args):
    np, fd, *_ = workflow.imports()
    from homogenization import strain_tensor
    from periodic import periodic_pairs
    from quadrature_audit import duffy_rule
    from scipy.stats import spearmanr
    from stress_recovery import stresses, gram
    from tet10_integration import audit_saved

    audit_saved(args.source, args.run / "parent-field-consistency.json")
    row = json.loads((args.source / "result.json").read_text())
    oracle_row = json.loads((args.oracle / "result.json").read_text())
    assert Path(oracle_row["source_parent"]).resolve() == args.source.resolve()
    assert row["order"] == 2 and row["qualified"] and row["elements"] <= workflow.CAP
    mesh = np.load(args.source / "mesh.npz")
    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    assert cells.shape[1] == 10 and len(points) <= 250000
    oracle = np.load(args.oracle / "positive125.npz")
    observed, volumes = oracle["energy"], oracle["volume"]
    assert observed.shape == volumes.shape == (len(cells),)
    assert (observed >= 0).all() and (volumes > 0).all() and observed.sum() > 0
    assert abs(observed.sum() / oracle_row["rules"][-1]["enrichment_energy"] - 1) < 1e-12
    _, classes, _ = periodic_pairs(points, surface)
    element = fd.lib_elements.element_list.get_element("tet10")(15)
    linear = fd.lib_elements.element_list.get_element("tet4")(4)
    strain = np.asarray(oracle_row["engineering_strain"])
    fields = np.einsum("l,lnc->nc", strain,
        np.load(args.source / "fields.npz")["displacement"])[None]
    started = time.perf_counter()
    samples = np.stack([stresses(points, cells, fields, derivative)[:, 0]
        for derivative in element.shape_function_derivative(element.xi_nd[:4])], axis=1)
    recovered = recover_vertices(samples, volumes, cells[:, :4], classes)
    affine_strain = strain_tensor(strain)
    affine_fields = (points @ affine_strain.T)[None]
    affine_samples = np.stack([stresses(points, cells, affine_fields, derivative)[:, 0]
        for derivative in element.shape_function_derivative(element.xi_nd[:4])], axis=1)
    expected = 2 * affine_strain + np.trace(affine_strain) * np.eye(3)
    affine_error = float(abs(affine_samples - expected).max())
    assert affine_error < 1e-10
    affine_recovered = recover_vertices(affine_samples, volumes, cells[:, :4], classes)
    assert abs(affine_recovered - expected).max() < 1e-10
    results = []
    previous = None
    for label, quadrature in (("own15", (element.xi_pg, element.w_pg)),
        ("positive64", duffy_rule(4)), ("positive125", duffy_rule(5))):
        locations, weights = quadrature
        estimated = np.zeros(len(cells))
        magnitude = np.zeros(len(cells))
        integrated_volume = np.zeros(len(cells))
        for location, derivative, weight in zip(locations,
            element.shape_function_derivative(locations), weights, strict=True):
            determinant = np.linalg.det(np.einsum("an,tnb->tab", derivative, points[cells]))
            assert determinant.min() > 0 and weight > 0
            stress = stresses(points, cells, fields, derivative)[:, 0]
            smooth = np.einsum("n,tnij->tij", linear.shape_function(location[None])[0], recovered)
            estimated += weight * determinant * gram((smooth - stress)[:, None])[:, 0, 0]
            magnitude += weight * determinant * gram(stress[:, None])[:, 0, 0]
            integrated_volume += weight * determinant
        assert estimated.min() >= 0 and estimated.sum() > 0
        np.testing.assert_allclose(integrated_volume, volumes, atol=1e-15, rtol=1e-10)
        count = int(np.ceil(.2 * len(cells)))
        captures = {}
        for name, indicator in (("recovery_density", estimated / volumes),
            ("recovery_mass", estimated), ("volume", volumes),
            ("stress_magnitude", magnitude), ("oracle", observed)):
            chosen = np.argsort(indicator)[-count:]
            captures[name] = float(observed[chosen].sum() / observed.sum())
        correlation = float(spearmanr(estimated / volumes, observed / volumes).statistic)
        results.append({"quadrature": label, "estimated_energy": float(estimated.sum()),
            "discrete_energy_effectivity_squared": float(estimated.sum() / observed.sum()),
            "density_spearman": correlation, "top_twenty_percent_enrichment_capture": captures,
            "relative_energy_change_from_previous_rule": None if previous is None else
                abs(float(estimated.sum() / previous - 1))})
        previous = float(estimated.sum())
        np.savez_compressed(args.run / (label + ".npz"), estimated=estimated,
            observed=observed, volumes=volumes, stress_magnitude=magnitude)
    measured = results[-1]
    captures = measured["top_twenty_percent_enrichment_capture"]
    passes = (measured["density_spearman"] >= .6 and captures["recovery_density"] >= .6
        and captures["recovery_density"] >= max(captures["volume"], captures["stress_magnitude"]) + .05
        and measured["relative_energy_change_from_previous_rule"] <= .01)
    workflow.save(args.run / "result.json", {"source": str(args.source), "oracle": str(args.oracle),
        "elements": len(cells), "nodes": len(points), "engineering_strain": strain.tolist(),
        "affine_patch_maximum_stress_error": affine_error, "recovery_and_diagnostic_s": time.perf_counter() - started,
        "rules": results, "predeclared_pilot_screen_passed": bool(passes),
        "scope": "Saved-field Tet10 recovery ranking against same-domain discrete cubic correction; fixed witness is a diagnostic input, not a free live estimator; no stiffness or local-stress certificate"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--oracle", required=True, type=Path)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), args.oracle / "manifest.json", args.oracle / "result.json",
        args.oracle / "positive125.npz"]
    paths.extend(args.source / name for name in ("result.json", "mesh.npz", "fields.npz"))
    paths.extend(workflow.HERE / name for name in ("stress_recovery.py", "tet10_integration.py",
        "quadrature_audit.py", "microgen_workflow.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    frozen = json.loads((args.oracle / "manifest.json").read_text())["inputs"]
    assert all(workflow.digest(Path(path)) == digest for path, digest in frozen.items())
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": 300, "threads": 1, "minimum_density_spearman": .6,
        "minimum_primary_capture": .6, "minimum_capture_gain_over_both_controls": .05,
        "maximum_final_quadrature_energy_change": .01})
    with Job(parts=["p2-recovery"]) as job:
        with job.part("p2-recovery") as part:
            record = workflow.bounded([sys.executable, __file__, *sys.argv[1:], "--worker"],
                args.run / "worker.log", part, timeout=300)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            result = json.loads((args.run / "result.json").read_text())
            job.result("p2-recovery", seconds=record["wall_s"],
                peak_gib=record["peak_rss_bytes"] / 1024**3,
                pilot_screen_passed=int(result["predeclared_pilot_screen_passed"]))
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
