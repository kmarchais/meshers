"""Locate saved-field disagreement without changing meshes or accuracy gates."""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[name] = "1"
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report
import numpy as np
from scipy.linalg import eigh
from scipy.spatial import cKDTree
from scipy.stats import qmc
from accuracy_metrics import field_score, sample_fields
from homogenization import MATERIAL
from run import digest, write_json


def canonical_coordinates(xyz, quantize=False):
    values = np.rint(xyz * 1e12).astype(np.int64) if quantize else xyz
    order = np.lexsort((values[:, :, 2], values[:, :, 1], values[:, :, 0]), axis=1)
    return np.take_along_axis(values, order[:, :, None], axis=1)


def geometry(path, cells):
    with np.load(path.parent / "mesh.npz") as mesh:
        xyz = mesh["points"][mesh["tetrahedra"][cells]]
        count = len(mesh["tetrahedra"])
    det = np.linalg.det(xyz[:, 1:] - xyz[:, :1])
    edges = [(a, b) for a in range(4) for b in range(a + 1, 4)]
    edge_sum = sum(np.sum((xyz[:, a] - xyz[:, b])**2, axis=1) for a, b in edges)
    return xyz, det / 6, np.sqrt(432 * det**2 / edge_sum**3), count


def main():
    output = Path(__file__).resolve().parents[1]
    manifest = json.loads((output / "manifest.json").read_text())
    report(0, "Verify frozen localization inputs", eta_s=120)
    for name, expected in manifest["source_sha256"].items():
        assert digest(output / "source" / name) == expected, name
    for name, expected in manifest["input_sha256"].items():
        assert digest(name) == expected, name
    tests = subprocess.run([sys.executable, "-m", "pytest", "-q", *map(str, sorted((output / "source").glob("test_*.py")))], capture_output=True, text=True)
    (output / "tests.log").write_text(tests.stdout + tests.stderr)
    if tests.returncode:
        raise RuntimeError("Localization regressions failed")
    prior = json.loads(Path(manifest["comparison"]).read_text())
    seeds = [20261013 + i for i in range(4)]
    probes = np.vstack([qmc.Sobol(3, scramble=True, seed=s).random_base2(manifest["power"] - 2) for s in seeds]) - .5
    compliance = np.linalg.inv(MATERIAL)
    document = {"state": "running", "references": [], "scope": "Common-probe stress-error localization in the native seed worst mixed load; wall distance uses periodic reference triangle centroids, not exact surface distance; tetrahedron sharing is a topology diagnostic, not an accuracy gate", "manifest_sha256": digest(output / "manifest.json")}
    started = time.monotonic()
    reference_count, candidate_count = len(manifest["references"]), len(manifest["candidates"])
    for ri, reference_name in enumerate(manifest["references"]):
        reference_path = Path(reference_name)
        report((ri + .02) / reference_count, "Load reference fields for localization", eta_s=120,
               probes=len(probes))
        reference = sample_fields(reference_path, probes, manifest["geometry_gate"], include_cells=True)
        sampled = []
        common = reference[4].copy()
        for ci, candidate_name in enumerate(manifest["candidates"]):
            report((ri + .05 + .40 * ci / candidate_count) / reference_count,
                   f"Locate candidate {ci + 1}/{candidate_count}", eta_s=120, probes=len(probes))
            candidate = sample_fields(Path(candidate_name), probes, manifest["geometry_gate"], include_cells=True)
            common &= candidate[4]
            sampled.append(candidate)
        reference_stress = reference[3][common]
        anchor = field_score(sampled[0][3][common], reference_stress)
        a, b = np.asarray(anchor["stress_error_energy_gram"]), np.asarray(anchor["reference_energy_gram"])
        _, vectors = eigh((a + a.T) / 2, (b + b.T) / 2)
        load = vectors[:, -1]
        normalization = float(load @ b @ load)
        assert abs(normalization - 1) < 1e-10
        with np.load(reference_path.parent / "mesh.npz") as mesh:
            triangles = mesh["points"][mesh["surface"]]
        cap = np.any(np.all(np.isclose(triangles, -.5, atol=1e-10, rtol=0), axis=1), axis=1)
        cap |= np.any(np.all(np.isclose(triangles, .5, atol=1e-10, rtol=0), axis=1), axis=1)
        tree = cKDTree((triangles[~cap].mean(axis=1) + .5) % 1, boxsize=1)
        distance, _ = tree.query((probes[common] + .5) % 1, workers=1)
        ref_xyz, ref_volume, _, _ = geometry(reference_path, reference[5][common])
        ref_exact, ref_quantized = canonical_coordinates(ref_xyz), canonical_coordinates(ref_xyz, True)
        result = {"reference": reference[2], "common_probes": int(common.sum()), "reference_material_coverage": float(common.sum() / reference[4].sum()), "anchor_engineering_strain": load.tolist(), "anchor_reference_energy": normalization, "rows": []}
        for ci, (candidate_name, candidate) in enumerate(zip(manifest["candidates"], sampled)):
            candidate_path = Path(candidate_name)
            score = field_score(candidate[3][common], reference_stress)
            report((ri + .50 + .49 * (ci + 1) / candidate_count) / reference_count,
                   f"Localize error {ci + 1}/{candidate_count}", eta_s=120,
                   coverage={"value": 100 * common.sum() / reference[4].sum(), "unit": "%"},
                   stress={"value": 100 * score["worst_relative_stress_energy_norm_difference"], "unit": "%"})
            expected = prior["comparisons"][ri]["rows"][ci]
            assert candidate_name == expected["path"]
            np.testing.assert_allclose(score["stress_error_energy_gram"], expected["stress_error_energy_gram"], rtol=1e-11, atol=1e-13)
            assert abs(score["worst_relative_stress_energy_norm_difference"] - expected["worst_relative_stress_energy_norm_difference"]) < 1e-11
            difference = (candidate[3][common] - reference_stress) @ load
            energy = np.einsum("ni,ij,nj->n", difference, compliance, difference)
            assert energy.min() >= -1e-15
            anchor_error = float(load @ np.asarray(score["stress_error_energy_gram"]) @ load)
            assert np.isclose(energy.mean(), anchor_error, rtol=1e-11, atol=1e-14)
            cells = candidate[5][common]
            xyz, volume, quality, count = geometry(candidate_path, cells)
            exact = np.all(canonical_coordinates(xyz) == ref_exact, axis=(1, 2))
            shared = np.all(canonical_coordinates(xyz, True) == ref_quantized, axis=(1, 2))
            total = float(energy.sum())
            def summary(mask):
                return {"probes": int(mask.sum()), "material_probe_share": float(mask.mean()), "anchor_error_share": float(energy[mask].sum() / total), "anchor_error_density": float(energy[mask].mean()) if mask.any() else None, "median_tet_volume": float(np.median(volume[mask])) if mask.any() else None, "median_reference_tet_volume": float(np.median(ref_volume[mask])) if mask.any() else None, "median_quality": float(np.median(quality[mask])) if mask.any() else None}
            edges = [0, .01, .02, .04, .08, float("inf")]
            bins = [{"lower_distance": lo, "upper_distance": hi if np.isfinite(hi) else None, **summary((distance >= lo) & (distance < hi))} for lo, hi in zip(edges[:-1], edges[1:])]
            np.testing.assert_allclose(sum(x["anchor_error_share"] for x in bins), 1, rtol=1e-12)
            error_by_cell = np.bincount(cells, weights=energy, minlength=count)
            probes_by_cell = np.bincount(cells, minlength=count)
            assert np.isclose(error_by_cell.sum(), total, rtol=1e-12)
            path = output / f"reference-{ri}-candidate-{ci}-cell-error.npz"
            np.savez_compressed(path, error_sum=error_by_cell, probe_count=probes_by_cell, engineering=load)
            result["rows"].append({"path": candidate_name, "worst_local_score": score["worst_relative_stress_energy_norm_difference"], "anchor_local_score": float(np.sqrt(anchor_error)), "anchor_error_energy": anchor_error, "wall_distance_bins": bins, "shared_reference_tet_exact": summary(exact), "shared_reference_tet_quantized_1e12": summary(shared), "nonshared_reference_tet_quantized_1e12": summary(~shared), "cell_error_artifact": path.name, "cell_error_sha256": digest(path)})
        document["references"].append(result)
        write_json(output / "localization.json", document)
    for name, expected in manifest["input_sha256"].items():
        assert digest(name) == expected, name
    for name, expected in manifest["source_sha256"].items():
        assert digest(output / "source" / name) == expected, name
    document.update(state="complete", inputs_unchanged=True, elapsed_s=time.monotonic() - started, checks=["Global Gram matrices and worst scores reproduce the retained candidate report", "Per-probe energy integrates to the anchor mixed-load Gram quadratic form", "Wall bins and per-cell localization conserve total sampled error"])
    write_json(output / "localization.json", document)
    report(1, "Stress-error localization complete; archived inputs unchanged", eta_s=0)


if __name__ == "__main__":
    main()
