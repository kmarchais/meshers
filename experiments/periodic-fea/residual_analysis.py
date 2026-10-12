"""Locate traction residuals in saved P1 fields; this is not an error bound.

Use max of six energy-normalized unit-load indicators. Periodic cap faces are
interior faces in the quotient mesh; physical surfaces have zero applied traction.
The constitutive law is the retained experiment's lambda = mu = 1.
"""

import argparse
import json
import os
import sys
import threading
import time
from pathlib import Path

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[name] = "1"

import numpy as np

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report
from run import digest, write_json


def faces(points, tets, classes):
    """Pair volume faces by periodic vertex classes, not by floating coordinates."""
    oriented = tets[:, [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]]].reshape(-1, 3)
    keys = np.sort(classes[oriented], axis=1)
    assert np.all(np.diff(keys, axis=1) > 0), "A face collapses in the periodic quotient"
    _, first, inverse, counts = np.unique(keys, axis=0, return_index=True,
                                          return_inverse=True, return_counts=True)
    assert counts.max() <= 2, "Non-manifold quotient faces"
    last = np.zeros(len(first), dtype=np.int64)
    np.maximum.at(last, inverse, np.arange(len(oriented)))
    left, right = first // 4, last // 4
    interior = counts == 2
    xyz = points[oriented[first]]
    normal = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    double_area = np.linalg.norm(normal, axis=1)
    assert np.all(double_area > 0)
    normal /= double_area[:, None]
    area = double_area / 2
    diameter = np.maximum.reduce([
        np.linalg.norm(xyz[:, i] - xyz[:, j], axis=1) for i, j in [(0, 1), (1, 2), (2, 0)]
    ])
    other = points[oriented[last]]
    other_normal = np.cross(other[:, 1] - other[:, 0], other[:, 2] - other[:, 0])
    other_normal /= np.linalg.norm(other_normal, axis=1)[:, None]
    if interior.any():
        assert np.max(np.linalg.norm(normal[interior] + other_normal[interior], axis=1)) < 1e-9
    boundary = np.zeros(len(tets), dtype=bool)
    boundary[left[~interior]] = True
    return left, right, interior, normal, diameter * area / 2, boundary


def stress_and_strain(displacement, tets, gradients):
    gradient = np.einsum("tia,tib->tab", displacement[tets], gradients)
    strain = (gradient + gradient.transpose(0, 2, 1)) / 2
    stress = 2 * strain + np.trace(strain, axis1=1, axis2=2)[:, None, None] * np.eye(3)
    return stress, strain


def indicator(stress, topology, count):
    left, right, interior, normal, weight, _ = topology
    jump = np.einsum("fij,fj->fi", stress[left], normal)
    # Sum outward tractions: sigma_L n + sigma_R (-n) = (sigma_L-sigma_R) n.
    jump[interior] -= np.einsum("fij,fj->fi", stress[right[interior]], normal[interior])
    contribution = weight * np.einsum("fi,fi->f", jump, jump)
    local = np.bincount(left, weights=contribution * np.where(interior, .5, 1), minlength=count)
    local += np.bincount(right[interior], weights=contribution[interior] / 2, minlength=count)
    assert np.isclose(local.sum(), contribution.sum(), rtol=1e-12, atol=1e-25)
    return local


def concentration(score, volume, boundary):
    total = score.sum()
    if total < 1e-20:
        return {"indicator_sum": float(total), "effectively_zero": True}
    ordered = np.argsort(-score)
    top = ordered[:max(1, int(np.ceil(.2 * len(score))))]
    cumulative = np.cumsum(score[ordered]) / total
    return {
        "indicator_sum": float(total),
        "top20_cells_indicator_share": float(score[top].sum() / total),
        "top20_cells_volume_share": float(volume[top].sum() / volume.sum()),
        "cells_fraction_for80_indicator": float((np.searchsorted(cumulative, .8) + 1) / len(score)),
        "boundary_adjacent_cells_fraction": float(boundary.mean()),
        "boundary_adjacent_indicator_share": float(score[boundary].sum() / total),
    }


def analyze(output, artifacts):
    started = time.monotonic()
    live = {"done": 0, "message": "Verify residual-analysis sources"}
    stop = threading.Event()

    def heartbeat():
        while not stop.is_set():
            fraction = live["done"] / (7 * len(artifacts))
            elapsed = time.monotonic() - started
            eta = elapsed * (1 - fraction) / fraction if fraction else 90
            report(fraction, live["message"], eta_s=eta, seconds=round(elapsed))
            stop.wait(2)

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    try:
        manifest = json.loads((output / "manifest.json").read_text())
        for relative, expected in manifest["source_sha256"].items():
            assert digest(output / "source" / relative) == expected, relative
        document = {"state": "running", "rows": [], "manifest_sha256": digest(output / "manifest.json"),
                    "scope": "Max of six normalized unit-load traction-jump indicators on the meshed domain; no geometry-error or mixed-strain error bound"}
        isolated = faces(np.vstack([np.zeros(3), np.eye(3)]), np.array([[0, 1, 2, 3]]), np.arange(4))
        assert indicator(np.eye(3)[None], isolated, 1).sum() > 0, "Free-surface traction omitted"
        document["checks"] = ["Nonzero constant stress on an isolated free tetrahedron has a residual"]
        for index, directory in enumerate(artifacts):
            live["message"] = f"Build periodic quotient faces: {directory.name}"
            row = json.loads((directory / "result.json").read_text())
            if index == 0:
                assert row["case"] == "solid-cell", "Analytic patch control must run first"
            assert all(row.get(k) for k in ["homogenization_audited", "geometry_eligible", "quality_eligible", "periodic_eligible"])
            assert digest(directory / "mesh.npz") == row["mesh_sha256"]
            mesh = np.load(directory / "mesh.npz")
            points, tets = mesh["points"], mesh["tetrahedra"]
            pairing = np.load(directory / "periodic.npz")
            topology = faces(points, tets, pairing["classes"])
            jacobian = points[tets[:, 1:]] - points[tets[:, :1]]
            volume = np.linalg.det(jacobian) / 6
            assert volume.min() > 0
            gradients = np.empty((len(tets), 4, 3))
            gradients[:, 1:] = np.linalg.inv(jacobian).transpose(0, 2, 1)
            gradients[:, 0] = -gradients[:, 1:].sum(axis=1)
            live["done"] += 1
            score = np.zeros(len(tets))
            columns, load_summaries = [], []
            fields = row["homogenization"]["displacement_fields"]
            assert len(fields) == 6
            for load, field_record in enumerate(fields):
                live["message"] = f"Traction residual load {load + 1}/6: {directory.name}"
                path = directory / field_record["artifact"]
                assert digest(path) == field_record["sha256"]
                field = np.load(path)
                assert np.array_equal(field["engineering"], np.eye(6)[load])
                displacement = field["displacement"]
                stress, strain = stress_and_strain(displacement, tets, gradients)
                # a(u,u) = integral sigma:epsilon = twice the physical strain energy.
                energy_norm_squared = np.einsum("t,tij,tij", volume, stress, strain)
                assert energy_norm_squared > 0
                local = indicator(stress, topology, len(tets)) / energy_norm_squared
                score = np.maximum(score, local)
                load_summaries.append(concentration(local, volume, topology[-1]))
                average = np.einsum("t,tij->ij", volume, stress)
                columns.append(average[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]])
                if row["case"] == "solid-cell":
                    assert local.sum() < 1e-20, "Affine periodic solid has nonzero traction residual"
                    if load == 0:
                        perturbed = displacement.copy()
                        perturbed[:, 0] += .01 * np.sin(2 * np.pi * points[:, 1]) * np.sin(2 * np.pi * points[:, 2])
                        probe, _ = stress_and_strain(perturbed, tets, gradients)
                        assert indicator(probe, topology, len(tets)).sum() > 1e-6
                        document["checks"].append("Periodic affine solid has zero residual; a periodic non-affine perturbation is detected")
                live["done"] += 1
            reported = np.asarray(row["homogenization"]["stiffness_engineering_voigt"])
            assert np.linalg.norm(np.asarray(columns).T - reported) / np.linalg.norm(reported) < 1e-10
            target = output / f"indicators-{index}.npz"
            np.savez_compressed(target, score=score, boundary_adjacent=topology[-1])
            summary = {"artifact": str(directory), "artifact_sha256": digest(directory / "result.json"),
                       "mesh_sha256": row["mesh_sha256"], "indicator_sha256": digest(target),
                       "case": row["case"], "tetrahedra": len(tets), "loads": load_summaries,
                       **concentration(score, volume, topology[-1])}
            document["rows"].append(summary)
            write_json(output / "residual-analysis.json", document)
            print(json.dumps({k: v for k, v in summary.items() if k != "loads"}), flush=True)
        assert artifacts and document["rows"][0]["case"] == "solid-cell", "Analytic patch control must run first"
        document.update(state="complete", seconds=time.monotonic() - started)
        write_json(output / "residual-analysis.json", document)
        report(1, "Saved-field residual concentration analysis complete", eta_s=0, artifacts=len(artifacts))
    finally:
        stop.set()
        thread.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--artifacts", type=Path, nargs="+", required=True)
    args = parser.parse_args()
    analyze(args.run.resolve(), [path.resolve() for path in args.artifacts])
