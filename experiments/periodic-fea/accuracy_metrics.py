"""Compare effective stiffness and local stress fields on independent probes.

The local score is a sampled complementary-energy norm on the common material
domain. It is reference disagreement, not a continuum error bound. Material
missing from either mesh is reported separately, never treated as zero stress.
All six loads participate in the Gram matrix, including their cross terms.
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[name] = "1"

import numpy as np
from scipy.linalg import eigvalsh
from scipy.stats import qmc

from homogenization import MATERIAL, element_geometry
from periodic import tensor_difference
from reference_sensitivity import read
from run import digest, write_json

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report


def field_score(candidate, reference):
    """Worst relative stress energy-norm difference over all mixed loads."""
    if candidate.shape != reference.shape or candidate.ndim != 3 or candidate.shape[1:] != (6, 6):
        raise ValueError("Expected matching sample x stress-component x load arrays")
    if not np.isfinite(candidate).all() or not np.isfinite(reference).all() or not len(reference):
        raise ValueError("Need nonempty finite fields")
    compliance = np.linalg.inv(MATERIAL)
    difference = candidate - reference
    gram = np.einsum("nsi,st,ntj->ij", reference, compliance, reference, optimize=True) / len(reference)
    error = np.einsum("nsi,st,ntj->ij", difference, compliance, difference, optimize=True) / len(reference)
    spectrum = eigvalsh((error + error.T) / 2, (gram + gram.T) / 2)
    return {
        "worst_relative_stress_energy_norm_difference": float(np.sqrt(max(0, spectrum[-1]))),
        "per_unit_load_relative_stress_energy_norm_difference": np.sqrt(np.diag(error) / np.diag(gram)).tolist(),
        "reference_energy_gram": gram.tolist(),
        "stress_error_energy_gram": error.tolist(),
    }


def probe_cells(points, tets, probes):
    """Locate probes with a cell-order-independent choice on shared faces."""
    import pyvista as pv
    from vtkmodules.vtkFiltersCore import vtkProbeFilter
    from vtkmodules.vtkCommonDataModel import vtkStaticCellLocator

    # VTK selects the first matching cell at a shared face/edge. Canonical
    # vertex-set ordering prevents equivalent meshes with permuted cell rows
    # from sampling different sides of a P1 stress jump. Keep original IDs.
    keys = np.sort(tets, axis=1)
    order = np.lexsort(keys.T[::-1])
    grid = pv.UnstructuredGrid(np.c_[np.full(len(tets), 4), tets[order]].ravel(),
                              np.full(len(tets), pv.CellType.TETRA), points)
    grid.cell_data["element_id"] = order.astype(np.int64)
    query = pv.PolyData(probes)
    probe = vtkProbeFilter()
    probe.SetInputData(query)
    probe.SetSourceData(grid)
    locator = vtkStaticCellLocator()
    locator.SetDataSet(grid)
    locator.BuildLocator()
    probe.SetCellLocator(locator)
    probe.SetComputeTolerance(False)
    probe.SetTolerance(1e-12)
    probe.Update()
    sampled = pv.wrap(probe.GetOutput())
    return np.asarray(sampled["element_id"]).copy(), np.asarray(sampled["vtkValidPointMask"], dtype=bool).copy()


def sample_fields(path, probes, geometry_gate=None, *, include_cells=False):
    row, tensor, provenance = read(path, geometry_gate)
    archive = np.load(path.parent / "mesh.npz")
    points, tets = archive["points"], archive["tetrahedra"]
    cells, valid = probe_cells(points, tets, probes)
    selected = tets[cells[valid]]
    _, gradients = element_geometry(SimpleNamespace(points=points, tetrahedra=selected))
    stresses = np.zeros((len(probes), 6, 6))
    for load, field in enumerate(row["homogenization"]["displacement_fields"]):
        field_path = path.parent / field["artifact"]
        if digest(field_path) != field["sha256"]:
            raise ValueError(f"Changed field: {field_path}")
        values = np.load(field_path)
        if not np.array_equal(values["engineering"], np.eye(6)[load]):
            raise ValueError("Unexpected load convention")
        gradient = np.einsum("nia,nib->nab", values["displacement"][selected], gradients)
        strain = (gradient + gradient.transpose(0, 2, 1)) / 2
        stress = 2 * strain + np.trace(strain, axis1=1, axis2=2)[:, None, None] * np.eye(3)
        stresses[valid, :, load] = stress[:, [0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]
    result = row, tensor, provenance, stresses, valid
    return (*result, cells) if include_cells else result


def engineering_constants(tensor):
    compliance = np.linalg.inv(tensor)
    return {
        "young_modulus_xyz": (1 / np.diag(compliance)[:3]).tolist(),
        "shear_modulus_xy_xz_yz": (1 / np.diag(compliance)[3:]).tolist(),
        "hydrostatic_bulk_modulus": float(1 / compliance[:3, :3].sum()),
        "poisson_xy_xz_yz": [-float(compliance[j, i] / compliance[i, i]) for i, j in [(0, 1), (0, 2), (1, 2)]],
    }


def wall_normal_error(path):
    """Area-weighted wall-normal discrepancy, excluding artificial periodic caps."""
    from fixtures import make_case

    archive = np.load(path.parent / "mesh.npz")
    triangles = archive["points"][archive["surface"]]
    cap = np.any(np.all(np.isclose(triangles, -.5, atol=1e-10, rtol=0), axis=1), axis=1)
    cap |= np.any(np.all(np.isclose(triangles, .5, atol=1e-10, rtol=0), axis=1), axis=1)
    triangles = triangles[~cap]
    definition = json.loads((path.parent.parent / "manifest.json").read_text())["specimens"][json.loads(path.read_text())["case"]]
    case = make_case("normal-audit", definition)
    center = triangles.mean(axis=1)
    normal = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    double_area = np.linalg.norm(normal, axis=1)
    normal /= double_area[:, None]
    gradient = np.empty_like(center)
    for axis in range(3):
        delta = np.eye(3)[axis] * 1e-5
        gradient[:, axis] = (case.field(*(center + delta).T) - case.field(*(center - delta).T)) / 2e-5
    gradient /= np.linalg.norm(gradient, axis=1)[:, None]
    # For two-sided sheets, only the plane-normal alignment matters here.
    cosine = np.clip(np.abs(np.einsum("ni,ni->n", normal, gradient)), 0, 1)
    angles = np.arccos(cosine)
    return {"area_weighted_rms_angle_degrees": float(np.rad2deg(np.sqrt(np.dot(double_area, angles**2) / double_area.sum()))),
            "maximum_centroid_angle_degrees": float(np.rad2deg(angles.max())),
            "scope": "Finite-difference analytic normal at wall triangle centroids; cap faces excluded"}


def compare(reference_path, candidates, output, power=17, geometry_gate=None):
    started = time.monotonic()
    report(0, "Load independently verified reference fields", eta_s=180)
    seeds = [20261013 + i for i in range(4)]
    probes = np.vstack([qmc.Sobol(3, scramble=True, seed=seed).random_base2(power - 2)
                       for seed in seeds]) - .5
    blocks = np.repeat(np.arange(4), 2 ** (power - 2))
    reference = sample_fields(reference_path, probes, geometry_gate)
    sampled_candidates = []
    common = reference[4].copy()
    for index, path in enumerate(candidates):
        report((index + 1) / (len(candidates) + 1), f"Probe saved fields: {path.parent.name}",
               eta_s=(time.monotonic() - started) * (len(candidates) - index))
        candidate = sample_fields(path, probes, geometry_gate)
        if candidate[0]["case"] != reference[0]["case"]:
            raise ValueError("Reference and candidate geometries differ")
        sampled_candidates.append(candidate)
        common &= candidate[4]
    result = {
        "scope": "Sampled local stress disagreement on common material; geometry omissions reported separately; no continuum certificate",
        "reference": reference[2], "probe_count": len(probes), "sobol_scramble_seeds": seeds,
        "source_sha256": digest(Path(__file__)), "rows": [],
        "point_location_policy": "VTK static cell locator with canonical vertex-set cell order; original element IDs and 1e-12 tolerance retained",
        "comparison_domain": "Identical intersection of reference and every candidate material domain",
        "common_material_probe_count": int(common.sum()),
        "reference_material_coverage": float(common.sum() / reference[4].sum()),
        "reference_wall_normal_error": wall_normal_error(reference_path),
    }
    for path, candidate in zip(candidates, sampled_candidates, strict=True):
        candidate_constants = engineering_constants(candidate[1])
        reference_constants = engineering_constants(reference[1])
        constants = {}
        for name in candidate_constants:
            a, b = np.asarray(candidate_constants[name]), np.asarray(reference_constants[name])
            constants[name] = {"candidate": a.tolist(), "reference": b.tolist(), "absolute_difference": (a - b).tolist()}
            if name != "poisson_xy_xz_yz":
                constants[name]["relative_difference"] = (a / b - 1).tolist()
        row = {
            **candidate[2], "workflow": candidate[0]["workflow"], "cells": candidate[0]["cells"],
            **{key: candidate[0].get(key) for key in ["nodes", "tetrahedra", "q_min", "inserted_body_nodes",
                                                      "body_marker_kind", "linear_solver", "generation_timing_scope", "estimated_live_recipe_s"]},
            "stiffness": tensor_difference(candidate[1], reference[1]),
            "compliance": tensor_difference(np.linalg.inv(candidate[1]), np.linalg.inv(reference[1])),
            "engineering_constants": constants,
            "wall_normal_error": wall_normal_error(path),
            "common_material_probe_count": int(common.sum()),
            "reference_material_coverage": float(common.sum() / reference[4].sum()),
            "candidate_material_coverage": float(common.sum() / candidate[4].sum()),
            "material_symmetric_difference_fraction_of_cell": float(np.mean(candidate[4] ^ reference[4])),
            **field_score(candidate[3][common], reference[3][common]),
            "four_independent_scramble_scores": [field_score(candidate[3][common & (blocks == b)], reference[3][common & (blocks == b)])["worst_relative_stress_energy_norm_difference"] for b in range(4)],
        }
        result["rows"].append(row)
        write_json(output, result)
        print(json.dumps({k: row[k] for k in ["workflow", "cells", "stiffness", "worst_relative_stress_energy_norm_difference", "reference_material_coverage"]}), flush=True)
    result["seconds"] = time.monotonic() - started
    result["state"] = "complete"
    write_json(output, result)
    report(1, "Local stress comparisons complete", eta_s=0, candidates=len(candidates))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--power", type=int, default=17)
    parser.add_argument("--common-geometry-gate", type=float, help="Recheck all artifacts against this explicit shared geometry gate")
    args = parser.parse_args()
    if not 8 <= args.power <= 21:
        parser.error("--power must be between 8 and 21")
    if args.common_geometry_gate is not None and (not np.isfinite(args.common_geometry_gate) or args.common_geometry_gate <= 0):
        parser.error("--common-geometry-gate must be finite and positive")
    compare(args.reference, args.candidates, args.output, args.power, args.common_geometry_gate)
