"""Shared-domain preflight and resource probe for saved parent25 h/p fields."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import microgen_workflow as workflow

RUNS = workflow.HERE / "runs"
ROOT = RUNS / "g3-lattice-cubic-parent25-fields-v1"
PREPARED = RUNS / "macbook-lattice-parent25-preferred-only-preflight-v4"
H_SOLVE = RUNS / "macbook-lattice-parent25-preferred-only-v1"
P_MESH = RUNS / "g3-lattice-cubic25-assembly-memory-v1"
P_SOLVE = RUNS / "g3-lattice-cubic25-reference-blocked-symmetry-v2"
P_AUDIT = RUNS / "g3-lattice-cubic25-assessment-v1"


def old_path(name):
    path = Path(name.replace("/home/kmarchais/dev/", str(workflow.ROOT.parent)+"/", 1))
    return path if path.is_absolute() else workflow.ROOT / path


def preflight(output):
    workflow.imports()
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4
    from cubic_reference import register_cubic, basis_checks
    from hp_coupled_audit import cubic_nodal_values
    from hp_shared_domain import periodic_jumps
    from nested_refinement import interpolation_check
    from hp_cross_energy import integrate

    paths = [Path(__file__)]
    paths.extend(workflow.HERE / name for name in ("microgen_workflow.py", "hp_shared_domain.py", "hp_cross_energy.py",
        "hp_coupled_audit.py", "cubic_reference.py", "nested_refinement.py", "quadrature_audit.py",
        "test_hp_cross_energy.py", "test_hp_shared_domain.py", "test_hp_coupled_audit.py"))
    for directory in (ROOT, PREPARED, H_SOLVE, P_MESH, P_SOLVE, P_AUDIT):
        paths.extend(directory / name for name in ("result.json", "manifest.json", "mesh.npz", "fields.npz", "prolongation.npz", "assessment.json")
            if (directory / name).exists())
    paths.append(H_SOLVE / (H_SOLVE.name+".json"))
    transfer_manifest = workflow.ROOT.parent / "meshers-handoff-20261024/data-manifest.json"
    paths.append(transfer_manifest)
    inputs = {str(path): workflow.digest(path) for path in paths}
    workflow.save(output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, output.parent / "analysis-driver.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_cross_energy.py",
        "experiments/adaptive/test_hp_shared_domain.py", "experiments/adaptive/test_hp_coupled_audit.py"], check=True, timeout=45)
    transferred = json.loads(transfer_manifest.read_text())
    frozen = {str((workflow.ROOT.parent / repo["name"] / item["path"]).resolve()): item["sha256"]
        for repo in transferred["repositories"] for item in repo["files"]}
    frozen_checked = 0
    for path in paths:
        if path.parent.name.startswith("g3-") and path.suffix == ".npz":
            assert inputs[str(path)] == frozen[str(path.resolve())]
            frozen_checked += 1
    for directory in (P_SOLVE, P_AUDIT):
        manifest = json.loads((directory / "manifest.json").read_text())
        for name, expected in manifest.get("inputs", manifest.get("sources")).items():
            path = old_path(name)
            if path.suffix == ".npz":
                assert workflow.digest(path) == expected
    h_row, p_row, p_audit = (json.loads((directory / "result.json").read_text()) for directory in (H_SOLVE, P_SOLVE, P_AUDIT))
    for row in (h_row, p_row):
        assert row["case"] == "smooth-lattice" and row["qualified"] and row["completed_loads"] == 6
        assert len(row["solver_records"]) == 6 and max(record["relative_residual"] for record in row["solver_records"]) < 2e-9
    assert p_audit["fresh_process_audit_passed"] and p_audit["retained_model_count"] == 24
    h_audit = json.loads((H_SOLVE / (H_SOLVE.name+".json")).read_text())
    assert h_audit["rules"][-1]["energy_vs_previous_rule"]["worst_relative_strain_energy_difference"] < 1e-4
    h_manifest = json.loads((H_SOLVE / "manifest.json").read_text())["inputs"]
    assert h_manifest[str(PREPARED / "mesh.npz")] == workflow.digest(PREPARED / "mesh.npz")
    root, mesh, cubic = (dict(np.load(directory / "mesh.npz")) for directory in (ROOT, H_SOLVE, P_MESH))
    for key in ("points", "tetrahedra", "surface"):
        np.testing.assert_array_equal(mesh[key], np.load(PREPARED / "mesh.npz")[key])
    h, p = (np.load(directory / "fields.npz")["displacement"] for directory in (H_SOLVE, P_SOLVE))
    for data, fields in ((mesh, h), (cubic, p)):
        assert fields.shape == (6, len(data["points"]), 3) and np.isfinite(fields).all()
        assert len(data["points"]) <= 250000 and len(data["tetrahedra"]) <= workflow.CAP
    maps = dict(np.load(PREPARED / "prolongation.npz"))
    reproduction = interpolation_check(root["points"], root["tetrahedra"], np.load(ROOT / "fields.npz")["displacement"], {**mesh, **maps})
    quadratic, linear, cubic_element = Tet10(), Tet4(), register_cubic()()
    basis_checks(cubic_element)
    assert len(cubic["tetrahedra"]) == len(root["tetrahedra"])
    expected = np.einsum("pn,cnd->cpd", quadratic.shape_function(cubic_element.xi_nd), root["points"][root["tetrahedra"]])
    cubic_error = float(abs(expected-cubic["points"][cubic["tetrahedra"]]).max())
    assert cubic_error < 1e-12
    vertices = maps["parent_barycentric_maps"] @ linear.xi_nd
    ratios = np.linalg.det(vertices[:, 1:]-vertices[:, :1])/np.linalg.det(linear.xi_nd[1:]-linear.xi_nd[:1])
    assert ratios.min() > 0
    np.testing.assert_allclose(np.bincount(maps["parents"], weights=ratios, minlength=len(root["tetrahedra"])), 1., atol=1e-13)
    periodic = dict(nested=periodic_jumps(mesh["points"], mesh["surface"], h),
        cubic=periodic_jumps(cubic["points"], cubic["surface"], p))
    nodal_p = cubic_nodal_values(mesh["points"], mesh["tetrahedra"], cubic["tetrahedra"], p, maps["parents"], maps["parent_barycentric_maps"])
    periodic["cubic_on_nested_nodes"] = periodic_jumps(mesh["points"], mesh["surface"], nodal_p)
    pin = int(np.linalg.norm(mesh["points"], axis=1).argmin())
    _, groups = np.unique(maps["parent_barycentric_maps"].reshape(-1, 16), axis=0, return_inverse=True)
    counts = np.bincount(groups)
    phase = time.perf_counter()
    element = Tet10(4)
    _, volume, minimum = integrate(mesh["points"], mesh["tetrahedra"], h, cubic["tetrahedra"], p,
        maps["parents"], maps["parent_barycentric_maps"], (element.xi_pg, element.w_pg))
    benchmark_s = time.perf_counter()-phase
    assert minimum > 0 and volume > 0
    workflow.save(output, dict(passed=True, parent25_shared_domain_preflight=True,
        root=str(ROOT), prepared=str(PREPARED), nested=str(H_SOLVE), cubic_mesh=str(P_MESH), cubic_solve=str(P_SOLVE), cubic_audit=str(P_AUDIT),
        frozen_npz_inputs_verified=frozen_checked, original_reproduction=reproduction, maximum_cubic_embedding_error=cubic_error,
        minimum_reference_subcell_ratio=float(ratios.min()), periodic=periodic, mode_pin_constants=(nodal_p[:, pin]-h[:, pin]).tolist(),
        map_count=len(counts), map_group_size_quantiles=np.quantile(counts, [0,.5,.9,1]).tolist(),
        physical_nested_nodes=len(mesh["points"]), physical_nested_tets=len(mesh["tetrahedra"]), cubic_nodes=len(cubic["points"]),
        four_point_benchmark_s=benchmark_s, conservative_positive_rule_time_estimate_s=2*benchmark_s*(125+216)/4,
        resource_probe_only=True, new_equilibrium_solve=False, reference_promoted=False,
        scope="Exact shared parent25 domain and field provenance/periodicity plus low-order resource timing. The four-point tensor is discarded and supports no stiffness, trial-gain, continuum or live-workflow claim."))
    assert all(workflow.digest(Path(path)) == expected for path, expected in inputs.items())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    preflight(args.output)


if __name__ == "__main__":
    main()
