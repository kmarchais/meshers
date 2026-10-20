"""Bounded cubic-domain sensitivity of fixed saved fields, without equilibrium claims."""

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


def prepare(args):
    np, *_ = workflow.imports()
    from cubic_reference import register_cubic, basis_checks
    from periodic import periodic_pairs
    from periodic_crease_audit import project_box_curve
    from quadrature_audit import duffy_rule

    row = json.loads((args.source / "result.json").read_text())
    assert row["qualified"] and row["completed_loads"] == 6 and row["order"] == 3
    assert row["quadrature_sensitivity_passed"]
    archive = np.load(args.prepared / "mesh.npz")
    points, cells, surface = (archive[name] for name in ("points", "tetrahedra", "surface"))
    assert cells.shape == (row["elements"], 20) and len(points) == row["nodes"] <= 250000
    assert len(cells) <= workflow.CAP and surface.shape[1] == 10
    element_class = register_cubic()
    checks = basis_checks(element_class())
    _, _, baseline_periodic = periodic_pairs(points, surface)
    baseline = workflow.geometry(args.run, "smooth-lattice", points, cells, surface,
        element_type=element_class.name, quadrature=duffy_rule(5))
    assert baseline["geometry_admitted"]
    assert abs(baseline["volume"] / row["geometry"]["volume"] - 1) < 1e-10
    assert abs(baseline["sampled_two_sided_max"] - row["geometry"]["sampled_two_sided_max"]) < 1e-10
    corners = points[surface[:, :3]]
    cap = np.any(np.all(workflow.on_box_plane(corners, -.5), axis=1) |
                 np.all(workflow.on_box_plane(corners, .5), axis=1), axis=1)
    wall_nodes = np.unique(surface[~cap])
    corner_nodes = np.unique(cells[:, :4])
    selected = np.setdiff1d(wall_nodes, corner_nodes)
    case, _, _ = workflow.specimen("smooth-lattice")
    modified = points.copy()
    modified[selected], projection = project_box_curve(case.function, points[selected])
    np.testing.assert_array_equal(modified[corner_nodes], points[corner_nodes])
    _, _, projected_periodic = periodic_pairs(modified, surface)
    projected = workflow.geometry(args.run, "smooth-lattice", modified, cells, surface,
        element_type=element_class.name, quadrature=duffy_rule(5))
    np.savez_compressed(args.run / "projected-mesh.npz", points=modified,
                        tetrahedra=cells, surface=surface)
    movement = np.linalg.norm(modified[selected] - points[selected], axis=1)
    workflow.save(args.run / "geometry.json", {
        "source": str(args.source), "prepared": str(args.prepared), "basis_checks": checks,
        "nodes": len(points), "elements": len(cells), "projected_wall_nodes": len(selected),
        "movement_quantiles": np.quantile(movement, [0, .1, .5, .9, 1]).tolist(),
        "projection": projection, "baseline": baseline, "projected": projected,
        "baseline_periodic": baseline_periodic, "projected_periodic": projected_periodic,
        "geometry_admitted": bool(baseline["geometry_admitted"] and projected["geometry_admitted"]),
        "corner_coordinates_unchanged": True, "connectivity_unchanged": True,
        "scope": "Explicit research cubic geometry, not a native Meshers feature or an equilibrium result"})


def reintegrate(args):
    np, *_ = workflow.imports()
    from cubic_reference import register_cubic
    from homogenization import MATERIAL, strain_tensor
    from periodic import tensor_difference
    from quadrature_audit import duffy_rule
    from stress_recovery import periodic_displacement_jump_error
    from tet10_integration import integrated_tensors

    geometry = json.loads((args.run / "geometry.json").read_text())
    assert geometry["geometry_admitted"]
    element_class = register_cubic()
    element = element_class()
    archives = {"baseline": np.load(args.prepared / "mesh.npz"),
                "projected": np.load(args.run / "projected-mesh.npz")}
    fields = np.load(args.source / "fields.npz")["displacement"]
    assert fields.shape == (6, geometry["nodes"], 3) and np.isfinite(fields).all()
    np.testing.assert_array_equal(archives["baseline"]["tetrahedra"], archives["projected"]["tetrahedra"])
    np.testing.assert_array_equal(archives["baseline"]["surface"], archives["projected"]["surface"])
    affine_fields = np.array([element.xi_nd @ strain_tensor(load).T for load in np.eye(6)])
    saved = np.asarray(json.loads((args.source / "result.json").read_text())["stiffness"])
    previous, rules = {}, []
    for order in (4, 5):
        quadrature = duffy_rule(order)
        patch_stress, patch_energy, patch_volume = integrated_tensors(element.xi_nd,
            np.arange(20)[None], affine_fields, quadrature, element_type=element_class.name)
        np.testing.assert_allclose(patch_stress, MATERIAL / 6, atol=1e-13, rtol=1e-13)
        np.testing.assert_allclose(patch_energy, MATERIAL / 6, atol=1e-13, rtol=1e-13)
        assert abs(patch_volume - 1 / 6) < 1e-14
        energies, domains = {}, {}
        for name, archive in archives.items():
            points, cells, surface = (archive[key] for key in ("points", "tetrahedra", "surface"))
            jump = periodic_displacement_jump_error(points, surface, fields)
            assert jump < 1e-9
            stress, energy, volume = integrated_tensors(points, cells, fields, quadrature,
                                                       element_type=element_class.name)
            np.testing.assert_allclose(energy, energy.T, atol=1e-12, rtol=0)
            assert np.linalg.eigvalsh(energy).min() > 0
            assert abs(volume / geometry[name]["volume"] - 1) < 1e-10
            domains[name] = {"fixed_field_energy": energy.tolist(), "integrated_stress": stress.tolist(),
                "volume": volume, "periodic_jump_max": jump,
                "stress_energy_relative_difference": float(np.linalg.norm(stress - energy) / np.linalg.norm(energy)),
                "energy_vs_previous_rule": None if name not in previous else tensor_difference(energy, previous[name])}
            energies[name] = energy
            if name == "baseline":
                domains[name]["energy_vs_saved_equilibrium"] = tensor_difference(energy, saved)
                assert domains[name]["energy_vs_saved_equilibrium"]["worst_relative_strain_energy_difference"] < 1e-4
        rules.append({"points": order**3, "domains": domains,
                      "projected_vs_baseline_fixed_energy": tensor_difference(energies["projected"], energies["baseline"])})
        previous = energies
    sensitivity = max(domain["energy_vs_previous_rule"]["worst_relative_strain_energy_difference"]
                      for domain in rules[-1]["domains"].values())
    shift = rules[-1]["projected_vs_baseline_fixed_energy"]["worst_relative_strain_energy_difference"]
    workflow.save(args.run / "sensitivity.json", {"rules": rules,
        "quadrature_sensitivity": sensitivity, "quadrature_sensitivity_passed": sensitivity <= 1e-4,
        "fixed_field_geometry_shift": shift,
        "predeclared_material_geometry_signal": bool(sensitivity <= 1e-4 and shift >= 5e-4),
        "scope": "Same nodal displacement coefficients on two admitted domains; projected fields are not equilibrated, and their energy is not a homogenized stiffness reference or an accuracy certificate"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=workflow.HERE / "runs/g3-lattice-cubic-reference-cache-release-v1")
    parser.add_argument("--prepared", type=Path, default=workflow.HERE / "runs/g3-lattice-cubic-assembly-memory-v1")
    parser.add_argument("--worker", choices=("geometry", "sensitivity"))
    args = parser.parse_args()
    if args.worker:
        (prepare if args.worker == "geometry" else reintegrate)(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    sources = [Path(__file__)] + [workflow.HERE / name for name in (
        "microgen_workflow.py", "cubic_reference.py", "periodic_crease_audit.py",
        "quadrature_audit.py", "stress_recovery.py", "tet10_integration.py", "test_cubic_reference.py")]
    inputs = sources + [args.source / "result.json", args.source / "fields.npz", args.prepared / "mesh.npz"]
    samples = workflow.HERE / "runs/g3-microgen-defaults-v4"
    for suffix in ("-specimen.json", "-samples.npz"):
        path = samples / ("smooth-lattice" + suffix)
        inputs.append(path)
        shutil.copyfile(path, args.run / path.name)
    hashes = {str(path): workflow.digest(path) for path in inputs}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 300, "tetrahedra_cap": workflow.CAP,
        "node_cap": 250000, "geometry_shift_screen": 5e-4, "quadrature_sensitivity_limit": 1e-4,
        "no_equilibrium_solve": True, "fixed_specimen_preserved": True})
    with Job(parts=["controls", "geometry", "sensitivity"]) as job:
        for stage in ("controls", "geometry", "sensitivity"):
            if stage == "sensitivity" and not json.loads((args.run / "geometry.json").read_text())["geometry_admitted"]:
                job.result("geometry-rejected", qualified=0)
                break
            with job.part(stage) as part:
                if stage == "controls":
                    command = [sys.executable, "-m", "pytest", "-q"] + [
                        "experiments/adaptive/test_cubic_reference.py::" + name for name in (
                            "test_cubic_basis_reproduces_quadratic_geometry",
                            "test_cubic_surface_basis_matches_reference_element_faces",
                            "test_nonstandard_geometry_rejects_implicit_linear_fallback")]
                else:
                    command = [sys.executable, __file__, "--run", str(args.run), "--source", str(args.source),
                               "--prepared", str(args.prepared), "--worker", stage]
                record = workflow.bounded(command, args.run / (stage + ".log"), part, timeout=300)
                workflow.save(args.run / (stage + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
