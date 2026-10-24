"""Ranked-space coupling and AMG setup gates with a non-equilibrium lift."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as workflow
import numpy as np

from hp_coupled_preflight import inherited_inputs
from hp_shared_domain import ROOT, CUBIC, CUBIC_SOLVE, compose, periodic_jumps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--resource", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--coupling", type=Path, help="Passed coupling record for a setup-only stage")
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from hp_cross_energy import integrate, combine
    from hp_coupled_audit import cubic_nodal_values
    from quadrature_audit import duffy_rule
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4

    prepared = json.loads((args.prepared / "result.json").read_text())
    assert prepared["passed"] and prepared["source_augmented_space_included"] and prepared["geometry"]["geometry_admitted"]
    resource = json.loads((args.resource / "result.json").read_text())
    assert resource["captured"] and resource["no_solve"]
    assert json.loads((args.resource / "worker.json").read_text())["peak_rss_bytes"] < workflow.MEMORY
    paths = [Path(__file__), args.prepared / "result.json", args.prepared / "mesh.npz", args.prepared / "prolongation.npz",
        args.resource / "result.json", args.resource / "worker.json"]
    paths.extend(workflow.HERE / name for name in ("hp_cross_energy.py", "test_hp_cross_energy.py", "hp_shared_domain.py",
        "hp_coupled_audit.py", "stress_recovery.py", "amg_reference.py", "capture_reference_system.py", "hp_coupled_preflight.py"))
    inputs = inherited_inputs(args.prepared)
    if args.coupling:
        inputs.update(inherited_inputs(args.coupling))
        paths.extend(args.coupling / name for name in ("result.json", "coupling.npz", "root-maps.npz"))
    inputs.update({str(path): workflow.digest(path) for path in paths})
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    for name in ("hp_cross_energy.py", "test_hp_cross_energy.py", "hp_coupled_preflight.py"):
        shutil.copyfile(workflow.HERE / name, args.output.parent / name)
    mesh = np.load(args.prepared / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    prolongation = np.load(args.prepared / "prolongation.npz")
    lift = prolongation["displacement"]
    assert int(np.linalg.norm(points, axis=1).argmin()) == prepared["preserved_source_pin_node"]
    if not args.coupling:
        subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_cross_energy.py",
            "experiments/adaptive/test_hp_shared_domain.py"], check=True, timeout=45)
        previous_solve = json.loads((Path(prepared["source_coupled_solve"]) / "result.json").read_text())
        previous_consistent = json.loads((Path(previous_solve["consistent_source"]) / "result.json").read_text())
        screen = json.loads((Path(previous_consistent["source"]) / "result.json").read_text())
        old_maps = np.load(Path(screen["preflight"]) / "composed-maps.npz")
        parents, maps = compose(old_maps["parents"], old_maps["parent_barycentric_maps"],
            prolongation["parents"], prolongation["parent_barycentric_maps"])
        root = np.load(ROOT / "mesh.npz")
        element, linear = Tet10(15), Tet4()
        barycentric = linear.shape_function(element.xi_nd)
        maximum_geometry_error = 0.
        for start in range(0, len(cells), 512):
            subset = slice(start, start + 512)
            locations = np.einsum("pa,cab->cpb", barycentric, maps[subset]) @ linear.xi_nd
            weights = element.shape_function(locations.reshape(-1, 3)).reshape(-1, 10, 10)
            expected = np.einsum("cpn,cnd->cpd", weights, root["points"][root["tetrahedra"][parents[subset]]])
            maximum_geometry_error = max(maximum_geometry_error, float(abs(expected - points[cells[subset]]).max()))
        assert maximum_geometry_error < 1e-12
        np.savez_compressed(args.output.parent / "root-maps.npz", parents=parents, parent_barycentric_maps=maps)
        cubic = np.load(CUBIC / "mesh.npz")
        p_fields = np.load(CUBIC_SOLVE / "fields.npz")["displacement"]
        nodal_p = cubic_nodal_values(points, cells, cubic["tetrahedra"], p_fields, parents, maps)
        periodic_jumps(points, surface, lift)
        periodic_jumps(points, surface, nodal_p)
        from periodic import periodic_pairs
        pairs, _, _ = periodic_pairs(points, surface)
        mode_values = nodal_p - lift
        jump = max(float(abs(mode_values[:, pair[:, 1]] - mode_values[:, pair[:, 0]]).max(initial=0)) for pair in pairs)
        assert jump < 1e-9
        records = []
        for count in (15, 216):
            quadrature = (element.xi_pg, element.w_pg) if count == 15 else duffy_rule(6)
            integrated = integrate(points, cells, lift, cubic["tetrahedra"], p_fields, parents, maps,
                quadrature, return_coupling=count == 15)
            energy, volume, minimum = integrated[:3]
            if count == 15:
                force = integrated[3]
                work = lift.transpose(0, 2, 1).reshape(6, -1) @ force
                np.testing.assert_allclose(work, energy[:6, 6:] - energy[:6, :6], atol=1e-12)
                assert np.linalg.norm(force.reshape(3, len(points), 6).sum(axis=1)) < 1e-10
                matrix = energy[:6, :6] + energy[6:, 6:] - energy[:6, 6:] - energy[6:, :6]
                np.savez_compressed(args.output.parent / "coupling.npz", forces=force, energy=energy, mode_matrix=matrix)
            record = {"points": count, "energy_matrix": energy.tolist(), "volume": volume,
                "minimum_sampled_jacobian": minimum, **combine(energy)}
            records.append(record)
            workflow.save(args.output.parent / f"rule-{count}.json", record)
            print(f"Ranked lift/mode rule {count} integrated", flush=True)
        energies = [np.asarray(row["energy_matrix"]) for row in records]
        grams = [e[:6, :6] + e[6:, 6:] - e[:6, 6:] - e[6:, :6] for e in energies]
        gram_change = float(np.linalg.norm(grams[0] - grams[1]) / np.linalg.norm(grams[1]))
        gains = [np.asarray(row["gain_tensor"]) for row in records]
        gain_change = float(np.linalg.norm(gains[0] - gains[1]) / np.linalg.norm(gains[1]))
        workflow.save(args.output, {"passed": max(gram_change, gain_change) < .01, "prepared": str(args.prepared),
            "physical_nodes": len(points), "physical_tets": len(cells), "rules": records,
            "mode_gram_relative_rule_change": gram_change, "fixed_lift_gain_relative_rule_change": gain_change,
            "precision_gate": .01, "maximum_root_geometry_error": maximum_geometry_error,
            "mode_jump_max": jump, "mode_pin_constants": mode_values[:, prepared["preserved_source_pin_node"]].tolist(),
            "coupling_bytes": force.nbytes, "coupling_sha256": workflow.digest(args.output.parent / "coupling.npz"),
            "root_maps_sha256": workflow.digest(args.output.parent / "root-maps.npz"), "lift_equilibrium_assumed": False,
            "new_equilibrium_solve": False, "reference_promoted": False})
    else:
        from stress_recovery import solve
        from amg_reference import AuditedAMG
        from capture_reference_system import CapturedSystem
        from unittest.mock import patch

        coupling_report = json.loads((args.coupling / "result.json").read_text())
        assert coupling_report["passed"] and coupling_report["coupling_sha256"] == workflow.digest(args.coupling / "coupling.npz")
        forces = np.load(args.coupling / "coupling.npz")["forces"]

        class Setup:
            def __init__(self):
                self.problem, self.records = None, []

            def __call__(self, matrix, rhs, **kwargs):
                matrix = matrix.tocsr()
                transform = self.problem._MatCB
                assert np.array_equal(self.problem.space.get_rank_vector("Disp"), np.arange(3))
                full = np.zeros((transform.shape[0], 6)); full[:len(forces)] = forces
                reduced = transform.T @ full
                base = np.zeros_like(full); base[:len(forces)] = lift.transpose(0, 2, 1).reshape(6, -1).T
                internal = self.problem.get_A() @ base
                residual = transform.T @ internal
                assert reduced.shape == residual.shape == (len(rhs), 6)
                amg = AuditedAMG(elastic_blocks=True); amg.problem = self.problem
                with patch("scipy.sparse.linalg.cg", side_effect=CapturedSystem) as intercepted:
                    try:
                        amg(matrix, reduced[:, 0])
                    except CapturedSystem:
                        assert intercepted.call_count == 1 and not amg.records
                    else:
                        raise AssertionError("Setup must stop before CG")
                workflow.save(args.output, {"passed": True, "prepared": str(args.prepared), "coupling": str(args.coupling),
                    "physical_nodes": len(points), "physical_tets": len(cells), "reduced_unknowns": len(rhs),
                    "augmented_unknowns": len(rhs) + 6, "mode_columns": 6, "lifted_residual_columns": 6,
                    "lift_reduced_residual_norms": np.linalg.norm(residual, axis=0).tolist(),
                    "lift_reduced_residual_relative_to_internal_force": (np.linalg.norm(residual, axis=0) / np.linalg.norm(internal, axis=0)).tolist(),
                    "full_coupling_bytes": forces.nbytes, "reduced_coupling_bytes": reduced.nbytes,
                    "lift_residual_bytes": residual.nbytes, "amg_levels": len(amg.hierarchy.levels),
                    "amg_operator_complexity": amg.hierarchy.operator_complexity(), "cg_started": False,
                    "new_equilibrium_solve": False, "reference_promoted": False})
                raise CapturedSystem()

        try:
            solve(points, cells, n_elm_gp=15, solver=Setup())
        except CapturedSystem:
            assert args.output.exists()
        else:
            raise AssertionError("Setup-only capture must halt")
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
