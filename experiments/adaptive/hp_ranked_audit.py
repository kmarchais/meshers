"""Full-field and fresh-operator audit of a generalized ranked correction."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as workflow
import numpy as np
from scipy.linalg import eigvalsh

from hp_coupled_preflight import inherited_inputs
from hp_shared_domain import CUBIC, CUBIC_SOLVE, periodic_jumps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--previous-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from hp_cross_energy import integrate
    from hp_coupled_audit import cubic_nodal_values
    from quadrature_audit import duffy_rule
    from fedoo.lib_elements.tetrahedron import Tet10
    from stress_recovery import solve
    from capture_reference_system import CapturedSystem
    from amg_reference import reduced_matrix_key
    from cubic_refinement import MODEL_REPORTS

    solved = json.loads((args.source / "result.json").read_text())
    assert solved["fresh_auxiliary_responses"] == 12
    assert solved["field_archive_sha256"] == workflow.digest(args.source / "hp-fields.npz")
    assert solved["auxiliary_archive_sha256"] == workflow.digest(args.source / "auxiliary-responses.npz")
    prepared, coupling_path = Path(solved["prepared"]), Path(solved["coupling"])
    previous = json.loads((args.previous_audit / "result.json").read_text())
    assert previous["passed"] and Path(previous["source"]).resolve() == Path(json.loads((prepared / "result.json").read_text())["source_coupled_solve"]).resolve()
    paths = [Path(__file__), args.source / "result.json", args.source / "hp-fields.npz", args.source / "auxiliary-responses.npz",
        args.previous_audit / "result.json"]
    paths.extend(workflow.HERE / name for name in ("hp_cross_energy.py", "hp_coupled_audit.py", "test_hp_coupled_audit.py",
        "stress_recovery.py", "amg_reference.py", "quadrature_audit.py", "cubic_refinement.py"))
    paths.extend(workflow.HERE / "runs" / name for name in MODEL_REPORTS)
    inputs = {**inherited_inputs(args.source), **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    for name in ("test_hp_coupled_audit.py", "hp_cross_energy.py", "hp_coupled_audit.py"):
        shutil.copyfile(workflow.HERE / name, args.output.parent / name)
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_coupled_audit.py"], check=True, timeout=30)
    mesh, cubic = np.load(prepared / "mesh.npz"), np.load(CUBIC / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    lift = np.load(prepared / "prolongation.npz")["displacement"]
    fields = np.load(args.source / "hp-fields.npz")
    base, coefficients = fields["base_displacement"], fields["mode_coefficients"]
    p_fields = np.load(CUBIC_SOLVE / "fields.npz")["displacement"]
    maps = np.load(coupling_path / "root-maps.npz")
    nodal_p = cubic_nodal_values(points, cells, cubic["tetrahedra"], p_fields, maps["parents"], maps["parent_barycentric_maps"])
    pin = int(np.linalg.norm(points, axis=1).argmin())
    gauge = nodal_p[:, pin] - lift[:, pin]
    np.testing.assert_allclose(gauge, json.loads((coupling_path / "result.json").read_text())["mode_pin_constants"], atol=1e-12)
    nodal = base + np.einsum("ij,ind->jnd", coefficients, nodal_p - lift - gauge[:, None])
    periodic = periodic_jumps(points, surface, nodal)
    assert abs(nodal[:, pin]).max() < 1e-12
    np.savez_compressed(args.output.parent / "nodal-values.npz", displacement=nodal, mode_pin_constants=gauge)
    basis = np.vstack((-coefficients, coefficients, np.eye(6)))
    predicted = np.asarray(solved["predicted_energy_tensor"])
    rules = []
    for count in (15, 125, 216):
        element = Tet10(15)
        quadrature = (element.xi_pg, element.w_pg) if count == 15 else duffy_rule(5 if count == 125 else 6)
        energy, volume, minimum, mean = integrate(points, cells, lift, cubic["tetrahedra"], p_fields,
            maps["parents"], maps["parent_barycentric_maps"], quadrature, additional_fields=base, return_stress=True)
        tensor, stress = basis.T @ energy @ basis, mean @ basis
        shift = float(abs(eigvalsh(tensor - predicted, predicted)).max())
        balance = float(np.linalg.norm(tensor - stress) / np.linalg.norm(tensor))
        assert minimum > 0 and shift < (1e-7 if count == 15 else 1e-4)
        if count == 15:
            assert balance < 1e-7
        row = {"points": count, "energy_tensor": tensor.tolist(), "stress_tensor": stress.tolist(),
            "energy_vs_predicted_maximum_difference": shift, "stress_energy_relative_difference": balance,
            "minimum_sampled_jacobian": minimum, "volume": volume}
        rules.append(row); workflow.save(args.output.parent / f"rule-{count}.json", row)
        print(f"Ranked complete fields at {count} points: energy shift {shift:.9g}", flush=True)
    last = np.asarray(rules[-1]["energy_tensor"])
    old_rule = np.asarray(rules[-2]["energy_tensor"])
    precision = float(abs(eigvalsh(last - old_rule, last)).max())
    assert precision < 1e-4
    old_tensor = np.asarray(previous["rules"][-1]["energy_tensor"])
    nested_change = eigvalsh(last - old_tensor, old_tensor)
    assert nested_change.max() <= 1e-4
    data = np.load(coupling_path / "coupling.npz")
    forces, gram = data["forces"], data["mode_matrix"]
    base_vectors = base.transpose(0, 2, 1).reshape(6, -1).T
    mode_residual = forces.T @ base_vectors + gram @ coefficients
    normalization = forces.T @ lift.transpose(0, 2, 1).reshape(6, -1).T
    mode_relative = float(np.linalg.norm(mode_residual) / np.linalg.norm(normalization))
    assert mode_relative < 2e-9
    auxiliary = np.load(args.source / "auxiliary-responses.npz")
    residual_report = {}

    class Replay:
        def __init__(self):
            self.problem, self.records = None, []

        def __call__(self, matrix, rhs, **kwargs):
            matrix = matrix.tocsr(); key = reduced_matrix_key(matrix)[1]
            assert all(row["reduced_matrix_sha256"] == key for row in solved["solver_records"])
            np.testing.assert_array_equal(auxiliary["free_dofs"], self.problem._dof_free)
            transform = self.problem._MatCB
            full = np.zeros((transform.shape[0], 6)); full[:len(forces)] = forces
            coupling = transform.T @ full
            padded_lift = np.zeros_like(full); padded_lift[:len(forces)] = lift.transpose(0, 2, 1).reshape(6, -1).T
            lifted_residual = transform.T @ (self.problem.get_A() @ padded_lift)
            errors = []
            for name, columns in (("mode_responses", coupling), ("lift_responses", lifted_residual)):
                error = np.linalg.norm(matrix @ auxiliary[name] - columns, axis=0) / np.linalg.norm(columns, axis=0)
                errors.extend(error.tolist())
            padded_base = np.zeros_like(full); padded_base[:len(forces)] = base_vectors
            internal = self.problem.get_A() @ padded_base
            complete = transform.T @ internal + coupling @ coefficients
            relative = np.linalg.norm(complete, axis=0) / np.linalg.norm(internal, axis=0)
            assert max(max(errors), relative.max()) < 2e-9
            residual_report.update(auxiliary_reduced_residuals=errors, complete_field_reduced_residuals=relative.tolist())
            raise CapturedSystem()

    try:
        solve(points, cells, n_elm_gp=15, solver=Replay())
    except CapturedSystem:
        assert residual_report
    else:
        raise AssertionError("Residual replay must not solve")
    models = []
    for name in MODEL_REPORTS:
        path = workflow.HERE / "runs" / name
        report = json.loads(path.read_text())
        assert report["case"] == "smooth-lattice" and report["stiffness_target"] == .01
        for index, model in enumerate(report["models"]):
            reference = np.asarray(model["stiffness"])
            values = eigvalsh(last - reference, reference)
            models.append({"source": str(path), "index": index, "signed_relative_energy_eigenvalues": values.tolist(),
                "worst_relative_difference": float(abs(values).max())})
    assert len(models) == 24
    worst = max(row["worst_relative_difference"] for row in models)
    workflow.save(args.output, {"passed": True, "source": str(args.source), "previous_audit": str(args.previous_audit),
        "rules": rules, "positive_rule_energy_difference": precision, "periodic": periodic,
        "pin_maximum_displacement": float(abs(nodal[:, pin]).max()), "mode_equilibrium_relative_residual": mode_relative,
        **residual_report, "nested_energy_change_eigenvalues": nested_change.tolist(),
        "retained_model_count": len(models), "models": models, "worst_retained_model_disagreement": worst,
        "passes_target_all_models": worst <= .01, "reference_promoted": False, "native_candidate_promoted": False,
        "scope": "Complete ranked modal fields and twelve saved auxiliary responses replayed on a fresh operator; all original models retained, no continuum or workflow-timing claim."})
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())
    print(f"Ranked complete audit passes; 24-model disagreement {worst:.9g}")


if __name__ == "__main__":
    main()
