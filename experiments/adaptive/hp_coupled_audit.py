"""Audit full modal fields, positive rules, reduced residuals and all models."""

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
from hp_shared_domain import SECOND_SOLVE, CUBIC_SOLVE, CUBIC, periodic_jumps


def cubic_nodal_values(points, cells, cubic_cells, fields, parents, maps):
    from cubic_reference import register_cubic
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4

    cubic, p2, p1 = register_cubic()(), Tet10(), Tet4()
    nodes, first = np.unique(cells.ravel(), return_index=True)
    np.testing.assert_array_equal(nodes, np.arange(len(points)))
    child, local = np.divmod(first, 10)
    barycentric = p1.shape_function(p2.xi_nd)
    result = np.empty((6, len(points), 3))
    for start in range(0, len(points), 2048):
        selection = slice(start, start + 2048)
        root = np.einsum("na,nab->nb", barycentric[local[selection]], maps[child[selection]]) @ p1.xi_nd
        values = cubic.shape_function(root)
        result[:, selection] = np.einsum("tn,ltnd->ltd", values, fields[:, cubic_cells[parents[child[selection]]]])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from hp_cross_energy import integrate
    from quadrature_audit import duffy_rule
    from fedoo.lib_elements.tetrahedron import Tet10
    from capture_reference_system import CapturedSystem
    from stress_recovery import solve
    from amg_reference import reduced_matrix_key
    from cubic_refinement import MODEL_REPORTS

    solved = json.loads((args.source / "result.json").read_text())
    assert solved["fresh_auxiliary_responses"] == 6
    assert solved["field_archive_sha256"] == workflow.digest(args.source / "hp-fields.npz")
    consistent_path = Path(solved["consistent_source"])
    consistent = json.loads((consistent_path / "result.json").read_text())
    screen = json.loads((Path(consistent["source"]) / "result.json").read_text())
    parent25 = solved.get("parent25", False)
    nested_source, cubic_source, cubic_mesh = SECOND_SOLVE, CUBIC_SOLVE, CUBIC / "mesh.npz"
    if parent25:
        from hp_parent25_probe import H_SOLVE, P_MESH, P_SOLVE, PREPARED
        assert consistent["parent25"] and consistent["nested_source"] == solved["nested_source"] == str(H_SOLVE)
        assert solved["cubic_source"] == str(P_SOLVE) and solved["cubic_mesh"] == str(P_MESH / "mesh.npz")
        assert solved["prepared"] == str(PREPARED)
        nested_source, cubic_source, cubic_mesh, shared = H_SOLVE, P_SOLVE, P_MESH / "mesh.npz", PREPARED
    else:
        shared = Path(screen["preflight"])
    paths = [Path(__file__), args.source / "result.json", args.source / "hp-fields.npz"]
    paths.append(workflow.HERE / "test_hp_coupled_audit.py")
    paths.extend(workflow.HERE / name for name in ("hp_cross_energy.py", "hp_coupled_solve.py", "stress_recovery.py",
        "amg_reference.py", "quadrature_audit.py"))
    if parent25:
        paths.extend(workflow.HERE / name for name in ("hp_batched_energy.py", "hp_parent25_probe.py", "test_hp_cross_energy.py"))
    paths.extend(workflow.HERE / "runs" / name for name in MODEL_REPORTS)
    inputs = {**inherited_inputs(args.source), **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    shutil.copyfile(workflow.HERE / "test_hp_coupled_audit.py", args.output.parent / "test_hp_coupled_audit.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_coupled_audit.py"], check=True, timeout=30)
    mesh, cubic = np.load(nested_source / "mesh.npz"), np.load(cubic_mesh)
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    h = np.load(nested_source / "fields.npz")["displacement"]
    p = np.load(cubic_source / "fields.npz")["displacement"]
    fields = np.load(args.source / "hp-fields.npz")
    responses, coefficients = fields["response_displacement"], fields["mode_coefficients"]
    maps = np.load(shared / ("prolongation.npz" if parent25 else "composed-maps.npz"))
    nodal_p = cubic_nodal_values(points, cells, cubic["tetrahedra"], p, maps["parents"], maps["parent_barycentric_maps"])
    pin = int(np.linalg.norm(points, axis=1).argmin())
    gauge = nodal_p[:, pin] - h[:, pin]
    nodal = h + np.einsum("ij,ind->jnd", coefficients, nodal_p - h - responses - gauge[:, None])
    periodic = periodic_jumps(points, surface, nodal)
    assert abs(nodal[:, pin]).max() < 1e-12
    np.savez_compressed(args.output.parent / "nodal-values.npz", displacement=nodal, mode_pin_constants=gauge)
    basis = np.vstack((np.eye(6) - coefficients, coefficients, -coefficients))
    rules = []
    predicted = np.asarray(solved["predicted_energy_tensor"])
    for count in (15, 125, 216):
        element = Tet10(15)
        quadrature = (element.xi_pg, element.w_pg) if count == 15 else duffy_rule(5 if count == 125 else 6)
        if parent25 and count != 15:
            from hp_batched_energy import integrate_batched
            quadratic_part = h - np.einsum("ij,ind->jnd", coefficients, h + responses)
            cubic_part = np.einsum("ij,ind->jnd", coefficients, p)
            energy, volume, minimum, _, _, mean = integrate_batched(points, cells, quadratic_part,
                cubic["tetrahedra"], cubic_part, maps["parents"], maps["parent_barycentric_maps"], quadrature,
                return_stress=True)
            combined_basis = np.vstack((np.eye(6), np.eye(6)))
            tensor, stress = combined_basis.T @ energy @ combined_basis, mean @ combined_basis
        else:
            energy, volume, minimum, mean = integrate(points, cells, h, cubic["tetrahedra"], p,
                maps["parents"], maps["parent_barycentric_maps"], quadrature,
                additional_fields=responses, return_stress=True)
            tensor, stress = basis.T @ energy @ basis, mean @ basis
        shift = float(abs(eigvalsh(tensor - predicted, predicted)).max())
        balance = float(np.linalg.norm(stress - tensor) / np.linalg.norm(tensor))
        assert minimum > 0 and shift < (1e-7 if count == 15 else 1e-4)
        if count == 15:
            assert balance < 1e-7
        row = {"points": count, "energy_tensor": tensor.tolist(), "stress_tensor": stress.tolist(),
            "energy_vs_predicted_maximum_difference": shift, "stress_energy_relative_difference": balance,
            "minimum_sampled_jacobian": minimum, "volume": volume}
        rules.append(row)
        workflow.save(args.output.parent / f"rule-{count}.json", row)
        print(f"Full modal fields at {count} points: energy shift {shift:.9g}", flush=True)
    previous, last = (np.asarray(row["energy_tensor"]) for row in rules[-2:])
    precision = float(abs(eigvalsh(previous - last, last)).max())
    assert precision < 1e-4
    data = np.load(consistent_path / "coupling.npz")
    coupling, old_energy, gram = data["forces"], data["energy"], data["mode_matrix"]
    response_vectors = responses.transpose(0, 2, 1).reshape(6, -1).T
    schur = gram - coupling.T @ response_vectors
    residual = coupling.T @ h.transpose(0, 2, 1).reshape(6, -1).T
    mode_residual = float(np.linalg.norm(schur @ coefficients + residual) / np.linalg.norm(residual))
    assert mode_residual < 2e-9
    coarse = h - np.einsum("ij,ind->jnd", coefficients, responses)
    residual_report = {}

    class Replay:
        def __init__(self):
            self.problem, self.records = None, []

        def __call__(self, matrix, rhs, **kwargs):
            matrix = matrix.tocsr()
            key = reduced_matrix_key(matrix)[1]
            assert all(record["reduced_matrix_sha256"] == key for record in solved["solver_records"])
            transform = self.problem._MatCB
            padded = np.zeros((transform.shape[0], 6)); padded[:len(coupling)] = coupling
            reduced = transform.T @ padded
            vectors = response_vectors[self.problem._dof_free]
            response_residual = np.linalg.norm(matrix @ vectors - reduced, axis=0) / np.linalg.norm(reduced, axis=0)
            full = np.zeros_like(padded); full[:len(coupling)] = coarse.transpose(0, 2, 1).reshape(6, -1).T
            internal = self.problem.get_A() @ full
            complete = transform.T @ internal + reduced @ coefficients
            complete_residual = np.linalg.norm(complete, axis=0) / np.linalg.norm(internal, axis=0)
            assert max(response_residual.max(), complete_residual.max()) < 2e-9
            residual_report.update({"auxiliary_reduced_residuals": response_residual.tolist(),
                "complete_field_reduced_residuals": complete_residual.tolist()})
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
        for index, model in enumerate(json.loads(path.read_text())["models"]):
            matrix = np.asarray(model["stiffness"])
            values = eigvalsh(last - matrix, matrix)
            models.append({"source": str(path), "index": index, "signed_relative_energy_eigenvalues": values.tolist(),
                "worst_relative_difference": float(abs(values).max())})
    assert len(models) == 24
    worst = max(row["worst_relative_difference"] for row in models)
    workflow.save(args.output, {"passed": True, "source": str(args.source), "parent25": parent25, "rules": rules,
        "positive_rule_energy_difference": precision, "periodic": periodic,
        "pin_maximum_displacement": float(abs(nodal[:, pin]).max()), "mode_equilibrium_relative_residual": mode_residual,
        **residual_report, "retained_model_count": len(models), "models": models,
        "worst_retained_model_disagreement": worst, "passes_target_all_models": worst <= .01,
        "reference_promoted": False, "native_candidate_promoted": False,
        "scope": "Full modal-field energy, stress, periodic, pin and freshly reassembled weak-equilibrium audit; retained empirical model check, not a continuum or workflow-timing claim."})
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())
    print(f"Complete coupled audit passes; 24-model disagreement {worst:.9g}")


if __name__ == "__main__":
    main()
