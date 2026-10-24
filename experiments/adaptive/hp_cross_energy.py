"""Cross-energy of saved conforming h/p fields on a verified shared domain."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as workflow
import numpy as np
from scipy.linalg import eigvalsh

from hp_shared_domain import SECOND_SOLVE, CUBIC_SOLVE, CUBIC, ROOT


def integrate(points, cells, h_fields, cubic_cells, p_fields, parents, maps, quadrature):
    """Integrate all 144 cross terms, grouping identical affine child maps."""
    from cubic_reference import register_cubic
    from fedoo.lib_elements.tetrahedron import Tet10, Tet4

    linear, quadratic, cubic = Tet4(), Tet10(), register_cubic()()
    locations, weights = quadrature
    child_derivatives = quadratic.shape_function_derivative(locations)
    barycentric = linear.shape_function(locations)
    unique, groups = np.unique(maps.reshape(-1, 16), axis=0, return_inverse=True)
    reference_derivative = linear.shape_function_derivative([[.25] * 3])[0]
    energy, volume, minimum = np.zeros((12, 12)), 0., np.inf
    for group, flat in enumerate(unique):
        mapping = flat.reshape(4, 4)
        selected = np.flatnonzero(groups == group)
        root_locations = barycentric @ mapping @ linear.xi_nd
        transform = reference_derivative @ mapping @ linear.xi_nd
        p_derivatives = np.einsum("ab,qbn->qan", transform, cubic.shape_function_derivative(root_locations))
        # Bound scratch arrays independently of the largest map group.
        for start in range(0, len(selected), 2048):
            subset = selected[start:start + 2048]
            geometry = points[cells[subset]]
            h_values = h_fields[:, cells[subset]]
            p_values = p_fields[:, cubic_cells[parents[subset]]]
            for derivative_h, derivative_p, weight in zip(child_derivatives, p_derivatives, weights, strict=True):
                jacobian = np.einsum("an,tnb->tab", derivative_h, geometry)
                determinant = np.linalg.det(jacobian)
                assert determinant.min() > 0
                minimum = min(minimum, float(determinant.min()))
                inverse = np.linalg.inv(jacobian)
                gradient_h = np.einsum("tba,an->tbn", inverse, derivative_h)
                gradient_p = np.einsum("tba,an->tbn", inverse, derivative_p)
                gradients = np.concatenate((np.einsum("ltni,tjn->ltij", h_values, gradient_h),
                    np.einsum("ltni,tjn->ltij", p_values, gradient_p)))
                strain = (gradients + gradients.swapaxes(-1, -2)) / 2
                stress = 2 * strain + np.trace(strain, axis1=-2, axis2=-1)[..., None, None] * np.eye(3)
                energy += np.einsum("t,itkl,jtkl->ij", weight * determinant, strain, stress, optimize=True)
                volume += float(weight * determinant.sum())
    return energy, volume, minimum


def combine(energy):
    """Optimize six periodic p-minus-h modes with the original fields fixed."""
    np.testing.assert_allclose(energy, energy.T, atol=1e-12, rtol=1e-12)
    h, cross, p = energy[:6, :6], energy[:6, 6:], energy[6:, 6:]
    matrix = h + p - cross - cross.T
    residual = cross.T - h
    matrix = (matrix + matrix.T) / 2
    assert eigvalsh(matrix).min() > 1e-12 * eigvalsh(matrix).max(), "Mode Gram matrix is not resolved SPD"
    coefficients = -np.linalg.solve(matrix, residual)
    gain = -residual.T @ coefficients
    trial = h - gain
    basis = np.vstack((np.eye(6) - coefficients, coefficients))
    direct = basis.T @ energy @ basis
    np.testing.assert_allclose(trial, direct, atol=1e-12, rtol=1e-10)
    assert eigvalsh(trial).min() > 0
    assert eigvalsh(gain).min() > -1e-12
    assert eigvalsh(p - trial).min() > -1e-12
    return {"coefficients": coefficients.tolist(), "mode_matrix_eigenvalues": eigvalsh(matrix).tolist(),
        "mode_solve_relative_residual": float(np.linalg.norm(matrix @ coefficients + residual) / np.linalg.norm(residual)),
        "trial_energy_tensor": trial.tolist(), "gain_tensor": gain.tolist(),
        "maximum_relative_gain": float(eigvalsh(gain, h).max()),
        "direct_energy_identity_relative_error": float(np.linalg.norm(direct - trial) / np.linalg.norm(h))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from periodic import tensor_difference
    from quadrature_audit import duffy_rule

    preflight = json.loads((args.preflight / "result.json").read_text())
    assert preflight["passed"] and preflight["maximum_composed_geometry_error"] < 1e-12
    assert preflight["composed_maps_sha256"] == workflow.digest(args.preflight / "composed-maps.npz")
    paths = [Path(__file__), workflow.HERE / "test_hp_cross_energy.py", workflow.HERE / "hp_shared_domain.py",
        workflow.HERE / "microgen_workflow.py", workflow.HERE / "cubic_reference.py", workflow.HERE / "quadrature_audit.py",
        args.preflight / "result.json", args.preflight / "composed-maps.npz", args.preflight / "analysis-inputs.json"]
    inherited = json.loads((args.preflight / "analysis-inputs.json").read_text())
    assert all(workflow.digest(Path(path)) == expected for path, expected in inherited.items())
    inputs = {**inherited, **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    shutil.copyfile(workflow.HERE / "test_hp_cross_energy.py", args.output.parent / "test_hp_cross_energy.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_cross_energy.py"], check=True, timeout=45)
    mesh = np.load(SECOND_SOLVE / "mesh.npz")
    cubic = np.load(CUBIC / "mesh.npz")
    h_fields = np.load(SECOND_SOLVE / "fields.npz")["displacement"]
    p_fields = np.load(CUBIC_SOLVE / "fields.npz")["displacement"]
    maps = np.load(args.preflight / "composed-maps.npz")
    records = []
    for order in (5, 6):
        energy, volume, minimum = integrate(mesh["points"], mesh["tetrahedra"], h_fields, cubic["tetrahedra"], p_fields,
            maps["parents"], maps["parent_barycentric_maps"], duffy_rule(order))
        record = {"positive_rule_points": order**3, "energy_matrix": energy.tolist(), "volume": volume,
            "minimum_sampled_jacobian": minimum, **combine(energy)}
        for name, block, source in (("nested", energy[:6, :6], SECOND_SOLVE), ("cubic", energy[6:, 6:], CUBIC_SOLVE)):
            tensor = np.asarray(json.loads((source / "result.json").read_text())["stiffness"])
            difference = tensor_difference(block, tensor)
            assert difference["worst_relative_strain_energy_difference"] < 1e-4
            record[name + "_own_energy_difference"] = difference
        records.append(record)
        workflow.save(args.output.parent / f"rule-{order**3}.json", record)
        print(f"Positive rule {order**3}: trial gain {record['maximum_relative_gain']:.9g}", flush=True)
    gains = [np.asarray(record["gain_tensor"]) for record in records]
    precision = float(np.linalg.norm(gains[-1] - gains[0]) / np.linalg.norm(gains[-1]))
    workflow.save(args.output, {"preflight": str(args.preflight), "root": str(ROOT), "rules": records,
        "gain_rule_relative_change": precision, "precision_gate": .01, "precision_passed": precision <= .01,
        "signal_screen": .0005, "signal_passed": records[-1]["maximum_relative_gain"] >= .0005,
        "all_sources_unchanged": True, "new_equilibrium_solve": False, "reference_promoted": False,
        "scope": "Six common periodic p-minus-h trial modes; fixed old coefficients and physical domain; no re-equilibration, continuum bound, native candidate or workflow timing claim."})
    assert all(workflow.digest(Path(path)) == expected for path, expected in inputs.items())


if __name__ == "__main__":
    main()
