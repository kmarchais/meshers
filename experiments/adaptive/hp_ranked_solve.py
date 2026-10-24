"""General lifted-base coupled correction: twelve fresh auxiliary solves."""

import argparse
import json
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

import microgen_workflow as workflow
import numpy as np
from scipy.linalg import eigvalsh

from hp_coupled_preflight import inherited_inputs


def lifted_condensation(energy, coupling, lift_residual, mode_responses, lift_responses):
    h, cross, p = energy[:6, :6], energy[:6, 6:], energy[6:, 6:]
    gram = h + p - cross - cross.T
    mode_residual = cross.T - h
    schur = gram - coupling.T @ mode_responses
    np.testing.assert_allclose(schur, schur.T, atol=1e-12, rtol=1e-9)
    schur = (schur + schur.T) / 2
    assert eigvalsh(schur).min() > 1e-12 * eigvalsh(schur).max()
    condensed_residual = mode_residual - coupling.T @ lift_responses
    coefficients = -np.linalg.solve(schur, condensed_residual)
    correction = -lift_responses - mode_responses @ coefficients
    gain = lift_residual.T @ lift_responses - condensed_residual.T @ coefficients
    direct = -(lift_residual.T @ correction + mode_residual.T @ coefficients)
    np.testing.assert_allclose(gain, direct, atol=1e-12, rtol=1e-9)
    np.testing.assert_allclose(gain, gain.T, atol=1e-12, rtol=1e-9)
    gain = (gain + gain.T) / 2
    tensor = h - gain
    assert eigvalsh(gain).min() >= -1e-12 and eigvalsh(tensor).min() > 0
    return correction, coefficients, tensor, gain, schur


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--setup", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from amg_reference import AuditedAMG
    from capture_reference_system import CapturedSystem
    from stress_recovery import solve
    from hp_shared_domain import periodic_jumps

    setup = json.loads((args.setup / "result.json").read_text())
    assert setup["passed"] and not setup["cg_started"]
    assert json.loads((args.setup / "worker.json").read_text())["peak_rss_bytes"] < workflow.MEMORY
    prepared, coupling_path = Path(setup["prepared"]), Path(setup["coupling"])
    coupling_report = json.loads((coupling_path / "result.json").read_text())
    assert coupling_report["passed"]
    paths = [Path(__file__), args.setup / "result.json", args.setup / "worker.json", coupling_path / "result.json",
        coupling_path / "coupling.npz", coupling_path / "root-maps.npz", prepared / "mesh.npz", prepared / "prolongation.npz"]
    paths.extend(workflow.HERE / name for name in ("test_hp_ranked_solve.py", "hp_coupled_preflight.py",
        "amg_reference.py", "stress_recovery.py", "hp_cross_energy.py"))
    inputs = {**inherited_inputs(args.setup), **{str(path): workflow.digest(path) for path in paths}}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    shutil.copyfile(workflow.HERE / "test_hp_ranked_solve.py", args.output.parent / "test_hp_ranked_solve.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_ranked_solve.py"], check=True, timeout=30)
    mesh = np.load(prepared / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    lift = np.load(prepared / "prolongation.npz")["displacement"]
    data = np.load(coupling_path / "coupling.npz")
    forces, energy, gram = data["forces"], data["energy"], data["mode_matrix"]

    class Correction:
        def __init__(self):
            self.problem, self.records = None, []

        def __call__(self, matrix, rhs, **kwargs):
            matrix = matrix.tocsr()
            transform = self.problem._MatCB
            full = np.zeros((transform.shape[0], 6)); full[:len(forces)] = forces
            coupling = transform.T @ full
            base = np.zeros_like(full); base[:len(forces)] = lift.transpose(0, 2, 1).reshape(6, -1).T
            internal = self.problem.get_A() @ base
            residual = transform.T @ internal
            np.testing.assert_allclose(forces.T @ base[:len(forces)], energy[6:, :6] - energy[:6, :6], atol=1e-12)
            amg = AuditedAMG(elastic_blocks=True); amg.problem = self.problem
            responses = [np.empty_like(coupling), np.empty_like(residual)]
            started = time.monotonic()
            def expired(*_):
                raise TimeoutError("Auxiliary response exceeded 240 seconds")
            previous = signal.signal(signal.SIGALRM, expired)
            try:
                for family, right_hand_sides in enumerate((coupling, residual)):
                    for column in range(6):
                        signal.setitimer(signal.ITIMER_REAL, 240)
                        try:
                            responses[family][:, column] = amg(matrix, right_hand_sides[:, column])
                        finally:
                            signal.setitimer(signal.ITIMER_REAL, 0)
                        workflow.save(args.output.parent / "checkpoint.json", {"completed_auxiliary_responses": 6 * family + column + 1,
                            "response_family": "mode" if family == 0 else "lift-residual", "solver_records": amg.records,
                            "elapsed_s": time.monotonic() - started})
                    if family == 0:
                        schur = gram - coupling.T @ responses[0]
                        assert eigvalsh((schur + schur.T) / 2).min() > 0
            finally:
                signal.signal(signal.SIGALRM, previous)
            correction, coefficients, tensor, gain, schur = lifted_condensation(energy, coupling, residual, *responses)
            final_residual = matrix @ correction + coupling @ coefficients + residual
            relative = np.linalg.norm(final_residual, axis=0) / np.linalg.norm(residual, axis=0)
            assert relative.max() < 2e-9
            mode_residual = coupling.T @ correction + gram @ coefficients + energy[6:, :6] - energy[:6, :6]
            mode_relative = float(np.linalg.norm(mode_residual) / np.linalg.norm(energy[6:, :6] - energy[:6, :6]))
            assert mode_relative < 2e-9
            physical = (transform @ correction)[:3 * len(points)].T.reshape(6, 3, len(points)).transpose(0, 2, 1)
            final_base = lift + physical
            periodic = periodic_jumps(points, surface, final_base)
            pin = int(np.linalg.norm(points, axis=1).argmin())
            assert abs(final_base[:, pin]).max() < 1e-12
            np.savez_compressed(args.output.parent / "hp-fields.npz", base_displacement=final_base, mode_coefficients=coefficients)
            np.savez_compressed(args.output.parent / "auxiliary-responses.npz", mode_responses=responses[0],
                lift_responses=responses[1], free_dofs=self.problem._dof_free)
            workflow.save(args.output, {"status": "ranked-coupled-solved-needs-independent-audit", "prepared": str(prepared),
                "coupling": str(coupling_path), "physical_nodes": len(points), "physical_tets": len(cells),
                "reduced_mesh_unknowns": len(rhs), "augmented_unknowns": len(rhs) + 6, "extra_scalar_mode_unknowns": 6,
                "fresh_auxiliary_responses": 12, "response_families": ["six mode forces", "six lifted-base residuals"],
                "solver_records": amg.records, "schur_eigenvalues": eigvalsh(schur).tolist(),
                "predicted_energy_tensor": tensor.tolist(), "energy_gain_from_lift": gain.tolist(),
                "maximum_relative_gain_from_lift": float(eigvalsh(gain, energy[:6, :6]).max()),
                "complete_field_reduced_residuals": relative.tolist(), "mode_equilibrium_relative_residual": mode_relative,
                "periodic_base": periodic, "base_pin_max": float(abs(final_base[:, pin]).max()),
                "field_archive_sha256": workflow.digest(args.output.parent / "hp-fields.npz"),
                "auxiliary_archive_sha256": workflow.digest(args.output.parent / "auxiliary-responses.npz"),
                "auxiliary_solve_s": time.monotonic() - started, "reference_promoted": False, "native_candidate_promoted": False,
                "field_representation": "Final quadratic base plus six fixed P3-minus-lift periodic modes; mode pin constants are in the coupling record; nodal Tet10 values alone are incomplete",
                "scope": "General non-equilibrium-lift Schur correction, with twelve fresh auxiliary solves; no live workflow or continuum claim."})
            raise CapturedSystem()

    try:
        solve(points, cells, n_elm_gp=15, solver=Correction())
    except CapturedSystem:
        assert args.output.exists()
    else:
        raise AssertionError("Coupled capture must not return a fake macro solution")
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
