"""Six fresh mode-response solves and a small Schur correction on a saved domain."""

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
from hp_shared_domain import SECOND_SOLVE, CUBIC_SOLVE, CUBIC


def condense(energy, coupling, responses):
    h, cross, p = energy[:6, :6], energy[:6, 6:], energy[6:, 6:]
    gram = h + p - cross - cross.T
    residual = cross.T - h
    schur = gram - coupling.T @ responses
    np.testing.assert_allclose(schur, schur.T, atol=1e-12, rtol=1e-9)
    schur = (schur + schur.T) / 2
    assert eigvalsh(schur).min() > 1e-12 * eigvalsh(schur).max(), "Unresolved/nonpositive Schur complement"
    coefficients = -np.linalg.solve(schur, residual)
    gain = -residual.T @ coefficients
    fixed_gain = residual.T @ np.linalg.solve(gram, residual)
    assert eigvalsh(gain - fixed_gain).min() > -1e-12
    tensor = h - gain
    assert eigvalsh(tensor).min() > 0
    return coefficients, tensor, gain, schur


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--consistent", type=Path, required=True)
    parser.add_argument("--setup", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from amg_reference import AuditedAMG
    from capture_reference_system import CapturedSystem
    from stress_recovery import solve
    from periodic import periodic_pairs

    consistent = json.loads((args.consistent / "result.json").read_text())
    setup = json.loads((args.setup / "result.json").read_text())
    assert consistent["passed"] and setup["passed"] and setup["amg_setup_performed"] and not setup["cg_started"]
    nested_source, cubic_source, cubic_mesh = SECOND_SOLVE, CUBIC_SOLVE, CUBIC
    parent25 = consistent.get("parent25", False)
    if parent25:
        from hp_parent25_probe import H_SOLVE, P_MESH, P_SOLVE
        assert consistent["nested_source"] == str(H_SOLVE) and consistent["cubic_source"] == str(P_SOLVE)
        assert consistent["cubic_mesh"] == str(P_MESH)
        assert setup["parent25"] and setup["nested_source"] == str(H_SOLVE)
        assert setup["consistent_source"] == str(args.consistent)
        nested_source, cubic_source, cubic_mesh = H_SOLVE, P_SOLVE, P_MESH
    assert json.loads((args.setup / "worker.json").read_text())["peak_rss_bytes"] < workflow.MEMORY
    inputs = {**inherited_inputs(args.consistent), **inherited_inputs(args.setup)}
    paths = [Path(__file__), args.consistent / "result.json", args.consistent / "coupling.npz",
        args.setup / "result.json", args.setup / "worker.json"]
    paths.extend(workflow.HERE / name for name in ("test_hp_coupled_solve.py", "test_hp_cross_energy.py",
        "hp_cross_energy.py", "amg_reference.py", "stress_recovery.py", "hp_coupled_preflight.py"))
    if parent25:
        paths.append(workflow.HERE / "hp_parent25_probe.py")
    inputs.update({str(path): workflow.digest(path) for path in paths})
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    for name in ("test_hp_coupled_solve.py", "hp_cross_energy.py", "test_hp_cross_energy.py"):
        shutil.copyfile(workflow.HERE / name, args.output.parent / name)
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_coupled_solve.py",
        "experiments/adaptive/test_hp_cross_energy.py"], check=True, timeout=45)
    mesh = np.load(nested_source / "mesh.npz")
    points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    source_fields = np.load(nested_source / "fields.npz")["displacement"]
    data = np.load(args.consistent / "coupling.npz")
    forces, energy = data["forces"], data["energy"]

    class Correction:
        def __init__(self):
            self.problem, self.records = None, []

        def __call__(self, matrix, rhs, **kwargs):
            matrix = matrix.tocsr()
            transform = self.problem._MatCB
            full_fields = np.zeros((transform.shape[0], 6))
            full_fields[:3 * len(points)] = source_fields.transpose(0, 2, 1).reshape(6, -1).T
            internal = self.problem.get_A() @ full_fields
            residuals = np.linalg.norm(transform.T @ internal, axis=0) / np.maximum(np.linalg.norm(internal, axis=0), 1e-30)
            assert residuals.max() < 2e-9, "Saved baseline does not equilibrate on the actual assembled operator"
            padded = np.zeros_like(full_fields)
            padded[:len(forces)] = forces
            reduced = transform.T @ padded
            responses = np.empty_like(reduced)
            amg = AuditedAMG(elastic_blocks=True)
            amg.problem = self.problem
            started = time.monotonic()
            def expired(*_):
                raise TimeoutError("Auxiliary response exceeded 240 seconds")
            previous = signal.signal(signal.SIGALRM, expired)
            try:
                for column in range(6):
                    signal.setitimer(signal.ITIMER_REAL, 240)
                    try:
                        responses[:, column] = amg(matrix, reduced[:, column])
                    finally:
                        signal.setitimer(signal.ITIMER_REAL, 0)
                    workflow.save(args.output.parent / "checkpoint.json", {"completed_auxiliary_responses": column + 1,
                        "solver_records": amg.records, "elapsed_s": time.monotonic() - started})
            finally:
                signal.signal(signal.SIGALRM, previous)
            coefficients, tensor, gain, schur = condense(energy, reduced, responses)
            physical = (transform @ responses)[:3 * len(points)].T.reshape(6, 3, len(points)).transpose(0, 2, 1)
            pairs, _, periodic = periodic_pairs(points, surface)
            jump = max(float(abs(physical[:, matching[:, 1]] - physical[:, matching[:, 0]]).max(initial=0)) for matching in pairs)
            pin = int(np.linalg.norm(points, axis=1).argmin())
            assert jump < 1e-9 and abs(physical[:, pin]).max() < 1e-12
            np.savez_compressed(args.output.parent / "hp-fields.npz", response_displacement=physical,
                mode_coefficients=coefficients)
            workflow.save(args.output, {"status": "coupled-solved-needs-independent-audit", "consistent_source": str(args.consistent),
                "nested_source": str(nested_source), "cubic_source": str(cubic_source), "cubic_mesh": str(cubic_mesh / "mesh.npz"),
                "parent25": parent25, "prepared": consistent.get("prepared"),
                "physical_nodes": len(points), "physical_tets": len(cells), "extra_scalar_mode_unknowns": 6,
                "reduced_mesh_unknowns": len(rhs), "augmented_unknowns": len(rhs) + 6,
                "fresh_auxiliary_responses": 6, "saved_baseline_reduced_residuals": residuals.tolist(),
                "solver_records": amg.records, "mode_coefficients": coefficients.tolist(),
                "schur_eigenvalues": eigvalsh(schur).tolist(), "predicted_energy_tensor": tensor.tolist(),
                "energy_gain_tensor": gain.tolist(), "maximum_relative_gain": float(eigvalsh(gain, energy[:6, :6]).max()),
                "zero_macro_response_jump_max": jump, "response_pin_max": float(abs(physical[:, pin]).max()),
                "periodic": periodic, "auxiliary_solve_s": time.monotonic() - started,
                "field_archive_sha256": workflow.digest(args.output.parent / "hp-fields.npz"),
                "field_representation": "Immutable h/p baseline fields plus six quadratic zero-macro response fields and a common 6x6 mode coefficient matrix; nodal interpolation alone is not the complete field",
                "pin_gauge": "Subtract each p-minus-h mode value at the fine pin when reconstructing displacements; constant shifts do not change energy",
                "reference_promoted": False, "native_candidate_promoted": False,
                "scope": "Six auxiliary zero-macro mode responses and Schur correction on the frozen domain, not fresh macro workflow timings or continuum certification."})
            raise CapturedSystem()

    try:
        solve(points, cells, n_elm_gp=15, solver=Correction())
    except CapturedSystem:
        assert args.output.exists()
    else:
        raise AssertionError("All coupled results must be captured without returning a fake macro solution")
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
