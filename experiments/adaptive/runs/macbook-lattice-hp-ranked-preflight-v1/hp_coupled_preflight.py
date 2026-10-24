"""Consistent-rule and resource gates before any coupled-space solve."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as workflow
import numpy as np
from scipy.linalg import eigvalsh

from hp_shared_domain import SECOND_SOLVE, CUBIC_SOLVE, CUBIC


def inherited_inputs(source):
    """Resolve measured Python snapshots; never relax a data hash."""
    inputs = {}
    measured = json.loads((source / "analysis-inputs.json").read_text())
    snapshot_directories = {source, *(Path(name).parent for name in measured
        if Path(name).parent.is_relative_to(workflow.HERE / "runs"))}
    for name, expected in measured.items():
        path = Path(name)
        if workflow.digest(path) != expected and path.suffix == ".py":
            matches = [candidate for directory in sorted(snapshot_directories) for candidate in directory.glob("*.py")
                if workflow.digest(candidate) == expected]
            assert matches, f"Missing exact measured source: {path}"
            path = matches[0]
        assert workflow.digest(path) == expected, str(path)
        inputs[str(path)] = expected
    return inputs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--assembly", type=Path, help="Read a passed consistent-rule result and measure reduced assembly")
    parser.add_argument("--amg-setup-only", action="store_true")
    parser.add_argument("--resource", type=Path, help="Required passed reduced-assembly record before AMG setup")
    args = parser.parse_args()
    assert not args.output.exists()
    assert not args.amg_setup_only or (args.assembly and args.resource)
    workflow.imports()
    from fedoo.lib_elements.tetrahedron import Tet10
    from hp_cross_energy import integrate, combine

    report = json.loads((args.source / "result.json").read_text())
    assert report["precision_passed"] and report["signal_passed"]
    shared = Path(report["preflight"])
    paths = [Path(__file__), workflow.HERE / "hp_cross_energy.py", workflow.HERE / "test_hp_cross_energy.py",
        workflow.HERE / "stress_recovery.py", workflow.HERE / "capture_reference_system.py",
        args.source / "result.json", args.source / "analysis-inputs.json"]
    inputs = inherited_inputs(args.source)
    if args.assembly:
        paths.extend((args.assembly / "result.json", args.assembly / "coupling.npz", args.assembly / "analysis-inputs.json"))
        inputs.update(inherited_inputs(args.assembly))
    if args.amg_setup_only:
        paths.extend((args.resource / "result.json", args.resource / "worker.json"))
        assert json.loads((args.resource / "result.json").read_text())["passed"]
        assert json.loads((args.resource / "worker.json").read_text())["peak_rss_bytes"] < workflow.MEMORY
    inputs.update({str(path): workflow.digest(path) for path in paths})
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    shutil.copyfile(workflow.HERE / "hp_cross_energy.py", args.output.parent / "hp_cross_energy.py")
    shutil.copyfile(workflow.HERE / "test_hp_cross_energy.py", args.output.parent / "test_hp_cross_energy.py")
    mesh = np.load(SECOND_SOLVE / "mesh.npz")
    points, cells = mesh["points"], mesh["tetrahedra"]
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    if not args.assembly:
        subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_hp_cross_energy.py"], check=True, timeout=45)
        cubic = np.load(CUBIC / "mesh.npz")
        h_fields = np.load(SECOND_SOLVE / "fields.npz")["displacement"]
        p_fields = np.load(CUBIC_SOLVE / "fields.npz")["displacement"]
        maps = np.load(shared / "composed-maps.npz")
        element = Tet10(15)
        energy, volume, minimum, coupling = integrate(points, cells, h_fields, cubic["tetrahedra"], p_fields,
            maps["parents"], maps["parent_barycentric_maps"], (element.xi_pg, element.w_pg), return_coupling=True)
        combined = combine(energy)
        higher = report["rules"][-1]
        old_energy = np.asarray(higher["energy_matrix"])
        matrix = energy[:6, :6] + energy[6:, 6:] - energy[:6, 6:] - energy[6:, :6]
        old_matrix = old_energy[:6, :6] + old_energy[6:, 6:] - old_energy[:6, 6:] - old_energy[6:, :6]
        gram_difference = float(np.linalg.norm(matrix - old_matrix) / np.linalg.norm(old_matrix))
        gain_difference = float(np.linalg.norm(np.asarray(combined["gain_tensor"]) - higher["gain_tensor"])
            / np.linalg.norm(higher["gain_tensor"]))
        work = h_fields.transpose(0, 2, 1).reshape(6, -1) @ coupling
        np.testing.assert_allclose(work, energy[:6, 6:] - energy[:6, :6], atol=1e-12, rtol=1e-9)
        assert np.linalg.norm(coupling.reshape(3, len(points), 6).sum(axis=1)) < 1e-10
        np.savez_compressed(args.output.parent / "coupling.npz", forces=coupling, energy=energy, mode_matrix=matrix)
        workflow.save(args.output, {"passed": gram_difference <= .01 and gain_difference <= .01,
            "source": str(args.source), "quadrature_points": 15, "positive_rule_points_compared": higher["positive_rule_points"],
            "mode_gram_relative_difference": gram_difference, "gain_relative_difference": gain_difference,
            "precision_gate": .01, "volume": volume, "minimum_sampled_jacobian": minimum,
            "full_coupling_shape": list(coupling.shape), "coupling_bytes": coupling.nbytes,
            "coupling_sha256": workflow.digest(args.output.parent / "coupling.npz"),
            "virtual_work_maximum_error": float(abs(work - (energy[:6, 6:] - energy[:6, :6])).max()),
            "mode_matrix_eigenvalues": eigvalsh(matrix).tolist(), **combined,
            "physical_nodes": len(points), "physical_tets": len(cells), "extra_scalar_mode_unknowns": 6,
            "new_equilibrium_solve": False, "reference_promoted": False})
        print(f"15-point mode Gram difference {gram_difference:.9g}, gain difference {gain_difference:.9g}")
    else:
        from capture_reference_system import CapturedSystem
        from stress_recovery import solve
        from amg_reference import verify_symmetric

        consistent = json.loads((args.assembly / "result.json").read_text())
        assert consistent["passed"] and consistent["coupling_sha256"] == workflow.digest(args.assembly / "coupling.npz")
        coupling = np.load(args.assembly / "coupling.npz")["forces"]

        class Capture:
            def __init__(self):
                self.problem, self.records = None, []

            def __call__(self, matrix, rhs, **kwargs):
                matrix = matrix.tocsr()
                assert np.isfinite(matrix.data).all() and (matrix.diagonal() > 0).all()
                if not args.amg_setup_only:
                    verify_symmetric(matrix)
                transform = self.problem._MatCB
                assert np.array_equal(self.problem.space.get_rank_vector("Disp"), np.arange(3))
                assert coupling.shape == (3 * len(points), 6)
                padded = np.zeros((transform.shape[0], 6))
                padded[:len(coupling)] = coupling
                reduced = transform.T @ padded
                assert reduced.shape == (len(rhs), 6) and np.isfinite(reduced).all()
                rng = np.random.default_rng(20261024)
                vector = rng.normal(size=len(rhs))
                np.testing.assert_allclose(vector @ reduced, (transform @ vector) @ padded, atol=1e-10)
                fields = np.load(SECOND_SOLVE / "fields.npz")["displacement"]
                first_field = fields[0].T.ravel()[self.problem._dof_free]
                residual = float(np.linalg.norm(matrix @ first_field - rhs) / np.linalg.norm(rhs))
                assert residual < 2e-9
                setup = {}
                if args.amg_setup_only:
                    from amg_reference import AuditedAMG
                    from unittest.mock import patch
                    amg = AuditedAMG(elastic_blocks=True)
                    amg.problem = self.problem
                    with patch("scipy.sparse.linalg.cg", side_effect=CapturedSystem) as intercepted:
                        try:
                            amg(matrix, reduced[:, 0])
                        except CapturedSystem:
                            assert intercepted.call_count == 1 and not amg.records
                        else:
                            raise AssertionError("AMG capture must stop at the CG boundary")
                    assert amg.symmetry_key == amg.key
                    setup = {"amg_setup_performed": True, "cg_started": False,
                        "amg_levels": len(amg.hierarchy.levels),
                        "amg_operator_complexity": amg.hierarchy.operator_complexity(),
                        "amg_level_shapes": [list(level.A.shape) for level in amg.hierarchy.levels]}
                workflow.save(args.output, {"passed": True, "no_solve": True, "physical_nodes": len(points),
                    "physical_tets": len(cells), "matrix_shape": list(matrix.shape), "matrix_nonzeros": matrix.nnz,
                    "matrix_csr_bytes": matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes,
                    "full_coupling_bytes": coupling.nbytes, "reduced_coupling_bytes": reduced.nbytes,
                    "extra_scalar_mode_unknowns": 6, "augmented_reduced_unknowns": len(rhs) + 6,
                    "complete_symmetry_check_passed": True, "constraint_virtual_work_passed": True,
                    "saved_first_load_reduced_relative_residual": residual,
                    "dense_augmented_matrix_built": False, "amg_setup_performed": False,
                    "new_equilibrium_solve": False, "reference_promoted": False, **setup})
                raise CapturedSystem()

        try:
            solve(points, cells, n_elm_gp=15, solver=Capture())
        except CapturedSystem:
            assert args.output.exists()
        else:
            raise AssertionError("Capture must stop before a solution")
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
