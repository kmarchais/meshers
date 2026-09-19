"""Periodic small-strain homogenization of Meshers cells with fedoo 1.0 beta.

Six prescribed engineering strains are solved with fedoo PeriodicBC. Independent
stress/energy integration and Meshers pair metadata audit those solutions. The
result is compared with fedoo's automatic stress-controlled homogenization API.
"""

import hashlib
import importlib
import json
import platform
import time
from importlib.metadata import version
from pathlib import Path

import fedoo as fd
import meshers
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import cg

ROOT = Path(__file__).resolve().parents[2]
COMPONENTS = ["E_xx", "E_yy", "E_zz", "E_xy", "E_xz", "E_yz"]
MATERIAL = np.eye(6)
MATERIAL[:3, :3] = 1
np.fill_diagonal(MATERIAL[:3, :3], 3)


class AuditedCG:
    """Fedoo solver callback that checks every reduced linear system."""

    def __init__(self):
        """Create a record of the six load-case solves."""
        self.records = []

    def __call__(self, matrix, rhs, **_kwargs):
        """Solve the fedoo-reduced system.

        Returns:
            Solution after checking CG status and the recomputed residual.
        """
        diagonal = matrix.diagonal()
        assert np.all(diagonal > 0)
        skew = matrix - matrix.T
        assert np.linalg.norm(skew.data) / np.linalg.norm(matrix.data) < 1e-12
        iterations = 0

        def count(_):
            nonlocal iterations
            iterations += 1

        solution, info = cg(
            matrix,
            rhs,
            M=diags(1 / diagonal),
            rtol=1e-12,
            atol=1e-14,
            maxiter=20000,
            callback=count,
        )
        residual = float(np.linalg.norm(matrix @ solution - rhs) / max(np.linalg.norm(rhs), 1e-30))
        assert info == 0, (info, residual)
        assert residual < 2e-9, residual
        self.records.append({
            "iterations": iterations,
            "relative_residual": residual,
            "free_dofs": len(rhs),
        })
        return solution


def strain_tensor(engineering):
    """Convert engineering strain [xx, yy, zz, xy, xz, yz] to a tensor.

    Returns:
        Symmetric 3 by 3 strain tensor.
    """
    xx, yy, zz, xy, xz, yz = engineering
    return np.array([[xx, xy / 2, xz / 2], [xy / 2, yy, yz / 2], [xz / 2, yz / 2, zz]])


def pair_classes(mesh):
    """Audit face-node and triangle matching.

    Returns:
        Representative node ID for each periodic vertex equivalence class.
    """
    parent = np.arange(len(mesh.points))

    def root(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for axis, pairs in enumerate(mesh.periodic_pairs):
        low = np.flatnonzero(np.isclose(mesh.points[:, axis], 0, atol=1e-12, rtol=0))
        high = np.flatnonzero(np.isclose(mesh.points[:, axis], 1, atol=1e-12, rtol=0))
        assert np.array_equal(np.sort(pairs[:, 0]), low)
        assert np.array_equal(np.sort(pairs[:, 1]), high)
        shift = np.eye(3)[axis]
        assert np.max(np.abs(mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]] - shift)) < 1e-12
        mapping = np.arange(len(mesh.points))
        mapping[pairs[:, 1]] = pairs[:, 0]
        sides = []
        for value in [0, 1]:
            triangles = mesh.surface[
                np.all(
                    np.isclose(mesh.points[mesh.surface, axis], value, atol=1e-12, rtol=0), axis=1
                )
            ]
            sides.append({tuple(sorted(row)) for row in mapping[triangles]})
        assert sides[0] == sides[1], axis
        for low_node, high_node in pairs:
            parent[root(high_node)] = root(low_node)
    return np.array([root(i) for i in range(len(parent))])


def element_geometry(mesh):
    """Return signed volumes and physical P1 shape-function gradients."""
    points = mesh.points[mesh.tetrahedra]
    jac = points[:, 1:] - points[:, :1]
    volumes = np.linalg.det(jac) / 6
    assert np.all(volumes > 0)
    gradients = np.empty((len(points), 4, 3))
    gradients[:, 1:] = np.swapaxes(np.linalg.inv(jac), 1, 2)
    gradients[:, 0] = -gradients[:, 1:].sum(axis=1)
    return volumes, gradients


def homogenize(mesh):  # ruff: ignore[too-many-locals, too-many-statements] - keep all six-load-case audits together.
    """Run strain- and stress-controlled fedoo homogenization and cross-check.

    Returns:
        Effective stiffness, directional constants and verification diagnostics.
    """
    start = time.perf_counter()
    classes = pair_classes(mesh)
    volumes, gradients = element_geometry(mesh)
    fd.ModelingSpace("3D")
    fem_mesh = fd.Mesh(mesh.points.copy(), mesh.tetrahedra.copy(), "tet4")
    assembly = fd.Assembly.create(
        fd.weakform.StressEquilibrium(fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)),
        fem_mesh,
        n_elm_gp=1,
    )
    center = int(np.linalg.norm(mesh.points - 0.5, axis=1).argmin())
    solver = AuditedCG()
    stiffness = np.zeros((6, 6))
    reactions = np.zeros((6, 6))
    all_strains, all_stresses, checks = [], [], []
    for case, engineering in enumerate(np.eye(6)):
        problem = fd.problem.Linear(assembly)
        problem.bc.add(fd.constraint.PeriodicBC("small_strain", meshperio=True))
        problem.bc.add("Dirichlet", [center], "Disp", 0)
        for variable, value in zip(COMPONENTS, engineering, strict=True):
            problem.bc.add("Dirichlet", variable, float(value))
        problem.set_solver(solver)
        problem.solve()
        displacement = problem.get_disp().T
        assert np.all(np.isfinite(displacement))
        assert np.max(np.abs(problem.get_dof_solution("MeanStrain").ravel() - engineering)) < 1e-12
        macro = strain_tensor(engineering)
        jump_error = max(
            float(
                np.max(
                    np.abs(
                        displacement[pairs[:, 1]]
                        - displacement[pairs[:, 0]]
                        - (mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]]) @ macro.T
                    )
                )
            )
            for pairs in mesh.periodic_pairs
        )
        assert jump_error < 1e-10, jump_error
        gradient = np.einsum("tia,tib->tab", displacement[mesh.tetrahedra], gradients)
        strain = (gradient + np.swapaxes(gradient, -1, -2)) / 2
        stress = 2 * strain + np.trace(strain, axis1=1, axis2=2)[:, None, None] * np.eye(3)
        # Integrate over the solid, normalize by the full unit-cell volume (1).
        mean_stress = np.einsum("t,tij->ij", volumes, stress)
        stiffness[:, case] = mean_stress[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]
        reactions[:, case] = problem.get_ext_forces("MeanStrain").ravel()
        all_strains.append(strain)
        all_stresses.append(stress)
        # Internal forces must balance within every periodic equivalence class,
        # including edge/corner classes; interior classes contain one node.
        forces = problem.get_ext_forces("Disp", include_mpc=False).T
        balanced = np.zeros_like(forces)
        np.add.at(balanced, classes, forces)
        force_error = float(np.linalg.norm(balanced) / np.linalg.norm(forces))
        assert force_error < 1e-7, force_error
        checks.append({
            "load": COMPONENTS[case],
            "periodic_jump_error": jump_error,
            "periodic_force_balance_relative_error": force_error,
            **solver.records[-1],
        })
    energy_matrix = np.einsum("t,itkl,jtkl->ij", volumes, all_strains, all_stresses)
    hill_mandel = float(np.linalg.norm(energy_matrix - stiffness) / np.linalg.norm(stiffness))
    reaction_error = float(np.linalg.norm(reactions - stiffness) / np.linalg.norm(stiffness))
    symmetry = float(np.linalg.norm(stiffness - stiffness.T) / np.linalg.norm(stiffness))
    assert max(hill_mandel, reaction_error, symmetry) < 1e-7
    eigenvalues = np.linalg.eigvalsh((stiffness + stiffness.T) / 2)
    assert np.all(eigenvalues > 0)
    # Uniform-strain (Voigt) upper bound using the measured solid volume.
    voigt_gap = np.linalg.eigvalsh(volumes.sum() * MATERIAL - (stiffness + stiffness.T) / 2)
    assert voigt_gap.min() > -1e-7
    automatic_solver = AuditedCG()
    automatic = fd.homogen.get_homogenized_stiffness(
        assembly,
        meshperio=True,
        solver=automatic_solver,
        rigid_body_constraint="pin",
    )
    automatic_error = float(np.linalg.norm(automatic - stiffness) / np.linalg.norm(stiffness))
    assert automatic_error < 1e-7, automatic_error
    compliance = np.linalg.inv(stiffness)
    return {
        "tetrahedra": len(mesh.tetrahedra),
        "points": len(mesh.points),
        "solid_volume_fraction": float(volumes.sum()),
        "minimum_quality": mesh.diagnostics["minimum_mmg_quality"],
        "periodic_pair_counts": [len(p) for p in mesh.periodic_pairs],
        "stiffness_engineering_voigt": stiffness.tolist(),
        "stiffness_eigenvalues": eigenvalues.tolist(),
        "young_modulus_xyz": (1 / np.diag(compliance)[:3]).tolist(),
        "poisson_xy_xz_yz": [
            -float(compliance[j, i] / compliance[i, i]) for i, j in [(0, 1), (0, 2), (1, 2)]
        ],
        "hill_mandel_relative_error": hill_mandel,
        "reaction_stress_relative_difference": reaction_error,
        "stiffness_symmetry_relative_error": symmetry,
        "automatic_homogenization_relative_difference": automatic_error,
        "strain_controlled_loads": checks,
        "automatic_stress_controlled_solves": automatic_solver.records,
        "seconds": time.perf_counter() - start,
    }


def specimen(name, cells):
    """Generate a full solid or gyroid sheet in a periodic unit cube.

    Returns:
        Meshers mesh with exact translational pairing on all three axes.
    """
    if name == "solid":
        return meshers.generate(
            lambda x, _y, _z: x * 0 - 1, cells=cells, periodic=(True, True, True), optimize_passes=0
        )
    return meshers.generate(
        "gyroid",
        band=(-0.55, 0.55),
        cells=cells,
        periodic=(True, True, True),
        geometry_tolerance=0.04,
        optimize_passes=4,
    )


def main():
    """Verify homogeneous recovery and measure TPMS refinement behavior."""
    extension = Path(importlib.import_module("meshers._meshers").__file__)
    report = {
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            **{n: version(n) for n in ["fedoo", "meshers", "numpy", "scipy"]},
        },
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "material": {"young_modulus": 2.5, "poisson_ratio": 0.25},
        "geometry": {
            "domain": [0, 1, 0, 1, 0, 1],
            "field": "Meshers built-in gyroid, period 1",
            "band": [-0.55, 0.55],
            "geometry_tolerance": 0.04,
            "optimize_passes": 4,
            "snap": 0.2,
        },
        "voigt_order": ["xx", "yy", "zz", "xy", "xz", "yz"],
        "shear_convention": "Engineering strain gamma_ij=2*epsilon_ij; physical shear stress",
        "normalization": "Full unit-cell volume, including void",
        "periodicity": (
            "Translational in x, y and z; traction-free internal surfaces; one translation pin"
        ),
        "solver": (
            "fedoo assembly/PeriodicBC/homogen API with audited SciPy CG, "
            "Jacobi, rtol=1e-12, atol=1e-14"
        ),
        "limitations": [
            "Linear elasticity only",
            "No rotational periodicity, nonlinear buckling or MMG comparison",
            "TPMS refinement is an observed stability study, not a certified continuum limit",
        ],
        "cases": [],
    }
    for name, levels in [("solid", [6]), ("gyroid_sheet", [12, 18, 24, 32])]:
        previous = None
        for cells in levels:
            result = {"geometry": name, "cells": cells, **homogenize(specimen(name, cells))}
            tensor = np.array(result["stiffness_engineering_voigt"])
            if name == "solid":
                error = float(np.linalg.norm(tensor - MATERIAL) / np.linalg.norm(MATERIAL))
                assert error < 1e-8, error
                result["analytic_stiffness_relative_error"] = error
            if previous is not None:
                result["stiffness_relative_change"] = float(
                    np.linalg.norm(tensor - previous) / np.linalg.norm(tensor)
                )
            previous = tensor
            report["cases"].append(result)
            print(
                json.dumps({
                    k: v
                    for k, v in result.items()
                    if k not in {"strain_controlled_loads", "automatic_stress_controlled_solves"}
                }),
                flush=True,
            )
            (ROOT / "docs/fedoo-homogenization-results.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )


if __name__ == "__main__":
    main()
