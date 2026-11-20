# Copyright (c) 2026 Kevin Marchais

"""Periodic small-strain homogenization of Meshers cells with fedoo 1.0 beta.

Six prescribed engineering strains are solved with fedoo PeriodicBC. Independent
stress/energy integration and Meshers pair metadata audit those solutions. The
result is compared with fedoo's automatic stress-controlled homogenization API.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import platform
import time
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import fedoo as fd
import meshers
import numpy as np
from reporting import write_json_line
from scipy.sparse import diags
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse.linalg import cg
from solver_checks import (
    MAX_RELATIVE_RESIDUAL,
    require,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping, Sequence

    from numpy.typing import NDArray
    from scipy.sparse import sparray, spmatrix


MATRIX_SYMMETRY_RTOL = 1e-12
PERIODIC_PAIR_ATOL = 1e-12
MEAN_STRAIN_ATOL = 1e-12
PERIODIC_JUMP_ATOL = 1e-10
MAX_PERIODIC_FORCE_ERROR = 1e-7
HOMOGENIZATION_TOLERANCE = 1e-7
SOLID_REFERENCE_RTOL = 1e-8
ROOT = Path(__file__).resolve().parents[2]
COMPONENTS = ["E_xx", "E_yy", "E_zz", "E_xy", "E_xz", "E_yz"]
MATERIAL = np.eye(6)
MATERIAL[:3, :3] = 1
np.fill_diagonal(MATERIAL[:3, :3], 3)


class AuditedCG:
    """Fedoo solver callback that checks every reduced linear system."""

    def __init__(self) -> None:
        """Create a record of the six load-case solves."""
        self.records = []

    def __call__(
        self,
        matrix: spmatrix | sparray,
        rhs: NDArray[np.float64],
        **_kwargs: object,
    ) -> NDArray[np.float64]:
        """Solve the fedoo-reduced system.

        Returns:
            Solution after checking CG status and the recomputed residual.
        """
        diagonal = matrix.diagonal()
        require(np.all(diagonal > 0), "CG requires a positive diagonal")
        skew = matrix - matrix.T
        require(
            np.linalg.norm(skew.data) / np.linalg.norm(matrix.data) < MATRIX_SYMMETRY_RTOL,
            "Reduced stiffness matrix is asymmetric",
        )
        iterations = 0

        def count(_: NDArray[np.float64]) -> None:
            nonlocal iterations
            iterations += 1

        solve_start = time.perf_counter()
        solution, info = cg(
            matrix,
            rhs,
            M=diags(1 / diagonal),
            rtol=1e-12,
            atol=1e-14,
            maxiter=20000,
            callback=count,
        )
        solve_seconds = time.perf_counter() - solve_start
        residual = float(np.linalg.norm(matrix @ solution - rhs) / max(np.linalg.norm(rhs), 1e-30))
        require(info == 0, f"CG returned status {info} with relative residual {residual}")
        require(
            residual < MAX_RELATIVE_RESIDUAL,
            f"CG relative residual {residual} exceeds the limit",
        )
        self.records.append({
            "iterations": iterations,
            "cg_seconds": solve_seconds,
            "relative_residual": residual,
            "free_dofs": len(rhs),
        })
        return solution


class ReorderedCG(AuditedCG):
    """Run the same Jacobi CG after symmetric reverse Cuthill-McKee ordering."""

    def __call__(
        self,
        matrix: spmatrix | sparray,
        rhs: NDArray[np.float64],
        **kwargs: object,
    ) -> NDArray[np.float64]:
        """Return the solution in the original degree-of-freedom order."""
        start = time.perf_counter()
        csr = matrix.tocsr()
        order = reverse_cuthill_mckee(csr, symmetric_mode=True)
        permuted = csr[order][:, order]
        ordering_seconds = time.perf_counter() - start
        permuted_solution = super().__call__(permuted, rhs[order], **kwargs)
        x = np.empty_like(permuted_solution)
        x[order] = permuted_solution
        residual = float(np.linalg.norm(matrix @ x - rhs) / max(np.linalg.norm(rhs), 1e-30))
        require(
            residual < MAX_RELATIVE_RESIDUAL,
            f"Original-system CG relative residual {residual} exceeds the limit",
        )
        self.records[-1].update(
            ordering_seconds=ordering_seconds,
            original_relative_residual=residual,
        )
        return x


class _PeriodicMesh(Protocol):
    """Tet4 mesh data used by both the native and reconstructed solver inputs."""

    @property
    def points(self) -> NDArray[np.float64]: ...

    @property
    def tetrahedra(self) -> NDArray[np.int64]: ...

    @property
    def surface(self) -> NDArray[np.int64]: ...

    @property
    def periodic_pairs(self) -> Sequence[NDArray[np.int64]]: ...

    @property
    def diagnostics(self) -> Mapping[str, object]: ...


def strain_tensor(engineering: NDArray[np.float64]) -> NDArray[np.float64]:
    """Convert engineering strain [xx, yy, zz, xy, xz, yz] to a tensor.

    Returns:
        Symmetric 3 by 3 strain tensor.
    """
    xx, yy, zz, xy, xz, yz = engineering
    return np.array([[xx, xy / 2, xz / 2], [xy / 2, yy, yz / 2], [xz / 2, yz / 2, zz]])


def pair_classes(mesh: _PeriodicMesh) -> NDArray[np.intp]:
    """Audit face-node and triangle matching.

    Returns:
        Representative node ID for each periodic vertex equivalence class.
    """
    parent = np.arange(len(mesh.points))

    def root(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for axis, pairs in enumerate(mesh.periodic_pairs):
        low = np.flatnonzero(np.isclose(mesh.points[:, axis], 0, atol=1e-12, rtol=0))
        high = np.flatnonzero(np.isclose(mesh.points[:, axis], 1, atol=1e-12, rtol=0))
        require(np.array_equal(np.sort(pairs[:, 0]), low), "Incomplete low-side periodic pairs")
        require(np.array_equal(np.sort(pairs[:, 1]), high), "Incomplete high-side periodic pairs")
        shift = np.eye(3)[axis]
        require(
            np.max(np.abs(mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]] - shift))
            < PERIODIC_PAIR_ATOL,
            "Periodic pair has incorrect translation",
        )
        mapping = np.arange(len(mesh.points))
        mapping[pairs[:, 1]] = pairs[:, 0]
        sides = []
        for value in [0, 1]:
            triangles = mesh.surface[
                np.all(
                    np.isclose(mesh.points[mesh.surface, axis], value, atol=1e-12, rtol=0),
                    axis=1,
                )
            ]
            sides.append({tuple(sorted(row)) for row in mapping[triangles]})
        require(sides[0] == sides[1], f"Opposite face triangles differ on axis {axis}")
        for low_node, high_node in pairs:
            parent[root(high_node)] = root(low_node)
    return np.array([root(i) for i in range(len(parent))])


def element_geometry(mesh: _PeriodicMesh) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return signed volumes and physical P1 shape-function gradients."""
    points = mesh.points[mesh.tetrahedra]
    jac = points[:, 1:] - points[:, :1]
    volumes = np.linalg.det(jac) / 6
    require(np.all(volumes > 0), "Nonpositive tetrahedron volume")
    gradients = np.empty((len(points), 4, 3))
    gradients[:, 1:] = np.swapaxes(np.linalg.inv(jac), 1, 2)
    gradients[:, 0] = -gradients[:, 1:].sum(axis=1)
    return volumes, gradients


@dataclass
class _StrainResponse:
    strain: NDArray[np.float64]
    stress: NDArray[np.float64]
    mean_stress: NDArray[np.float64]
    reaction: NDArray[np.float64]
    checks: dict[str, object]


def _prescribed_strain(
    assembly: fd.Assembly,
    center: int,
    engineering: NDArray[np.float64],
    solver: AuditedCG,
) -> tuple[fd.problem.Linear, NDArray[np.float64]]:
    problem = fd.problem.Linear(assembly)
    problem.bc.add(fd.constraint.PeriodicBC("small_strain", meshperio=True))
    problem.bc.add("Dirichlet", [center], "Disp", 0)
    for variable, value in zip(COMPONENTS, engineering, strict=True):
        problem.bc.add("Dirichlet", variable, float(value))
    problem.set_solver(solver)
    problem.solve()
    displacement = problem.get_disp().T
    require(np.all(np.isfinite(displacement)), "Nonfinite displacement")
    require(
        np.max(np.abs(problem.get_dof_solution("MeanStrain").ravel() - engineering))
        < MEAN_STRAIN_ATOL,
        "Prescribed mean strain was lost",
    )
    return problem, displacement


def _periodic_jump(
    mesh: _PeriodicMesh,
    displacement: NDArray[np.float64],
    engineering: NDArray[np.float64],
) -> float:
    macro = strain_tensor(engineering)
    jump_error = max(
        float(
            np.max(
                np.abs(
                    displacement[pairs[:, 1]]
                    - displacement[pairs[:, 0]]
                    - (mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]]) @ macro.T,
                ),
            ),
        )
        for pairs in mesh.periodic_pairs
    )
    require(jump_error < PERIODIC_JUMP_ATOL, f"Periodic displacement jump error {jump_error}")
    return jump_error


def _force_balance(problem: fd.problem.Linear, classes: NDArray[np.intp]) -> float:
    # Internal forces must balance within every periodic equivalence class,
    # including edge/corner classes; interior classes contain one node.
    forces = problem.get_ext_forces("Disp", include_mpc=False).T
    balanced = np.zeros_like(forces)
    np.add.at(balanced, classes, forces)
    force_error = float(np.linalg.norm(balanced) / np.linalg.norm(forces))
    require(force_error < MAX_PERIODIC_FORCE_ERROR, f"Periodic force error {force_error}")
    return force_error


@dataclass
class _PeriodicCell:
    mesh: _PeriodicMesh
    volumes: NDArray[np.float64]
    gradients: NDArray[np.float64]
    classes: NDArray[np.intp]
    assembly: fd.Assembly
    center: int

    @classmethod
    def from_mesh(cls, mesh: _PeriodicMesh) -> _PeriodicCell:
        classes = pair_classes(mesh)
        volumes, gradients = element_geometry(mesh)
        fd.ModelingSpace("3D")
        fem_mesh = fd.Mesh(mesh.points.copy(), mesh.tetrahedra.copy(), "tet4")
        assembly = fd.Assembly.create(
            fd.weakform.StressEquilibrium(fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)),
            fem_mesh,
            n_elm_gp=1,
        )
        center = int(np.argmin(np.linalg.norm(mesh.points - 0.5, axis=1)))
        return cls(mesh, volumes, gradients, classes, assembly, center)

    def strain_response(
        self,
        case: int,
        engineering: NDArray[np.float64],
        solver: AuditedCG,
        stress_probe_cell: int | None,
    ) -> _StrainResponse:
        problem, displacement = _prescribed_strain(self.assembly, self.center, engineering, solver)
        jump_error = _periodic_jump(self.mesh, displacement, engineering)
        gradient = np.einsum("tia,tib->tab", displacement[self.mesh.tetrahedra], self.gradients)
        strain = (gradient + np.swapaxes(gradient, -1, -2)) / 2
        stress = 2 * strain + np.trace(strain, axis1=1, axis2=2)[:, None, None] * np.eye(3)
        # The reference cell has unit volume. Integrate stress over its solid part.
        mean_stress = np.einsum("t,tij->ij", self.volumes, stress)
        reaction = problem.get_ext_forces("MeanStrain").ravel().copy()
        checks = {
            "load": COMPONENTS[case],
            "periodic_jump_error": jump_error,
            "periodic_force_balance_relative_error": _force_balance(problem, self.classes),
            **solver.records[-1],
        }
        if stress_probe_cell is not None:
            checks["probe_sigma_xx"] = float(stress[stress_probe_cell, 0, 0])
        return _StrainResponse(
            strain,
            stress,
            mean_stress[[0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]],
            reaction,
            checks,
        )


def _stiffness_checks(
    volumes: NDArray[np.float64],
    responses: list[_StrainResponse],
    stiffness: NDArray[np.float64],
) -> dict[str, float | list[float]]:
    reactions = np.column_stack([r.reaction for r in responses])
    energy_matrix = np.einsum(
        "t,itkl,jtkl->ij",
        volumes,
        [r.strain for r in responses],
        [r.stress for r in responses],
    )
    hill_mandel = float(np.linalg.norm(energy_matrix - stiffness) / np.linalg.norm(stiffness))
    reaction_error = float(np.linalg.norm(reactions - stiffness) / np.linalg.norm(stiffness))
    symmetry = float(np.linalg.norm(stiffness - stiffness.T) / np.linalg.norm(stiffness))
    require(
        max(hill_mandel, reaction_error, symmetry) < HOMOGENIZATION_TOLERANCE,
        "Energy, reaction or symmetry check failed",
    )
    eigenvalues = np.linalg.eigvalsh((stiffness + stiffness.T) / 2)
    require(np.all(eigenvalues > 0), "Effective stiffness must be positive definite")
    # Uniform-strain (Voigt) upper bound using the measured solid volume.
    voigt_gap = np.linalg.eigvalsh(volumes.sum() * MATERIAL - (stiffness + stiffness.T) / 2)
    require(
        voigt_gap.min() > -HOMOGENIZATION_TOLERANCE,
        "Effective stiffness exceeds the Voigt bound",
    )
    return {
        "stiffness_eigenvalues": eigenvalues.tolist(),
        "hill_mandel_relative_error": hill_mandel,
        "reaction_stress_relative_difference": reaction_error,
        "stiffness_symmetry_relative_error": symmetry,
    }


def homogenize(
    mesh: _PeriodicMesh,
    *,
    stress_probe_cell: int | None = None,
    solver_factory: Callable[[], AuditedCG] = AuditedCG,
) -> dict[str, object]:
    """Run strain- and stress-controlled fedoo homogenization and cross-check.

    Returns:
        Effective stiffness, directional constants and verification diagnostics.
    """
    start = time.perf_counter()
    cell = _PeriodicCell.from_mesh(mesh)
    solver = solver_factory()
    responses = [
        cell.strain_response(case, engineering, solver, stress_probe_cell)
        for case, engineering in enumerate(np.eye(6))
    ]
    stiffness = np.column_stack([r.mean_stress for r in responses])
    checks = _stiffness_checks(cell.volumes, responses, stiffness)
    automatic_solver = solver_factory()
    automatic = fd.homogen.get_homogenized_stiffness(
        cell.assembly,
        meshperio=True,
        solver=automatic_solver,
        rigid_body_constraint="pin",
    )
    automatic_error = float(np.linalg.norm(automatic - stiffness) / np.linalg.norm(stiffness))
    require(
        automatic_error < HOMOGENIZATION_TOLERANCE,
        f"Stress-control discrepancy {automatic_error}",
    )
    compliance = np.linalg.inv(stiffness)
    return {
        "tetrahedra": len(mesh.tetrahedra),
        "points": len(mesh.points),
        "solid_volume_fraction": float(cell.volumes.sum()),
        "minimum_quality": mesh.diagnostics["minimum_mmg_quality"],
        "periodic_pair_counts": [len(p) for p in mesh.periodic_pairs],
        "stiffness_engineering_voigt": stiffness.tolist(),
        **checks,
        "young_modulus_xyz": (1 / np.diag(compliance)[:3]).tolist(),
        "poisson_xy_xz_yz": [
            -float(compliance[j, i] / compliance[i, i]) for i, j in [(0, 1), (0, 2), (1, 2)]
        ],
        "automatic_homogenization_relative_difference": automatic_error,
        "strain_controlled_loads": [r.checks for r in responses],
        "automatic_stress_controlled_solves": automatic_solver.records,
        "seconds": time.perf_counter() - start,
    }


def specimen(name: str, cells: int) -> meshers.Mesh:
    """Generate a full solid or gyroid sheet in a periodic unit cube.

    Returns:
        Meshers mesh with exact translational pairing on all three axes.
    """
    if name == "solid":
        return meshers.generate(
            lambda x, _y, _z: x * 0 - 1,
            cells=cells,
            periodic=(True, True, True),
            optimize_passes=0,
        )
    return meshers.generate(
        "gyroid",
        band=(-0.55, 0.55),
        cells=cells,
        periodic=(True, True, True),
        geometry_tolerance=0.04,
        optimize_passes=4,
    )


def main() -> None:
    """Verify homogeneous recovery and measure TPMS refinement behavior."""
    extension = Path(importlib.import_module("meshers._meshers").__file__)
    report = {
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            **{n: version(n) for n in ["fedoo", "meshers", "numpy", "scipy"]},
        },
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks_sha256": hashlib.sha256(
            Path(__file__).with_name("solver_checks.py").read_bytes(),
        ).hexdigest(),
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
                require(error < SOLID_REFERENCE_RTOL, f"Solid reference discrepancy {error}")
                result["analytic_stiffness_relative_error"] = error
            if previous is not None:
                result["stiffness_relative_change"] = float(
                    np.linalg.norm(tensor - previous) / np.linalg.norm(tensor),
                )
            previous = tensor
            report["cases"].append(result)
            write_json_line(
                {
                    k: v
                    for k, v in result.items()
                    if k not in {"strain_controlled_loads", "automatic_stress_controlled_solves"}
                },
            )
            (ROOT / "docs/fedoo-homogenization-results.json").write_text(
                json.dumps(report, indent=2) + "\n",
            )


if __name__ == "__main__":
    main()
