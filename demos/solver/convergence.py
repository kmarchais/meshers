# Copyright (c) 2026 Kevin Marchais

"""P1 tetrahedral elasticity verification on actual Meshers meshes.

Manufactured displacement (sin(pi*x)sin(pi*y)sin(pi*z), 0, 0), lambda=mu=1.
Dirichlet data on every boundary node. This tests discretization on each
polyhedral mesh, not geometry-error convergence or periodic homogenization.
"""

import hashlib
import importlib
import json
import platform
import time
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import meshers
import numpy as np
from numpy.typing import NDArray
from reporting import write_json_line
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import cg
from solver_checks import (
    AFFINE_DISPLACEMENT_ATOL,
    AFFINE_ENERGY_ATOL,
    MATRIX_SYMMETRY_ATOL,
    MAX_RELATIVE_RESIDUAL,
    require,
)

ROOT = Path(__file__).resolve().parents[2]
QUADRATURE = np.full((4, 4), (5 - np.sqrt(5)) / 20)
np.fill_diagonal(QUADRATURE, (5 + 3 * np.sqrt(5)) / 20)


def exact(
    points: NDArray[np.float64],
    *,
    affine: bool = False,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Return displacement, displacement gradient and body force."""
    if affine:
        matrix = np.array([[0.2, 0.1, -0.05], [0.03, -0.1, 0.04], [0.02, 0.06, 0.12]])
        return (
            points @ matrix.T + [0.1, -0.2, 0.3],
            np.broadcast_to(matrix, (*points.shape[:-1], 3, 3)),
            np.zeros_like(points),
        )
    s, c = np.sin(np.pi * points), np.cos(np.pi * points)
    value = np.prod(s, axis=-1)
    u = np.zeros_like(points)
    u[..., 0] = value
    gradient = np.zeros((*points.shape[:-1], 3, 3))
    for a in range(3):
        other = [i for i in range(3) if i != a]
        gradient[..., 0, a] = np.pi * c[..., a] * s[..., other[0]] * s[..., other[1]]
    f = np.zeros_like(points)
    f[..., 0] = 5 * np.pi**2 * value
    f[..., 1] = -2 * np.pi**2 * c[..., 0] * c[..., 1] * s[..., 2]
    f[..., 2] = -2 * np.pi**2 * c[..., 0] * s[..., 1] * c[..., 2]
    return u, gradient, f


@dataclass
class ElementQuadrature:
    """P1 shape gradients and manufactured fields at the four integration points."""

    mesh: meshers.Mesh
    volumes: NDArray[np.float64]
    gradients: NDArray[np.float64]
    values: NDArray[np.float64]
    exact_gradients: NDArray[np.float64]
    force: NDArray[np.float64]
    max_edge: float

    @classmethod
    def from_mesh(cls, mesh: meshers.Mesh, *, affine: bool) -> "ElementQuadrature":
        """Evaluate the manufactured solution on the physical tetrahedra.

        Returns:
            Element volumes, shape gradients and reference fields at quadrature points.
        """
        p = mesh.points[mesh.tetrahedra]
        jac = p[:, 1:] - p[:, :1]
        volumes = np.linalg.det(jac) / 6
        require(np.all(volumes > 0), "Mesh contains a nonpositive tetrahedron volume")
        gradients = np.empty((len(p), 4, 3))
        gradients[:, 1:] = np.swapaxes(np.linalg.inv(jac), 1, 2)
        gradients[:, 0] = -gradients[:, 1:].sum(axis=1)
        q = np.einsum("qi,tia->tqa", QUADRATURE, p)
        values, exact_gradients, force = exact(q, affine=affine)
        max_edge = max(
            np.max(np.linalg.norm(p[:, i] - p[:, j], axis=1))
            for i in range(4)
            for j in range(i + 1, 4)
        )
        return cls(mesh, volumes, gradients, values, exact_gradients, force, float(max_edge))

    def assemble(self) -> tuple[csr_matrix, NDArray[np.float64]]:
        """Assemble the reference system in node-major order.

        Returns:
            Sparse stiffness matrix and integrated body force.
        """
        cells = self.mesh.tetrahedra
        dot = np.einsum("tia,tja->tij", self.gradients, self.gradients)
        # Node/component indices (i,a,j,b), lambda=mu=1.
        stiffness = (
            np.einsum("tia,tjb->tiajb", self.gradients, self.gradients)
            + np.einsum("tib,tja->tiajb", self.gradients, self.gradients)
            + np.einsum("tij,ab->tiajb", dot, np.eye(3))
        ) * self.volumes[:, None, None, None, None]
        dofs = (cells[..., None] * 3 + np.arange(3)).reshape(-1, 12)
        rows = np.broadcast_to(dofs[:, :, None], (len(cells), 12, 12)).ravel()
        cols = np.broadcast_to(dofs[:, None, :], (len(cells), 12, 12)).ravel()
        size = 3 * len(self.mesh.points)
        matrix = coo_matrix((stiffness.reshape(-1), (rows, cols)), shape=(size, size)).tocsr()
        local_load = np.einsum("qi,tqa,t->tia", QUADRATURE, self.force, self.volumes) / 4
        rhs = np.zeros(size)
        np.add.at(rhs, dofs.ravel(), local_load.ravel())
        return matrix, rhs

    def errors(self, displacement: NDArray[np.float64]) -> dict[str, float]:
        """Integrate errors from the nodal displacement field.

        Returns:
            Displacement and energy norms, minimum volume and maximum edge length.
        """
        local = displacement[self.mesh.tetrahedra]
        discrete_values = np.einsum("qi,tia->tqa", QUADRATURE, local)
        discrete_gradient = np.einsum("tia,tib->tab", local, self.gradients)
        delta = discrete_gradient[:, None] - self.exact_gradients
        strain = (delta + np.swapaxes(delta, -1, -2)) / 2
        energy = np.trace(strain, axis1=-2, axis2=-1) ** 2 + 2 * np.sum(strain**2, axis=(-2, -1))
        l2 = np.sqrt(
            np.sum(self.volumes[:, None] * np.sum((discrete_values - self.values) ** 2, axis=-1))
            / 4,
        )
        energy_error = np.sqrt(np.sum(self.volumes[:, None] * energy) / 4)
        return {
            "minimum_signed_volume": float(np.min(self.volumes)),
            "maximum_edge_length": self.max_edge,
            "l2_displacement_error": float(l2),
            "energy_error": float(energy_error),
        }


def _solve_dirichlet(
    mesh: meshers.Mesh,
    matrix: csr_matrix,
    rhs: NDArray[np.float64],
    *,
    affine: bool,
) -> tuple[NDArray[np.float64], int, int, float]:
    boundary = np.unique(mesh.surface)
    locked = (boundary[:, None] * 3 + np.arange(3)).ravel()
    free = np.ones(len(rhs), dtype=bool)
    free[locked] = False
    require(np.any(free), "Elasticity system has no free degrees of freedom")
    displacement = exact(mesh.points, affine=affine)[0].ravel().copy()
    reduced = matrix[free][:, free]
    load = rhs[free] - matrix[free][:, locked] @ displacement[locked]
    solution, iterations, residual = _conjugate_gradient(reduced, load)
    displacement[free] = solution
    return displacement.reshape(-1, 3), int(free.sum()), iterations, residual


def _conjugate_gradient(
    reduced: csr_matrix,
    load: NDArray[np.float64],
) -> tuple[NDArray[np.float64], int, float]:
    symmetry = reduced - reduced.T
    require(
        np.max(np.abs(symmetry.data), initial=0) < MATRIX_SYMMETRY_ATOL,
        "Reduced stiffness matrix is asymmetric",
    )
    require(np.all(reduced.diagonal() > 0), "Reduced stiffness matrix has a nonpositive diagonal")
    iterations = 0

    def count(_: NDArray[np.float64]) -> None:
        nonlocal iterations
        iterations += 1

    solution, info = cg(
        reduced,
        load,
        M=diags(1 / reduced.diagonal()),
        rtol=1e-10,
        atol=1e-12,
        maxiter=10000,
        callback=count,
    )
    residual = float(np.linalg.norm(reduced @ solution - load) / max(np.linalg.norm(load), 1e-30))
    require(info == 0, f"CG returned status {info}")
    require(residual < MAX_RELATIVE_RESIDUAL, f"Relative residual {residual} exceeds the limit")
    return solution, iterations, residual


def solve(mesh: meshers.Mesh, *, affine: bool = False) -> dict[str, float]:
    """Assemble independent linear elasticity and measure solution errors.

    Returns:
        Mesh size, quality, displacement/energy errors and solver measurements.
    """
    start = time.perf_counter()
    quadrature = ElementQuadrature.from_mesh(mesh, affine=affine)
    matrix, rhs = quadrature.assemble()
    displacement, free_dofs, iterations, residual = _solve_dirichlet(
        mesh,
        matrix,
        rhs,
        affine=affine,
    )
    errors = quadrature.errors(displacement)
    l2, energy_error = errors["l2_displacement_error"], errors["energy_error"]
    if affine:
        require(l2 < AFFINE_DISPLACEMENT_ATOL, f"Affine displacement error {l2} exceeds the limit")
        require(
            energy_error < AFFINE_ENERGY_ATOL,
            f"Affine energy error {energy_error} exceeds the limit",
        )
    return {
        "tetrahedra": len(mesh.tetrahedra),
        "points": len(mesh.points),
        "free_dofs": free_dofs,
        "minimum_quality": mesh.diagnostics["minimum_mmg_quality"],
        "maximum_edge_length": quadrature.max_edge,
        "l2_displacement_error": float(l2),
        "energy_error": float(energy_error),
        "cg_iterations": iterations,
        "relative_residual": residual,
        "seconds": time.perf_counter() - start,
    }


def specimen(name: str, n: int) -> meshers.Mesh:
    """Build a fixed geometry at the requested resolution through the public API.

    Returns:
        Physical tetrahedral mesh of the requested fixed specimen.
    """
    if name == "cube":
        return meshers.generate(lambda x, _y, _z: x * 0 - 1, cells=n, optimize_passes=0)

    def gyroid(
        x: NDArray[np.float64],
        y: NDArray[np.float64],
        z: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        a, b, c = (
            2 * np.pi * (x / 0.75 + 0.125),
            2 * np.pi * (y / 0.75 + 0.125),
            2 * np.pi * (z + 0.125),
        )
        return np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)

    constraints = {
        "upper": lambda x, y, z: gyroid(x, y, z) - 0.55,
        "lower": lambda x, y, z: -gyroid(x, y, z) - 0.55,
        "wall": lambda x, y, _z: np.sqrt(x * x + y * y) - 0.75,
    }
    mapping = None
    if name == "torus":

        def mapping(
            x: NDArray[np.float64],
            y: NDArray[np.float64],
            z: NDArray[np.float64],
        ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
            return (3 + x) * np.cos(z / 3), (3 + x) * np.sin(z / 3), -y

    return meshers.generate_intersection(
        constraints,
        bounds=(-0.8, 0.8, -0.8, 0.8, 0, 1),
        cells=n,
        coordinate_map=mapping,
        snap=0.2,
        optimize_passes=8,
    )


def main() -> None:
    """Run affine patch tests and three refinement levels for each geometry."""
    extension = Path(importlib.import_module("meshers._meshers").__file__)
    report = {
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "meshers": version("meshers"),
            "numpy": np.__version__,
            "scipy": version("scipy"),
        },
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks_sha256": hashlib.sha256(
            Path(__file__).with_name("solver_checks.py").read_bytes(),
        ).hexdigest(),
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "method": "P1 linear elasticity, lambda=mu=1, 4-point tetrahedron quadrature",
        "boundary_conditions": "Manufactured displacement on every boundary node",
        "solver": "SciPy CG, diagonal preconditioner, rtol=1e-10, atol=1e-12",
        "limitations": [
            "No comparison against MMG or commercial FE software",
            "No periodic homogenization, nonlinear response or geometry-error convergence",
        ],
        "cases": [],
    }
    for name in ["cube", "cylinder", "torus"]:
        levels = [6, 9, 12] if name == "cube" else [12, 18, 24]
        results = []
        patch = None
        for n in levels:
            mesh = specimen(name, n)
            if patch is None:
                patch = solve(mesh, affine=True)
            result = {"cells": n, **solve(mesh)}
            if results:
                ratio = result["maximum_edge_length"] / results[-1]["maximum_edge_length"]
                for metric in ["l2_displacement_error", "energy_error"]:
                    result[metric + "_rate"] = float(
                        np.log(result[metric] / results[-1][metric]) / np.log(ratio),
                    )
            results.append(result)
            write_json_line({"geometry": name, **result})
        require(
            results[-1]["l2_displacement_error"] < results[0]["l2_displacement_error"],
            f"{name}: l2_displacement_error did not decrease under refinement",
        )
        require(
            results[-1]["energy_error"] < results[0]["energy_error"],
            f"{name}: energy_error did not decrease under refinement",
        )
        report["cases"].append({"geometry": name, "affine_patch": patch, "refinement": results})
        (ROOT / "docs/solver-convergence-results.json").write_text(
            json.dumps(report, indent=2) + "\n",
        )


if __name__ == "__main__":
    main()
