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
from importlib.metadata import version
from pathlib import Path

import meshers
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import cg

ROOT = Path(__file__).resolve().parents[2]
QUADRATURE = np.full((4, 4), (5 - np.sqrt(5)) / 20)
np.fill_diagonal(QUADRATURE, (5 + 3 * np.sqrt(5)) / 20)


def exact(points, *, affine=False):
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


def solve(mesh, *, affine=False):  # ruff: ignore[too-many-locals, too-many-statements] - independent element assembly and error integration.
    """Assemble independent linear elasticity and measure solution errors.

    Returns:
        Mesh size, quality, displacement/energy errors and solver measurements.
    """
    start = time.perf_counter()
    cells = mesh.tetrahedra
    p = mesh.points[cells]
    jac = p[:, 1:] - p[:, :1]
    volumes = np.linalg.det(jac) / 6
    assert np.all(volumes > 0)
    gradients = np.empty((len(cells), 4, 3))
    gradients[:, 1:] = np.swapaxes(np.linalg.inv(jac), 1, 2)
    gradients[:, 0] = -gradients[:, 1:].sum(axis=1)
    dot = np.einsum("tia,tja->tij", gradients, gradients)
    # Node/component indices (i,a,j,b), lambda=mu=1.
    stiffness = (
        np.einsum("tia,tjb->tiajb", gradients, gradients)
        + np.einsum("tib,tja->tiajb", gradients, gradients)
        + np.einsum("tij,ab->tiajb", dot, np.eye(3))
    ) * volumes[:, None, None, None, None]
    dofs = (cells[..., None] * 3 + np.arange(3)).reshape(-1, 12)
    rows = np.broadcast_to(dofs[:, :, None], (len(cells), 12, 12)).ravel()
    cols = np.broadcast_to(dofs[:, None, :], (len(cells), 12, 12)).ravel()
    size = 3 * len(mesh.points)
    matrix = coo_matrix((stiffness.reshape(-1), (rows, cols)), shape=(size, size)).tocsr()
    q = np.einsum("qi,tia->tqa", QUADRATURE, p)
    values, analytic_gradient, force = exact(q, affine=affine)
    local_load = np.einsum("qi,tqa,t->tia", QUADRATURE, force, volumes) / 4
    rhs = np.zeros(size)
    np.add.at(rhs, dofs.ravel(), local_load.ravel())
    boundary = np.unique(mesh.surface)
    locked = (boundary[:, None] * 3 + np.arange(3)).ravel()
    free = np.ones(size, dtype=bool)
    free[locked] = False
    assert np.any(free)
    displacement = exact(mesh.points, affine=affine)[0].ravel().copy()
    reduced = matrix[free][:, free]
    load = rhs[free] - matrix[free][:, locked] @ displacement[locked]
    symmetry = reduced - reduced.T
    assert np.max(np.abs(symmetry.data), initial=0) < 1e-10
    assert np.all(reduced.diagonal() > 0)
    iterations = 0

    def count(_):
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
    assert info == 0, info
    assert residual < 2e-9, residual
    displacement[free] = solution
    local = displacement.reshape(-1, 3)[cells]
    discrete_values = np.einsum("qi,tia->tqa", QUADRATURE, local)
    discrete_gradient = np.einsum("tia,tib->tab", local, gradients)
    delta = discrete_gradient[:, None] - analytic_gradient
    strain = (delta + np.swapaxes(delta, -1, -2)) / 2
    energy = np.trace(strain, axis1=-2, axis2=-1) ** 2 + 2 * np.sum(strain**2, axis=(-2, -1))
    l2 = np.sqrt(np.sum(volumes[:, None] * np.sum((discrete_values - values) ** 2, axis=-1)) / 4)
    energy_error = np.sqrt(np.sum(volumes[:, None] * energy) / 4)
    if affine:
        assert l2 < 1e-8, l2
        assert energy_error < 1e-7, energy_error
    edges = np.concatenate([
        np.linalg.norm(p[:, i] - p[:, j], axis=1) for i in range(4) for j in range(i + 1, 4)
    ])
    return {
        "tetrahedra": len(cells),
        "points": len(mesh.points),
        "free_dofs": int(free.sum()),
        "minimum_quality": mesh.diagnostics["minimum_mmg_quality"],
        "maximum_edge_length": float(edges.max()),
        "l2_displacement_error": float(l2),
        "energy_error": float(energy_error),
        "cg_iterations": iterations,
        "relative_residual": residual,
        "seconds": time.perf_counter() - start,
    }


def specimen(name, n):
    """Build a fixed geometry at the requested resolution through the public API.

    Returns:
        Physical tetrahedral mesh of the requested fixed specimen.
    """
    if name == "cube":
        return meshers.generate(lambda x, _y, _z: x * 0 - 1, cells=n, optimize_passes=0)

    def gyroid(x, y, z):
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

        def mapping(x, y, z):
            return (3 + x) * np.cos(z / 3), (3 + x) * np.sin(z / 3), -y

    return meshers.generate_intersection(
        constraints,
        bounds=(-0.8, 0.8, -0.8, 0.8, 0, 1),
        cells=n,
        coordinate_map=mapping,
        snap=0.2,
        optimize_passes=8,
    )


def main():
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
                        np.log(result[metric] / results[-1][metric]) / np.log(ratio)
                    )
            results.append(result)
            print(json.dumps({"geometry": name, **result}), flush=True)
        assert results[-1]["l2_displacement_error"] < results[0]["l2_displacement_error"]
        assert results[-1]["energy_error"] < results[0]["energy_error"]
        report["cases"].append({"geometry": name, "affine_patch": patch, "refinement": results})
        (ROOT / "docs/solver-convergence-results.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
