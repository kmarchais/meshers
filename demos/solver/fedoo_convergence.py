"""Verify Meshers tetrahedra with fedoo's assembly, loads, BCs and linear solve.

Uses the same manufactured solution and specimens as convergence.py. Only the
geometry and analytic error integration are shared; fedoo assembles the system.
"""

import hashlib
import importlib
import json
import platform
import time
from importlib.metadata import version
from pathlib import Path

import fedoo as fd
import numpy as np
from convergence import QUADRATURE, ROOT, exact, specimen
from convergence import solve as reference_solve


def errors(mesh, displacement, *, affine=False):  # ruff: ignore[too-many-locals] - explicit error integration.
    """Integrate displacement and energy errors independently of fedoo.

    Returns:
        Positive volumes, maximum edge length and the two error norms.
    """
    points = mesh.points[mesh.tetrahedra]
    jac = points[:, 1:] - points[:, :1]
    volumes = np.linalg.det(jac) / 6
    assert np.all(volumes > 0)
    gradients = np.empty((len(points), 4, 3))
    gradients[:, 1:] = np.swapaxes(np.linalg.inv(jac), 1, 2)
    gradients[:, 0] = -gradients[:, 1:].sum(axis=1)
    q = np.einsum("qi,tia->tqa", QUADRATURE, points)
    values, analytic_gradient, _ = exact(q, affine=affine)
    local = displacement[mesh.tetrahedra]
    discrete_values = np.einsum("qi,tia->tqa", QUADRATURE, local)
    discrete_gradient = np.einsum("tia,tib->tab", local, gradients)
    delta = discrete_gradient[:, None] - analytic_gradient
    strain = (delta + np.swapaxes(delta, -1, -2)) / 2
    energy = np.trace(strain, axis1=-2, axis2=-1) ** 2 + 2 * np.sum(strain**2, axis=(-2, -1))
    l2 = np.sqrt(np.sum(volumes[:, None] * np.sum((discrete_values - values) ** 2, axis=-1)) / 4)
    energy_error = np.sqrt(np.sum(volumes[:, None] * energy) / 4)
    h = max(
        np.linalg.norm(points[:, i] - points[:, j], axis=1).max()
        for i in range(4)
        for j in range(i + 1, 4)
    )
    return {
        "minimum_signed_volume": float(volumes.min()),
        "maximum_edge_length": float(h),
        "l2_displacement_error": float(l2),
        "energy_error": float(energy_error),
    }


def solve(mesh, *, affine=False):  # ruff: ignore[too-many-locals] - retain the complete fedoo setup in one example.
    """Assemble and solve exclusively through fedoo, then audit its result.

    Returns:
        Error norms, mesh statistics, boundary and algebraic residual checks.
    """
    start = time.perf_counter()
    fd.ModelingSpace("3D")
    fem_mesh = fd.Mesh(mesh.points.copy(), mesh.tetrahedra.copy(), "tet4")
    # E=2.5, nu=0.25 corresponds to lambda=mu=1.
    material = fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)
    equilibrium = fd.weakform.StressEquilibrium(material)
    force = exact(fem_mesh.gausspoint_coordinates(4), affine=affine)[2]
    load = fd.weakform.DistributedLoad(list(force.T))
    assembly = fd.Assembly.create(equilibrium + load, fem_mesh, n_elm_gp=4)
    problem = fd.problem.Linear(assembly)
    boundary = np.unique(mesh.surface)
    boundary_values = exact(mesh.points[boundary], affine=affine)[0]
    for component, variable in enumerate(["DispX", "DispY", "DispZ"]):
        problem.bc.add("Dirichlet", boundary, variable, boundary_values[:, component])
    iterations = 0

    def count(_):
        nonlocal iterations
        iterations += 1

    problem.set_solver("cg", precond=True, rtol=1e-10, atol=1e-12, maxiter=10000, callback=count)
    problem.solve()
    displacement = problem.get_disp().T
    assert displacement.shape == mesh.points.shape
    assert np.all(np.isfinite(displacement))
    boundary_error = float(np.max(np.abs(displacement[boundary] - boundary_values)))
    assert boundary_error < 1e-12, boundary_error
    # fedoo stores component-major DOFs; audit only unconstrained equations.
    locked = (boundary[None, :] + len(mesh.points) * np.arange(3)[:, None]).ravel()
    free = np.ones(3 * len(mesh.points), dtype=bool)
    free[locked] = False
    assert np.any(free)
    matrix = problem.get_A()
    vector = problem.get_B() + problem.get_D()
    reduced = matrix[free][:, free]
    rhs = vector[free] - matrix[free][:, locked] @ displacement.T.ravel()[locked]
    symmetry = reduced - reduced.T
    assert np.max(np.abs(symmetry.data), initial=0) < 1e-10
    assert np.all(reduced.diagonal() > 0)
    residual = float(
        np.linalg.norm(reduced @ displacement.T.ravel()[free] - rhs)
        / max(np.linalg.norm(rhs), 1e-30)
    )
    assert residual < 2e-9, residual
    result = errors(mesh, displacement, affine=affine)
    if affine:
        assert result["l2_displacement_error"] < 1e-8, result
        assert result["energy_error"] < 1e-7, result
    return {
        "tetrahedra": len(mesh.tetrahedra),
        "points": len(mesh.points),
        "free_dofs": int(free.sum()),
        "minimum_quality": mesh.diagnostics["minimum_mmg_quality"],
        **result,
        "cg_iterations": iterations,
        "relative_residual": residual,
        "maximum_boundary_error": boundary_error,
        "seconds": time.perf_counter() - start,
    }


def main():
    """Run fedoo patch tests and convergence, checking the independent assembly."""
    extension = Path(importlib.import_module("meshers._meshers").__file__)
    report = {
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            **{name: version(name) for name in ["meshers", "fedoo", "numpy", "scipy"]},
        },
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "shared_source_sha256": hashlib.sha256(
            Path(__file__).with_name("convergence.py").read_bytes()
        ).hexdigest(),
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "method": (
            "fedoo tet4, ElasticIsotrop(E=2.5, nu=0.25), "
            "StressEquilibrium + DistributedLoad, 4 Gauss points"
        ),
        "boundary_conditions": "fedoo Dirichlet: manufactured displacement on every boundary node",
        "solver": (
            "fedoo Linear with SciPy CG backend, diagonal preconditioner, rtol=1e-10, atol=1e-12"
        ),
        "limitations": [
            "No periodic homogenization, nonlinear response or geometry-error convergence",
            "Same analytic error integration as the independent study; same SciPy CG backend",
            "No MMG comparison; timings are not controlled performance benchmarks",
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
            reference = reference_solve(mesh)
            for metric in ["l2_displacement_error", "energy_error"]:
                delta = abs(result[metric] / reference[metric] - 1)
                assert delta < 1e-6, (name, n, metric, delta)
                result[metric + "_reference_relative_difference"] = delta
                if results:
                    ratio = result["maximum_edge_length"] / results[-1]["maximum_edge_length"]
                    result[metric + "_rate"] = float(
                        np.log(result[metric] / results[-1][metric]) / np.log(ratio)
                    )
            results.append(result)
            print(json.dumps({"geometry": name, **result}), flush=True)
        assert results[-1]["l2_displacement_error_rate"] > 1.7
        assert results[-1]["energy_error_rate"] > 0.85
        report["cases"].append({"geometry": name, "affine_patch": patch, "refinement": results})
        (ROOT / "docs/fedoo-convergence-results.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
