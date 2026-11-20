# Copyright (c) 2026 Kevin Marchais

"""Verify Meshers tetrahedra with fedoo's assembly, loads, BCs and linear solve.

Uses the same manufactured solution and specimens as convergence.py. Only the
geometry and analytic error integration are shared; fedoo assembles the system.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import platform
import time
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING

import fedoo as fd
import numpy as np
from convergence import ROOT, ElementQuadrature, exact, specimen
from convergence import solve as reference_solve
from reporting import write_json_line
from solver_checks import (
    AFFINE_DISPLACEMENT_ATOL,
    AFFINE_ENERGY_ATOL,
    MATRIX_SYMMETRY_ATOL,
    MAX_RELATIVE_RESIDUAL,
    require,
)

if TYPE_CHECKING:
    import meshers
    from numpy.typing import NDArray

PRESCRIBED_DISPLACEMENT_ATOL = 1e-12
MAX_REFERENCE_DISAGREEMENT = 1e-6
MIN_DISPLACEMENT_RATE = 1.7
MIN_ENERGY_RATE = 0.85


def errors(
    mesh: meshers.Mesh,
    displacement: NDArray[np.float64],
    *,
    affine: bool = False,
) -> dict[str, float]:
    """Integrate solution errors independently of fedoo.

    Returns:
        Mesh size and displacement and energy error norms.
    """
    return ElementQuadrature.from_mesh(mesh, affine=affine).errors(displacement)


def _elasticity_problem(mesh: meshers.Mesh, *, affine: bool) -> fd.problem.Linear:
    fd.ModelingSpace("3D")
    fem_mesh = fd.Mesh(mesh.points.copy(), mesh.tetrahedra.copy(), "tet4")
    # E=2.5, nu=0.25 corresponds to lambda=mu=1.
    material = fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)
    equilibrium = fd.weakform.StressEquilibrium(material)
    force = exact(fem_mesh.gausspoint_coordinates(4), affine=affine)[2]
    load = fd.weakform.DistributedLoad(list(force.T))
    assembly = fd.Assembly.create(equilibrium + load, fem_mesh, n_elm_gp=4)
    return fd.problem.Linear(assembly)


def _audit_free_equations(
    problem: fd.problem.Linear,
    displacement: NDArray[np.float64],
    boundary: NDArray[np.intp],
) -> tuple[int, float]:
    # fedoo stores component-major DOFs; audit only unconstrained equations.
    locked = (boundary[None, :] + len(displacement) * np.arange(3)[:, None]).ravel()
    free = np.ones(3 * len(displacement), dtype=bool)
    free[locked] = False
    require(np.any(free), "Elasticity system has no free degrees of freedom")
    matrix = problem.get_A()
    vector = problem.get_B() + problem.get_D()
    reduced = matrix[free][:, free]
    rhs = vector[free] - matrix[free][:, locked] @ displacement.T.ravel()[locked]
    symmetry = reduced - reduced.T
    require(
        np.max(np.abs(symmetry.data), initial=0) < MATRIX_SYMMETRY_ATOL,
        "Reduced stiffness matrix is asymmetric",
    )
    require(np.all(reduced.diagonal() > 0), "Reduced stiffness matrix has a nonpositive diagonal")
    residual = float(
        np.linalg.norm(reduced @ displacement.T.ravel()[free] - rhs)
        / max(np.linalg.norm(rhs), 1e-30),
    )
    require(residual < MAX_RELATIVE_RESIDUAL, f"Relative residual {residual} exceeds the limit")
    return int(free.sum()), residual


def solve(mesh: meshers.Mesh, *, affine: bool = False) -> dict[str, float]:
    """Assemble and solve exclusively through fedoo, then audit its result.

    Returns:
        Error norms, mesh statistics, boundary and algebraic residual checks.
    """
    start = time.perf_counter()
    problem = _elasticity_problem(mesh, affine=affine)
    boundary = np.unique(mesh.surface)
    boundary_values = exact(mesh.points[boundary], affine=affine)[0]
    for component, variable in enumerate(["DispX", "DispY", "DispZ"]):
        problem.bc.add("Dirichlet", boundary, variable, boundary_values[:, component])
    iterations = 0

    def count(_: NDArray[np.float64]) -> None:
        nonlocal iterations
        iterations += 1

    problem.set_solver("cg", precond=True, rtol=1e-10, atol=1e-12, maxiter=10000, callback=count)
    problem.solve()
    displacement = problem.get_disp().T
    require(
        displacement.shape == mesh.points.shape,
        "Displacement shape differs from the mesh coordinates",
    )
    require(np.all(np.isfinite(displacement)), "Displacement contains nonfinite values")
    boundary_error = float(np.max(np.abs(displacement[boundary] - boundary_values)))
    require(
        boundary_error < PRESCRIBED_DISPLACEMENT_ATOL,
        f"Boundary displacement error {boundary_error} exceeds the limit",
    )
    free_dofs, residual = _audit_free_equations(problem, displacement, boundary)
    result = errors(mesh, displacement, affine=affine)
    if affine:
        require(
            result["l2_displacement_error"] < AFFINE_DISPLACEMENT_ATOL,
            f"Affine displacement error exceeds the limit: {result}",
        )
        require(
            result["energy_error"] < AFFINE_ENERGY_ATOL,
            f"Affine energy error exceeds the limit: {result}",
        )
    return {
        "tetrahedra": len(mesh.tetrahedra),
        "points": len(mesh.points),
        "free_dofs": free_dofs,
        "minimum_quality": mesh.diagnostics["minimum_mmg_quality"],
        **result,
        "cg_iterations": iterations,
        "relative_residual": residual,
        "maximum_boundary_error": boundary_error,
        "seconds": time.perf_counter() - start,
    }


def main() -> None:
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
            Path(__file__).with_name("convergence.py").read_bytes(),
        ).hexdigest(),
        "checks_sha256": hashlib.sha256(
            Path(__file__).with_name("solver_checks.py").read_bytes(),
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
                require(
                    delta < MAX_REFERENCE_DISAGREEMENT,
                    f"{name} at {n} cells: {metric} differs from the reference by {delta:.6g}",
                )
                result[metric + "_reference_relative_difference"] = delta
                if results:
                    ratio = result["maximum_edge_length"] / results[-1]["maximum_edge_length"]
                    result[metric + "_rate"] = float(
                        np.log(result[metric] / results[-1][metric]) / np.log(ratio),
                    )
            results.append(result)
            write_json_line({"geometry": name, **result})
        require(
            results[-1]["l2_displacement_error_rate"] > MIN_DISPLACEMENT_RATE,
            f"{name}: l2_displacement_error convergence rate must exceed 1.7",
        )
        require(
            results[-1]["energy_error_rate"] > MIN_ENERGY_RATE,
            f"{name}: energy_error convergence rate must exceed 0.85",
        )
        report["cases"].append({"geometry": name, "affine_patch": patch, "refinement": results})
        (ROOT / "docs/fedoo-convergence-results.json").write_text(
            json.dumps(report, indent=2) + "\n",
        )


if __name__ == "__main__":
    main()
