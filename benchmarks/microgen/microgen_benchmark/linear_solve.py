"""One common elasticity linear solve, with an independent direct option."""

import time

import numpy as np
from scipy.sparse.linalg import cg, splu


def solve_system(matrix, rhs, coordinates, free, method, ordering, checkpoint=None):
    if method == "direct":
        factor = splu(matrix, permc_spec=ordering)
        return factor.solve(rhs), {
            "linear_solver": method,
            "solver_ordering": ordering,
            "factor_nonzeros": factor.L.nnz + factor.U.nnz,
        }
    if method != "amg":
        raise ValueError(f"Unknown linear solver {method}")
    import pyamg

    modes = np.zeros((len(coordinates), 3, 6))
    modes[:, :, :3] = np.eye(3)
    for axis in range(3):
        modes[:, :, 3 + axis] = np.cross(np.eye(3)[axis], coordinates - coordinates.mean(axis=0))
    hierarchy = pyamg.smoothed_aggregation_solver(
        matrix.tocsr(),
        B=modes.reshape(-1, 6)[free],
        smooth="energy",
        symmetry="symmetric",
        max_coarse=200,
        coarse_solver="splu",
    )
    preconditioner = hierarchy.aspreconditioner()
    result = np.empty_like(rhs)
    counts = []
    residuals = []
    for column in range(rhs.shape[1]):
        count = 0
        last_report = time.monotonic()

        def iteration(approximate, load=column):
            nonlocal count, last_report
            count += 1
            if checkpoint and time.monotonic() - last_report >= 2:
                residual = np.linalg.norm(matrix @ approximate - rhs[:, load]) / max(
                    np.linalg.norm(rhs[:, load]), 1e-30
                )
                checkpoint(
                    fea_rhs=load + 1,
                    fea_cg_iterations=count,
                    fea_cg_relative_residual=float(residual),
                )
                last_report = time.monotonic()

        result[:, column], info = cg(
            matrix,
            rhs[:, column],
            M=preconditioner,
            rtol=1e-11,
            atol=0,
            maxiter=1000,
            callback=iteration,
        )
        if info:
            raise ValueError(f"AMG load {column} did not converge in 1000 iterations: {info}")
        counts.append(count)
        residual = float(
            np.linalg.norm(matrix @ result[:, column] - rhs[:, column])
            / max(np.linalg.norm(rhs[:, column]), 1e-30)
        )
        residuals.append(residual)
        if checkpoint:
            checkpoint(
                fea_rhs=column + 1, fea_cg_iterations=count, fea_cg_relative_residual=residual
            )
    return result, {
        "linear_solver": method,
        "amg_iterations": counts,
        "amg_final_relative_residuals": residuals,
        "amg_operator_complexity": hierarchy.operator_complexity(),
        "amg_relative_tolerance": 1e-11,
        "amg_max_iterations": 1000,
    }
