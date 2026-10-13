"""Optional direct-solver experiment; the retained CG baseline stays unchanged."""

import hashlib
import time

import numpy as np


class AuditedPardiso:
    """Reuse an exact reduced-matrix factorization and audit each solution."""

    def __init__(self):
        from pypardiso import PyPardisoSolver

        self.engine = PyPardisoSolver(mtype=11)
        self.records = []
        self.matrix_signature = None
        self.factorizations = 0

    def __call__(self, matrix, rhs, **_kwargs):
        started = time.perf_counter()
        assert np.all(matrix.diagonal() > 0)
        skew = matrix - matrix.T
        assert np.linalg.norm(skew.data) / np.linalg.norm(matrix.data) < 1e-12
        sparse = matrix.tocsr(copy=True)
        sparse.sum_duplicates()
        sparse.sort_indices()
        signature = hashlib.sha256(str(sparse.shape).encode())
        for array in [sparse.indptr, sparse.indices, sparse.data]:
            signature.update(array.dtype.str.encode())
            signature.update(memoryview(np.ascontiguousarray(array)))
        signature = signature.hexdigest()
        factor_seconds = 0.
        if signature != self.matrix_signature:
            factor_started = time.perf_counter()
            self.engine.factorize(sparse)
            factor_seconds = time.perf_counter() - factor_started
            self.matrix_signature = signature
            self.factorizations += 1
        solution = self.engine.solve(sparse, rhs)
        residual = float(np.linalg.norm(matrix @ solution - rhs) / max(np.linalg.norm(rhs), 1e-30))
        assert np.all(np.isfinite(solution))
        assert residual < 2e-9, residual
        self.records.append({
            "linear_solver": "pypardiso", "matrix_type": 11,
            "iterations": None, "relative_residual": residual, "free_dofs": len(rhs),
            "seconds": time.perf_counter() - started,
            "factorization_seconds": factor_seconds, "factorizations": self.factorizations,
            "matrix_sha256": signature,
        })
        return solution

    def free_memory(self):
        """Release the factorization before assembling the next formulation."""
        if self.matrix_signature is not None:
            self.engine.free_memory(everything=True)
            self.matrix_signature = None
