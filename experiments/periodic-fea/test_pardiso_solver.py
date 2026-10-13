"""Check the real optional solver, including exact reuse and invalidation."""

import numpy as np
import pytest
from scipy.sparse import csr_matrix

pytest.importorskip("pypardiso")
from pardiso_solver import AuditedPardiso


def test_direct_factorization_reuses_only_an_identical_matrix():
    dense = np.array([[4., -1., 0.], [-1., 4., -1.], [0., -1., 3.]])
    matrix = csr_matrix(dense)
    solver = AuditedPardiso()
    try:
        for rhs in [np.array([1., 2., 3.]), np.array([-2., 1., 4.])]:
            np.testing.assert_allclose(solver(matrix.copy(), rhs), np.linalg.solve(dense, rhs), rtol=1e-12)
        assert solver.factorizations == 1
        changed = matrix.copy()
        changed[0, 0] += 1
        rhs = np.array([3., -1., 2.])
        np.testing.assert_allclose(solver(changed, rhs), np.linalg.solve(changed.toarray(), rhs), rtol=1e-12)
        assert solver.factorizations == 2
        assert max(record["relative_residual"] for record in solver.records) < 1e-12
        solver.free_memory()
        np.testing.assert_allclose(solver(matrix, rhs), np.linalg.solve(dense, rhs), rtol=1e-12)
        assert solver.factorizations == 3
    finally:
        solver.free_memory()
