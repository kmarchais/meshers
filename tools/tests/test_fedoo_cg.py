# Copyright (c) 2026 Kevin Marchais

"""Check solver admission with a deliberately failed CG result."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
pytest.importorskip("fedoo")
import fedoo_homogenization as fea
from scipy.sparse import csr_matrix, eye


@pytest.mark.parametrize(("info", "message"), [(7, "status 7"), (0, "relative residual")])
def test_failed_solve_is_rejected(monkeypatch, info, message):
    monkeypatch.setattr(fea, "cg", lambda *_args, **_kwargs: (np.zeros(2), info))
    solver = fea.AuditedCG()
    with pytest.raises(ValueError, match=message):
        solver(eye(2, format="csr"), np.ones(2))
    assert solver.records == []


def test_reordered_cg_returns_original_degrees_of_freedom():
    matrix = csr_matrix([
        [4.0, 0.0, -1.0, 0.0],
        [0.0, 2.0, 0.0, -1.0],
        [-1.0, 0.0, 5.0, -1.0],
        [0.0, -1.0, -1.0, 3.0],
    ])
    rhs = np.array([2.0, -3.0, 4.0, 1.0])
    solver = fea.ReorderedCG()
    np.testing.assert_allclose(matrix @ solver(matrix, rhs), rhs, atol=1e-12)
    assert solver.records[-1]["original_relative_residual"] < 1e-12


def test_reordered_solver_rejects_bad_original_residual(monkeypatch):
    monkeypatch.setattr(fea.AuditedCG, "__call__", lambda _self, _matrix, rhs: np.zeros_like(rhs))
    with pytest.raises(ValueError, match="Original-system CG relative residual"):
        fea.ReorderedCG()(eye(3, format="csr"), np.ones(3))
