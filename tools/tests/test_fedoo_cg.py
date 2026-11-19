"""Check solver admission with a deliberately failed CG result."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
pytest.importorskip("fedoo")
import fedoo_homogenization as fea
from scipy.sparse import eye


@pytest.mark.parametrize(("info", "message"), [(7, "status 7"), (0, "relative residual")])
def test_failed_solve_is_rejected(monkeypatch, info, message):
    monkeypatch.setattr(fea, "cg", lambda *_args, **_kwargs: (np.zeros(2), info))
    solver = fea.AuditedCG()
    with pytest.raises(ValueError, match=message):
        solver(eye(2, format="csr"), np.ones(2))
    assert solver.records == []
