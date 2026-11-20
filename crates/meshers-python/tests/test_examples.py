# Copyright (c) 2026 Kevin Marchais

"""Run user examples unchanged and validate the meshes they export."""

from pathlib import Path

import numpy as np
import pyfive
import pytest
from python_process import run_python

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.mark.parametrize(
    "name",
    ["basic", "evaluators", "graded", "periodic", "intersection", "cancellation"],
)
def test_user_example(name, tmp_path):
    run_python(str(EXAMPLES / f"{name}.py"), cwd=tmp_path)
    for path in tmp_path.glob("*.vtkhdf"):
        with pyfive.File(path) as file:
            group = file["VTKHDF"]
            points = group["Points"][:]
            cells = group["Connectivity"][:]
            assert len(points) > 0
            assert len(cells) > 0
            assert np.isfinite(points).all()
            assert cells.min() >= 0
            assert cells.max() < len(points)
            if "Volume" in group["CellData"]:
                assert np.all(group["CellData/Volume"][:] > 0)
