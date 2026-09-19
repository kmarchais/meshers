"""Run user examples unchanged and validate the meshes they export."""

import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.mark.parametrize(
    "name", ["basic", "evaluators", "graded", "periodic", "intersection", "cancellation"]
)
def test_user_example(name, tmp_path):
    result = subprocess.run(
        [sys.executable, str(EXAMPLES / f"{name}.py")],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": os.pathsep.join(sys.path)},
        text=True,
        capture_output=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for path in tmp_path.glob("*.vtkhdf"):
        with h5py.File(path) as file:
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
