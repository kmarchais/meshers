# Copyright (c) 2026 Kevin Marchais

"""Check distance bounds for the frozen spinodoid witness."""

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("scipy")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
from spinodoid_witness import distance_bounds


@pytest.mark.parametrize(
    ("free", "distance"),
    [([True] * 3, 0.08 / np.sqrt(2)), ([False, True, True], 0.08)],
)
def test_oblique_zero_plane_at_box_cap(free, distance):
    spectrum = SimpleNamespace(
        modes=np.array([[1, 1, 0]]),
        amplitudes=np.ones(1),
        phases=np.zeros(1),
        cell_size=1,
        threshold=0,
    )
    lower, upper = distance_bounds(spectrum, np.array([-0.5, -0.33, 0]), np.array(free))
    assert 0.05 < lower <= distance
    assert upper == pytest.approx(distance)
