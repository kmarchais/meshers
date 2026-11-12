"""Experimental diagnostics must preserve the public generation result."""

from contextlib import contextmanager
from unittest.mock import patch

import meshers
import numpy as np
import pytest
from meshers import _meshers

pytestmark = pytest.mark.skipif(
    not getattr(_meshers, "_LINEAR_OBSERVATIONS", False),
    reason="requires experimental-surfaces",
)


@contextmanager
def observe(sink):
    original = _meshers.generate

    def generate(*args, **kwargs):
        return original(*args, **kwargs, _observations=sink)

    # Test-local, single-thread dispatch. The public parser and compiler run normally.
    with patch.object(_meshers, "generate", generate):
        yield


@pytest.mark.parametrize("adaptive", [False, True])
@pytest.mark.parametrize("compiled", [False, True])
def test_scalar_observations_preserve_output(adaptive, compiled):
    def field(x, y, z):
        return -1 + 0.05 * np.cos(2 * np.pi * x)

    options = {
        "cells": 4,
        "periodic": (True,) * 3,
        "order": 2,
        "optimize_passes": 1,
        "max_tetrahedra": 10000,
        "adapt": adaptive,
        "adapt_iterations": 0,
        "compile": compiled,
    }
    ordinary = meshers.generate(field, **options)
    sink = []
    with observe(sink):
        observed = meshers.generate(field, **options)
    for key in ("points", "tetrahedra", "surface", "boundary_tags"):
        np.testing.assert_array_equal(getattr(ordinary, key), getattr(observed, key))
    for a, b in zip(ordinary.periodic_pairs, observed.periodic_pairs, strict=True):
        np.testing.assert_array_equal(a, b)
    for key in (
        "volume",
        "minimum_mmg_quality",
        "sampled_surface_error",
        "elements_below_quality_01",
        "callback_calls",
        "evaluator",
    ):
        assert ordinary.diagnostics[key] == observed.diagnostics[key]
    assert [row["stage"] for row in sink] == (
        ["clipped", "adapted", "optimized", "geometry-closed", "final-linear"]
        if adaptive
        else ["clipped", "optimized", "geometry-refined", "final-linear"]
    )
    final = sink[-1]
    assert final["admission"] == "unchecked-linear-snapshot"
    np.testing.assert_array_equal(final["points"], observed.points[: len(final["points"])])
    np.testing.assert_array_equal(final["tetrahedra"], observed.tetrahedra[:, :4])
    assert final["surface"].shape[1] == 3


def test_scalar_observations_preserve_exact_failure():
    options = {
        "bounds": (0, 1, 0, 1, 0, 0.001),
        "cells": 4,
        "optimize_passes": 0,
        "minimum_quality": 0.05,
        "max_tetrahedra": 10000,
    }

    def field(x, y, z):
        return x * 0 - 1

    with pytest.raises(meshers.MeshingError, match="minimum MMG quality") as ordinary:
        meshers.generate(field, **options)
    sink = []
    with observe(sink), pytest.raises(meshers.MeshingError) as observed:
        meshers.generate(field, **options)
    assert type(observed.value) is type(ordinary.value)
    assert str(observed.value) == str(ordinary.value)
    assert sink[-1]["stage"] == "final-linear"
    assert len(sink[-1]["tetrahedra"]) > 0


def test_scalar_observations_reject_non_list_sink():
    with observe({}), pytest.raises(TypeError, match="_observations"):
        meshers.generate("gyroid", cells=4)
