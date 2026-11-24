# Copyright (c) 2026 Kevin Marchais

"""Check that `strict`, `fallback` and `defeaturing` refuse degraded meshes."""

import warnings

import meshers
import numpy as np
import pytest

BOX = (-0.5, 0.5) * 3


def ball(x, y, z):
    return np.sqrt(x * x + y * y + z * z) - 0.35


def fibre(_x, y, z):
    # A periodic fibre along x whose wall passes 0.004 from four box edges.
    centre = 0.5 - 0.2 / np.sqrt(2) + 0.004
    distance = np.full_like(y, np.inf)
    for dy in (-1.0, 0.0, 1.0):
        for dz in (-1.0, 0.0, 1.0):
            distance = np.minimum(distance, np.hypot(y - centre - dy, z - centre - dz))
    return distance - 0.2


# Adaptation aims at a looser tolerance than the final checks accept.
REFUSED_ADAPTATION = {
    "bounds": BOX,
    "cells": 16,
    "periodic": (True, True, True),
    "geometry_tolerance": 0.004,
    "minimum_quality": 0.1,
    "adapt": True,
    "adapt_tolerance": 0.05,
    "min_length": 0.1,
    "max_length": 0.25,
}
GRAZING_FIBRE = {
    "compile_expression": False,
    "bounds": BOX,
    "cells": 16,
    "periodic": (True, True, True),
    "geometry_tolerance": 0.02,
    "minimum_quality": 0.1,
}


def test_fallback_is_returned_with_a_warning_and_refused_when_strict():
    with pytest.warns(UserWarning, match="uniform lattice mesh"):
        mesh = meshers.generate(ball, **REFUSED_ADAPTATION)
    assert mesh.diagnostics["adaptation_abandoned"]

    for refusal in ({"strict": True}, {"fallback": False}):
        with pytest.raises(meshers.MeshingError, match="sampled geometry error"):
            meshers.generate(ball, **REFUSED_ADAPTATION, **refusal)

    # Refusing only defeaturing, or re-allowing the fallback under strict, keeps it.
    for kept in ({"defeaturing": False}, {"strict": True, "fallback": True}):
        with pytest.warns(UserWarning, match="uniform lattice mesh"):
            mesh = meshers.generate(ball, **REFUSED_ADAPTATION, **kept)
        assert mesh.diagnostics["adaptation_abandoned"]


def test_defeaturing_is_returned_with_a_warning_and_refused_when_strict():
    with pytest.warns(UserWarning, match="off the implicit surface"):
        mesh = meshers.generate(fibre, **GRAZING_FIBRE)
    assert mesh.diagnostics["defeaturing_moves"] > 0

    for refusal in ({"strict": True}, {"defeaturing": False}):
        with pytest.raises(meshers.MeshingError, match="minimum MMG quality"):
            meshers.generate(fibre, **GRAZING_FIBRE, **refusal)

    for kept in ({"fallback": False}, {"strict": True, "defeaturing": True}):
        with pytest.warns(UserWarning, match="off the implicit surface"):
            meshers.generate(fibre, **GRAZING_FIBRE, **kept)


def test_strict_changes_nothing_for_a_mesh_that_is_not_degraded():
    options = {"bounds": BOX, "cells": 12, "geometry_tolerance": 0.02}
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        plain = meshers.generate(ball, **options)
        strict = meshers.generate(ball, **options, strict=True)
    np.testing.assert_array_equal(plain.points, strict.points)
    np.testing.assert_array_equal(plain.tetrahedra, strict.tetrahedra)


@pytest.mark.parametrize("name", ["strict", "fallback", "defeaturing"])
def test_degradation_options_must_be_bools(name):
    with pytest.raises(TypeError, match=f"{name} must be a bool"):
        meshers.generate(ball, bounds=BOX, cells=8, **{name: "fallback"})
