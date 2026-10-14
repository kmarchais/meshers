import meshers
import numpy as np
import pytest


def test_optional_interior_refinement_preserves_periodic_boundary():
    options = {"cells": 4, "optimize_passes": 0, "minimum_quality": 0.1, "periodic": (True,) * 3}

    def field(x, y, z):
        return x * 0 - 1

    seed = meshers.generate(field, **options)
    refined = meshers.generate(field, interior_refine_nodes=20, **options)
    assert 0 < refined.diagnostics["inserted_interior_nodes"] <= 20
    np.testing.assert_array_equal(refined.points[: len(seed.points)], seed.points)
    np.testing.assert_array_equal(refined.surface, seed.surface)
    np.testing.assert_array_equal(refined.boundary_tags, seed.boundary_tags)
    for old, new in zip(seed.periodic_pairs, refined.periodic_pairs, strict=True):
        np.testing.assert_array_equal(old, new)
    assert abs(refined.diagnostics["volume"] - seed.diagnostics["volume"]) < 1e-12
    assert refined.diagnostics["minimum_mmg_quality"] >= 0.1
    bounded = meshers.generate(
        field, interior_refine_nodes=20, max_tetrahedra=len(seed.tetrahedra), **options
    )
    assert bounded.diagnostics["inserted_interior_nodes"] == 0
    np.testing.assert_array_equal(bounded.points, seed.points)
    np.testing.assert_array_equal(bounded.tetrahedra, seed.tetrahedra)


def test_minimum_quality_gate_and_diagnostics():
    args = {
        "bounds": (0, 1, 0, 1, 0, 0.001),
        "cells": 4,
        "optimize_passes": 0,
    }
    mesh = meshers.generate(lambda x, y, z: x * 0 - 1, **args)
    assert mesh.diagnostics["minimum_mmg_quality"] < 0.01
    assert mesh.diagnostics["elements_below_quality_01"] > 0
    assert mesh.diagnostics["reverted_snap_vertices"] == 0
    with pytest.raises(meshers.MeshingError, match="minimum MMG quality"):
        meshers.generate(lambda x, y, z: x * 0 - 1, minimum_quality=0.05, **args)
    accepted = meshers.generate(
        lambda x, y, z: x * 0 - 1, cells=4, optimize_passes=0, minimum_quality=0.05
    )
    assert accepted.diagnostics["minimum_mmg_quality"] >= 0.05
    assert meshers.generate(lambda x, y, z: x * 0 + 1, cells=4, minimum_quality=1.0).diagnostics[
        "empty"
    ]


@pytest.mark.parametrize("value", [-0.1, 1.1, np.nan, np.inf])
def test_invalid_minimum_quality(value):
    with pytest.raises(ValueError, match="minimum_quality"):
        meshers.generate("gyroid", cells=4, minimum_quality=value)


@pytest.mark.parametrize("compiled", [True, False])
def test_safe_snapping_preserves_periodic_pairs(compiled):
    def sheet(x, y, z):
        a, b, c = [2 * np.pi * (2 * v + 0.125) for v in (x, y, z)]
        return abs(np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)) - 0.45

    mesh = meshers.generate(
        sheet,
        cells=16,
        periodic=(True,) * 3,
        snap=0.2,
        optimize_passes=0,
        geometry_tolerance=0.04,
        compile=compiled,
    )
    assert mesh.diagnostics["reverted_snap_vertices"] > 0
    assert mesh.diagnostics["sampled_surface_error"] <= 0.04
    assert mesh.diagnostics["minimum_mmg_quality"] > 0
    for axis, pairs in enumerate(mesh.periodic_pairs):
        assert len(pairs) > 0
        shift = np.zeros(3)
        shift[axis] = 1
        np.testing.assert_allclose(
            mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]],
            np.broadcast_to(shift, (len(pairs), 3)),
            atol=1e-10,
        )
