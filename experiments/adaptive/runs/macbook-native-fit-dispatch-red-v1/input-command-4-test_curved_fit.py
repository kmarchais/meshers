"""Public bounded fitting contract and exact coupled movement of periodic nodes."""

import meshers
import numpy as np
import pytest


def test_curved_fit_preserves_counts_corners_and_periodic_associations():
    options = {
        "cells": 9,
        "bounds": (-0.5, 0.5) * 3,
        "band": (-0.55, 0.55),
        "order": 2,
        "periodic": (True,) * 3,
        "geometry_tolerance": 1.0,
        "optimize_passes": 1,
    }
    before = meshers.generate("gyroid", **options)
    after = meshers.generate("gyroid", curved_fit=True, **options)
    np.testing.assert_array_equal(after.tetrahedra, before.tetrahedra)
    np.testing.assert_array_equal(after.surface, before.surface)
    corners = np.unique(after.tetrahedra[:, :4])
    np.testing.assert_array_equal(after.points[corners], before.points[corners])
    assert after.points.shape == before.points.shape
    fit = after.diagnostics["quadratic"]["fit"]
    assert fit["moved_nodes"] > 0
    assert fit["final_distance_estimate"] <= fit["initial_distance_estimate"]
    assert after.diagnostics["quadratic"]["minimum_jacobian_ratio"] > 0.001
    for axis, (a, b) in enumerate(zip(before.periodic_pairs, after.periodic_pairs, strict=True)):
        np.testing.assert_array_equal(a, b)
        difference = after.points[b[:, 1]] - after.points[b[:, 0]] - np.eye(3)[axis]
        assert abs(difference).max() <= 4 * np.finfo(float).eps


@pytest.mark.parametrize("analytic", [True, False])
def test_curved_fit_requires_quadratic_elements_and_leaves_a_plane_unchanged(analytic):
    with pytest.raises(ValueError, match="curved_fit requires order 2"):
        meshers.generate("gyroid", curved_fit=True, cells=4)

    def field(x, _y, _z):
        return x - 0.15

    compiled = meshers.compile_field(field, gradient=None if analytic else False)

    options = {
        "cells": 4,
        "order": 2,
        "geometry_tolerance": 1.0,
        "optimize_passes": 0,
        "bounds": (-0.5, 0.5) * 3,
    }
    before = meshers.generate(compiled, **options)
    after = meshers.generate(compiled, curved_fit=True, **options)
    np.testing.assert_array_equal(before.points, after.points)
    assert after.diagnostics["quadratic"]["fit"]["moved_nodes"] == 0


def test_curved_fit_rejects_linear_retained_intersection_before_dispatch():
    def intersection(x, y, _z):
        return np.maximum(x - 0.15, y - 0.2)

    with pytest.raises(ValueError, match="curved_fit requires order 2"):
        meshers.generate(intersection, curved_fit=True, cells=4)
