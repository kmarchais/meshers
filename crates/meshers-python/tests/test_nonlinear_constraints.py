"""Retained nonlinear intersections through the scalar generation interface."""

import meshers
import numpy as np
import pytest


def bored_prism(x, y, z):
    box = np.maximum(np.maximum(abs(x) - 0.35, abs(y) - 0.35), abs(z) - 0.3)
    return np.maximum(box, 0.12 - np.sqrt((x - 0.03) ** 2 + (y + 0.02) ** 2))


def test_retained_branches_preserve_values_gradients_and_cached_owners():
    field = meshers.compile_field(bored_prism)
    branches = field._intersection_fields()
    assert len(branches) == 4
    assert all(a is b for a, b in zip(branches, field._intersection_fields(), strict=True))
    points = np.random.default_rng(42).uniform(-0.5, 0.5, (1024, 3))
    np.testing.assert_allclose(
        np.max([f.values(points) for f in branches], axis=0), field.values(points), atol=2e-16
    )
    np.testing.assert_allclose(branches[-1].gradient(0.15, -0.02, 0), (-1, 0, 0), atol=1e-14)
    assert all(f.has_gradient for f in branches)


def test_scalar_nonlinear_intersection_keeps_curves_geometry_and_quality():
    mesh = meshers.generate(
        bored_prism,
        bounds=(-0.5, 0.5) * 3,
        cells=34,
        geometry_tolerance=0.002,
        minimum_quality=0.1,
        max_tetrahedra=250000,
        optimize_passes=2,
    )
    assert mesh.diagnostics["retained_constraints"] == 7
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    exact = 0.7 * 0.7 * 0.6 - np.pi * 0.12**2 * 0.6
    assert abs(mesh.diagnostics["volume"] / exact - 1) <= 0.01
    assert len(mesh.tetrahedra) <= 250000
    assert len(mesh.points) <= 250000
    edges = mesh.feature_edges
    masks = mesh.constraint_masks[edges]
    on_bore = (masks & (1 << 6) != 0).all(axis=1)
    on_cap = ((masks & ((1 << 4) | (1 << 5))) != 0).all(axis=1)
    rim = mesh.points[np.unique(edges[on_bore & on_cap])]
    assert len(rim) > 0
    assert np.max(abs(abs(rim[:, 2]) - 0.3)) < 1e-12
    assert np.max(abs(np.linalg.norm(rim[:, :2] - (0.03, -0.02), axis=1) - 0.12)) <= 0.002


def test_custom_gradient_policies_do_not_get_replaced_by_branch_gradients():
    for gradient in (False, lambda x, y, z: (0, 0, 0)):
        assert meshers.compile_field(bored_prism, gradient=gradient)._intersection_fields() is None


def test_native_retained_roots_cannot_describe_a_different_field():
    from meshers import _meshers

    nodes = [("input", [0], 0.0), ("input", [1], 0.0), ("maximum", [0, 1], 0.0)]
    with pytest.raises(ValueError, match="outer maximum"):
        _meshers._compile_expression(nodes, [2], [[1], [0]])


def test_retained_evaluators_outlive_their_parent():
    import gc

    from meshers import _meshers

    nodes = [
        ("input", [0], 0.0),
        ("input", [1], 0.0),
        ("mul", [0, 0], 0.0),
        ("maximum", [2, 1], 0.0),
    ]
    parent = _meshers._compile_expression(nodes, [3], [[2], [1]])
    branches = parent._intersection_fields()
    del parent
    gc.collect()
    points = np.array(((0.2, 0.3, 0.4), (-0.7, -0.1, 0.2)))
    np.testing.assert_array_equal(branches[0].values(points), points[:, 0] ** 2)
    np.testing.assert_array_equal(branches[1].values(points), points[:, 1])
