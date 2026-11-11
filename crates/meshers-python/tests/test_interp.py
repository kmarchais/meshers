"""Native constant interpolation tables retain NumPy values and frozen captures."""

import meshers
import numpy as np
import pytest


@pytest.mark.parametrize(
    ("knots", "values", "left", "right"),
    [
        ([-0.5, 0.0, 0.2, 0.5], [0.3, 0.4, 0.35, 0.45], None, None),
        ([-2.0, -0.5, 0.0, 3.0], [-2.0, 1.0, 0.2, 2.0], -0.6, 3.0),
        ([0.2], [1.7], -0.2, 2.0),
    ],
)
def test_native_interp_values_at_knots_adjacent_and_exterior(knots, values, left, right):
    def field(x, y, z):
        return np.interp(z, knots, values, left=left, right=right)

    compiled = meshers.compile_field(field)
    knots = np.asarray(knots)
    positions = np.r_[
        np.linspace(-4.0, 4.0, 1001),
        knots,
        np.nextafter(knots, -np.inf),
        np.nextafter(knots, np.inf),
    ]
    points = np.zeros((len(positions), 3))
    points[:, 2] = positions
    np.testing.assert_allclose(compiled.values(points), field(*points.T), atol=1e-14, rtol=1e-14)


def test_interp_chain_rule_and_knot_policy():
    field = meshers.compile_field(
        lambda x, y, z: np.interp(2 * x - y, [0.0, 0.4, 1.0], [0.0, 1.0, 0.0])
    )
    np.testing.assert_allclose(field.gradient(0.1, 0.05, 0), [5.0, -2.5, 0.0])
    np.testing.assert_allclose(field.gradient(0.5, 0.2, 0), [-10 / 3, 5 / 3, 0.0])
    np.testing.assert_allclose(field.gradient(0.2, 0, 0), [-10 / 3, 5 / 3, 0.0])
    for x in (-1.0, 0.0, 0.5, 1.0):
        np.testing.assert_array_equal(field.gradient(x, 0, 0), [0.0, 0.0, 0.0])


def test_interp_frozen_capture_and_explicit_finite_difference_mode():
    knots = np.array([0.0, 0.5, 1.0])
    values = np.array([0.0, 0.4, 0.8])

    def field(x, y, z):
        return np.interp(z, knots, values)

    first = meshers.compile_field(field, gradient=False)
    assert not first.has_gradient
    np.testing.assert_equal(first(0, 0, 0.5), 0.4)
    values[1] = 0.6
    second = meshers.compile_field(field, gradient=False)
    np.testing.assert_equal(second(0, 0, 0.5), 0.6)
    np.testing.assert_equal(first(0, 0, 0.5), 0.4)


@pytest.mark.parametrize(
    ("knots", "values", "options"),
    [
        ([0.0, 1.0], [0.0, 1.0], {"period": 1.0}),
        ([0.0, 0.0, 1.0], [0.0, 0.5, 1.0], {}),
        ([1.0, 0.0], [0.0, 1.0], {}),
        ([0.0, np.nan], [0.0, 1.0], {}),
        ([0.0, 1.0], [0.0, 1j], {}),
        ([], [], {}),
        ([[0.0, 1.0]], [[0.0, 1.0]], {}),
    ],
)
def test_unsupported_interpolation_stays_on_callback_contract(knots, values, options):
    with pytest.raises(meshers.TraceError, match="interp"):
        meshers.compile_field(lambda x, y, z: np.interp(z, knots, values, **options))


def test_interp_automatic_generation_uses_native_evaluator():
    def field(x, y, z):
        return x * x + y * y + z * z - np.interp(z, [-0.5, 0.0, 0.5], [0.08, 0.1, 0.13])

    mesh = meshers.generate(
        field, bounds=(-0.5, 0.5) * 3, cells=8, order=2, optimize_passes=1, geometry_tolerance=0.05
    )
    assert len(mesh.tetrahedra) > 0
    assert mesh.diagnostics["evaluator"] == "compiled"
    assert mesh.diagnostics["callback_calls"] == 0
    assert mesh.diagnostics["quadratic"]["minimum_jacobian_ratio"] > 1e-3
