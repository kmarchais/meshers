import meshers
import numpy as np
import pytest


@pytest.mark.parametrize("op", [np.arccos, np.arcsin, np.arctan])
def test_inverse_trig_value_and_derivative(op):
    compiled = meshers.compile_field(lambda x, y, z: op(x + 0.2 * y))
    for x in [-0.7, 0.2, 0.7]:
        point = np.array([x, 0.1, 0.0])
        assert compiled(*point) == pytest.approx(op(x + 0.02))
        value = point[0] + 0.2 * point[1]
        derivative = (op(value + 1e-6) - op(value - 1e-6)) / 2e-6
        np.testing.assert_allclose(
            compiled.gradient(*point), [derivative, 0.2 * derivative, 0], rtol=1e-6
        )


def test_atan2_and_clip():
    angle = meshers.compile_field(lambda x, y, z: np.arctan2(y, x))
    for x, y in [(1.0, 2.0), (-1.0, 2.0), (-1.0, -2.0), (1.0, -2.0)]:
        assert angle(x, y, 0) == pytest.approx(np.arctan2(y, x))
        np.testing.assert_allclose(
            angle.gradient(x, y, 0), [-y / (x * x + y * y), x / (x * x + y * y), 0]
        )
    np.testing.assert_allclose(angle.gradient(0, 0, 0), [0, 0, 0])
    clipped = meshers.compile_field(lambda x, y, z: np.clip(x, -0.5, 0.5))
    for x, grad in [(-1, 0), (0, 1), (1, 0)]:
        assert clipped(x, 0, 0) == pytest.approx(np.clip(x, -0.5, 0.5))
        np.testing.assert_allclose(clipped.gradient(x, 0, 0), [grad, 0, 0])


def test_norm_origin_subgradient_and_singular_sqrt():
    norm = meshers.compile_field(lambda x, y, z: np.sqrt(x * x + y * y + z * z))
    np.testing.assert_allclose(norm.gradient(0, 0, 0), [0, 0, 0])
    np.testing.assert_allclose(norm.gradient(3, 4, 0), [0.6, 0.8, 0])
    singular = meshers.compile_field(lambda x, y, z: np.sqrt(x))
    with pytest.raises(ValueError, match="nonfinite compiled gradient"):
        singular.gradient(0, 0, 0)


def test_clipped_norm_inside_box():
    field = meshers.compile_field(
        lambda x, y, z: np.sqrt(
            np.maximum(abs(x) - 1, 0) ** 2
            + np.maximum(abs(y) - 1, 0) ** 2
            + np.maximum(abs(z) - 1, 0) ** 2
        )
    )
    np.testing.assert_allclose(field.gradient(0, 0, 0), [0, 0, 0])
    np.testing.assert_allclose(field.gradient(2, 0, 0), [1, 0, 0])
