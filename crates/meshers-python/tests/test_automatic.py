import meshers
import numpy as np
import pytest
from test_rust_compile import grad, gyroid


def test_automatic_derivatives_and_grading():
    f = meshers.compile_field(gyroid)
    for point in np.random.default_rng(42).uniform(-1, 1, (30, 3)):
        np.testing.assert_allclose(f.gradient(*point), grad(*point), atol=2e-14)

    def graded(x, y, z):
        return ((1 - x) * gyroid(x, y, z) + x * (np.cos(x) + np.cos(y) + np.cos(z))) / (
            0.5 + 0.1 * z
        )

    f = meshers.compile_field(graded)
    for point in np.random.default_rng(42).uniform(-0.5, 0.5, (20, 3)):
        numeric = []
        for axis in range(3):
            delta = np.eye(3)[axis] * 1e-6
            numeric.append((graded(*(point + delta)) - graded(*(point - delta))) / 2e-6)
        np.testing.assert_allclose(f.gradient(*point), numeric, rtol=1e-7, atol=1e-8)


def test_cache_captured_parameters_and_automatic_generate():
    radius = [0.3]

    def field(x, y, z):
        return (x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 - radius[0] ** 2

    first = meshers.compile_field(field)
    assert meshers.compile_field(field) is first
    radius[0] = 0.4
    second = meshers.compile_field(field)
    assert first is not second
    assert abs(second(0.5, 0.5, 0.5) + 0.16) < 1e-14
    m = meshers.generate(field, cells=8, geometry_tolerance=0.1)
    assert m.diagnostics["evaluator"] == "compiled"
    assert m.diagnostics["callback_calls"] == 0


def test_piecewise_derivatives_and_integer_power():
    f = meshers.compile_field(
        lambda x, y, z: np.where(x > 0, x * x, -x) + np.abs(y) + np.maximum(z, 0)
    )
    np.testing.assert_allclose(f.gradient(-2, 0, -1), [-1, 0, 0])
    np.testing.assert_allclose(f.gradient(2, 1, 1), [4, 1, 1])
    f = meshers.compile_field(lambda x, y, z: x**3)
    np.testing.assert_allclose(f.gradient(-2, 0, 0), [12, 0, 0])


def test_unsupported_numpy_falls_back_visibly():
    def field(x, y, z):
        return np.full_like(x, -1.0)

    with pytest.warns(RuntimeWarning, match="NumPy callbacks"):
        m = meshers.generate(field, cells=4, optimize_passes=0)
    assert m.diagnostics["evaluator"] == "python"
    assert m.diagnostics["callback_calls"] > 0


def test_invalid_graph_rejected():
    from meshers import _meshers

    for nodes, outputs in [
        ([("add", [0, 0], 0.0)], [0]),
        ([("input", [3], 0.0)], [0]),
        ([("constant", [], 1.0)], [4]),
        ([("invalid", [], 0.0)], [0]),
    ]:
        with pytest.raises(ValueError):
            _meshers._compile_expression(nodes, outputs)


def test_array_gradient_automatically_falls_back():
    def field(x, y, z):
        return (x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 - 0.1

    def gradient(x, y, z):
        return 2 * np.column_stack((x - 0.5, y - 0.5, z - 0.5))

    with pytest.warns(RuntimeWarning, match="NumPy callbacks"):
        result = meshers.generate(field, gradient=gradient, cells=8, geometry_tolerance=0.1)
    assert result.diagnostics["evaluator"] == "python"


def test_explicit_auto_threads_and_conservative_default(monkeypatch):
    from meshers import _meshers

    assert 1 <= _meshers.available_threads() <= 256
    monkeypatch.setattr(_meshers, "available_threads", lambda: 2)
    mesh = meshers.generate(lambda x, y, z: -1.0, cells=4, optimize_passes=0, threads=None)
    assert mesh.diagnostics["threads"] == 2
    explicit = meshers.generate(lambda x, y, z: -1.0, cells=4, optimize_passes=0)
    assert explicit.diagnostics["threads"] == 1


def test_compatibility_compile_import():
    from meshers.compiler import compile_field

    assert compile_field is meshers.compile_field
