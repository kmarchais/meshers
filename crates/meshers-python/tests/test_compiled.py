import gc
import subprocess
import sys
import threading

import meshers
import numpy as np
import pytest


def gyroid(x, y, z):
    x = x * 2 * np.pi
    y = y * 2 * np.pi
    z = z * 2 * np.pi
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def sphere(x, y, z):
    return np.sqrt((x - 0.47) ** 2 + (y - 0.51) ** 2 + (z - 0.49) ** 2) - 0.29


def sphere_gradient(x, y, z):
    d = np.sqrt((x - 0.47) ** 2 + (y - 0.51) ** 2 + (z - 0.49) ** 2)
    return ((x - 0.47) / d, (y - 0.51) / d, (z - 0.49) / d)


def test_compiled_lifetime_and_native_match():
    f = meshers.compile_field(gyroid)
    gc.collect()
    assert abs(f(0.1, 0.2, 0.3) - gyroid(0.1, 0.2, 0.3)) < 1e-14
    args = {
        "cells": 12,
        "band": (-0.5, 0.5),
        "periodic": (True,) * 3,
        "optimize_passes": 1,
        "geometry_tolerance": 0.1,
        "threads": 4,
    }
    a = meshers.generate(f, **args)
    b = meshers.generate("gyroid", **args)
    assert a.diagnostics["callback_calls"] == 0
    np.testing.assert_allclose(a.points, b.points, rtol=0, atol=5e-9)
    np.testing.assert_array_equal(a.tetrahedra, b.tetrahedra)
    assert abs(a.diagnostics["volume"] - b.diagnostics["volume"]) < 1e-10
    assert abs(a.diagnostics["minimum_mmg_quality"] - b.diagnostics["minimum_mmg_quality"]) < 1e-8
    assert all(len(p) for p in a.periodic_pairs)


def test_compiled_gradient_and_captured_grading():
    width = 0.29

    def graded(x, y, z):
        return np.sqrt((x - 0.47) ** 2 + (y - 0.51) ** 2 + (z - 0.49) ** 2) - (width + 0.01 * x)

    f = meshers.compile_field(graded)
    a = meshers.generate(f, cells=10, geometry_tolerance=0.1, optimize_passes=1)
    b = meshers.generate(graded, cells=10, geometry_tolerance=0.1, optimize_passes=1)
    assert abs(a.diagnostics["volume"] - b.diagnostics["volume"]) < 1e-10
    f = meshers.compile_field(sphere, gradient=sphere_gradient)
    m = meshers.generate(f, cells=10, geometry_tolerance=0.1, optimize_passes=1)
    assert f.has_gradient
    assert m.diagnostics["callback_calls"] == 0


def test_compiled_cancellation():
    token = meshers.CancellationToken()
    f = meshers.compile_field(gyroid)
    timer = threading.Timer(0.05, token.cancel)
    timer.start()
    try:
        with pytest.raises(meshers.CancelledError):
            meshers.generate(f, cells=64, band=(-0.5, 0.5), cancel=token)
    finally:
        timer.join()


def test_no_external_compiler_or_optional_dependencies():
    import os

    code = """
import importlib.abc, subprocess, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('numba','sympy'):
            raise ImportError('deliberately unavailable')
sys.meta_path.insert(0,Block())
def forbidden(*a, **kw): raise AssertionError('external process requested')
subprocess.Popen=forbidden
import meshers
m=meshers.generate(lambda x,y,z:x*x+y*y+z*z-.5,cells=8,geometry_tolerance=.1)
assert m.diagnostics['evaluator']=='compiled'
assert m.diagnostics['callback_calls']==0
"""
    subprocess.run([sys.executable, "-c", code], env={**os.environ, "PATH": ""}, check=True)


@pytest.mark.parametrize("name", ["gyroid", "schwarz_p", "schwarz_d"])
def test_builtin_and_analytic_compiled_paths(name):
    def field(x, y, z):
        a, b, c = 2 * np.pi * x, 2 * np.pi * y, 2 * np.pi * z
        if name == "gyroid":
            return np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)
        if name == "schwarz_p":
            return np.cos(a) + np.cos(b) + np.cos(c)
        return (
            np.sin(a) * np.sin(b) * np.sin(c)
            + np.sin(a) * np.cos(b) * np.cos(c)
            + np.cos(a) * np.sin(b) * np.cos(c)
            + np.cos(a) * np.cos(b) * np.sin(c)
        )

    options = {
        "cells": 8,
        "band": (-0.5, 0.5),
        "periodic": (True, True, True),
        "geometry_tolerance": 0.1,
        "optimize_passes": 1,
    }
    native = meshers.generate(name, **options)
    compiled = meshers.generate(field, **options)
    np.testing.assert_array_equal(native.tetrahedra, compiled.tetrahedra)
    np.testing.assert_allclose(native.points, compiled.points, atol=5e-9, rtol=0)


@pytest.mark.parametrize(
    ("field", "expected", "gradient"),
    [
        (lambda x, y, z: z * z, 9.0, [0.0, 0.0, 6.0]),
        (lambda x, y, z: np.sin(y), float(np.sin(2.0)), [0.0, float(np.cos(2.0)), 0.0]),
        (lambda x, y, z: x * 0 + 4.0, 4.0, [0.0, 0.0, 0.0]),
    ],
)
def test_pruned_kernels_keep_input_coordinate_indices(field, expected, gradient):
    compiled = meshers.compile_field(field)
    assert compiled(1.0, 2.0, 3.0) == expected
    np.testing.assert_allclose(compiled.gradient(1.0, 2.0, 3.0), gradient, rtol=0, atol=1e-15)


def test_value_kernel_does_not_require_finite_derivative():
    compiled = meshers.compile_field(lambda x, y, z: np.sqrt(x))
    np.testing.assert_equal(compiled(0.0, 1.0, 2.0), 0.0)
    with pytest.raises(ValueError, match="gradient"):
        compiled.gradient(0.0, 1.0, 2.0)


def test_invalid_unused_nodes_are_still_rejected():
    with pytest.raises(ValueError, match="operation"):
        meshers._meshers._compile_expression([("constant", [], 1.0), ("unknown", [], 0.0)], [0])


def test_compiled_batch_matches_scalar_for_strided_points():
    f = meshers.compile_field(gyroid)
    storage = np.random.default_rng(43).normal(size=(9000, 6))
    points = storage[::2, ::2]
    expected = np.array([f(*p) for p in points])
    np.testing.assert_array_equal(f.values(points), expected)
    assert f.values(np.empty((0, 3))).shape == (0,)
    with pytest.raises(ValueError, match="shape"):
        f.values(np.zeros((3, 2)))


def test_compiled_batch_rejects_nonfinite_values():
    f = meshers.compile_field(lambda x, y, z: np.sqrt(x))
    with pytest.raises(ValueError, match="nonfinite"):
        f.values(np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]))
