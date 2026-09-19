import gc

import meshers
import numpy as np
import pytest


def gyroid(x, y, z):
    x = x * 2 * np.pi
    y = y * 2 * np.pi
    z = z * 2 * np.pi
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def grad(x, y, z):
    x = x * 2 * np.pi
    y = y * 2 * np.pi
    z = z * 2 * np.pi
    return (
        2 * np.pi * (np.cos(x) * np.cos(y) - np.sin(z) * np.sin(x)),
        2 * np.pi * (np.cos(y) * np.cos(z) - np.sin(x) * np.sin(y)),
        2 * np.pi * (np.cos(z) * np.cos(x) - np.sin(y) * np.sin(z)),
    )


def test_rust_scalar_and_gradient_match():
    for derivative in (None, grad):
        f = meshers.compile_field(gyroid, gradient=derivative, backend="rust")
        gc.collect()
        assert abs(f(0.1, 0.2, 0.3) - gyroid(0.1, 0.2, 0.3)) < 1e-14
        a = meshers.generate(
            f,
            cells=12,
            band=(-0.5, 0.5),
            periodic=(True,) * 3,
            geometry_tolerance=0.1,
            optimize_passes=1,
            threads=4,
        )
        assert a.diagnostics["callback_calls"] == 0
        assert len(a.tetrahedra) == 9422
        assert abs(a.diagnostics["volume"] - 0.326240000731) < 1e-10


def test_rust_grading_and_unsupported_branch():
    def graded(x, y, z):
        return gyroid(x, y, z) / (0.35 + 0.15 * x)

    f = meshers.compile_field(graded, backend="rust")
    assert abs(f(0.2, 0.3, 0.4) - graded(0.2, 0.3, 0.4)) < 1e-14

    def unsupported(x, y, z):
        return x if x == 0 else y

    with pytest.raises(TypeError):
        meshers.compile_field(unsupported, backend="rust")


def test_rust_compiler_does_not_require_numba():
    import subprocess
    import sys

    code = """
import sys,importlib.abc
class BlockNumba(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='numba' or fullname.startswith('numba.'):
            raise ImportError('Numba deliberately unavailable')
sys.meta_path.insert(0,BlockNumba())
import meshers
f=meshers.compile_field(lambda x,y,z:x*x+y*y+z*z-.5,backend='rust')
assert abs(f(.1,.2,.3)+.36)<1e-14
assert 'numba' not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)
