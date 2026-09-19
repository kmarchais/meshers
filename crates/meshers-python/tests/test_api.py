import gc
import threading
import time
from functools import partial

import meshers
import numpy as np
import pytest

generate_callback = partial(meshers.generate, compile=False)


def sphere(x, y, z):
    return np.sqrt((x - 0.47) ** 2 + (y - 0.51) ** 2 + (z - 0.49) ** 2) - 0.29


def test_custom_arrays_and_batching():
    sizes = []

    def field(x, y, z):
        sizes.append(len(x))
        return sphere(x, y, z)

    result = generate_callback(field, cells=8, geometry_tolerance=0.1, optimize_passes=0)
    gc.collect()
    assert result.points.dtype == np.float64
    assert result.tetrahedra.dtype == np.int64
    assert result.tetrahedra.shape[1] == 4
    assert result.points[result.tetrahedra].shape[2] == 3
    assert max(sizes) > 6
    assert max(sizes) <= 4096
    assert result.diagnostics["volume"] > 0.08
    assert all(p.shape == (0, 2) for p in result.periodic_pairs)


def test_batch_sizes_do_not_change_geometry():
    args = {"cells": 8, "geometry_tolerance": 0.1, "optimize_passes": 0}
    a = generate_callback(sphere, batch_size=7, **args)
    b = generate_callback(sphere, batch_size=4096, **args)
    np.testing.assert_array_equal(a.tetrahedra, b.tetrahedra)
    np.testing.assert_array_equal(a.points, b.points)


def test_native_and_numpy_periodic_field_agree():
    def gyroid(x, y, z):
        x, y, z = [2 * np.pi * v for v in (x, y, z)]
        return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)

    args = {
        "cells": 12,
        "band": (-0.5, 0.5),
        "periodic": (True,) * 3,
        "geometry_tolerance": 0.1,
        "optimize_passes": 0,
    }
    a = generate_callback("gyroid", **args)
    b = generate_callback(gyroid, **args)
    assert all(len(p) for p in a.periodic_pairs)
    assert abs(a.diagnostics["volume"] - b.diagnostics["volume"]) < 1e-8
    assert a.diagnostics["callback_calls"] == 0


def test_gradient_and_original_exception():
    def gradient(x, y, z):
        d = np.column_stack([x - 0.47, y - 0.51, z - 0.49])
        return d / np.linalg.norm(d, axis=1)[:, None]

    assert len(
        generate_callback(
            sphere,
            gradient=gradient,
            cells=8,
            geometry_tolerance=0.1,
            optimize_passes=0,
        ).points
    )
    marker = LookupError("original callback failure")

    def bad(x, y, z):
        raise marker

    with pytest.raises(LookupError) as caught:
        generate_callback(bad)
    assert caught.value is marker


@pytest.mark.parametrize(
    "field",
    [
        lambda x, y, z: 0.0,
        lambda x, y, z: np.zeros((len(x), 1)),
        lambda x, y, z: np.zeros(len(x), dtype=np.float32),
        lambda x, y, z: np.full(len(x), np.nan),
    ],
)
def test_invalid_callback_results(field):
    with pytest.raises((ValueError, TypeError)):
        generate_callback(field, cells=8)


def test_options_and_mesh_failure():
    with pytest.raises(ValueError):
        generate_callback(sphere, cells=0)
    with pytest.raises(ValueError):
        generate_callback(sphere, periodic=(1, 0, 0))
    with pytest.raises(meshers.MeshingError):
        generate_callback(sphere, cells=8, geometry_tolerance=1e-12)


def test_cancellation_from_another_python_thread():
    token = meshers.CancellationToken()
    worker = threading.Thread(target=lambda: (time.sleep(0.05), token.cancel()))
    worker.start()
    try:
        with pytest.raises(meshers.CancelledError):
            meshers.generate("gyroid", cells=60, band=(-0.5, 0.5), cancel=token)
    finally:
        worker.join(timeout=5)
    assert not worker.is_alive()


def test_vtkhdf_round_trip_and_physical_shifts(tmp_path):
    import pyvista as pv

    result = generate_callback(
        lambda x, y, z: np.full_like(x, -1.0),
        cells=4,
        bounds=(2, 4, -1, 2, 3, 7),
        periodic=(True, True, True),
        optimize_passes=0,
    )
    path = tmp_path / "mesh.vtkhdf"
    result.write_vtkhdf(path)
    restored = pv.read(path)
    np.testing.assert_array_equal(restored.points, result.points)
    np.testing.assert_array_equal(restored.cells.reshape(-1, 5)[:, 1:], result.tetrahedra)
    np.testing.assert_allclose(restored.point_data["PeriodicShift"].max(axis=0), [2, 3, 4])
    assert restored.cell_data["MMGQuality"].min() > 0
    surface = tmp_path / "surface.vtkhdf"
    result.write_vtkhdf(surface, surface=True)
    assert pv.read(surface).n_cells == len(result.surface)


def test_empty_mesh_shapes():
    result = generate_callback(lambda x, y, z: np.ones_like(x), cells=4)
    assert result.points.shape == (0, 3)
    assert result.tetrahedra.shape == (0, 4)
    assert result.surface.shape == (0, 3)
    assert result.diagnostics["empty"]


def test_keyboard_interrupt_in_native_work():
    import os
    import subprocess
    import sys

    if os.name != "posix":
        pytest.skip("POSIX signal test")
    code = """
import os, signal, threading
import meshers
timer = threading.Timer(.1, lambda: os.kill(os.getpid(), signal.SIGINT))
timer.start()
try:
    meshers.generate('gyroid', cells=60, band=(-.5,.5))
except KeyboardInterrupt:
    print('interrupted')
else:
    raise AssertionError('interrupt was ignored')
finally:
    timer.join()
"""
    run = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, check=False, text=True, timeout=15
    )
    assert run.returncode == 0, run.stderr
    assert "interrupted" in run.stdout


def test_optimizer_batch_size_and_gradient_batches():
    sizes = []

    def grad(x, y, z):
        sizes.append(len(x))
        d = np.column_stack((x - 0.47, y - 0.51, z - 0.49))
        return d / np.linalg.norm(d, axis=1)[:, None]

    args = {
        "cells": 8,
        "optimize_passes": 1,
        "geometry_tolerance": 0.1,
        "threads": 1,
        "gradient": grad,
    }
    a = generate_callback(sphere, batch_size=7, **args)
    assert max(sizes) <= 7
    sizes.clear()
    b = generate_callback(sphere, batch_size=4096, **args)
    assert max(sizes) > 7
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(a.tetrahedra, b.tetrahedra)


def test_error_and_cancellation_during_batched_gradient():
    problem = RuntimeError("batch gradient failed")

    def grad(x, y, z):
        if len(x) > 6:
            raise problem
        return np.ones((len(x), 3), dtype=np.float64)

    with pytest.raises(RuntimeError) as caught:
        generate_callback(sphere, gradient=grad, cells=8)
    assert caught.value is problem
    token = meshers.CancellationToken()

    def cancel_gradient(x, y, z):
        token.cancel()
        return np.ones((len(x), 3), dtype=np.float64)

    with pytest.raises(meshers.CancelledError):
        generate_callback(sphere, gradient=cancel_gradient, cancel=token, cells=8)


def test_native_field_with_python_gradient_uses_calling_thread():
    caller = threading.get_ident()
    sizes = []

    def grad(x, y, z):
        assert threading.get_ident() == caller
        sizes.append(len(x))
        x, y, z = (2 * np.pi * v for v in (x, y, z))
        return (
            2
            * np.pi
            * np.column_stack((
                np.cos(x) * np.cos(y) - np.sin(z) * np.sin(x),
                np.cos(y) * np.cos(z) - np.sin(x) * np.sin(y),
                np.cos(z) * np.cos(x) - np.sin(y) * np.sin(z),
            ))
        )

    m = generate_callback(
        "gyroid",
        gradient=grad,
        cells=8,
        band=(-0.5, 0.5),
        geometry_tolerance=0.1,
        optimize_passes=1,
        threads=4,
    )
    assert len(m.tetrahedra) > 0
    assert max(sizes) > 6
