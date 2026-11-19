import gc
import os
import re
import threading
import time
from functools import partial
from pathlib import Path

import meshers
import numpy as np
import pytest
from python_process import run_python

generate_callback = partial(meshers.generate, compile=False)


@pytest.mark.parametrize("compiled", [True, False])
def test_geometry_rejection_locates_a_sample_with_the_reported_error(compiled):
    def schwarz(x, y, z):
        return (np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y) + np.cos(2 * np.pi * z)) / 0.3

    with pytest.raises(meshers.MeshingError, match="sampled geometry error") as failure:
        meshers.generate(
            schwarz,
            cells=8,
            bounds=(-0.5, 0.5) * 3,
            band=(-1, 1),
            geometry_tolerance=0.01,
            optimize_passes=0,
            compile=compiled,
        )
    match = re.search(
        r"geometry error (\S+) exceeds.*worst sample at \[([^\]]+)\]", str(failure.value)
    )
    assert match
    point = np.fromstring(match[2], sep=",")
    assert point.shape == (3,)
    assert np.all(abs(point) <= 0.5)
    gradient = -2 * np.pi / 0.3 * np.sin(2 * np.pi * point)
    distance = abs(abs(schwarz(*point)) - 1) / np.linalg.norm(gradient)
    assert distance == pytest.approx(float(match[1]), rel=1e-8)
    assert distance > 0.01


@pytest.mark.parametrize("compiled", [True, False])
def test_adaptive_schwarz_d_corner_distance_is_finite(compiled):
    def field(x, y, z):
        x, y, z = [2 * np.pi * v for v in (x, y, z)]
        return (
            np.sin(x) * np.sin(y) * np.sin(z)
            + np.sin(x) * np.cos(y) * np.cos(z)
            + np.cos(x) * np.sin(y) * np.cos(z)
            + np.cos(x) * np.cos(y) * np.sin(z)
        )

    result = meshers.generate(
        field,
        bounds=(-0.5, 0.5) * 3,
        cells=12,
        periodic=(True,) * 3,
        adapt=True,
        adapt_iterations=1,
        adapt_tolerance=0.1,
        min_length=0.01,
        max_length=0.1,
        geometry_tolerance=0.02,
        minimum_quality=0.1,
        optimize_passes=1,
        compile=compiled,
    )
    assert np.isfinite(result.diagnostics["sampled_surface_error"])
    assert result.diagnostics["sampled_surface_error"] < 0.02
    assert result.diagnostics["volume"] == pytest.approx(0.5, abs=0.001)
    assert all(len(pairs) for pairs in result.periodic_pairs)


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

    def bad(_x, _y, _z):
        raise marker

    with pytest.raises(LookupError) as caught:
        generate_callback(bad)
    assert caught.value is marker


@pytest.mark.parametrize(
    "field",
    [
        lambda _x, _y, _z: 0.0,
        lambda x, _y, _z: np.zeros((len(x), 1)),
        lambda x, _y, _z: np.zeros(len(x), dtype=np.float32),
        lambda x, _y, _z: np.full(len(x), np.nan),
    ],
)
def test_invalid_callback_results(field):
    with pytest.raises((ValueError, TypeError)):
        generate_callback(field, cells=8)


def test_options_and_mesh_failure():
    with pytest.raises(ValueError, match=r"cells 4\.\.128"):
        generate_callback(sphere, cells=0)
    with pytest.raises(ValueError, match="periodic must contain three booleans"):
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


def test_keyboard_interrupt_in_native_work():
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
    assert "interrupted" in run_python("-c", code, timeout=15)


@pytest.mark.parametrize("delay", [0.15, 2.0])
@pytest.mark.parametrize("threads", [1, 4])
def test_native_interruption_on_windows_and_parallel_workers(delay, threads):
    # interrupt_main sets Python's pending SIGINT on Windows too. Exercise early
    # and later work; a wall-clock delay does not identify the native phase.
    # Explicit worker counts cover both the default and actual parallel paths.
    code = f"""
import sys, _thread, threading, time
sys.path.insert(0, {str(Path(meshers.__file__).parent.parent)!r})
import meshers
sent = []
def interrupt():
    sent.append(time.perf_counter())
    _thread.interrupt_main()
timer = threading.Timer({delay}, interrupt)
timer.start()
try:
    meshers.generate('gyroid', cells=96, band=(-.5,.5), max_tetrahedra=4_000_000, threads={threads})
except KeyboardInterrupt:
    elapsed = time.perf_counter() - sent[0]
    assert elapsed < 2, dict(interrupt_latency_s=elapsed, threads={threads}, delay={delay})
else:
    raise AssertionError('interrupt was ignored')
finally:
    timer.join()
"""
    run_python("-c", code, env={**os.environ, "RAYON_NUM_THREADS": "4"}, timeout=30)


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

    def grad(x, _y, _z):
        if len(x) > 6:
            raise problem
        return np.ones((len(x), 3), dtype=np.float64)

    with pytest.raises(RuntimeError) as caught:
        generate_callback(sphere, gradient=grad, cells=8)
    assert caught.value is problem
    token = meshers.CancellationToken()

    def cancel_gradient(x, _y, _z):
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


@pytest.mark.parametrize("order", [1, 2])
def test_expected_components_exposes_missing_island_in_nonempty_mesh(order):
    def solid(x, y, z):
        sphere = (x - 0.125) ** 2 + (y - 0.125) ** 2 + (z - 0.25) ** 2 - 0.06**2
        return np.minimum(0.75 - x, sphere)

    opts = {"order": order, "optimize_passes": 0, "snap": 0, "geometry_tolerance": 0.02}
    with pytest.raises(meshers.MeshingError, match="Expected 2 face-connected components, found 1"):
        meshers.generate(solid, cells=4, expected_components=2, **opts)
    fine = meshers.generate(solid, cells=16, expected_components=2, **opts)
    assert fine.diagnostics["face_connected_components"] == 2


def test_expected_components_exposes_lost_thin_slab():
    def slab(x, _y, _z):
        return (x - 0.06) ** 2 - 0.025**2

    opts = {"expected_components": 1, "optimize_passes": 0, "snap": 0, "geometry_tolerance": 0.02}
    with pytest.raises(meshers.MeshingError, match="Expected 1 face-connected components, found 0"):
        meshers.generate(slab, cells=4, **opts)
    fine = meshers.generate(slab, cells=24, **opts)
    assert fine.diagnostics["volume"] == pytest.approx(0.05)


def test_expected_components_on_retained_planes():
    with pytest.raises(meshers.MeshingError, match="found 1"):
        meshers.generate(
            lambda x, y, _z: np.maximum(x - 0.4, y - 0.6), cells=4, expected_components=2
        )


@pytest.mark.parametrize("invalid", [-1, True, 1.5])
def test_invalid_expected_components(invalid):
    with pytest.raises(ValueError, match="expected_components"):
        meshers.generate("gyroid", cells=4, expected_components=invalid)


def test_expected_components_count_before_periodic_identification():
    mesh = meshers.generate(
        lambda x, _y, _z: -np.cos(2 * np.pi * x),
        cells=8,
        periodic=(True, False, False),
        expected_components=2,
    )
    assert len(mesh.periodic_pairs[0]) > 0
    assert mesh.diagnostics["face_connected_components"] == 2


@pytest.mark.parametrize("height", [0.125, 0.25])
def test_translated_component_witness(height):
    center = np.array([0.125, 0.125, height])

    def sphere(x, y, z):
        return (x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2 - 0.06**2

    # Require the interior witness at resolved spacing. Check the oracle itself
    # with a deleted component, independent of the producer's own diagnostics.
    mesh = meshers.generate(sphere, cells=16, optimize_passes=0, snap=0, geometry_tolerance=0.02)
    vertices = mesh.points[mesh.tetrahedra]

    def covers_center(cells):
        coords = np.linalg.solve(
            (cells[:, 1:] - cells[:, :1]).transpose(0, 2, 1), (center - cells[:, 0])[..., None]
        )[..., 0]
        return np.any(np.all(coords >= -1e-10, axis=1) & (coords.sum(axis=1) <= 1 + 1e-10))

    assert covers_center(vertices)
    assert not covers_center(vertices[:0])
