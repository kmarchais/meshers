"""Exact feature recovery from bounded convex piecewise-affine expression graphs."""

from itertools import product

import meshers
import numpy as np
import pytest


def octahedron(x, y, z):
    return (abs(x) + abs(y) + abs(z) - 0.35) / np.sqrt(3)


def test_default_scalar_route_preserves_polyhedral_creases():
    mesh = meshers.generate(
        octahedron,
        bounds=(-0.5, 0.5) * 3,
        cells=24,
        geometry_tolerance=0.002,
        minimum_quality=0.1,
        max_tetrahedra=250000,
        optimize_passes=1,
    )
    assert mesh.diagnostics["polyhedral_constraints"] == 8
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    assert abs(mesh.diagnostics["volume"] - 4 * 0.35**3 / 3) < 1e-12
    assert len(mesh.feature_edges) > 0
    probes = np.concatenate((
        mesh.points[np.unique(mesh.surface)],
        mesh.points[mesh.surface].mean(axis=1),
    ))
    assert np.max(abs(octahedron(*probes.T))) < 1e-12
    for face in mesh.surface:
        points = mesh.points[face]
        assert any(
            np.max(abs(points @ normal - 0.35)) < 1e-12 for normal in product((1, -1), repeat=3)
        )


def test_halfspaces_match_the_compiled_field_and_reject_nonconvex_graphs():
    field = meshers.compile_field(octahedron)
    planes = np.asarray(field._halfspaces)
    assert planes.shape == (8, 4)
    points = np.random.default_rng(42).uniform(-0.5, 0.5, (256, 3))
    np.testing.assert_allclose(
        (points @ planes[:, :3].T + planes[:, 3]).max(axis=1), field.values(points), atol=3e-16
    )
    for function in (
        lambda x, y, z: np.minimum(x, y),
        lambda x, y, z: -abs(x),
        lambda x, y, z: x * x + y * y + z * z - 0.1,
        lambda x, y, z: abs(x) * (y + 2),
        lambda x, y, z: sum(abs(x + i * y + i * i * z) for i in range(7)) - 0.1,
        lambda x, y, z: np.maximum(0, x),
    ):
        assert meshers.compile_field(function)._halfspaces is None


@pytest.mark.parametrize("scale", [0.5, 2.0])
def test_affine_scaling_translation_and_maximum_remain_exact(scale):
    def field(x, y, z):
        return scale * np.maximum(abs(x - 0.02) + abs(y + 0.03), 2 * abs(z)) - 0.3

    compiled = meshers.compile_field(field)
    planes = np.asarray(compiled._halfspaces)
    assert planes.shape == (6, 4)
    points = np.random.default_rng(12).uniform(-0.5, 0.5, (256, 3))
    np.testing.assert_allclose(
        (points @ planes[:, :3].T + planes[:, 3]).max(axis=1), compiled.values(points), atol=4e-16
    )


@pytest.mark.parametrize(
    "options",
    [
        {"threads": 2},
        {"order": 2},
        {"periodic": (True,) * 3},
        {
            "gradient": lambda x, y, z: tuple(
                np.where(v >= 0, 1.0, -1.0) / np.sqrt(3) for v in (x, y, z)
            )
        },
    ],
)
def test_other_modes_keep_the_scalar_pipeline(options):
    mesh = meshers.generate(
        octahedron,
        bounds=(-0.5, 0.5) * 3,
        cells=12,
        geometry_tolerance=1.0,
        optimize_passes=0,
        **options,
    )
    assert "polyhedral_constraints" not in mesh.diagnostics
    if "threads" in options:
        assert mesh.diagnostics["threads"] == options["threads"]
    if "order" in options:
        assert mesh.tetrahedra.shape[1] == 10


def test_polyhedral_route_keeps_cancellation_and_budget_gates():
    cancel = meshers.CancellationToken()
    cancel.cancel()
    with pytest.raises(meshers.CancelledError):
        meshers.generate(octahedron, bounds=(-0.5, 0.5) * 3, cells=12, cancel=cancel)
    with pytest.raises(meshers.MeshingError, match="budget"):
        meshers.generate(octahedron, bounds=(-0.5, 0.5) * 3, cells=12, max_tetrahedra=1)
    with pytest.raises(meshers.MeshingError, match="minimum MMG quality"):
        meshers.generate(octahedron, bounds=(-0.5, 0.5) * 3, cells=12, minimum_quality=1.0)


@pytest.mark.parametrize("adapt", [False, True])
def test_rotated_translated_polyhedron_preserves_its_exact_faces(adapt):
    cosine, sine = np.cos(0.37), np.sin(0.37)

    def field(x, y, z):
        u = cosine * (x - 0.02) + sine * (y + 0.03)
        v = -sine * (x - 0.02) + cosine * (y + 0.03)
        return (abs(u) + abs(v) + abs(z - 0.01) - 0.3) / np.sqrt(3)

    mesh = meshers.generate(
        field,
        bounds=(-0.5, 0.5) * 3,
        cells=24,
        geometry_tolerance=0.002,
        minimum_quality=0.1,
        max_tetrahedra=250000,
        optimize_passes=1,
        adapt=adapt,
        min_length=0.002,
        max_length=0.1,
        adapt_tolerance=0.002,
    )
    assert mesh.diagnostics["polyhedral_constraints"] == 8
    assert abs(mesh.diagnostics["volume"] - 4 * 0.3**3 / 3) < 1e-12
    probes = np.concatenate((
        mesh.points[np.unique(mesh.surface)],
        mesh.points[mesh.surface].mean(axis=1),
    ))
    assert np.max(abs(field(*probes.T))) < 1e-12
    if adapt:
        assert "adapt" in mesh.diagnostics
    else:
        planes = np.asarray(meshers.compile_field(field)._halfspaces)
        for index, plane in enumerate(planes):
            faces = mesh.surface_for(f"plane{index}")
            assert len(faces) > 0
            assert np.max(abs(mesh.points[faces] @ plane[:3] + plane[3])) < 1e-12


@pytest.mark.parametrize("options", [{"threads": 1.0}, {"interior_refine_nodes": 0.0}])
def test_polyhedral_dispatch_keeps_integer_option_validation(options):
    with pytest.raises(TypeError):
        meshers.generate(octahedron, bounds=(-0.5, 0.5) * 3, cells=12, **options)


def test_polyhedral_dispatch_keeps_zero_budget_validation():
    with pytest.raises(ValueError):
        meshers.generate(octahedron, bounds=(-0.5, 0.5) * 3, cells=12, max_tetrahedra=0)


@pytest.mark.parametrize("option", ["threads", "interior_refine_nodes", "optimize_passes"])
def test_polyhedral_dispatch_accepts_numpy_integer_options(option):
    options = {"threads": 1, "interior_refine_nodes": 0, "optimize_passes": 0}
    options[option] = np.int64(options[option])
    mesh = meshers.generate(octahedron, bounds=(-0.5, 0.5) * 3, cells=12, **options)
    assert mesh.diagnostics["polyhedral_constraints"] == 8


def test_adaptation_preserves_exact_polyhedral_faces_and_performs_work():
    mesh = meshers.generate(
        octahedron,
        bounds=(-0.5, 0.5) * 3,
        cells=24,
        adapt=True,
        geometry_tolerance=0.002,
        adapt_tolerance=0.002,
        min_length=0.002,
        max_length=0.1,
        gradation=0.3,
        minimum_quality=0.1,
        max_tetrahedra=250000,
        optimize_passes=1,
    )
    probes = np.concatenate((
        mesh.points[np.unique(mesh.surface)],
        mesh.points[mesh.surface].mean(axis=1),
    ))
    assert np.max(abs(octahedron(*probes.T))) < 1e-12
    assert abs(mesh.diagnostics["volume"] - 4 * 0.35**3 / 3) < 1e-12
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    assert mesh.diagnostics["polyhedral_constraints"] == 8
    assert mesh.diagnostics["adapt"]["collapses"] > 0


def test_polyhedral_adaptation_honors_the_spatial_size_callback():
    calls = []

    def size(x, y, z):
        calls.append(len(x))
        return np.full_like(x, 0.035)

    options = {
        "bounds": (-0.5, 0.5) * 3,
        "cells": 24,
        "adapt": True,
        "geometry_tolerance": 0.002,
        "adapt_tolerance": 0.002,
        "min_length": 0.002,
        "max_length": 0.1,
        "minimum_quality": 0.1,
        "max_tetrahedra": 250000,
        "optimize_passes": 1,
    }
    coarse = meshers.generate(octahedron, **options)
    fine = meshers.generate(octahedron, size_field=size, **options)
    assert calls
    assert max(calls) <= 4096
    assert len(fine.tetrahedra) > len(coarse.tetrahedra)
    assert fine.diagnostics["adapt"]["splits"] > 0
    assert fine.diagnostics["minimum_mmg_quality"] >= 0.1
    assert abs(fine.diagnostics["volume"] - 4 * 0.35**3 / 3) < 1e-12
    assert np.max(abs(octahedron(*fine.points[fine.surface].mean(axis=1).T))) < 1e-12


def test_polyhedral_adaptation_preserves_callback_exceptions():
    marker = LookupError("spatial size callback marker")

    def size(x, y, z):
        raise marker

    with pytest.raises(LookupError) as caught:
        meshers.generate(
            octahedron,
            bounds=(-0.5, 0.5) * 3,
            cells=24,
            adapt=True,
            min_length=0.002,
            max_length=0.1,
            size_field=size,
        )
    assert caught.value is marker


def test_polyhedral_refinement_keeps_the_hard_budget():
    with pytest.raises(meshers.MeshingError, match="max_tetrahedra"):
        meshers.generate(
            octahedron,
            bounds=(-0.5, 0.5) * 3,
            cells=24,
            adapt=True,
            geometry_tolerance=0.002,
            min_length=0.002,
            max_length=0.1,
            minimum_quality=0.1,
            max_tetrahedra=20000,
            size_field=lambda x, y, z: np.full_like(x, 0.02),
        )


def test_polyhedral_adaptation_keeps_invalid_sizing_rejection():
    with pytest.raises(meshers.MeshingError, match="invalid adaptation options"):
        meshers.generate(
            octahedron, bounds=(-0.5, 0.5) * 3, cells=24, adapt=True, min_length=0.2, max_length=0.1
        )
