"""Nonconvex sharp-feature regressions through the public scalar interface."""

import meshers
import numpy as np
import pytest


def notched_prism(x, y, z):
    return np.maximum(
        np.maximum(np.maximum(abs(x) - 0.35, abs(y) - 0.35), abs(z) - 0.3),
        np.minimum(x - 0.05, y - 0.05),
    )


def test_region_metadata_reproduces_the_actual_field():
    compiled = meshers.compile_field(notched_prism)
    planes, regions = compiled._polyhedral_regions
    assert compiled._halfspaces is None
    assert len(planes) == 8
    assert len(regions) == 2
    planes = np.asarray(planes)
    points = np.random.default_rng(42).uniform(-0.5, 0.5, (1024, 3))
    values = points @ planes[:, :3].T + planes[:, 3]
    branches = [
        values[:, [bool(mask & (1 << i)) for i in range(8)]].max(axis=1) for mask in regions
    ]
    np.testing.assert_allclose(np.min(branches, axis=0), compiled.values(points), atol=2e-16)


@pytest.mark.parametrize("adapt", [False, True])
@pytest.mark.parametrize("snap", [0.0, 0.2])
@pytest.mark.parametrize("angle", [0.0, 0.37])
def test_nonconvex_sharp_prism_preserves_geometry_and_quality(adapt, snap, angle):
    cosine, sine = np.cos(angle), np.sin(angle)
    cx, cy, cz = (0.02, -0.03, 0.01) if angle else (0.0, 0.0, 0.0)

    def field(x, y, z):
        return notched_prism(
            cosine * (x - cx) + sine * (y - cy),
            -sine * (x - cx) + cosine * (y - cy),
            z - cz,
        )

    # The rotated, unsnapped nonadaptive recipe needs two explicit quality
    # rounds. Its one-round rejection remains frozen in the experiment record.
    passes = 2 if angle and not adapt and snap == 0 else 1
    mesh = meshers.generate(
        field,
        bounds=(-0.5, 0.5) * 3,
        cells=24,
        adapt=adapt,
        min_length=0.002,
        max_length=0.1,
        adapt_tolerance=0.002,
        geometry_tolerance=0.002,
        minimum_quality=0.1,
        max_tetrahedra=250000,
        optimize_passes=passes,
        snap=snap,
    )
    assert mesh.diagnostics["polyhedral_regions"] == 2
    assert mesh.diagnostics["polyhedral_constraints"] == 8
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    assert abs(mesh.diagnostics["volume"] - 0.24) < 1e-12
    probes = np.concatenate((
        mesh.points[np.unique(mesh.surface)],
        mesh.points[mesh.surface].mean(axis=1),
    ))
    assert np.max(abs(field(*probes.T))) < 1e-12
    assert np.max(field(*mesh.points[mesh.tetrahedra].mean(axis=1).T)) <= 1e-12
    if adapt:
        assert mesh.diagnostics["adapt"]["collapses"] > 0


def test_unsupported_nonconvex_arithmetic_keeps_scalar_behavior():
    for field in (
        lambda x, y, z: np.minimum(x, y) + 0.1,
        lambda x, y, z: abs(np.minimum(x, y)),
        lambda x, y, z: np.minimum(x * x + y * y + z * z - 0.1, x),
    ):
        assert meshers.compile_field(field)._polyhedral_regions is None


def test_nonconvex_route_keeps_budget_quality_and_cancellation_gates():
    options = {"bounds": (-0.5, 0.5) * 3, "cells": 24, "optimize_passes": 1}
    with pytest.raises(meshers.MeshingError, match="budget"):
        meshers.generate(notched_prism, max_tetrahedra=100, **options)
    with pytest.raises(meshers.MeshingError, match="quality"):
        meshers.generate(notched_prism, minimum_quality=0.99, **options)
    cancel = meshers.CancellationToken()
    cancel.cancel()
    with pytest.raises(meshers.CancelledError):
        meshers.generate(notched_prism, cancel=cancel, **options)


def test_callback_metadata_cannot_override_the_requested_geometry():
    def sphere(x, y, z):
        return x * x + y * y + z * z - 0.04

    sphere._polyhedral_regions = meshers.compile_field(notched_prism)._polyhedral_regions
    mesh = meshers.generate(
        sphere,
        compile=False,
        cells=16,
        bounds=(-0.5, 0.5) * 3,
        geometry_tolerance=0.005,
        minimum_quality=0,
        max_tetrahedra=250000,
    )
    assert "polyhedral_regions" not in mesh.diagnostics
    assert np.max(np.linalg.norm(mesh.points, axis=1)) < 0.21


def test_periodic_reentrant_normal_jump_uses_verified_region_traces():
    def field(x, y, z):
        return np.maximum(np.maximum(abs(x) + 0.3 * abs(z) - 0.35, abs(y) - 0.2),
            np.minimum(x - 0.05, y - 0.05))

    mesh = meshers.generate(
        field, bounds=(-0.5, 0.5) * 3, cells=24, periodic=(False, False, True),
        adapt=True, min_length=0.002, max_length=0.1, gradation=0.3,
        layers=0, adapt_iterations=4, quality_guard=0.3,
        optimize_passes=1, geometry_tolerance=0.002, minimum_quality=0.1,
        max_tetrahedra=250000,
    )
    assert mesh.diagnostics["polyhedral_regions"] == 2
    assert mesh.diagnostics["polyhedral_constraints"] == 8
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    assert abs(mesh.diagnostics["volume"] - 0.18625) < 1e-12
    corners = mesh.points[mesh.tetrahedra]
    assert np.all(np.linalg.det(corners[:, 1:] - corners[:, :1]) > 0)
    boundary = mesh.points[mesh.surface]
    wall = boundary[np.ptp(boundary[:, :, 2], axis=1) > 0]
    probes = np.concatenate((wall.reshape(-1, 3), wall.mean(axis=1)))
    assert np.max(abs(field(*probes.T))) < 1e-12
    pairs = mesh.periodic_pairs[2]
    assert len(pairs) > 0
    np.testing.assert_array_equal(mesh.points[pairs[:, 0], :2], mesh.points[pairs[:, 1], :2])
    mapping = np.arange(len(mesh.points))
    mapping[pairs[:, 1]] = pairs[:, 0]
    low = mesh.surface[np.all(np.equal(boundary[:, :, 2], -0.5), axis=1)]
    high = mesh.surface[np.all(np.equal(boundary[:, :, 2], 0.5), axis=1)]
    assert len(low) > 0
    assert {tuple(sorted(f)) for f in low} == {tuple(sorted(f)) for f in mapping[high]}
