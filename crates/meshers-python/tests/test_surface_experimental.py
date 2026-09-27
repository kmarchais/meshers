"""Checks for the feature-gated native surface-only generator."""

import meshers
import numpy as np
import pytest
import pyvista as pv


def gyroid(x, y, z):
    x, y, z = (2 * np.pi * v for v in (x, y, z))
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


@pytest.mark.skipif(
    not hasattr(meshers._meshers, "generate_surface"),
    reason="requires experimental-surfaces Rust feature",
)
def test_generates_closed_surface_without_tetrahedra():
    surface = meshers.generate_surface(
        gyroid,
        bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
        cells=12,
        band=(-0.25, 0.25),
        smoothing_iterations=0,
        improvement_rounds=0,
        polish_passes=0,
    )
    assert surface.points.shape[1] == 3
    assert surface.triangles.shape[1] == 3
    assert len(surface.labels) == len(surface.triangles)
    assert set(surface.labels) == set(range(8))
    assert surface.diagnostics["minimum_angle_degrees"] > 0
    mesh = pv.PolyData(
        surface.points,
        np.column_stack((np.full(len(surface.triangles), 3), surface.triangles)),
    )
    assert mesh.n_open_edges == 0
    assert mesh.volume > 0


def test_rejects_periodic_optimization_until_pairs_are_preserved():
    with pytest.raises(NotImplementedError, match="periodic"):
        meshers.generate_surface(gyroid, band=(-0.25, 0.25), periodic=(True,) * 3)
