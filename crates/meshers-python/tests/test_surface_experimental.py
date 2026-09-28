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


@pytest.mark.skipif(
    not hasattr(meshers._meshers, "generate_surface"),
    reason="requires experimental-surfaces Rust feature",
)
def test_graded_surface_default_is_closed_and_respects_implicit_walls():
    def graded_gyroid(x, y, z):
        thickness = 0.6 + (0.1 / 3) * (x + y + z)
        return gyroid(x, y, z) / (0.5 * thickness)

    surface = meshers.generate_surface(
        graded_gyroid,
        bounds=(-1, 1, -1, 1, -1, 1),
        cells=31,
        band=(-1, 1),
    )
    assert surface.diagnostics["minimum_angle_degrees"] > 10
    mesh = pv.PolyData(
        surface.points,
        np.column_stack((np.full(len(surface.triangles), 3), surface.triangles)),
    )
    assert mesh.n_open_edges == 0
    assert mesh.volume > 0
    for label, level in ((0, 1), (1, -1)):
        wall_nodes = np.unique(surface.triangles[surface.labels == label])
        x, y, z = surface.points[wall_nodes].T
        assert np.max(np.abs(graded_gyroid(x, y, z) - level)) < 1e-8


def split_p(x, y, z):
    x, y, z = (2 * np.pi * v for v in (x, y, z))
    return (
        1.1
        * (
            np.sin(2 * x) * np.cos(y) * np.sin(z)
            + np.sin(2 * y) * np.cos(z) * np.sin(x)
            + np.sin(2 * z) * np.cos(x) * np.sin(y)
        )
        - 0.2
        * (
            np.cos(2 * x) * np.cos(2 * y)
            + np.cos(2 * y) * np.cos(2 * z)
            + np.cos(2 * z) * np.cos(2 * x)
        )
        - 0.4 * (np.cos(2 * x) + np.cos(2 * y) + np.cos(2 * z))
    )


@pytest.mark.skipif(
    not hasattr(meshers._meshers, "generate_surface"),
    reason="requires experimental-surfaces Rust feature",
)
def test_split_p_periodic_surface_pairs_cap_triangles():
    surface = meshers.generate_surface(
        split_p,
        bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
        cells=12,
        band=(-0.25, 0.25),
        periodic=(True,) * 3,
    )
    assert surface.diagnostics["minimum_angle_degrees"] > 5
    mesh = pv.PolyData(
        surface.points,
        np.column_stack((np.full(len(surface.triangles), 3), surface.triangles)),
    )
    assert mesh.n_open_edges == 0
    for label, level in ((0, 0.25), (1, -0.25)):
        wall_nodes = np.unique(surface.triangles[surface.labels == label])
        x, y, z = surface.points[wall_nodes].T
        assert np.max(np.abs(split_p(x, y, z) - level)) < 1e-8
    for axis in range(3):
        caps = []
        for side in range(2):
            triangles = surface.triangles[surface.labels == 2 + 2 * axis + side]
            points = surface.points[triangles].copy()
            points[:, :, axis] = 0
            rounded = np.rint(points * 1e9).astype(np.int64)
            caps.append({tuple(sorted(map(tuple, face))) for face in rounded})
        assert caps[0]
        assert caps[0] == caps[1]


def test_periodic_surface_rejects_unpaired_topology_edits():
    with pytest.raises(ValueError, match="zero smoothing and improvement"):
        meshers.generate_surface(
            gyroid,
            band=(-0.25, 0.25),
            periodic=(True,) * 3,
            improvement_rounds=1,
        )


@pytest.mark.skipif(
    not hasattr(meshers._meshers, "generate_surface"),
    reason="requires experimental-surfaces Rust feature",
)
def test_periodic_gyroid_retries_degenerate_grid_alignment():
    surface = meshers.generate_surface(
        gyroid,
        bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
        cells=24,
        band=(-0.25, 0.25),
        periodic=(True,) * 3,
    )
    assert surface.diagnostics["background_cells"] == 26
    assert surface.diagnostics["resolution_retries"] == 2
    assert surface.diagnostics["minimum_angle_degrees"] > 15
