"""Repeated cells share seam nodes and retain only outer boundary caps."""

import meshers
import numpy as np
import pytest
import pyvista as pv


def gyroid(x, y, z):
    x, y, z = (2 * np.pi * v for v in (x, y, z))
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


BOUNDS = (-0.5, 0.5) * 3


@pytest.mark.skipif(
    not hasattr(meshers._meshers, "generate_surface"),
    reason="requires experimental-surfaces Rust feature",
)
def test_surface_tiles_close_without_internal_caps():
    cell = meshers.generate_surface(
        gyroid,
        bounds=BOUNDS,
        cells=12,
        band=(-0.25, 0.25),
        periodic=(True, True, True),
    )
    tiled = meshers.tile_periodic(cell, bounds=BOUNDS, repeats=(2, 2, 2))
    assert len(tiled.triangles) < 8 * len(cell.triangles)
    assert tiled.diagnostics["repeated_cells"] == 8
    poly = pv.PolyData(
        tiled.points,
        np.column_stack((np.full(len(tiled.triangles), 3), tiled.triangles)),
    )
    assert poly.n_open_edges == 0
    p = tiled.points[tiled.triangles]
    angles = []
    for axis in range(3):
        a = p[:, (axis + 1) % 3] - p[:, axis]
        b = p[:, (axis + 2) % 3] - p[:, axis]
        angles.append(
            np.degrees(
                np.arctan2(
                    np.linalg.norm(np.cross(a, b), axis=1),
                    np.einsum("ij,ij->i", a, b),
                )
            )
        )
    assert min(np.min(a) for a in angles) == pytest.approx(
        cell.diagnostics["minimum_angle_degrees"], abs=1e-6
    )
    for axis in range(3):
        caps = []
        for side in range(2):
            faces = tiled.triangles[tiled.labels == 2 + 2 * axis + side]
            points = tiled.points[faces].copy()
            points[:, :, axis] = 0
            rounded = np.rint(points * 1e9).astype(np.int64)
            caps.append({tuple(sorted(map(tuple, face))) for face in rounded})
        assert caps[0] == caps[1]


def test_volume_tiles_share_nodes_and_periodic_outer_faces():
    cell = meshers.generate(
        gyroid,
        bounds=BOUNDS,
        cells=12,
        band=(-0.25, 0.25),
        periodic=(True, True, True),
        geometry_tolerance=0.02,
        optimize_passes=2,
    )
    tiled = meshers.tile_periodic(cell, bounds=BOUNDS, repeats=(2, 2, 2))
    assert len(tiled.tetrahedra) == 8 * len(cell.tetrahedra)
    assert len(tiled.points) < 8 * len(cell.points)
    assert tiled.diagnostics["volume"] == pytest.approx(8 * cell.diagnostics["volume"])
    assert (
        tiled.diagnostics["minimum_mmg_quality"]
        == cell.diagnostics["minimum_mmg_quality"]
    )
    tetrahedra = tiled.points[tiled.tetrahedra]
    determinants = np.einsum(
        "ij,ij->i",
        tetrahedra[:, 1] - tetrahedra[:, 0],
        np.cross(
            tetrahedra[:, 2] - tetrahedra[:, 0],
            tetrahedra[:, 3] - tetrahedra[:, 0],
        ),
    )
    assert np.min(determinants) > 0
    edge_sum = sum(
        np.sum((tetrahedra[:, i] - tetrahedra[:, j]) ** 2, axis=1)
        for i in range(4)
        for j in range(i + 1, 4)
    )
    mmg_quality = np.sqrt(432 * determinants**2 / edge_sum**3)
    assert np.min(mmg_quality) == pytest.approx(
        cell.diagnostics["minimum_mmg_quality"], abs=1e-8
    )
    for axis, pairs in enumerate(tiled.periodic_pairs):
        assert len(pairs) > 0
        delta = tiled.points[pairs[:, 1]] - tiled.points[pairs[:, 0]]
        expected = np.tile(np.eye(3)[axis] * 2, (len(pairs), 1))
        np.testing.assert_allclose(delta, expected, rtol=0, atol=1e-9)
    poly = pv.PolyData(
        tiled.points,
        np.column_stack((np.full(len(tiled.surface), 3), tiled.surface)),
    )
    assert poly.n_open_edges == 0


def test_axis_graded_volume_tiles_only_in_periodic_directions():
    bounds = (-1.0, 1.0, -0.5, 0.5, -0.5, 0.5)

    def graded(x, y, z):
        return gyroid(x, y, z) / (0.3 + 0.05 * x)

    slab = meshers.generate(
        graded,
        bounds=bounds,
        cells=(23, 11, 11),
        band=(-1, 1),
        periodic=(False, True, True),
        geometry_tolerance=0.02,
        minimum_quality=0,
        optimize_passes=2,
    )
    tiled = meshers.tile_periodic(slab, bounds=bounds, repeats=(1, 2, 2))

    assert len(tiled.tetrahedra) == 4 * len(slab.tetrahedra)
    assert len(tiled.periodic_pairs[0]) == 0
    for axis in (1, 2):
        pairs = tiled.periodic_pairs[axis]
        assert len(pairs) > 0
        displacement = tiled.points[pairs[:, 1]] - tiled.points[pairs[:, 0]]
        np.testing.assert_allclose(displacement[:, axis], 2.0, rtol=0, atol=1e-9)
        np.testing.assert_allclose(
            displacement[:, [i for i in range(3) if i != axis]],
            0.0,
            rtol=0,
            atol=1e-9,
        )
    assert tiled.diagnostics["minimum_mmg_quality"] == pytest.approx(
        slab.diagnostics["minimum_mmg_quality"]
    )
