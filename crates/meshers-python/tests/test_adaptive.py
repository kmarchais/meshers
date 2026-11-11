import meshers
import numpy as np
import pytest


def sphere(x, y, z):
    return np.sqrt(x * x + y * y + z * z) - 0.35


def gyroid(x, y, z):
    a, b, c = 2 * np.pi * x, 2 * np.pi * y, 2 * np.pi * z
    return (np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)) / 0.55


BOX = (-0.5, 0.5) * 3


@pytest.mark.skipif(not hasattr(meshers, "generate_surface"), reason="surface build feature")
def test_adaptive_surface_meets_tolerance_with_fewer_triangles():
    tolerance = 1e-3
    fixed = meshers.generate_surface(sphere, bounds=BOX, cells=32, quality=True)
    adapted = meshers.generate_surface(
        sphere, bounds=BOX, cells=32, geometry_tolerance=tolerance, adaptive=True
    )
    assert adapted.diagnostics["maximum_sampled_surface_error"] <= tolerance
    assert len(adapted.triangles) < len(fixed.triangles)
    assert adapted.diagnostics["minimum_triangle_shape"] > 0.3
    exact = 4 / 3 * np.pi * 0.35**3
    assert abs(adapted.diagnostics["volume"] / exact - 1) < 0.01


@pytest.mark.skipif(not hasattr(meshers, "generate_surface"), reason="surface build feature")
@pytest.mark.parametrize("balanced", [False, True])
def test_fast_surface_meets_tolerance_with_fewer_triangles(balanced):
    tolerance = 1e-3
    fixed = meshers.generate_surface(gyroid, bounds=BOX, cells=24, band=(-1.0, 1.0))
    fast = meshers.generate_surface(
        gyroid,
        bounds=BOX,
        cells=24,
        band=(-1.0, 1.0),
        geometry_tolerance=tolerance,
        adaptive=True,
        iterations=2,
        balanced=balanced,
        fast=True,
    )
    assert fast.diagnostics["maximum_sampled_surface_error"] <= tolerance
    assert len(fast.triangles) < 0.8 * len(fixed.triangles)
    assert fast.diagnostics["minimum_triangle_shape"] > 0


def test_adapted_periodic_sheet_stays_periodic_and_shaped():
    mesh = meshers.generate(
        gyroid,
        bounds=BOX,
        cells=16,
        band=(-1.0, 1.0),
        periodic=(True, True, True),
        geometry_tolerance=4e-3,
        minimum_quality=0.1,
        adapt=True,
        max_length=0.1,
    )
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    assert mesh.diagnostics["sampled_surface_error"] <= 4e-3
    for axis, pairs in enumerate(mesh.periodic_pairs):
        assert len(pairs) > 0
        shift = mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]]
        np.testing.assert_allclose(shift, np.eye(3)[np.full(len(pairs), axis)], atol=1e-12)


def test_quadratic_wall_nodes_lie_on_the_surface(tmp_path):
    mesh = meshers.generate(sphere, bounds=BOX, cells=10, geometry_tolerance=0.1, order=2)
    assert mesh.tetrahedra.shape[1] == 10
    assert mesh.surface.shape[1] == 6
    assert mesh.diagnostics["quadratic"]["minimum_jacobian_ratio"] > 0
    mid = mesh.points[mesh.surface[:, 3:]].reshape(-1, 3)
    assert np.max(np.abs(sphere(*mid.T))) < 1e-12
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "quadratic.vtkhdf"
    mesh.write_vtkhdf(path)
    with h5py.File(path) as file:
        assert set(file["VTKHDF/Types"][()]) == {24}


def test_adapted_smooth_lattice_preserves_quality_and_periodic_pairs():
    corners = (
        (-0.5, -0.5, -0.5), (-0.5, -0.5, 0.5),
        (-0.5, 0.5, -0.5), (-0.5, 0.5, 0.5),
        (0.5, -0.5, -0.5), (0.5, -0.5, 0.5),
        (0.5, 0.5, -0.5), (0.5, 0.5, 0.5),
    )

    def field(coordinates_x, coordinates_y, coordinates_z):
        distance = None
        for corner in corners:
            parameter = np.clip(
                (coordinates_x * corner[0] + coordinates_y * corner[1]
                 + coordinates_z * corner[2]) / 0.75, 0, 1
            )
            radius = np.sqrt(
                (coordinates_x - parameter * corner[0]) ** 2
                + (coordinates_y - parameter * corner[1]) ** 2
                + (coordinates_z - parameter * corner[2]) ** 2
            )
            if distance is None:
                distance = radius
            else:
                difference = np.abs(distance - radius)
                blend = np.maximum(0.08 - difference, 0.0) / 0.08
                distance = (distance + radius - difference) / 2 - blend ** 3 * 0.08 / 6
        return np.minimum(distance - 0.12, 0.02)

    mesh = meshers.generate(
        field, bounds=BOX, cells=20, periodic=(True, True, True), threads=1,
        optimize_passes=1, minimum_quality=0.1, geometry_tolerance=1.0,
        max_tetrahedra=250_000, adapt=True, min_length=0.002, max_length=0.1,
        gradation=0.3, adapt_iterations=4, quality_guard=0.3,
        adapt_tolerance=0.004, order=2,
    )
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    assert mesh.diagnostics["elements_below_quality_01"] == 0
    assert mesh.diagnostics["quadratic"]["minimum_jacobian_ratio"] > 1e-3
    for axis, pairs in enumerate(mesh.periodic_pairs):
        assert len(pairs) > 0
        shift = mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]]
        np.testing.assert_allclose(shift, np.eye(3)[np.full(len(pairs), axis)], atol=1e-12)
