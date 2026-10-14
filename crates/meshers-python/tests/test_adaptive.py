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
