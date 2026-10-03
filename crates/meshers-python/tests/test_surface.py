"""Exercise the optional surface API and a real STL export."""

import meshers
import numpy as np
import pytest
import pyvista as pv

pytestmark = pytest.mark.skipif(
    not hasattr(meshers, "generate_surface"), reason="optional surface feature"
)


@pytest.mark.parametrize("quality", [False, True])
def test_closed_sphere_can_be_exported_to_stl(tmp_path, quality):
    def sphere(x, y, z):
        return (x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 - 0.35**2

    mesh = meshers.generate_surface(sphere, cells=12, quality=quality)
    assert mesh.triangles.ndim == 2 and mesh.triangles.shape[1] == 3
    faces = np.column_stack((np.full(len(mesh.triangles), 3), mesh.triangles))
    surface = pv.PolyData(mesh.points, faces)
    assert surface.n_open_edges == 0
    assert surface.volume > 0
    assert np.max(abs(sphere(*mesh.points.T))) < 1e-10
    filename = tmp_path / "sphere.stl"
    surface.save(filename)
    restored = pv.read(filename).clean()
    assert restored.n_open_edges == 0
    assert np.isclose(restored.volume, surface.volume, rtol=1e-6)
    assert mesh.diagnostics["experimental"]


def test_surface_rejects_insufficient_geometry_accuracy():
    with pytest.raises(meshers.MeshingError, match="tolerance"):
        meshers.generate_surface(
            lambda x, y, z: (x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 - 0.35**2,
            cells=8,
            geometry_tolerance=1e-8,
        )
