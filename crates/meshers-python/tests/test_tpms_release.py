import h5py
import meshers
import numpy as np
import pytest


@pytest.mark.parametrize("field", ["gyroid", "schwarz_p", "schwarz_d"])
def test_tpms_wheel_round_trip(field, tmp_path):
    mesh = meshers.generate(
        field,
        cells=12,
        band=(-0.5, 0.5),
        periodic=(True, True, True),
        geometry_tolerance=0.1,
        optimize_passes=1,
        threads=2,
    )
    assert len(mesh.tetrahedra) > 0
    assert np.isfinite(mesh.points).all()
    vertices = mesh.points[mesh.tetrahedra]
    volumes = np.linalg.det(vertices[:, 1:] - vertices[:, :1]) / 6
    assert np.all(volumes > 0)
    for axis, pairs in enumerate(mesh.periodic_pairs):
        assert len(pairs) > 0
        shifts = mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]]
        np.testing.assert_allclose(
            np.abs(shifts), np.broadcast_to(np.eye(3)[axis], shifts.shape), atol=1e-10
        )
    path = tmp_path / f"{field}.vtkhdf"
    mesh.write_vtkhdf(path)
    with h5py.File(path) as exported:
        grid = exported["VTKHDF"]
        np.testing.assert_array_equal(grid["Points"][:], mesh.points)
        np.testing.assert_array_equal(grid["Connectivity"][:], mesh.tetrahedra.ravel())
        np.testing.assert_allclose(grid["CellData/Volume"][:], volumes)
        assert np.all(grid["CellData/MMGQuality"][:] > 0)
