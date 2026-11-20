# Copyright (c) 2026 Kevin Marchais

"""Exercise VTKHDF export with independent PyVista and HDF readers."""

import meshers
import numpy as np
import pyfive
import pytest
import pyvista as pv


def test_vtkhdf_unwritable_path_raises(tmp_path):
    mesh = meshers.generate("gyroid", cells=4, geometry_tolerance=1)
    with pytest.raises(OSError, match="could not write"):
        mesh.write_vtkhdf(tmp_path / "missing-directory" / "mesh.vtkhdf")


def test_vtkhdf_round_trip_and_physical_shifts(tmp_path):
    result = meshers.generate(
        lambda x, _y, _z: np.full_like(x, -1.0),
        compile_expression=False,
        cells=4,
        bounds=(2, 4, -1, 2, 3, 7),
        periodic=(True, True, True),
        optimize_passes=0,
    )
    path = tmp_path / "mesh.vtkhdf"
    result.write_vtkhdf(path)
    restored = pv.read(path)
    np.testing.assert_array_equal(restored.points, result.points)
    np.testing.assert_array_equal(restored.cells.reshape(-1, 5)[:, 1:], result.tetrahedra)
    np.testing.assert_allclose(restored.point_data["PeriodicShift"].max(axis=0), [2, 3, 4])
    assert restored.cell_data["MMGQuality"].min() > 0
    surface = tmp_path / "surface.vtkhdf"
    result.write_vtkhdf(surface, surface=True)
    assert pv.read(surface).n_cells == len(result.surface)


@pytest.mark.parametrize(("order", "width", "face_width"), [(1, 4, 3), (2, 10, 6)])
def test_empty_mesh_shapes(order, width, face_width, tmp_path):
    result = meshers.generate(
        lambda x, _y, _z: np.ones_like(x),
        compile_expression=False,
        cells=4,
        order=order,
        expected_components=0,
        periodic=(True,) * 3,
        curved_fit=order == 2,
    )
    assert result.points.shape == (0, 3)
    assert result.tetrahedra.shape == (0, width)
    assert result.surface.shape == (0, face_width)
    assert result.diagnostics["empty"]
    assert result.diagnostics["face_connected_components"] == 0
    assert all(p.shape == (0, 2) for p in result.periodic_pairs)
    for surface in (False, True):
        path = tmp_path / "empty.vtkhdf"
        result.write_vtkhdf(path, surface=surface)
        with pyfive.File(path) as f:
            group = f["VTKHDF"]
            np.testing.assert_array_equal(group["Offsets"][:], [0])
            assert group["Connectivity"].shape == group["Types"].shape == (0,)
            if order == 2 and not surface:
                assert group["CellData/CornerVolume"].shape == (0,)
