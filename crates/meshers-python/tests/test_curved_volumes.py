"""Curved-cell volumes use the final mapping rather than its four corners."""

import h5py
import meshers
import numpy as np
import pytest

EDGES = np.array([(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)])


def reference_nodes():
    corners = np.vstack((np.zeros(3), np.eye(3)))
    return np.vstack((corners, corners[EDGES].mean(axis=1)))


def one_cell(points):
    return meshers.Mesh(
        points=points,
        tetrahedra=np.arange(10, dtype=np.int64)[None, :],
        surface=np.empty((0, 6), dtype=np.int64),
        boundary_tags=np.empty(0, dtype=np.uint8),
        periodic_pairs=tuple(np.empty((0, 2), dtype=np.int64) for _ in range(3)),
        diagnostics={},
    )


def test_export_integrates_current_curved_nodes(tmp_path):
    points = reference_nodes()
    points[4, 2] = 0.1
    mesh = one_cell(points)
    path = tmp_path / "curved.vtkhdf"
    mesh.write_vtkhdf(path)
    # z=t+0.4*L0*L1 gives detJ=1-0.4*r, which stays positive.
    with h5py.File(path) as file:
        data = file["VTKHDF/CellData"]
        assert data["Volume"][0] == pytest.approx(0.15, abs=1e-14)
        assert data["CornerVolume"][0] == pytest.approx(1 / 6, abs=1e-14)
        assert data["CornerMMGQuality"][0] > 0
        assert "MMGQuality" not in data
        assert file["VTKHDF/Types"][0] == 24
    points[4, 2] = 0.05
    mesh.write_vtkhdf(path)
    with h5py.File(path) as file:
        assert file["VTKHDF/CellData/Volume"][0] == pytest.approx(19 / 120, abs=1e-14)


def test_export_integrates_cubic_jacobian(tmp_path):
    points = reference_nodes()
    a, b, c = (0.2, 0.3, 0.4)
    points += np.array([a, b, c]) * points**2
    # Rotate, scale and translate the same known mapping.
    points = 2 * points[:, [1, 2, 0]] + [1e3, -2e3, 3e3]
    exact = 8 * (1 / 6 + (a + b + c) / 12 + (a * b + a * c + b * c) / 30 + a * b * c / 90)
    path = tmp_path / "cubic.vtkhdf"
    one_cell(points).write_vtkhdf(path)
    with h5py.File(path) as file:
        assert file["VTKHDF/CellData/Volume"][0] == pytest.approx(exact, abs=1e-11)


def test_generated_tet10_diagnostics_match_integrated_export(tmp_path):
    mesh = meshers.generate(
        lambda x, y, z: x * x + y * y + z * z - 0.16,
        bounds=(-0.5, 0.5) * 3,
        cells=6,
        order=2,
        geometry_tolerance=1,
        optimize_passes=0,
    )
    path = tmp_path / "sphere.vtkhdf"
    mesh.write_vtkhdf(path)
    with h5py.File(path) as file:
        curved = file["VTKHDF/CellData/Volume"][:].sum()
        linear = file["VTKHDF/CellData/CornerVolume"][:].sum()
    assert mesh.diagnostics["volume"] == pytest.approx(curved, abs=1e-12)
    assert mesh.diagnostics["quadratic"]["linear_volume"] == pytest.approx(linear, abs=1e-12)
    assert curved > 1.03 * linear
    assert mesh.diagnostics["minimum_mmg_quality_basis"] == "linear_corners"
