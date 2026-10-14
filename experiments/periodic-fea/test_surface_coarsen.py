"""Independent topology, partition and cap checks for surface contraction."""
import numpy as np

from microgen_benchmark.geometry import boundary, canonical_mesh
from surface_coarsen import coarsen
from test_body_refine import cube


def split_edge(points, tets, a, b):
    mask = np.any(tets == a, axis=1) & np.any(tets == b, axis=1)
    children = np.tile(tets[mask], (2, 1))
    positions = np.r_[np.argmax(tets[mask] == a, axis=1), np.argmax(tets[mask] == b, axis=1)]
    children[np.arange(len(children)), positions] = len(points)
    return canonical_mesh(np.vstack([points, (points[a] + points[b]) / 2]),
                          np.vstack([tets[~mask], children]))[:2]


def cells(tets):
    ordered = np.sort(tets, axis=1)
    return ordered[np.lexsort(ordered.T[::-1])]


def test_flat_boundary_midpoint_contraction_restores_original_volume_mesh():
    original_points, original_tets = cube()
    original_points *= .8
    points, tets = split_edge(original_points, original_tets, 0, 1)
    result, result_tets, retained, info = coarsen(points, tets, boundary(points, tets), 1, .45)
    assert info["removed_nodes"] == 1
    np.testing.assert_array_equal(result, original_points)
    np.testing.assert_array_equal(retained, np.arange(len(original_points)))
    np.testing.assert_array_equal(cells(result_tets), cells(original_tets))
    np.testing.assert_array_equal(cells(boundary(result, result_tets)),
                                  cells(boundary(original_points, original_tets)))
    volumes = np.linalg.det(result[result_tets[:, 1:]] - result[result_tets[:, :1]]) / 6
    assert volumes.min() > 0
    np.testing.assert_allclose(volumes.sum(), .8**3, atol=1e-14)


def test_periodic_cap_midpoint_is_protected():
    points, tets = cube()
    points, tets = split_edge(points, tets, 0, 1)
    result, result_tets, retained, info = coarsen(points, tets, boundary(points, tets), 1, .6)
    assert info["removed_nodes"] == 0
    np.testing.assert_array_equal(result, points)
    np.testing.assert_array_equal(result_tets, tets)
    np.testing.assert_array_equal(retained, np.arange(len(points)))


def test_link_condition_prevents_deleting_a_material_component():
    points = np.array([[0, 0, 0], [.2, 0, 0], [0, .2, 0], [0, 0, .2]])
    tets = np.array([[0, 1, 2, 3]])
    result, result_tets, _, info = coarsen(points, tets, boundary(points, tets), 1, .3)
    assert info["removed_nodes"] == 0
    assert info["rejections"]["link"] > 0
    np.testing.assert_array_equal(result, points)
    np.testing.assert_array_equal(result_tets, tets)


def test_geometry_guard_can_reject_an_otherwise_valid_boundary_collapse():
    points, tets = cube()
    points *= .8
    points, tets = split_edge(points, tets, 0, 1)
    result, result_tets, _, info = coarsen(points, tets, boundary(points, tets), 1, .45,
                                       geometry_check=lambda triangles: False)
    assert info["removed_nodes"] == 0
    assert info["rejections"]["geometry"] > 0
    np.testing.assert_array_equal(result, points)
    np.testing.assert_array_equal(result_tets, tets)


def test_analytic_plane_guard_rejects_reversed_and_distant_wall_triangles():
    from types import SimpleNamespace
    from run_surface_coarsening import wall_guard
    case = SimpleNamespace(field=lambda x, y, z: 1 + x)
    check = wall_guard(case, .001)
    triangles = np.array([[[0, 0, 0], [0, .2, 0], [0, 0, .2]]])
    assert check(triangles)
    assert not check(triangles[:, ::-1])
    moved = triangles.copy()
    moved[:, :, 0] += .003
    assert not check(moved)


def test_memory_monitor_includes_an_actual_allocating_child():
    import subprocess
    import sys
    import psutil
    from run_surface_coarsening import owned_process_rss
    child = "import time; data=bytearray(64*1024*1024); data[::4096]=b'x'*len(data[::4096]); print('ready', flush=True); time.sleep(20)"
    parent = "import subprocess,sys; p=subprocess.Popen([sys.executable,'-c',sys.argv[1]], stdout=subprocess.PIPE, text=True); print(p.stdout.readline().strip(), flush=True); p.wait()"
    process = subprocess.Popen([sys.executable, "-c", parent, child], stdout=subprocess.PIPE, text=True)
    try:
        assert process.stdout.readline().strip() == "ready"
        assert owned_process_rss(psutil.Process(process.pid)) >= 64*1024*1024
    finally:
        try:
            children = psutil.Process(process.pid).children(recursive=True)
        except psutil.NoSuchProcess:
            children = []
        for member in children:
            try:
                member.kill()
            except psutil.NoSuchProcess:
                pass
        if process.poll() is None:
            process.kill()
        process.wait(timeout=10)
