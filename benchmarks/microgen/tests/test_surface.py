import numpy as np
import pytest
from test_elasticity import cube

from microgen_benchmark.backends import generate_surface
from microgen_benchmark.cases import CASES
from microgen_benchmark.geometry import boundary, canonical_surface
from microgen_benchmark.reference import boundary_samples
from microgen_benchmark.surface import validate_surface
from microgen_benchmark.surface_cleanup import conform_cap_edges
from microgen_benchmark.surface_intersections import check


def test_closed_cube_surface_volume_and_winding():
    p, t = cube(2)
    faces = boundary(p, t)
    metric = validate_surface(p, faces, CASES[0], target_points=boundary_samples(CASES[0], 4))
    assert np.isclose(metric["volume"], 1)
    assert metric["sampled_reverse_distance"] < 1e-12
    with pytest.raises(ValueError, match="wound"):
        altered = faces.copy()
        altered[0] = altered[0, ::-1]
        validate_surface(p, altered, CASES[0])
    with pytest.raises(ValueError, match="signed"):
        validate_surface(p, faces[:, ::-1], CASES[0])


def test_nonadjacent_coplanar_and_transverse_intersections():
    a = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    for b in [a + [0.1, 0.1, 0], np.array([[0.2, 0.2, -1], [0.2, 0.2, 1], [0.8, 0.2, 0]])]:
        assert check(np.vstack([a, b]), np.array([[0, 1, 2], [3, 4, 5]]))["intersecting_pairs"] == 1
    with pytest.raises(RuntimeError, match="incomplete"):
        check(np.vstack([a, a + [0.1, 0.1, 0]]), np.array([[0, 1, 2], [3, 4, 5]]), max_pairs=0)


@pytest.mark.parametrize("workflow", ["vtk", "vtk-flying-edges"])
def test_real_vtk_surface_avoids_tetrahedra(workflow):
    config = {"threads": 1}
    points, faces, timing = generate_surface(CASES[0], workflow, 4, config)
    points, faces, flips, _ = canonical_surface(points, faces)
    assert flips in {0, len(faces)}
    assert timing["generated_tetrahedra"] == 0
    assert timing["mesh_s"] > 0
    assert np.isclose(validate_surface(points, faces, CASES[0])["volume"], 1)


def test_cap_stitching_preserves_geometry_and_leaves_real_holes_open():
    p, t = cube(1)
    faces = boundary(p, t)
    a, b, c = faces[0]
    p = np.vstack([p, (p[a] + p[b]) / 2])
    midpoint = len(p) - 1
    mismatched = np.vstack([[[a, midpoint, c], [midpoint, b, c]], faces[1:]])
    with pytest.raises(ValueError, match="open"):
        validate_surface(p, mismatched, CASES[0])
    fixed, diagnostics = conform_cap_edges(p, mismatched, 1)
    assert diagnostics["cap_edge_splits"] == 1
    assert np.isclose(validate_surface(p, fixed, CASES[0])["volume"], 1)
    missing_face, _ = conform_cap_edges(p, fixed[1:], 1)
    with pytest.raises(ValueError, match="open"):
        validate_surface(p, missing_face, CASES[0])
