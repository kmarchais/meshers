"""Background conversion must preserve exact clipping inputs and staged caps."""

import numpy as np
import pyvista as pv
import pytest

from preclip_probe import sampled_inputs, triangulate_background


def grid():
    background = pv.ImageData(dimensions=(3, 3, 3), spacing=(.31, .31, .31), origin=(-.5, -.5, -.5))
    points = background.points
    background["surface"] = points[:, 0] + 2 * points[:, 1] - points[:, 2]
    background["lower_surface"] = background["surface"] + .55
    background["upper_surface"] = background["surface"] - .55
    return background


@pytest.mark.parametrize("tetra_per_cell", [5, 6])
def test_triangulation_preserves_exact_field_and_point_arrays(tetra_per_cell):
    background = grid()
    before = sampled_inputs(background)
    converted = triangulate_background(background, tetra_per_cell=tetra_per_cell)
    assert sampled_inputs(converted) == before and converted.n_cells > background.n_cells
    np.testing.assert_array_equal(converted.points, background.points)


def test_five_tet_conversion_rejects_nonrectilinear_points():
    background = grid().cast_to_structured_grid()
    background.points[0, 0] += .01
    with pytest.raises(AssertionError):
        triangulate_background(background, tetra_per_cell=5)


def test_five_tet_conversion_preserves_permuted_point_order():
    background = grid().cast_to_structured_grid()
    background.points = background.points[:, [1, 0, 2]].copy()
    before = sampled_inputs(background)
    converted = triangulate_background(background, tetra_per_cell=5)
    assert sampled_inputs(converted) == before
    assert converted.n_cells == 5 * background.n_cells


def test_background_cap_is_checked_before_conversion():
    with pytest.raises(ValueError, match="Background tetrahedron upper bound"):
        triangulate_background(grid(), cap=1)
