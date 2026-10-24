"""Corner bounds must dominate actual triangle and adjacent-tet quality."""

import numpy as np

from cap_corner_probe import decimal_wedge_bound, polygon_corners, wedge_bound
from cap_quality_bound import bounds


def test_acute_corner_bound_matches_optimal_triangle_and_dominates_other_lengths():
    center = np.zeros(2)
    for degrees in (0.01, 1., 3., 10., 59.):
        radians = np.radians(degrees)
        previous, following = np.array([np.cos(radians), np.sin(radians)]), np.array([1., 0.])
        maximum, convex, angle = wedge_bound(previous, center, following)
        assert convex and abs(angle - degrees) < 1e-12
        assert abs(float(decimal_wedge_bound(previous, center, following)) - maximum) < 1e-14
        for ratio in np.geomspace(.01, 100, 31):
            triangle = np.column_stack((np.array([center, previous, following * ratio]), np.zeros(3)))
            actual = bounds(triangle[None])[0].item()
            assert actual <= maximum + 1e-14
        triangle = np.column_stack((np.array([center, previous, following]), np.zeros(3)))
        assert abs(bounds(triangle[None])[0].item() - maximum) < 1e-14
        squared_face = sum(((triangle[a] - triangle[b])**2).sum() for a, b in ((0, 1), (1, 2), (2, 0)))
        apex = triangle.mean(axis=0)
        apex[2] = np.sqrt(2 * squared_face / 9)
        tetrahedron = np.vstack((triangle, apex))
        determinant = abs(np.linalg.det(tetrahedron[1:] - tetrahedron[:1]))
        squared_tet = sum(((tetrahedron[a] - tetrahedron[b])**2).sum()
            for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)))
        quality = (12 * (determinant / 2)**(2 / 3) / squared_tet)**1.5
        assert abs(quality - maximum) < 1e-14


def test_square_is_orientation_independent_and_not_obstructed():
    points = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    cells = np.array([[0, 1, 2], [0, 2, 3]])
    rows, ambiguous = polygon_corners(points, cells)
    reverse, _ = polygon_corners(points, cells[:, ::-1])
    assert rows == reverse and not ambiguous
    assert len(rows) == 4 and all(row["maximum_incident_cap_quality"] == 1 for row in rows)


def test_reflex_corner_does_not_inherit_the_small_complementary_angle_bound():
    points = np.array([[0., 0.], [2., 0.], [2., 2.], [1., .02], [0., 2.]])
    rows, ambiguous = polygon_corners(points, np.array([[0, 1, 3], [1, 2, 3], [0, 3, 4]]))
    row = next(row for row in rows if row["node"] == 3)
    assert not ambiguous and not row["convex"] and row["interior_angle_deg"] > 300
    assert row["maximum_incident_cap_quality"] == 1
