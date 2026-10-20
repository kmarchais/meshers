"""Positive tet determinants alone must not hide same-side face overlaps."""

from collections import Counter
import itertools

import numpy as np

from facet_incidence_audit import audit


def surface(cells):
    counts = Counter(tuple(sorted(face)) for cell in cells for face in itertools.combinations(cell, 3))
    return np.array([face for face, count in counts.items() if count == 1])


def test_valid_bipyramid_has_opposite_face_orientations():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [.2, .2, 1.], [.2, .2, -1.]])
    cells = np.array([[0, 1, 2, 3], [0, 2, 1, 4]])
    row = audit(points, cells, surface(cells))
    assert row["oriented_facet_consistency_passed"] and row["same_orientation_shared_faces"] == 0


def test_positive_overlapping_tets_are_rejected():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [.2, .2, 1.], [.2, .2, .5]])
    cells = np.array([[0, 1, 2, 3], [0, 1, 2, 4]])
    row = audit(points, cells, surface(cells))
    assert row["nonpositive_corner_determinants"] == 0
    assert row["same_orientation_shared_faces"] == 1 and not row["oriented_facet_consistency_passed"]
    assert row["conflict_examples"][0]["same_physical_side"]


def test_missing_boundary_face_is_rejected():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    cells = np.arange(4)[None]
    row = audit(points, cells, surface(cells)[:-1])
    assert not row["provided_boundary_matches_extraction"] and not row["oriented_facet_consistency_passed"]
