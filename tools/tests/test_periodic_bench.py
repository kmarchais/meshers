# Copyright (c) 2026 Kevin Marchais

"""Known-answer controls for the independent measurements of the periodic bench."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "demos" / "solver"))
pytest.importorskip("pyvista")
pytest.importorskip("meshers")
import periodic_bench as bench


def test_regular_tetrahedron_shape():
    points = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    volume, quality, dihedral = bench.tet_shapes(points, np.array([[0, 2, 1, 3]]))
    assert volume[0] == pytest.approx(8 / 3)
    assert quality[0] == pytest.approx(1)
    assert dihedral == pytest.approx(np.degrees(np.arccos(1 / 3)))
    assert bench.tet_shapes(points, np.array([[0, 1, 2, 3]]))[0][0] < 0


def test_sliver_has_poor_quality_and_extreme_dihedrals():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0.01]], dtype=float)
    _, quality, dihedral = bench.tet_shapes(points, np.array([[0, 1, 2, 3]]))
    assert quality[0] < 0.02
    assert dihedral.min() < 2
    assert dihedral.max() > 178


@pytest.mark.parametrize("scale", [1.0, 7.0])
def test_mesh_to_wall_measures_a_chord_whatever_the_field_scaling(scale):
    radius = 0.3

    def sphere(p):
        return scale * ((p**2).sum(axis=1) - radius**2)

    points = radius * np.eye(3)
    distance = bench.mesh_to_wall(sphere, points, np.array([[0, 1, 2]]))
    assert distance[0] == pytest.approx(radius * (1 - 1 / np.sqrt(3)), rel=1e-6)


def test_wall_reference_lies_on_the_wall_inside_the_cell():
    field = bench.CASES["sphere"]
    reference = bench.wall_reference(field, grid=24)
    assert len(reference) > 1000
    assert abs(reference).max() <= 0.5
    assert abs(field(reference)).max() < 1e-8


def test_wall_to_mesh_sees_a_missing_feature():
    field = bench.spheres(np.array([[-0.2, 0, 0], [0.2, 0, 0]]), 0.1)
    reference = bench.wall_reference(field, grid=24)
    points = np.array([[-0.3, 0, 0], [-0.2, 0.1, 0], [-0.2, 0, 0.1]])
    distance = bench.wall_to_mesh(reference, points, np.array([[0, 1, 2]]))
    assert distance.max() > 0.3


def test_reference_volume_of_a_sphere():
    volume = bench.reference_volume(bench.CASES["sphere"], grid=64)
    assert volume == pytest.approx(4 / 3 * np.pi * 0.35**3, rel=2e-3)


def test_cap_mismatch_separates_inexact_and_missing_images():
    low = np.array([[-0.5, 0.1, 0.2], [-0.5, 0.3, 0.2], [-0.5, 0.1, 0.4]])
    points = np.concatenate([low, low + np.array([1.0, 0, 0])])
    surface = np.array([[0, 1, 2], [3, 5, 4]])
    assert bench.cap_mismatch(points, surface) == (0, 0, 0.0)
    points[4, 2] = np.nextafter(0.2, 1)
    unmatched, inexact, offset = bench.cap_mismatch(points, surface)
    assert (unmatched, inexact) == (0, 2)
    assert 0 < offset < 1e-15
    points[4, 2] = 0.25
    assert bench.cap_mismatch(points, surface)[:2] == (2, 0)
    # A wall facet that only touches the face is not a cap.
    assert bench.cap_mismatch(points, np.array([[0, 1, 5]])) == (0, 0, 0.0)


def test_every_case_is_periodic_and_has_solid_and_void():
    probes = np.random.default_rng(5).uniform(-0.5, 0.5, (4000, 3))
    for name, field in bench.CASES.items():
        values = field(probes)
        assert (values < 0).any(), name
        assert (values > 0).any(), name
        shifted = field(probes + np.array([1.0, -1.0, 1.0]))
        assert np.allclose(values, shifted, atol=1e-9), name


def test_measure_records_a_good_mesh_and_a_rejection():
    good = bench.measure("sphere", 12, adapt=False, tolerance=0.02, quality=0.1)
    assert good["missed_limits"] == ["volume"]
    assert good["unmatched_cap_nodes"] == 0
    assert good["wall_to_mesh"] < 0.02
    assert good["field_points"] > 0
    rejected = bench.measure("sphere", 12, adapt=False, tolerance=1e-5, quality=0.1)
    assert "sampled geometry error" in rejected["rejected"]
