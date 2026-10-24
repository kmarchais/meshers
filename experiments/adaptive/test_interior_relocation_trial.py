"""Whole-mesh repair cannot move boundary vertices or accept new bad tets."""

import itertools

import numpy as np
import pytest

from interior_relocation_trial import repair, validate_repair_source
from test_facet_incidence_audit import surface
from test_interior_star_probe import star


def test_lattice_repair_requires_declared_actual_route_and_geometry():
    import copy

    source = {"case": "smooth-lattice", "route": "microgen-vtk-mmg-first-pass-override", "resolution": 51,
        "mmg_calls": [{"nofem": True, "hausd": .0015, "hmax": .06, "hgrad": 1.3},
            {"hausd": .001, "hgrad": 1.3, "hmax": .06, "hmin": None, "hsiz": None, "ls": True, "nr": True}],
        "periodic": {}, "geometry": {"sampled_two_sided_max": .0017, "relative_volume_error": -.0049}}
    validate_repair_source(source)
    for key, value in (("resolution", 52), ("route", "meshers-seed"),
        ("geometry", {"sampled_two_sided_max": .0021, "relative_volume_error": -.0049})):
        changed = copy.deepcopy(source); changed[key] = value
        with pytest.raises(AssertionError):
            validate_repair_source(changed)
    changed = copy.deepcopy(source); changed["mmg_calls"][1]["hausd"] = .002
    with pytest.raises(AssertionError):
        validate_repair_source(changed)


def test_complete_interior_star_repairs_with_fixed_boundary():
    points, cells = star()
    moved, report = repair(points, cells, cells[:, 1:], maximum_attempts=4)
    assert report["initial_bad_elements"] > 0 and report["remaining_bad_elements"] == 0
    assert report["accepted_moves"] == 1
    np.testing.assert_array_equal(moved[1:], points[1:])


def test_all_boundary_bad_tet_stays_rejected():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [.5, .5, .001]])
    cells = np.arange(4)[None]
    surface = np.array(list(itertools.combinations(range(4), 3)))
    moved, report = repair(points, cells, surface)
    assert report["remaining_bad_elements"] == 1 and report["accepted_moves"] == 0
    assert report["attempts"][0]["status"] == "all-boundary"
    np.testing.assert_array_equal(moved, points)


def test_incoherent_source_is_rejected_before_relocation():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [.2, .2, 1.], [.2, .2, .5]])
    cells = np.array([[0, 1, 2, 3], [0, 1, 2, 4]])
    with pytest.raises(ValueError, match="oriented facet"):
        repair(points, cells, surface(cells))
