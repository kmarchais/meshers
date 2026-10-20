"""Numerical geometry gates cannot admit an incoherent oriented tet complex."""

import pytest

from microgen_workflow import geometry_is_admitted


def row():
    return {"sampled_two_sided_max": .001, "relative_volume_error": .001, "q_min_corner": .2,
        "projection_failed": 0, "minimum_sampled_jacobian": 1., "oriented_facet_consistency_passed": True}


def test_incoherent_complex_fails_even_when_all_numerical_gates_pass():
    record = row()
    assert geometry_is_admitted(record)
    record["oriented_facet_consistency_passed"] = False
    assert not geometry_is_admitted(record)


@pytest.mark.parametrize("change", [{"sampled_two_sided_max": .002001}, {"relative_volume_error": .010001},
    {"q_min_corner": .099999}, {"projection_failed": 1}, {"minimum_sampled_jacobian": 0.}])
def test_original_numerical_gates_stay_strict(change):
    record = row()
    record.update(change)
    assert not geometry_is_admitted(record)
