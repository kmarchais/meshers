"""Target levels and box locks remain hard constraints for the fallback."""

import numpy as np

from box_projection import constrained_foot


def test_locked_box_coordinate_and_target_level_are_preserved():
    probe = np.array([.5, .1, 0.])
    foot, result = constrained_foot(lambda x,y,z:x+y-.8, probe, 0., np.array([True, False, False]), np.array([.6,.2,0.]))
    np.testing.assert_allclose(foot, [.5,.3,0.], atol=1e-8)
    assert result["feasible_target"] and result["normalized_field_residual"] <= 1e-9


def test_unreachable_level_cannot_be_admitted():
    probe = np.zeros(3)
    _, result = constrained_foot(lambda x,y,z:x+y-2., probe, 0., np.zeros(3, dtype=bool), np.ones(3))
    assert not result["feasible_target"]


def test_fixed_corner_is_checked_without_optimizer_or_waiver():
    probe = np.full(3, .5)
    for level, expected in ((1., True), (.8, False)):
        foot, result = constrained_foot(lambda x,y,z:x+y, probe, level, np.ones(3, dtype=bool), probe)
        np.testing.assert_array_equal(foot, probe)
        assert result["feasible_target"] == expected
