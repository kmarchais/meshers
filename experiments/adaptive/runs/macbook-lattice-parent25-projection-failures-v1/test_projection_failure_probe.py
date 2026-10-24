"""A trimmed plane has a valid target despite its unbounded normal foot."""

import numpy as np

from projection_failure_probe import constrained_foot


def test_box_constrained_plane_projection():
    probe = np.array([.49, .1, 0.])
    unconstrained = probe + .105*np.array([1., 1., 0.])
    assert unconstrained[0] > .5
    candidate, result = constrained_foot(lambda x,y,z:x+y-.8, probe, 0., np.zeros(3, dtype=bool), unconstrained)
    np.testing.assert_allclose(candidate, [.5, .3, 0.], atol=1e-8)
    assert result["feasible_target"]
