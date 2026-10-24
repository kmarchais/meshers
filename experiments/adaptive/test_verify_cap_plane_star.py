"""Independent cap-only replay reproduces a known infeasible four-triangle fan."""

from decimal import Decimal

import numpy as np

from verify_cap_plane_star import replay


def test_symmetric_cap_fan_certificate_has_the_analytic_upper_bound():
    points = np.array([[-.5, 0, 0], [-.5, .5, 0], [-.5, 0, .5],
        [-.5, -.5, 0], [-.5, 0, -.5], [.5, 0, 0]])
    row = {"nodes": [0, 5], "axis": 0, "scale": 1., "normalized_expansion_point": [0., 0.],
        "weights": [.25] * 4, "incident_cap_triangles": [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]]}
    upper = replay(points, row, Decimal(".99"))
    assert upper < 0 and abs(float(upper) - (np.sqrt(3) / 2 - .99)) < 1e-14
