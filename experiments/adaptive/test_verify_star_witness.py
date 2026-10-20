"""Independent determinant arithmetic must reproduce both witness outcomes."""

import pytest

from star_feasibility import diagnose
from test_star_feasibility import paired_tets
from verify_star_witness import tangent_upper


@pytest.mark.parametrize("target", [.5, .99])
def test_extended_precision_replays_known_feasible_and_obstructed_stars(target):
    points, cells = paired_tets()
    row = diagnose(points, cells, 0, target)
    upper = tangent_upper(points, cells, 0, target, row["normalized_expansion_point"], row["weights"],
        row["physical_origin"], row["scale"])
    assert abs(upper - row["checked_tangent_upper_bound"]) < 1e-12
    assert (upper < -row["required_negative_margin"]) == row["fixed_neighbor_obstruction"]
