"""Curved boundary, volume, role, count and overlap checks of edge removal."""

import microgen_workflow as w


def test_disjoint_curved_quality_pass_preserves_its_boundary_and_handles_noop():
    w.imports()
    from quadratic_quality_preflight import check_application
    check_application()
