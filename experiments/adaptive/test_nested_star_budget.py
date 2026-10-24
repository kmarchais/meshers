"""Re-audited checkpoints can retain the cumulative count without legacy fields."""

from nested_star_budget import checkpoint_rounds


def test_reaudited_checkpoint_needs_no_legacy_completed_count():
    assert checkpoint_rounds({"total_rounds": 13}) == 13
    assert checkpoint_rounds({"completed_rounds": 8}) == 8
    assert checkpoint_rounds({"total_rounds": 16, "completed_rounds": 8}) == 16
