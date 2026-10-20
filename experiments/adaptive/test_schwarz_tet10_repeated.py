"""Ordering, warmup exclusion, and admission controls for repeated timings."""

from copy import deepcopy

import pytest

from schwarz_tet10_live import STAGE_TIMEOUTS, stage_names
from schwarz_tet10_repeated import schedule, summarize


def pairs():
    result = []
    for phase, native_first in schedule():
        times = {"microgen": 1000. if phase == "warmup" else 20., "meshers": 1000. if phase == "warmup" else 10.}
        result.append({"phase": phase, "native_first": native_first, "report": {
            "both_pass_declared_reference_models": True, "timing_clock_consistent": True,
            "stiffness_target": .001, "science_audited_workflow_s": times,
            "generation_lift_solve_s": {arm: value / 2 for arm, value in times.items()}}})
    return result


def test_stage_order_preserves_fresh_parent_dependency_and_alternates_arms():
    for native_first in (False, True):
        names = stage_names(native_first)
        assert set(names) == set(STAGE_TIMEOUTS) and len(names) == len(STAGE_TIMEOUTS)
        assert names[0] == "controls"
        assert names.index("microgen-parent") < names.index("microgen-tet10")
        assert (names.index("meshers-tet10") < names.index("microgen-parent")) == native_first
    assert [native_first for _, native_first in schedule()[1:]] == [True, False, True, False, True]


def test_repeat_summary_excludes_warmup_and_retains_both_costs():
    report = summarize(pairs())
    assert report["statistics"]["science_audited_workflow_s"]["microgen"]["median"] == 20.
    assert report["statistics"]["generation_lift_solve_s"]["meshers"]["median"] == 5.
    assert report["paired_ratios_microgen_over_meshers"]["science_audited_workflow_s"] == [2.] * 5
    assert report["meshers_faster_in_every_measured_pair"]


@pytest.mark.parametrize("key", ["both_pass_declared_reference_models", "timing_clock_consistent"])
def test_failed_or_clock_interrupted_pair_cannot_support_repeated_claim(key):
    rows = deepcopy(pairs())
    rows[3]["report"][key] = False
    with pytest.raises(AssertionError):
        summarize(rows)


def test_one_slower_pair_is_not_hidden_by_median():
    rows = pairs()
    rows[-1]["report"]["generation_lift_solve_s"]["meshers"] = 11.
    assert not summarize(rows)["meshers_faster_in_every_measured_pair"]
