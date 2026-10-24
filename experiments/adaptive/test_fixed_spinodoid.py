"""The comparison specimen is the frozen Fourier field, across FFT backends."""

import json
import numpy as np
import pytest

import microgen_workflow as workflow


def test_specimen_uses_exact_frozen_coefficients_before_sampling():
    _, shape, _ = workflow.specimen("spinodoid")
    frep = shape._frep
    frozen = json.loads((workflow.HERE / "runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json").read_text())
    actual = np.column_stack((frep.modes / frep.cell_size, frep.amplitudes, frep.phases))
    np.testing.assert_array_equal(actual, frozen["coefficients"])
    assert frep.threshold == frozen["threshold"] and shape.resolution == 16
    assert "grid" not in shape.__dict__


def test_restoration_rejects_a_different_realization():
    _, shape, _ = workflow.specimen("spinodoid")
    shape._frep.modes[0, 0] += 1
    with pytest.raises(AssertionError):
        workflow.restore_fixed_spinodoid(shape)


def test_grid_sampling_cannot_change_the_frozen_field():
    case, shape, _ = workflow.specimen("spinodoid")
    probes = np.random.default_rng(83).uniform(-.5, .5, (128, 3))
    before = case.function(*probes.T)
    shape.resolution = 6
    assert shape.grid.n_cells == 6**3
    np.testing.assert_array_equal(case.function(*probes.T), before)
    frozen = json.loads((workflow.HERE / "runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json").read_text())
    np.testing.assert_array_equal(np.column_stack((shape._frep.modes / shape._frep.cell_size,
        shape._frep.amplitudes, shape._frep.phases)), frozen["coefficients"])
