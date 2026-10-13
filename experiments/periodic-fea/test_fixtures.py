"""Shared-field checks before using the expanded fixture catalogue."""

import numpy as np
import pytest
import fixtures
from microgen.shape import surface_functions


@pytest.mark.parametrize("name", fixtures.TPMS)
def test_tpms_formula_and_all_periodic_axes(name):
    case = fixtures.make_case(name, {"function": name, "offset": 1.3, "microgen_callable": True})
    points = np.random.default_rng(42).uniform(-.5, .5, size=(25, 3))
    expected = getattr(surface_functions, name)(*(2 * np.pi * points.T)) / .65
    assert np.allclose(case.field(*points.T), expected, rtol=1e-12, atol=1e-12)
    for axis in range(3):
        assert np.allclose(case.field(*(points + np.eye(3)[axis]).T), expected, rtol=1e-11, atol=1e-11)


def test_spectral_solid_sign_periodicity_and_fixed_spectrum(monkeypatch):
    spectrum = {"modes": np.array([[1, 0, 0], [1, 2, -1]]), "amplitudes": np.array([1., .2]),
                "phases": np.array([.3, -.7]), "threshold": np.array(.1)}
    monkeypatch.setattr(fixtures, "spectral_data", lambda: spectrum)
    case = fixtures.make_case("spinodoid", {"function": "spinodoid_spectral", "offset": 1,
                                           "microgen_callable": True})
    points = np.random.default_rng(42).uniform(-.5, .5, size=(25, 3))
    raw = np.cos(2 * np.pi * ((points + .5) @ spectrum["modes"].T) + spectrum["phases"]) @ spectrum["amplitudes"] - spectrum["threshold"]
    assert np.allclose(case.field(*points.T), np.tanh(-raw))
    assert np.array_equal(case.material(*points.T), raw >= 0)
    assert case.band == (-1, 0)
    for axis in range(3):
        assert np.allclose(case.field(*(points + np.eye(3)[axis]).T), case.field(*points.T), atol=1e-12)
