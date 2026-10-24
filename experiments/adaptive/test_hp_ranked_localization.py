"""Independent physical-field algebra and an exactly displaced planar wall."""

import numpy as np

from hp_ranked_localization import increment_bases, wall_residuals


def test_increment_bases_match_direct_complete_fields():
    rng = np.random.default_rng(54)
    h, p, lift, base = rng.normal(size=(4, 6, 13, 3))
    old, new = rng.normal(size=(2, 6, 6))
    previous, current, delta = increment_bases(old, new)
    fields = np.concatenate((h, p, lift, base))
    old_direct = lift + np.einsum("ij,ind->jnd", old, p-h)
    new_direct = base + np.einsum("ij,ind->jnd", new, p-lift)
    np.testing.assert_allclose(np.einsum("ij,ind->jnd", previous, fields), old_direct)
    np.testing.assert_allclose(np.einsum("ij,ind->jnd", current, fields), new_direct)
    np.testing.assert_allclose(np.einsum("ij,ind->jnd", delta, fields), old_direct-new_direct)


def test_surface_residual_reproduces_exact_planar_area_and_distance():
    corners = np.array([[-.2, -.2, .02], [.2, -.2, .02], [-.2, .2, .02]])
    points = np.vstack((corners, (corners[[0, 1, 2]]+corners[[1, 2, 0]])/2))
    face = np.arange(6)[None]
    for order in (3, 5, 6):
        result = wall_residuals(points, face, lambda x,y,z:z, order)
        np.testing.assert_allclose(result["area"], .08, atol=1e-14)
        np.testing.assert_allclose(result["squared_distance_integral"], .08*.02**2, atol=1e-14)
        np.testing.assert_allclose(result["signed_normal_integral"], -.08*.02, atol=1e-14)
        np.testing.assert_allclose(result["maximum_sampled_projection"], .02, atol=1e-14)
        assert not result["failed_projections"].any()


def test_projection_outside_box_is_reported_without_admitting_it():
    corners = np.array([[-.2, -.2, .02], [.2, -.2, .02], [-.2, .2, .02]])
    points = np.vstack((corners, (corners[[0, 1, 2]]+corners[[1, 2, 0]])/2))
    result = wall_residuals(points, np.arange(6)[None], lambda x,y,z:z-.6, 3)
    np.testing.assert_array_equal(result["failed_projections"], [9])
    np.testing.assert_allclose(result["maximum_box_overshoot"], .1, atol=1e-14)
