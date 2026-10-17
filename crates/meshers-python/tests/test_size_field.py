import meshers
import numpy as np
import pytest


def solid(coordinates_x, coordinates_y, coordinates_z):
    return coordinates_x * 0 - 1


def generate(**options):
    return meshers.generate(
        solid,
        bounds=(-0.5, 0.5) * 3,
        cells=4,
        periodic=(True, True, True),
        adapt=True,
        min_length=0.03,
        max_length=0.2,
        adapt_iterations=3,
        minimum_quality=0.1,
        threads=1,
        **options,
    )


def test_constant_size_preserves_default_adaptation():
    original = generate()
    explicit = generate(size_field=lambda coordinates_x, *_: np.full_like(coordinates_x, 0.2))
    np.testing.assert_array_equal(original.points, explicit.points)
    np.testing.assert_array_equal(original.tetrahedra, explicit.tetrahedra)


def test_spatial_size_refines_target_region_and_preserves_periodicity():
    mesh = generate(
        size_field=lambda coordinates_x, *_: 0.06 + 0.14 * np.sin(np.pi * coordinates_x) ** 2
    )
    local_edges = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
    edges = np.unique(np.sort(mesh.tetrahedra[:, local_edges].reshape(-1, 2), axis=1), axis=0)
    centers = mesh.points[edges].mean(axis=1)
    lengths = np.linalg.norm(mesh.points[edges[:, 0]] - mesh.points[edges[:, 1]], axis=1)
    assert np.median(lengths[abs(centers[:, 0]) < 0.15]) < 0.85 * np.median(
        lengths[abs(centers[:, 0]) > 0.35]
    )
    for axis, pairs in enumerate(mesh.periodic_pairs):
        expected = np.broadcast_to(np.eye(3)[axis], (len(pairs), 3))
        np.testing.assert_allclose(mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]], expected)
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1


@pytest.mark.parametrize("value", [0.0, -0.1, np.nan, np.inf])
def test_invalid_lengths_are_rejected(value):
    with pytest.raises(ValueError, match=r"positive|nonfinite"):
        generate(size_field=lambda coordinates_x, *_: np.full_like(coordinates_x, value))


def test_nonperiodic_size_is_rejected():
    with pytest.raises(ValueError, match="size_field is incompatible with periodic"):
        generate(size_field=lambda coordinates_x, *_: coordinates_x + 0.6)


def test_size_requires_adaptation():
    with pytest.raises(ValueError, match="requires adapt=True"):
        meshers.generate(solid, size_field=lambda coordinates_x, *_: coordinates_x * 0 + 0.1)


def test_size_callback_exception_is_preserved():
    def fail(*coordinates):
        raise RuntimeError("size callback failed")

    with pytest.raises(RuntimeError, match="size callback failed"):
        generate(size_field=fail)
