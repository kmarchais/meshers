"""Public constrained meshing: geometry, quality, mapping, exports and failures."""

import json

import h5py
import meshers
import numpy as np
import pytest
import pyvista as pv


def constraints():
    return {"right": lambda x, y, z: x - 0.501, "back": lambda x, y, z: y - 0.751}


def volume(mesh):
    p = mesh.points[mesh.tetrahedra]
    v = np.linalg.det(p[:, 1:] - p[:, :1]) / 6
    assert np.all(v > 0)
    return v.sum()


def test_named_intersection_quality_and_export(tmp_path):
    args = {"cells": 8, "snap": 0, "periodic": (False, False, True)}
    original = meshers.generate_intersection(constraints(), optimize_passes=0, **args)
    mesh = meshers.generate_intersection(
        constraints(), optimize_passes=4, geometry_tolerance=1e-8, **args
    )
    assert (
        mesh.diagnostics["minimum_mmg_quality"] > original.diagnostics["minimum_mmg_quality"] * 10
    )
    assert abs(volume(mesh) - 0.501 * 0.751) < 1e-12
    np.testing.assert_allclose(mesh.points[mesh.surface_for("right"), 0], 0.501, atol=1e-12)
    assert len(mesh.feature_edges) > 0
    assert mesh.constraint_names == ("right", "back")
    assert len(mesh.periodic_pairs[2]) > 0
    path = tmp_path / "intersection.vtkhdf"
    mesh.write_vtkhdf(path)
    with h5py.File(path) as f:
        root = f["VTKHDF"]
        np.testing.assert_array_equal(root["PointData/ParameterCoordinates"], mesh.parameters)
        np.testing.assert_array_equal(root["PointData/SurfaceConstraints"], mesh.constraint_masks)
        np.testing.assert_array_equal(root["FieldData/FeatureEdges"], mesh.feature_edges)
        np.testing.assert_array_equal(root["FieldData/PeriodicPairs2"], mesh.periodic_pairs[2])
        np.testing.assert_array_equal(
            root["FieldData/PeriodicTransforms"], mesh.periodic_transforms.reshape(3, 16)
        )
        assert json.loads(bytes(root["FieldData/ConstraintNamesUTF8"][:])) == list(
            mesh.constraint_names
        )

    restored = pv.read(path)
    np.testing.assert_array_equal(restored.field_data["FeatureEdges"], mesh.feature_edges)
    assert json.loads(restored.field_data["ConstraintNamesUTF8"].tobytes()) == list(
        mesh.constraint_names
    )


def test_quality_gate_can_use_one_extra_optimization_pass():
    mesh = meshers.generate_intersection(
        constraints(),
        cells=8,
        snap=0,
        periodic=(False, False, True),
        optimize_passes=2,
        minimum_quality=0.1,
    )
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    assert mesh.diagnostics["quality_optimization_passes"] == 3


def test_rotational_periodicity():
    radius = 3.0
    angle = 1 / radius
    transform = np.eye(4)
    transform[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]

    def mapping(x, y, z):
        return ((radius + x) * np.cos(z / radius), (radius + x) * np.sin(z / radius), -y)

    mesh = meshers.generate_intersection(
        constraints(),
        cells=8,
        snap=0.1,
        coordinate_map=mapping,
        periodic=(False, False, True),
        periodic_transforms={2: transform},
        geometry_tolerance=0.02,
        optimize_passes=2,
    )
    assert volume(mesh) > 0
    a, b = mesh.periodic_pairs[2].T
    np.testing.assert_allclose(mesh.points[a] @ transform[:3, :3].T, mesh.points[b], atol=1e-12)
    assert mesh.diagnostics["sampled_surface_error"] < 0.02
    with pytest.raises(meshers.MeshingError, match="physical periodic"):
        meshers.generate_intersection(
            constraints(),
            cells=4,
            coordinate_map=mapping,
            periodic=(False, False, True),
            periodic_transforms={2: np.eye(4)},
            optimize_passes=0,
        )
    with pytest.raises(ValueError, match="explicit"):
        meshers.generate_intersection(
            constraints(), coordinate_map=mapping, periodic=(False, False, True)
        )


def test_callback_and_compiled_agree():
    args = {"cells": 4, "optimize_passes": 0, "snap": 0}
    a = meshers.generate_intersection(constraints(), **args)
    b = meshers.generate_intersection(constraints(), compile=False, **args)
    np.testing.assert_allclose(a.points, b.points, atol=1e-12)
    np.testing.assert_array_equal(a.tetrahedra, b.tetrahedra)
    assert a.diagnostics["callback_calls"] == 0
    assert b.diagnostics["callback_calls"] > 0


def test_failures_are_explicit():
    with pytest.raises(meshers.MeshingError, match="minimum MMG quality"):
        meshers.generate_intersection(constraints(), cells=4, minimum_quality=1)
    with pytest.raises(meshers.MeshingError, match="inverted"):
        meshers.generate_intersection(
            constraints(), cells=4, coordinate_map=lambda x, y, z: (-x, y, z)
        )
    with pytest.raises(meshers.MeshingError, match="periodic"):
        meshers.generate_intersection(constraints(), cells=4, periodic=(True, False, False))
    with pytest.raises(meshers.MeshingError, match="budget"):
        meshers.generate_intersection(constraints(), cells=4, max_tetrahedra=1)
    token = meshers.CancellationToken()
    token.cancel()
    with pytest.raises(meshers.CancelledError):
        meshers.generate_intersection(constraints(), cells=4, cancel=token)

    def fails(x, y, z):
        raise RuntimeError("user callback failed")

    with pytest.raises(RuntimeError, match="user callback failed"):
        meshers.generate_intersection(
            {"bad": fails, "back": lambda x, y, z: y - 0.5}, cells=4, compile=False
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"cells": 3},
        {"optimize_passes": 21},
        {"minimum_quality": np.nan},
        {"geometry_tolerance": -1},
        {"snap": 0.3},
        {"bounds": [0, 1]},
        {"periodic": (1, 0, 0)},
        {"periodic_transforms": {2: np.eye(4)}},
    ],
)
def test_invalid_options(kwargs):
    with pytest.raises(ValueError):
        meshers.generate_intersection(constraints(), **kwargs)


def test_empty_and_geometry_rejection():
    mesh = meshers.generate_intersection(
        {"outside": lambda x, y, z: x + 2, "back": lambda x, y, z: y - 0.5}, cells=4
    )
    assert mesh.points.shape == (0, 3)
    assert mesh.tetrahedra.shape == (0, 4)
    assert mesh.feature_edges.shape == (0, 2)
    with pytest.raises(meshers.MeshingError, match="sampled surface error"):
        meshers.generate_intersection(
            {
                "sphere": lambda x, y, z: (x - 0.5) ** 2 + (y - 0.5) ** 2 + (z - 0.5) ** 2 - 0.3**2,
                "cut": lambda x, y, z: z - 0.55,
            },
            cells=8,
            geometry_tolerance=1e-8,
            optimize_passes=0,
        )


def test_snapping_rollback_keeps_periodic_tpms_geometry():
    def gyroid(x, y, z):
        a, b, c = (
            2 * np.pi * (x / 0.75 + 0.125),
            2 * np.pi * (y / 0.75 + 0.125),
            2 * np.pi * (z + 0.125),
        )
        return np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)

    mesh = meshers.generate_intersection(
        {
            "upper": lambda x, y, z: gyroid(x, y, z) - 0.55,
            "lower": lambda x, y, z: -gyroid(x, y, z) - 0.55,
            "wall": lambda x, y, z: np.sqrt(x * x + y * y) - 0.75,
        },
        bounds=(-0.8, 0.8, -0.8, 0.8, 0, 3),
        cells=(24, 24, 40),
        periodic=(False, False, True),
        snap=0.2,
        optimize_passes=0,
    )
    assert mesh.diagnostics["reverted_snap_vertices"] > 0
    assert volume(mesh) > 0
    a, b = mesh.periodic_pairs[2].T
    np.testing.assert_allclose(mesh.points[a] + [0, 0, 3], mesh.points[b], atol=1e-12)
    vertices, degree = np.unique(mesh.feature_edges.ravel(), return_counts=True)
    end = np.isin(mesh.parameters[vertices, 2], [0, 3])
    assert np.all((degree == 2) | ((degree == 1) & end))


def test_cancellation_during_python_surface_verification():
    token = meshers.CancellationToken()

    def field(x, y, z):
        # Native construction uses scalar callbacks; verification samples batches.
        if x.size > 1:
            token.cancel()
        return x - 0.6

    with pytest.raises(meshers.CancelledError):
        meshers.generate_intersection(
            {"x": field, "y": lambda x, y, z: y - 0.7},
            cells=4,
            compile=False,
            optimize_passes=0,
            geometry_tolerance=0.1,
            cancel=token,
        )
    assert token.cancelled


@pytest.mark.parametrize("scale", [1e-14, 1e-9, 1.0, 1e9, 1e14])
@pytest.mark.parametrize("passes", [0, 2])
def test_constraint_projection_preserves_geometry_under_rescaling(scale, passes):
    mesh = meshers.generate_intersection(
        {"x": lambda x, y, z: scale * (x - 0.51), "y": lambda x, y, z: (y - 0.74) / scale},
        cells=4,
        snap=0.2,
        optimize_passes=passes,
    )
    assert mesh.constraint_masks is not None
    for bit, axis, position in [(1, 0, 0.51), (2, 1, 0.74)]:
        boundary = mesh.points[(mesh.constraint_masks & bit) != 0, axis]
        assert len(boundary)
        np.testing.assert_allclose(boundary, position, rtol=0, atol=1e-11)
        assert np.all(mesh.points[:, axis] <= position + 1e-11)
    assert mesh.diagnostics["volume"] == pytest.approx(0.51 * 0.74, abs=1e-11)


@pytest.mark.parametrize("batch_size", [1, 7, 4096])
@pytest.mark.parametrize("mapped", [False, True])
def test_surface_verification_respects_callback_batch_size(batch_size, mapped):
    sizes = []

    def record(x):
        sizes.append(x.size)
        assert x.size <= batch_size

    def field(x, y, z):
        record(x)
        return x - 0.6

    def mapping(x, y, z):
        record(x)
        return x + 0.1 * y * y, y, z

    mesh = meshers.generate_intersection(
        {"x": field, "y": lambda x, y, z: y - 0.7},
        cells=4,
        compile=False,
        optimize_passes=0,
        coordinate_map=mapping if mapped else None,
        batch_size=batch_size,
        geometry_tolerance=0.1,
    )
    assert len(mesh.tetrahedra)
    assert mesh.diagnostics["sampled_surface_error"] < 0.1
    assert sizes
    assert max(sizes) <= batch_size


@pytest.mark.parametrize("shift", [0.0, 1e6, -1e6])
def test_rotational_periodicity_tolerates_coordinate_roundoff(shift):
    radius = 3.0
    angle = 1 / radius
    transform = np.eye(4)
    transform[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    offset = np.full(3, shift)
    transform[:3, 3] = offset - transform[:3, :3] @ offset

    def mapping(x, y, z):
        return (
            (radius + x) * np.cos(z / radius) + shift,
            (radius + x) * np.sin(z / radius) + shift,
            -y + shift,
        )

    args = {
        "cells": 4,
        "coordinate_map": mapping,
        "periodic": (False, False, True),
        "geometry_tolerance": 0.1,
        "optimize_passes": 0,
    }
    mesh = meshers.generate_intersection(constraints(), periodic_transforms={2: transform}, **args)
    assert len(mesh.periodic_pairs[2])
    assert volume(mesh) > 0
    a, b = mesh.periodic_pairs[2].T
    mapped = mesh.points[a] @ transform[:3, :3].T + transform[:3, 3]
    np.testing.assert_allclose(mapped, mesh.points[b], rtol=0, atol=1e-9)
    incorrect = transform.copy()
    incorrect[0, 3] += 1e-6
    with pytest.raises(meshers.MeshingError, match="physical periodic"):
        meshers.generate_intersection(constraints(), periodic_transforms={2: incorrect}, **args)


@pytest.mark.parametrize("scale", [1e-9, 1.0, 1e9])
def test_mapped_surface_tolerance_tracks_physical_scale(scale):
    radius = 3.0

    def mapping(x, y, z):
        return (
            scale * (radius + x) * np.cos(z / radius),
            scale * (radius + x) * np.sin(z / radius),
            -scale * y,
        )

    args = {"cells": 4, "coordinate_map": mapping, "optimize_passes": 0}
    mesh = meshers.generate_intersection(constraints(), geometry_tolerance=0.01 * scale, **args)
    centers = mesh.points[mesh.surface_for("right")].mean(axis=1)
    radial_error = np.max(np.abs(np.linalg.norm(centers[:, :2], axis=1) - scale * (radius + 0.501)))
    assert radial_error > 1e-4 * scale
    assert mesh.diagnostics["sampled_surface_error"] >= radial_error * 0.99
    with pytest.raises(meshers.MeshingError, match="sampled surface error"):
        meshers.generate_intersection(constraints(), geometry_tolerance=1e-6 * scale, **args)
