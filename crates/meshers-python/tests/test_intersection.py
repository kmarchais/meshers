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


def test_octahedron_preserves_four_plane_corners_and_exact_volume():
    from itertools import product

    planes = {}
    signs = tuple(product((1, -1), repeat=3))
    for index, normal in enumerate(signs):

        def plane(x, y, z, normal=normal):
            return (normal[0] * x + normal[1] * y + normal[2] * z - 0.35) / np.sqrt(3)

        planes[f"plane{index}"] = plane
    mesh = meshers.generate_intersection(
        planes,
        bounds=(-0.5, 0.5) * 3,
        cells=24,
        geometry_tolerance=0.002,
        minimum_quality=0.1,
        max_tetrahedra=250000,
        optimize_passes=1,
    )
    assert abs(volume(mesh) - 4 * 0.35**3 / 3) < 1e-12
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    assert len(mesh.feature_edges) > 0
    for corner in 0.35 * np.concatenate((np.eye(3), -np.eye(3))):
        index = np.argmin(np.linalg.norm(mesh.points - corner, axis=1))
        assert np.linalg.norm(mesh.points[index] - corner) < 1e-12
        assert int(mesh.constraint_masks[index]).bit_count() == 4
    for index, normal in enumerate(signs):
        faces = mesh.surface_for(f"plane{index}")
        assert len(faces) > 0
        assert np.max(np.abs(mesh.points[faces] @ normal - 0.35)) < 1e-12


@pytest.mark.skipif(
    getattr(meshers._meshers, "MAX_CONSTRAINTS", 8) < 32, reason="optional wide-constraints feature"
)
def test_high_constraint_bits_and_planar_intersections():
    fields = {f"inactive{i}": lambda x, y, z: -1 + 0 * x for i in range(29)}
    fields.update({
        "x": lambda x, y, z: x - 0.63,
        "y": lambda x, y, z: y - 0.72,
        "z": lambda x, y, z: z - 0.81,
    })
    mesh = meshers.generate_intersection(fields, cells=8, geometry_tolerance=1e-10)
    assert mesh.constraint_masks.dtype == np.uint32
    assert np.isclose(volume(mesh), 0.63 * 0.72 * 0.81, atol=1e-12)
    for axis, (name, level) in enumerate(zip("xyz", (0.63, 0.72, 0.81), strict=True)):
        faces = mesh.surface_for(name)
        assert len(faces) > 0
        assert np.max(abs(mesh.points[faces, axis] - level)) < 1e-12


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


def _cylinder_sdf(x, y, z):
    radial = np.sqrt(y * y + z * z) - 0.23
    axial = np.abs(x) - 0.3
    return np.sqrt(np.maximum(radial, 0) ** 2 + np.maximum(axial, 0) ** 2) + np.minimum(
        np.maximum(radial, axial), 0
    )


def _box_sdf(x, y, z, scale=1.0):
    a, b, c = (scale * (np.abs(v) - 0.25) for v in (x, y, z))
    return np.sqrt(
        np.maximum(a, 0) ** 2 + np.maximum(b, 0) ** 2 + np.maximum(c, 0) ** 2
    ) + np.minimum(np.maximum(a, np.maximum(b, c)), 0)


def test_sdf_values_and_zero_sublevel_metadata_are_preserved():
    field = meshers.compile_field(_cylinder_sdf)
    points = np.random.default_rng(20261112).uniform(-0.5, 0.5, (4096, 3))
    before = field.values(points)
    branches = field._intersection_fields()
    assert branches is not None
    assert len(branches) == 2
    np.testing.assert_array_equal(before, field.values(points))
    np.testing.assert_allclose(before, _cylinder_sdf(*points.T), rtol=1e-14, atol=1e-14)
    np.testing.assert_array_equal(
        before <= 0, np.maximum.reduce([f.values(points) for f in branches]) <= 0
    )
    assert field._halfspaces is None
    assert field._polyhedral_regions is None


@pytest.mark.parametrize(
    "mode",
    [
        "extra-positive",
        "extra-negative",
        "missing",
        "unclamped",
        "other-root",
        "inside-minimum",
        "whole-offset",
    ],
)
def test_altered_sdf_graphs_keep_the_scalar_path(mode):
    def field(x, y, z):
        a, b = np.sqrt(y * y + z * z) - 0.23, np.abs(x) - 0.3
        square = np.maximum(a, 0) ** 2
        if mode != "missing":
            square += b * b if mode == "unclamped" else np.maximum(b, 0) ** 2
        if mode in {"extra-positive", "extra-negative"}:
            square += 0.0001 if mode == "extra-positive" else -0.0001
        if mode == "other-root":
            square += np.maximum(z, 0) ** 2
        inside = np.minimum(a, b) if mode == "inside-minimum" else np.maximum(a, b)
        value = np.sqrt(square) + np.minimum(inside, 0)
        return value + 0.01 if mode == "whole-offset" else value

    assert meshers.compile_field(field)._intersection_fields() is None


@pytest.mark.parametrize("gradient", [False, lambda x, y, z: (0 * x, 0 * y, 0 * z)])
def test_sdf_explicit_gradient_policy_does_not_retain_branches(gradient):
    assert meshers.compile_field(_cylinder_sdf, gradient=gradient)._intersection_fields() is None


@pytest.mark.parametrize("angle", [0.0, 0.37])
def test_adaptive_sdf_cylinder_preserves_the_rim_and_pose(angle):
    cosine, sine = np.cos(angle), np.sin(angle)
    center = np.array([0.03, -0.02, 0.01]) if angle else np.zeros(3)

    def field(x, y, z):
        x, y, z = x - center[0], y - center[1], z - center[2]
        return _cylinder_sdf(cosine * x + sine * y, -sine * x + cosine * y, z)

    mesh = meshers.generate(
        field,
        bounds=(-0.5, 0.5) * 3,
        cells=24,
        adapt=True,
        min_length=0.002,
        max_length=0.1,
        gradation=0.6,
        adapt_tolerance=0.0018,
        adapt_iterations=4,
        geometry_tolerance=0.002,
        minimum_quality=0.1,
        optimize_passes=1,
        max_tetrahedra=250000,
    )
    assert mesh.diagnostics["retained_constraints"] == 3
    assert mesh.diagnostics["minimum_mmg_quality"] >= 0.1
    corners = mesh.points[mesh.tetrahedra]
    determinants = np.linalg.det(corners[:, 1:] - corners[:, :1])
    assert determinants.min() > 0
    volume = determinants.sum() / 6
    assert abs(volume / (np.pi * 0.23**2 * 0.6) - 1) < 0.01
    theta = np.arange(512) * 2 * np.pi / 512
    local = np.concatenate([
        np.column_stack([np.full(512, end), 0.23 * np.cos(theta), 0.23 * np.sin(theta)])
        for end in (-0.3, 0.3)
    ])
    rim = (
        np.column_stack([
            cosine * local[:, 0] - sine * local[:, 1],
            sine * local[:, 0] + cosine * local[:, 1],
            local[:, 2],
        ])
        + center
    )
    surface = pv.PolyData(
        mesh.points, np.column_stack([np.full(len(mesh.surface), 3), mesh.surface])
    )
    _, closest = surface.find_closest_cell(rim, return_closest_point=True)
    assert np.linalg.norm(rim - closest, axis=1).max() <= 0.002


@pytest.mark.parametrize("adaptive", [False, True])
def test_sdf_original_parent_overflow_is_not_bypassed(adaptive):
    with pytest.raises((meshers.MeshingError, ValueError), match="nonfinite"):
        meshers.generate(
            lambda x, y, z: _box_sdf(x, y, z, 1e200),
            bounds=(-0.5, 0.5) * 3,
            cells=8,
            adapt=adaptive,
            geometry_tolerance=0.1,
            minimum_quality=0.1,
            optimize_passes=1,
        )


def test_sdf_underflow_does_not_enable_a_different_zero_set():
    from meshers.intersection import _retained_constraints

    field = meshers.compile_field(lambda x, y, z: _box_sdf(x, y, z, 1e-200))
    assert _retained_constraints(field, np.array((-0.5, 0.5) * 3)) is None


@pytest.mark.parametrize("band", [None, (-0.05, 0.05)])
def test_nonadaptive_sdf_and_bands_keep_the_original_scalar_path(band):
    mesh = meshers.generate(
        _box_sdf,
        bounds=(-0.5, 0.5) * 3,
        cells=8,
        band=band,
        geometry_tolerance=0.1,
        minimum_quality=0.1,
        optimize_passes=1,
    )
    assert mesh.diagnostics.get("retained_constraints", 0) == 0
    assert len(mesh.tetrahedra) > 0


@pytest.mark.parametrize(
    ("nodes", "output"),
    [([("input", [3], 0.0)], [0]), ([("constant", [], 1.0), ("add", [0, 2], 0.0)], [1])],
)
def test_sdf_recognition_still_validates_raw_graph_indices(nodes, output):
    with pytest.raises(ValueError, match="expression"):
        meshers._meshers._compile_expression(nodes, output, retain_sdf=True)


def test_localized_underflow_cannot_enable_a_different_zero_set():
    def field(x, y, z):
        scale = np.where(np.abs(x - 0.4) < 0.03, 1e-200, 1.0)
        return _box_sdf(x, y, z, scale)

    compiled = meshers.compile_field(field)
    assert compiled.values(np.array([[0.4, 0.0, 0.0]]))[0] == 0
    assert compiled._intersection_fields() is None


@pytest.mark.parametrize("kind", ["uniform-scale", "sheared-radius", "repeated-radius-basis"])
def test_nonunit_sdf_branches_keep_the_scalar_path(kind):
    def field(x, y, z):
        if kind == "uniform-scale":
            return _box_sdf(x, y, z, 2.0)
        square = y * y + ((y + z) ** 2 if kind == "sheared-radius" else y * y + z * z)
        radial, axial = np.sqrt(square) - 0.23, np.abs(x) - 0.3
        return np.sqrt(np.maximum(radial, 0) ** 2 + np.maximum(axial, 0) ** 2) + np.minimum(
            np.maximum(radial, axial), 0
        )

    assert meshers.compile_field(field)._intersection_fields() is None
