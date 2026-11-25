# Copyright (c) 2026 Kevin Marchais

"""Damage actual cap constraints while keeping valid corner incidence and translations."""

from dataclasses import replace

import meshers
import numpy as np
import pytest

from tools.performance_worker import check_periodic_axis


def dropped_pair(mesh):
    node = mesh.surface[mesh.boundary_tags == 1][0, -1]
    pairs = mesh.periodic_pairs[0]
    return replace(mesh, periodic_pairs=(pairs[pairs[:, 0] != node], *mesh.periodic_pairs[1:]))


def swapped_midpoints(mesh):
    surface = mesh.surface.copy()
    face = np.flatnonzero(mesh.boundary_tags == 2)[0]
    surface[face, [3, 4]] = surface[face, [4, 3]]
    return replace(mesh, surface=surface)


def flipped_diagonal(mesh):
    # Find two x+ triangles sharing an edge and a tet apex. Replacing the
    # diagonal in this pyramid changes only the triangulation of its base.
    cells = mesh.tetrahedra.copy()
    surface = mesh.surface.copy()
    faces = {tuple(sorted(surface[i])): i for i in np.flatnonzero(mesh.boundary_tags == 2)}
    on_plane = mesh.points[cells, 0] == 1
    edges = {}
    for cell in np.flatnonzero(on_plane.sum(axis=1) == 3):
        apex = cells[cell, ~on_plane[cell]][0]
        corners = cells[cell, on_plane[cell]]
        face = faces[tuple(sorted(corners))]
        for a, b, v in (np.roll(corners, k) for k in range(3)):
            edge = (apex, min(a, b), max(a, b))
            if edge in edges:
                other_cell, other_face, w = edges[edge]
                triangles = np.array([[v, w, a], [w, v, b]])
                # Both new faces must have nonzero area and consistent winding.
                uv = mesh.points[triangles, 1:]
                areas = np.linalg.det(uv[:, 1:] - uv[:, :1])
                if np.prod(areas) <= 0:
                    continue
                if areas[0] < 0:
                    triangles = triangles[:, [0, 2, 1]]
                surface[[face, other_face]] = triangles
                cells[[cell, other_cell]] = np.column_stack(([apex, apex], triangles))
                return replace(mesh, surface=surface, tetrahedra=cells)
            edges[edge] = (cell, face, v)
    pytest.fail("No strictly convex flippable cap pyramid in this fixture")


def unit_shift(axis):
    transform = np.eye(4)
    transform[axis, 3] = 1
    return transform


@pytest.fixture(scope="module", params=[1, 2])
def gyroid(request):
    return meshers.generate(
        "gyroid",
        cells=8,
        band=(-0.5, 0.5),
        periodic=(True,) * 3,
        geometry_tolerance=1,
        optimize_passes=1,
        threads=1,
        order=request.param,
    )


def test_accepts_reordered_cap_then_rejects_missing_pair(gyroid):
    # Rotate the local face numbering including the associated edge slots.
    perm = [1, 2, 0] if gyroid.surface.shape[1] == 3 else [1, 2, 0, 4, 5, 3]
    shuffled = replace(
        gyroid,
        surface=gyroid.surface[::-1, perm],
        boundary_tags=gyroid.boundary_tags[::-1],
        periodic_pairs=tuple(pairs[::-1] for pairs in gyroid.periodic_pairs),
    )
    check_periodic_axis(shuffled, 0, unit_shift(0))
    with pytest.raises(RuntimeError, match="Missing periodic cap pair"):
        check_periodic_axis(dropped_pair(shuffled), 0, unit_shift(0))


def test_rejects_damaged_cap(gyroid):
    if gyroid.surface.shape[1] == 3:
        damaged = flipped_diagonal(gyroid)
        reason = "Opposite cap traces differ"
    else:
        damaged = swapped_midpoints(gyroid)
        reason = "Tagged cap differs from volume trace"
    with pytest.raises(RuntimeError, match=reason):
        check_periodic_axis(damaged, 0, unit_shift(0))


def test_duplicate_pair_still_translates_but_repeats_a_constraint(gyroid):
    pairs = gyroid.periodic_pairs[0]
    repeated = replace(
        gyroid,
        periodic_pairs=(np.vstack((pairs, pairs[:1])), *gyroid.periodic_pairs[1:]),
    )
    with pytest.raises(RuntimeError, match="Invalid periodic correspondence"):
        check_periodic_axis(repeated, 0, unit_shift(0))


def rotational_fixture():
    angle = 1 / 3
    transform = np.eye(4)
    transform[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]

    def mapping(x, y, z):
        return ((3 + x) * np.cos(z / 3), (3 + x) * np.sin(z / 3), -y)

    mesh = meshers.generate_intersection(
        {"right": lambda x, _y, _z: x - 0.501, "back": lambda _x, y, _z: y - 0.751},
        cells=8,
        snap=0.1,
        coordinate_map=mapping,
        periodic=(False, False, True),
        periodic_transforms={2: transform},
        geometry_tolerance=0.02,
        optimize_passes=2,
    )
    return mesh, transform


def test_rotational_cap_requires_physical_map():
    mesh, transform = rotational_fixture()
    np.testing.assert_array_equal(mesh.periodic_transforms[2], transform)
    check_periodic_axis(mesh, 2, transform)
    with pytest.raises(RuntimeError, match="Invalid periodic correspondence"):
        check_periodic_axis(mesh, 2, np.eye(4))
