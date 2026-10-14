"""Verify exact partition and boundary preservation independently of the marker."""

import numpy as np
import pytest
from body_refine import split_interior_edges, verify_preserved_boundary
from microgen_benchmark.geometry import boundary, canonical_mesh


def cube():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
                       [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], float) - .5
    tets = np.array([[0, 1, 3, 7], [0, 3, 2, 7], [0, 2, 6, 7],
                     [0, 6, 4, 7], [0, 4, 5, 7], [0, 5, 1, 7]])
    return canonical_mesh(points, tets)[:2]


def test_cube_body_diagonal_split_preserves_boundary_volume_and_affine_field():
    points, tets = cube()
    original_boundary = boundary(points, tets)
    edges = np.unique(np.sort(tets[:, [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]].reshape(-1, 2), axis=1), axis=0)
    refined_points, children, info = split_interior_edges(points, tets, original_boundary, np.ones(len(edges)), 1)
    assert info["inserted_body_nodes"] == 1
    assert len(children) == 2 * len(tets)
    np.testing.assert_array_equal(refined_points[:len(points)], points)
    assert {tuple(sorted(t)) for t in boundary(refined_points, children)} == {tuple(sorted(t)) for t in original_boundary}
    refined_boundary = boundary(refined_points, children)
    verify_preserved_boundary(points, original_boundary, refined_points, refined_boundary)
    # Triangle order and cyclic permutations preserve the same oriented surface.
    verify_preserved_boundary(points, original_boundary, refined_points, np.roll(refined_boundary[::-1], 1, axis=1))
    moved = refined_points.copy()
    moved[0, 0] += .001
    with pytest.raises(ValueError, match="original vertex"):
        verify_preserved_boundary(points, original_boundary, moved, refined_boundary)
    with pytest.raises(ValueError, match="boundary triangle"):
        verify_preserved_boundary(points, original_boundary, refined_points, refined_boundary[:, ::-1])
    volumes = np.linalg.det(refined_points[children[:, 1:]] - refined_points[children[:, :1]]) / 6
    assert volumes.min() > 0
    np.testing.assert_allclose(volumes.sum(), 1, atol=1e-14)
    # Midpoint interpolation retains every affine displacement exactly.
    matrix = np.array([[1, .2, .3], [.4, 2, .5], [.6, .7, 3]])
    displacement = points @ matrix.T
    np.testing.assert_allclose((displacement[0] + displacement[7]) / 2,
                               refined_points[-1] @ matrix.T, atol=1e-14)


def test_zero_indicator_does_not_refine():
    points = np.array([[0, 0, 0], [.3, 0, 0], [0, .3, 0], [0, 0, .3]], float)
    tets = np.array([[0, 1, 2, 3]])
    p, t, info = split_interior_edges(points, tets, boundary(points, tets), np.zeros(6), 2)
    assert info["inserted_body_nodes"] == 0
    np.testing.assert_array_equal(p, points)
    np.testing.assert_array_equal(t, tets)


@pytest.mark.parametrize("marker", ["length", "volume"])
def test_geometry_control_preserves_boundary_without_reading_a_guide(tmp_path, monkeypatch, marker):
    import json
    from types import SimpleNamespace
    import body_refine
    from homogenization import MATERIAL
    from run import digest

    points, tets = cube()
    surface = boundary(points, tets)
    seed = tmp_path / "body-seed"
    seed.mkdir()
    np.savez_compressed(seed / "mesh.npz", points=points, tetrahedra=tets, surface=surface)
    row = {"case": "solid-cell", "cells": 1, "mesh_sha256": digest(seed / "mesh.npz"),
           "generation_s": 0., "solid_fraction": 1.,
           "homogenization": {"stiffness_engineering_voigt": MATERIAL.tolist()}}
    (seed / "result.json").write_text(json.dumps(row))
    (seed / "input-provenance.json").write_text(json.dumps({"specimen": {"function": "solid", "offset": 1}}))
    monkeypatch.setattr(body_refine, "__file__", str(tmp_path / "body_refine.py"))

    def forbidden_guide(_):
        raise AssertionError("Geometry-only refinement must not use the FEA guide")

    monkeypatch.setattr(body_refine, "guide_score", forbidden_guide)
    case = SimpleNamespace(id="solid-cell", function="solid", offset=1)
    p, t, info = body_refine.generate(case, 1, {"body_marker": marker, "body_refine_nodes": 1,
                                             "minimum_quality": .1, "max_tetrahedra": 100},
                                     tmp_path, lambda **_: None)
    assert info["inserted_body_nodes"] == 1
    assert info["body_marker_kind"] == marker
    assert info["archived_guide_fea_s"] == info["archived_guide_generation_s"] == 0
    verify_preserved_boundary(points, surface, p, boundary(p, t))


def test_volume_marker_targets_larger_cells_instead_of_the_longer_thin_block():
    points, tets = cube()
    thin_center, large_center = np.array([-.22, 0, 0]), np.array([.29, 0, 0])
    points = np.vstack([points * [.5, .2, .2] + thin_center,
                        points * .3 + large_center])
    tets = np.vstack([tets, tets + 8])
    surface = boundary(points, tets)
    edges = np.unique(np.sort(tets[:, [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]].reshape(-1, 2), axis=1), axis=0)
    length_points, _, _ = split_interior_edges(points, tets, surface, np.ones(len(edges)), 1, marker="length")
    np.testing.assert_allclose(length_points[-1], thin_center, atol=1e-14)
    refined, children, info = split_interior_edges(points, tets, surface, np.ones(len(edges)), 1, marker="volume")
    assert info["inserted_body_nodes"] == 1
    np.testing.assert_allclose(refined[-1], large_center, atol=1e-14)
    verify_preserved_boundary(points, surface, refined, boundary(refined, children))
    volumes = np.linalg.det(refined[children[:, 1:]] - refined[children[:, :1]]) / 6
    assert volumes.min() > 0
    np.testing.assert_allclose(volumes.sum(), .5 * .2 * .2 + .3**3, atol=1e-14)
    repeated, repeated_tets, _ = split_interior_edges(points, tets, surface, np.ones(len(edges)), 1, marker="volume")
    np.testing.assert_array_equal(repeated, refined)
    np.testing.assert_array_equal(repeated_tets, children)
