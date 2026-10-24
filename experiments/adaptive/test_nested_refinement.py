"""Checks that refinement retains curved geometry and the full old FE space."""

import numpy as np
import pytest
import json

import microgen_workflow as workflow
from nested_refinement import EDGES, FACES, child_maps, interpolation_check, positive_quality, refine


def test_fixture_provenance_rejects_changed_samples(tmp_path):
    from nested_preflight import fixture_inputs

    paths = [tmp_path / ("smooth-lattice" + suffix) for suffix in ("-samples.npz", "-specimen.json")]
    for path in paths:
        path.write_bytes(b"fixed fixture")
    (tmp_path / "manifest.json").write_text(json.dumps({"inputs": {
        "/home/kmarchais/dev/meshers/" + path.name: workflow.digest(path) for path in paths}}))
    assert fixture_inputs(tmp_path) == paths
    paths[0].write_bytes(b"different fixture")
    with pytest.raises(AssertionError, match="provenance mismatch"):
        fixture_inputs(tmp_path)


@pytest.mark.parametrize("offset", range(6))
def test_bisection_partitions_parent_reference_tetrahedron(offset):
    points = np.vstack((np.zeros(3), np.eye(3)))
    mapped = child_maps(offset) @ points
    determinants = np.linalg.det(mapped[:, 1:] - mapped[:, :1])
    np.testing.assert_allclose(determinants, .5, atol=1e-15)


def test_preferred_short_edge_is_bisected_with_exact_parent_map():
    workflow.imports()
    from fedoo.lib_elements.tetrahedron import Tet10

    points = Tet10().xi_nd.copy()
    cells = np.arange(10)[None]
    edges = {tuple(sorted(pair)): 4 + offset for offset, pair in enumerate(EDGES)}
    surface = np.array([list(face) + [edges[tuple(sorted((face[first], face[second])))]
        for first, second in ((0, 1), (1, 2), (2, 0))] for face in FACES])
    fields = np.random.default_rng(19).normal(size=(6, 10, 3))
    result = refine(points, cells, surface, fields, periodic=False, preferred_edges=[(0, 3)])
    np.testing.assert_array_equal(result["parent_barycentric_maps"], child_maps(3))
    assert result["diagnostics"]["accepted_preferred_edge_orbits"] == 1
    assert interpolation_check(points, cells, fields, result)["maximum_field_error"] < 1e-12
    with pytest.raises(AssertionError, match="must exist"):
        refine(points, cells, surface, fields, periodic=False, preferred_edges=[(0, 99)])


def test_curved_geometry_and_arbitrary_quadratic_fields_are_preserved():
    workflow.imports()
    from fedoo.lib_elements.tetrahedron import Tet10

    element = Tet10()
    points = element.xi_nd.copy()
    points[4:] += np.random.default_rng(73).normal(size=(6, 3)) * .015
    cells = np.arange(10)[None]
    edges = {tuple(sorted(pair)): 4 + offset for offset, pair in enumerate(EDGES)}
    surface = np.array([list(face) + [edges[tuple(sorted((face[first], face[second])))]
        for first, second in ((0, 1), (1, 2), (2, 0))] for face in FACES])
    fields = np.random.default_rng(29).normal(size=(6, 10, 3))
    result = refine(points, cells, surface, fields, periodic=False)
    assert len(result["tetrahedra"]) == 2
    assert interpolation_check(points, cells, fields, result)["maximum_field_error"] < 1e-12
    np.testing.assert_array_equal(result["points"][:10], points)
    np.testing.assert_array_equal(result["displacement"][:, :10], fields)


def test_periodic_cube_refinement_is_conforming_and_budgeted():
    _, _, meshers, *_ = workflow.imports()
    from homogenization import strain_tensor

    mesh = meshers.generate(lambda coordinates_x, *_: coordinates_x * 0 - 1,
        bounds=(-.5, .5) * 3, cells=4, periodic=(True, True, True), order=2, threads=1)
    fields = np.array([mesh.points @ strain_tensor(load).T for load in np.eye(6)])
    blocked = refine(mesh.points, mesh.tetrahedra, mesh.surface, fields, maximum_nodes=len(mesh.points))
    np.testing.assert_array_equal(blocked["tetrahedra"], mesh.tetrahedra)
    current_points, current_cells, current_surface = mesh.points, mesh.tetrahedra, mesh.surface
    for _ in range(2):
        result = refine(current_points, current_cells, current_surface, fields)
        assert result["diagnostics"]["accepted_edge_orbits"] > 0
        assert result["diagnostics"]["periodic"]["periodic_triangles_match"]
        assert positive_quality(result["points"][result["tetrahedra"][:, :4]]).min() >= .1
        assert interpolation_check(current_points, current_cells, fields, result)["maximum_geometry_error"] < 1e-12
        current_points, current_cells, current_surface = (result[key] for key in ("points", "tetrahedra", "surface"))
        fields = result["displacement"]
