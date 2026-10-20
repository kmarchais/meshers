"""Small checks for the experimental reference element and field checkpoints."""

from types import SimpleNamespace
import weakref

import numpy as np
import pytest

import microgen_workflow as workflow

pytest.importorskip("fedoo")


def test_cubic_basis_reproduces_quadratic_geometry():
    from cubic_reference import register_cubic, basis_checks

    result = basis_checks(register_cubic()())
    assert result["nodal_and_polynomial_identities"]
    assert result["map_derivative_embedding_error"] < 1e-12


def test_cubic_surface_basis_matches_reference_element_faces():
    from cubic_reference import register_cubic

    element = register_cubic()()
    selected = np.flatnonzero(element.multi_indices[:, 3] == 0)
    indices = element.multi_indices[selected, :3]
    nodal = indices / 3
    np.testing.assert_allclose(workflow.tri_values(nodal, degree=3), np.eye(10), atol=1e-14)
    barycentric = np.random.default_rng(20261020).dirichlet(np.ones(3), 128)
    locations = barycentric @ element.xi_nd[:3]
    actual = workflow.tri_values(barycentric, degree=3)
    np.testing.assert_allclose(actual, element.shape_function(locations)[:, selected], atol=1e-14)
    np.testing.assert_allclose(actual.sum(axis=1), 1, atol=1e-14)
    for first in range(4):
        for second in range(4 - first):
            expected = barycentric[:, 0]**first * barycentric[:, 1]**second
            values = nodal[:, 0]**first * nodal[:, 1]**second
            np.testing.assert_allclose(actual @ values, expected, atol=1e-14)


def test_cubic_node_budget_counts_shared_edges_and_faces():
    from cubic_reference import cubic_node_count

    assert cubic_node_count(np.array([[0, 1, 2, 3]])) == 20
    assert cubic_node_count(np.array([[0, 1, 2, 3], [0, 1, 2, 4]])) == 30


@pytest.mark.parametrize("block_rows", [1, 3, 32])
@pytest.mark.parametrize("perturbation", [0., 1e-14, 1e-8])
def test_blocked_symmetry_matches_full_frobenius_test(block_rows, perturbation):
    from scipy.sparse import csr_matrix
    from amg_reference import verify_symmetric

    matrix = np.random.default_rng(20261020).normal(size=(11, 11))
    matrix += matrix.T + 20 * np.eye(11)
    matrix[1, 9] += perturbation
    matrix = csr_matrix(matrix)
    skew = matrix - matrix.T
    expected = np.linalg.norm(skew.data) / np.linalg.norm(matrix.data) < 1e-12
    if expected:
        verify_symmetric(matrix, block_rows=block_rows)
    else:
        with pytest.raises(AssertionError):
            verify_symmetric(matrix, block_rows=block_rows)


def test_blocked_symmetry_never_subtracts_whole_matrix(monkeypatch):
    from scipy.sparse import csr_matrix
    from amg_reference import verify_symmetric

    original = csr_matrix.__sub__
    shapes = []
    def limited_subtraction(matrix, other):
        shapes.append(matrix.shape)
        assert matrix.shape[0] <= 3
        return original(matrix, other)
    monkeypatch.setattr(csr_matrix, "__sub__", limited_subtraction)
    verify_symmetric(csr_matrix(np.eye(11)), block_rows=3)
    assert shapes == [(3, 11), (3, 11), (3, 11), (2, 11)]


def test_nonstandard_geometry_rejects_implicit_linear_fallback(tmp_path):
    with pytest.raises(ValueError, match="explicit registered element"):
        workflow.geometry(tmp_path, "smooth-lattice", np.zeros((20, 3)),
                          np.arange(20)[None], np.arange(10)[None])


def test_solved_fields_are_checkpointed_in_load_order():
    _, _, meshers, *_ = workflow.imports()
    from stress_recovery import solve
    from homogenization import MATERIAL

    mesh = meshers.generate(lambda coordinates_x, *_: coordinates_x * 0 - 1,
        bounds=(-.5, .5) * 3, cells=4, minimum_quality=.1,
        periodic=(True, True, True), threads=1)
    saved = []
    fields, tensor, records = solve(mesh.points, mesh.tetrahedra,
        field_checkpoint=lambda load, field: saved.append((load, field.copy())))
    assert [load for load, _ in saved] == list(range(1, 7))
    np.testing.assert_array_equal(np.array([field for _, field in saved]), fields)
    np.testing.assert_allclose(tensor, MATERIAL, atol=1e-9, rtol=0)
    assert len(records) == 6


def test_nonstandard_elements_require_explicit_registration_choice():
    from stress_recovery import solve

    with pytest.raises(KeyError, match="20"):
        solve(np.zeros((20, 3)), np.arange(20)[None])


def test_four_point_affine_tet10_block_amg_matches_original_solver():
    _, _, meshers, *_ = workflow.imports()
    from amg_reference import AuditedAMG
    from homogenization import MATERIAL
    from stress_recovery import solve, vector_first_reactions
    from tet10_integration import integrated_tensors

    mesh = meshers.generate(lambda coordinates_x, *_: coordinates_x * 0 - 1,
        bounds=(-.5, .5) * 3, cells=4, order=2, minimum_quality=.1,
        periodic=(True, True, True), threads=1)
    original_fields, original_tensor, _ = solve(mesh.points, mesh.tetrahedra, n_elm_gp=4)
    solver = AuditedAMG(elastic_blocks=True)
    saved = []
    fields, tensor, records = solve(mesh.points, mesh.tetrahedra, n_elm_gp=4, solver=solver,
        field_checkpoint=lambda load, field: saved.append((load, field.copy())),
        reaction_reader=vector_first_reactions, release_reduced_cache=True)
    np.testing.assert_allclose(tensor, original_tensor, atol=1e-9, rtol=0)
    np.testing.assert_allclose(tensor, MATERIAL, atol=1e-9, rtol=0)
    np.testing.assert_allclose(fields, original_fields, atol=1e-9, rtol=0)
    assert [load for load, _ in saved] == list(range(1, 7))
    np.testing.assert_array_equal(np.array([field for _, field in saved]), fields)
    stress, energy, volume = integrated_tensors(mesh.points, mesh.tetrahedra, fields)
    np.testing.assert_allclose(stress, tensor, atol=1e-9, rtol=0)
    np.testing.assert_allclose(energy, tensor, atol=1e-9, rtol=0)
    assert abs(volume - 1) < 1e-12
    assert len(records) == 6 and sum(record["hierarchy_reused"] for record in records) == 5


@pytest.mark.parametrize("changed", ["rhs", "matrix"])
def test_saved_prefix_rejects_changed_equations(changed):
    from scipy.sparse import eye
    from amg_reference import VerifiedPrefixAMG

    solver = VerifiedPrefixAMG([np.ones((1, 3)), np.ones((1, 3))])
    solver.problem = SimpleNamespace(mesh=SimpleNamespace(n_nodes=1),
        n_node_dof=3, _dof_free=np.arange(3))
    matrix = eye(3, format="csr")
    if changed == "rhs":
        with pytest.raises(AssertionError, match="Saved prefix"):
            solver(matrix, np.ones(3) * 2)
    else:
        solver(matrix, np.ones(3))
        with pytest.raises(AssertionError):
            solver(matrix * 2, np.ones(3) * 2)


def test_prefix_problem_binding_does_not_retain_previous_problem():
    from amg_reference import VerifiedPrefixAMG

    class Problem:
        pass

    solver = VerifiedPrefixAMG([np.ones((1, 3))])
    previous = Problem()
    reference = weakref.ref(previous)
    solver.problem = previous
    assert solver.delegate.problem is previous
    replacement = Problem()
    solver.problem = replacement
    del previous
    assert reference() is None
    assert solver.problem is replacement and solver.delegate.problem is replacement


def test_completed_reduced_matrices_are_released():
    _, _, meshers, *_ = workflow.imports()
    from amg_reference import AuditedAMG
    from stress_recovery import solve
    from homogenization import MATERIAL

    class ObservedAMG(AuditedAMG):
        def __init__(self):
            super().__init__(elastic_blocks=True)
            self.matrices = []

        def __call__(self, matrix, rhs, **kwargs):
            self.matrices.append(weakref.ref(matrix))
            return super().__call__(matrix, rhs, **kwargs)

    mesh = meshers.generate(lambda coordinates_x, *_: coordinates_x * 0 - 1,
        bounds=(-.5, .5) * 3, cells=4, minimum_quality=.1,
        periodic=(True, True, True), threads=1)
    solver = ObservedAMG()
    _, tensor, records = solve(mesh.points, mesh.tetrahedra, solver=solver,
        release_reduced_cache=True)
    assert len(records) == 6 and all(reference() is None for reference in solver.matrices)
    np.testing.assert_allclose(tensor, MATERIAL, atol=1e-9, rtol=0)
