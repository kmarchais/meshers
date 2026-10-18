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
