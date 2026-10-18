"""Small checks for the experimental reference element and field checkpoints."""

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
