import itertools

import numpy as np

from microgen_benchmark.cases import CASES
from microgen_benchmark.elasticity import solve, worst_error
from microgen_benchmark.geometry import boundary, canonical_mesh


def cube(n):
    axis = np.linspace(-0.5, 0.5, n + 1)
    points = np.array(list(itertools.product(axis, repeat=3)))
    ids = np.arange(len(points)).reshape((n + 1,) * 3)
    tets = []
    for i, j, k in itertools.product(range(n), repeat=3):
        origin = np.array([i, j, k])
        for order in itertools.permutations(range(3)):
            nodes = [ids[tuple(origin)]]
            cursor = origin.copy()
            for a in order:
                cursor[a] += 1
                nodes.append(ids[tuple(cursor)])
            tets.append(nodes)
    p, t, _, _ = canonical_mesh(points, tets)
    return p, t


def test_exact_patch_and_refinement():
    errors = []
    for n in [2, 4]:
        p, t = cube(n)
        metric, fields = solve(p, t, boundary(p, t), 1)
        assert metric["patch_max_error"] < 1e-10
        assert metric["free_residual"] < 1e-10
        assert fields["displacements"].shape == (len(p), 24, 3)
        errors.append(metric["worst_energy_error"])
    assert 0.4 < errors[1] / errors[0] < 0.8


def test_units_do_not_change_fea():
    p, t = cube(2)
    surface = boundary(p, t)
    a, _ = solve(p, t, surface, 1)
    b, _ = solve(p * 0.001, t, surface, 0.001)
    assert np.isclose(a["worst_energy_error"], b["worst_energy_error"], rtol=1e-12)


def test_sparse_ordering_does_not_change_solution_or_energy():
    p, t = cube(4)
    surface = boundary(p, t)
    a, fa = solve(p, t, surface, 1, ordering="COLAMD")
    b, fb = solve(p, t, surface, 1, ordering="MMD_AT_PLUS_A")
    assert np.isclose(a["worst_energy_error"], b["worst_energy_error"], rtol=1e-10)
    assert np.allclose(fa["displacements"], fb["displacements"], atol=1e-10)
    assert np.allclose(fa["error_energy"], fb["error_energy"], atol=1e-10)


def test_worst_combination_is_basis_invariant():
    reference = np.diag([1.0, 2.0, 3.0])
    error = np.diag([0.01, 0.08, 0.03])
    transform = np.array([[1.0, 2.0, 0], [0, 1.0, 0.5], [0.2, 0, 1.0]])
    first, _ = worst_error(error, reference)
    second, _ = worst_error(transform.T @ error @ transform, transform.T @ reference @ transform)
    assert np.isclose(first, 0.2)
    assert np.isclose(first, second)


def test_grading_and_scale_are_actual_distinct_fields():
    uniform, graded, scaled = [
        next(c for c in CASES if c.id == key)
        for key in ["gyroid", "graded-gyroid", "scaled-gyroid"]
    ]
    x = np.array([-0.4, 0.1, 0.4])
    y = np.array([0.12, 0.13, 0.14])
    z = np.array([0.22, 0.23, 0.24])
    assert not np.allclose(uniform.field(x, y, z), graded.field(x, y, z))
    assert np.allclose(uniform.field(x, y, z), scaled.field(x * 0.001, y * 0.001, z * 0.001))
