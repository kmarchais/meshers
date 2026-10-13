"""The automatic API must retain six audited solves with a direct backend present."""

from types import SimpleNamespace

import numpy as np

import homogenization
from microgen_benchmark.geometry import boundary, canonical_mesh
from periodic import periodic_pairs


def test_both_periodic_formulations_use_the_supplied_callback(monkeypatch):
    points = np.array(np.meshgrid(*([np.arange(3) / 2 - .5] * 3), indexing="ij")).reshape(3, -1).T
    pattern = np.array([[0, 1, 3, 7], [0, 3, 2, 7], [0, 2, 6, 7],
                        [0, 6, 4, 7], [0, 4, 5, 7], [0, 5, 1, 7]])
    corners = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
                        [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]])
    tets = []
    for index in np.ndindex(2, 2, 2):
        ids = np.ravel_multi_index((corners + index).T, (3, 3, 3))
        tets.extend(ids[pattern])
    points, tets, _, _ = canonical_mesh(points, np.array(tets))
    surface = boundary(points, tets)
    pairs, _, _ = periodic_pairs(points, surface)
    solvers = []

    class Tracked(homogenization.AuditedCG):
        def __init__(self):
            super().__init__()
            solvers.append(self)

    monkeypatch.setattr(homogenization, "AuditedCG", Tracked)
    mesh = SimpleNamespace(points=points + .5, tetrahedra=tets, surface=surface,
                           periodic_pairs=pairs, diagnostics={"minimum_mmg_quality": .1})
    result = homogenization.homogenize(mesh)
    assert [len(solver.records) for solver in solvers] == [6, 6]
    np.testing.assert_allclose(result["stiffness_engineering_voigt"], homogenization.MATERIAL, atol=1e-10)
