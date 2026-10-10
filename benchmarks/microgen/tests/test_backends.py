from pathlib import Path

import numpy as np
import pytest

from microgen_benchmark.backends import generate, reopenable_temporary_file
from microgen_benchmark.cases import CASES, WORKFLOWS
from microgen_benchmark.elasticity import solve
from microgen_benchmark.geometry import canonical_mesh, validate


@pytest.mark.parametrize("workflow", WORKFLOWS)
def test_real_backend_cube_pipeline(workflow):
    config = {
        "threads": 1,
        "minimum_quality": 0.1,
        "geometry_tolerance": 0.01,
        "max_tetrahedra": 250_000,
        "mmg_boundary": "caps",
    }
    points, tets, timing = generate(CASES[0], workflow, 4, config)
    points, tets, _, _ = canonical_mesh(points, tets)
    geometry, surface = validate(points, tets, CASES[0])
    fea, _ = solve(points, tets, surface, 1)
    assert np.isclose(geometry["volume"], 1)
    assert fea["patch_max_error"] < 1e-8
    assert timing["mesh_s"] > 0


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.id)
def test_actual_microgen_fields_match_case(case):
    shape = case.shape(4)
    assert case.verify_field(shape, 4) == 125
    if case.implicit:
        return
    sample = np.asarray(shape.grid.points)
    raw = case.field(*sample.T) * 0.5 * case.width(*sample.T)
    assert np.allclose(raw, shape.grid.point_data["surface"], rtol=1e-12, atol=1e-12)
    offset = np.broadcast_to(np.asarray(shape.offset).ravel(order="F"), (len(sample),))
    assert np.allclose(offset, case.width(*sample.T), rtol=1e-12, atol=1e-12)


def test_implicit_lattice_known_material_points():
    case = next(c for c in CASES if c.id == "implicit-lattice")
    assert case.material(np.array([0, 0.25]), np.array([0, 0.25]), np.array([0, 0.25])).all()
    assert not case.material(np.array([0.5]), np.array([0]), np.array([0]))[0]


def test_spinodoid_resolution_does_not_change_field():
    from microgen_benchmark.cases import spinodoid_source

    case = next(c for c in CASES if c.id == "spinodoid")
    x = np.random.default_rng(71).uniform(-0.5, 0.5, (3, 91))
    low, high = case.shape(4), case.shape(12)
    assert np.array_equal(low.evaluate(*x), high.evaluate(*x))
    original = spinodoid_source(1).evaluate(*(a + 0.5 for a in x))
    assert np.array_equal(case.field(*x), np.tanh(original))
    assert np.array_equal(case.material(*x), original <= 0)


def test_microgen_temporary_file_reopens_and_cleans_up():
    with reopenable_temporary_file(suffix=".mesh", delete=True) as temporary:
        path = Path(temporary.name)
        path.write_bytes(b"reopened")
        assert path.read_bytes() == b"reopened"
    assert not path.exists()
