import hashlib
from dataclasses import replace
from itertools import product
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
    if workflow == "mmg-microgen":
        executable = Path(timing["mmg_cli_executable"])
        assert executable.parent.name == "bin"
        assert executable.parent.parent.name == "mmgpy"
        assert hashlib.sha256(executable.read_bytes()).hexdigest() == timing["mmg_cli_sha256"]
    if workflow == "mmg-tetra-tagged":
        assert set(timing["mmg_boundary_patch_counts"]) == set(map(str, range(3, 9)))
        assert timing["mmg_required_cap_triangles"] == 0
        assert geometry["q_min"] >= 0.1
    if workflow in {"mmg", "mmg-tetra", "mmg-tetra-clean"}:
        assert timing["mmg_required_cap_triangles"] > 0


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


@pytest.mark.parametrize("scale", [0.001, 1, 1000])
def test_traceable_bcc_matches_vector_field_and_has_finite_compiled_gradients(scale):
    import meshers

    case = replace(next(c for c in CASES if c.id == "implicit-lattice"), scale=scale)
    samples = np.vstack(
        [
            np.random.default_rng(20261012).uniform(-0.5, 0.5, (128, 3)),
            list(product((-0.5, 0, 0.5), repeat=3)),
        ]
    )
    points = samples * scale
    # Retain the original vector expression as an independent regression oracle.
    distance = np.full(len(samples), np.inf)
    for corner in product((-0.5, 0.5), repeat=3):
        end = np.asarray(corner)
        t = np.clip(np.sum(samples * end, axis=1) / np.dot(end, end), 0, 1)
        distance = np.minimum(distance, np.linalg.norm(samples - t[:, None] * end, axis=1))
    reference = np.tanh(distance / case.offset - 1)
    compiled = meshers.compile_field(case.field)
    assert np.allclose(case.field(*points.T), reference, atol=2e-14, rtol=0)
    assert np.allclose([compiled(*p) for p in points], reference, atol=2e-14, rtol=0)
    gradients = np.asarray([compiled.gradient(*p) for p in points])
    assert np.isfinite(gradients).all()
    # Random points avoid the capsule junctions and axes with nonunique gradients.
    numeric = []
    for axis in range(3):
        delta = np.zeros(3)
        delta[axis] = scale * 1e-6
        numeric.append(
            (case.field(*(points[:128] + delta).T) - case.field(*(points[:128] - delta).T))
            / (2 * delta[axis])
        )
    assert np.allclose(gradients[:128], np.asarray(numeric).T, rtol=1e-5, atol=1e-6 / scale)


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
