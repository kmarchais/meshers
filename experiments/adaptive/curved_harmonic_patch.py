"""Bounded geometry experiment: extend wall-midpoint projection into nearby edge nodes.

Corners and connectivity stay fixed. The scalar quadratic Laplacian is only a
geometry extension operator, not a material or effective-stiffness calculation.
"""

import argparse
from functools import partial
import json
from pathlib import Path
import shutil
import time

import microgen_workflow as w
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def laplacian(corners):
    from fedoo.lib_elements.tetrahedron import Tet4, Tet10

    a, b = 0.5854101966249685, 0.1381966011250105
    bary = np.full((4, 4), b)
    np.fill_diagonal(bary, a)
    derivative_linear = Tet4(1).shape_function_derivative(np.array([[0.25, 0.25, 0.25]]))[0]
    edge = np.einsum("an,tnb->tab", derivative_linear, corners)
    inverse = np.linalg.inv(edge)
    determinant = np.linalg.det(edge)
    assert np.all(determinant > 0)
    matrix = np.zeros((len(corners), 10, 10))
    for derivative in Tet10(15).shape_function_derivative(bary[:, 1:]):
        gradient = np.einsum("tab,bn->tan", inverse, derivative)
        matrix += np.einsum("tai,taj->tij", gradient, gradient) * determinant[:, None, None] / 24
    return matrix


def controls():
    from nested_refinement import EDGES

    corners = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    matrix = laplacian(corners[None])[0]
    np.testing.assert_allclose(matrix, matrix.T, atol=1e-14)
    np.testing.assert_allclose(matrix @ np.ones(10), 0, atol=1e-14)
    points = np.vstack((corners, corners[EDGES].mean(axis=1)))
    for axis in range(3):
        np.testing.assert_allclose(points[:, axis] @ matrix @ points[:, axis], 1 / 6, atol=1e-14)
    assert np.linalg.eigvalsh(matrix)[1] > 0
    return dict(
        symmetric=True,
        constants_null=True,
        affine_energy_exact=True,
        positive_on_nonconstant_space=True,
    )


def experiment(source, fixture_root, output):
    from comparison_mesh_audit import observe
    from periodic_trace import periodic_pairs
    from quadratic_edge_removal_probe import measure
    import tpms_family_fixtures as fixtures
    import scorecard

    saved = json.loads((source / "result.json").read_text())
    name = saved["case"]
    target = fixtures.case(name, fixture_root)
    hashes = {
        str(f): w.digest(f)
        for f in (
            Path(__file__),
            source / "result.json",
            source / "mesh.npz",
            w.HERE / "periodic_trace.py",
            w.HERE / "quadratic_edge_removal_probe.py",
            w.HERE / "comparison_mesh_audit.py",
            w.HERE / "scorecard.py",
            w.HERE / "microgen_workflow.py",
            fixture_root / name / (name + "-definition.json"),
        )
    }
    for suffix in ("-samples.npz", "-specimen.json"):
        path = source / (name + suffix)
        hashes[str(path)] = w.digest(path)
        shutil.copyfile(path, output / path.name)
    with np.load(source / "mesh.npz") as mesh:
        points, cells, surface = (mesh[k].copy() for k in ("points", "tetrahedra", "surface"))
    assert max(len(points), len(cells)) <= w.CAP
    checked = partial(periodic_pairs, declared_bounds=([-0.5] * 3, [0.5] * 3))
    _, classes, _ = checked(points, surface)
    corners = points[surface[:, :3]]
    cap = np.any(
        np.all(w.on_box_plane(corners, -0.5), axis=1)
        | np.all(w.on_box_plane(corners, 0.5), axis=1),
        axis=1,
    )
    wall = surface[~cap]
    mids = np.unique(wall[:, 3:])
    probe = points[mids]
    locked = np.isclose(abs(probe), 0.5, atol=1e-12, rtol=0)
    levels = np.array([v for v in target.levels if np.isfinite(v)])
    values = target.function(*probe.T)
    level = levels[np.argmin(abs(values[:, None] - levels), axis=1)]
    gradient = scorecard.gradient(target.function, probe)
    gradient[locked] = 0
    estimate = abs(values - level) / np.maximum(np.linalg.norm(gradient, axis=1), 1e-30)
    marked = np.unique(classes[mids[estimate > 0.001]])
    selected = mids[np.isin(classes[mids], marked)]
    assert len(selected) > 0
    selected_mask = np.isin(mids, selected)
    foot = probe[selected_mask].copy()
    selected_locks, selected_levels = locked[selected_mask], level[selected_mask]
    for _ in range(40):
        gradient = scorecard.gradient(target.function, foot)
        gradient[selected_locks] = 0
        residual = target.function(*foot.T) - selected_levels
        step = (
            residual[:, None]
            * gradient
            / np.maximum((gradient * gradient).sum(axis=1), 1e-30)[:, None]
        )
        foot -= step
        if np.linalg.norm(step, axis=1).max() < 1e-12:
            break
    residual = abs(target.function(*foot.T) - selected_levels)
    assert np.isfinite(foot).all() and residual.max() < 1e-9 and abs(foot).max() <= 0.5 + 1e-10
    np.testing.assert_array_equal(foot[selected_locks], points[selected][selected_locks])
    started = time.perf_counter()
    mid_classes = classes[cells[:, 4:]]
    seed_tets = np.any(np.isin(mid_classes, marked), axis=1)
    free = np.setdiff1d(np.unique(mid_classes[seed_tets]), np.unique(classes[mids]))
    moving = np.union1d(free, marked)
    support = np.any(np.isin(mid_classes, moving), axis=1)
    local_classes = np.unique(mid_classes[support])
    assert len(local_classes) <= 12000, "Local resource preflight: too many geometry unknowns"
    indices = np.searchsorted(local_classes, mid_classes[support])
    element_matrices = laplacian(points[cells[support, :4]])[:, 4:, 4:]
    rows = np.broadcast_to(indices[:, :, None], element_matrices.shape).ravel()
    cols = np.broadcast_to(indices[:, None, :], element_matrices.shape).ravel()
    matrix = coo_matrix(
        (element_matrices.ravel(), (rows, cols)), shape=(len(local_classes),) * 2
    ).tocsr()
    target_move = np.zeros((len(local_classes), 3))
    for group in marked:
        mask = classes[selected] == group
        target_move[np.searchsorted(local_classes, group)] = (
            foot[mask] - points[selected[mask]]
        ).mean(axis=0)
    harmonic = target_move.copy()
    residuals = []
    for axis in range(3):
        active = np.isin(local_classes, free) & (abs(points[local_classes, axis]) != 0.5)
        rhs = -(matrix @ target_move[:, axis])[active]
        system = matrix[active][:, active].tocsc()
        harmonic[active, axis] = spsolve(system, rhs)
        error = float(
            np.linalg.norm(system @ harmonic[active, axis] - rhs) / max(np.linalg.norm(rhs), 1e-30)
        )
        assert error < 1e-9
        residuals.append(error)
    move = np.zeros_like(points)
    members = np.flatnonzero(np.isin(classes, moving))
    move[members] = harmonic[np.searchsorted(local_classes, classes[members])]
    move[abs(points) == 0.5] = 0
    vertex_ids = np.unique(cells[:, :4])
    np.testing.assert_array_equal(move[vertex_ids], 0)
    limits = np.full(len(points), np.inf)
    from nested_refinement import EDGES

    length = np.linalg.norm(points[cells[:, EDGES[:, 1]]] - points[cells[:, EDGES[:, 0]]], axis=2)
    np.minimum.at(limits, cells[:, 4:].ravel(), 0.25 * length.ravel())
    trials = []
    result = points.copy()
    chosen = 0.0
    for fraction in (1.0, 0.5, 0.25, 0.125):
        candidate = points + fraction * move
        minimum, _ = measure(candidate[cells[support]])
        allowed = bool(
            minimum.min() > 0
            and abs(candidate).max() <= 0.5 + 1e-12
            and np.all(np.linalg.norm(fraction * move, axis=1) <= limits + 1e-12)
        )
        trials.append(
            dict(fraction=fraction, minimum_sampled_jacobian=float(minimum.min()), allowed=allowed)
        )
        if allowed:
            result, chosen = candidate, fraction
            break
    seconds = time.perf_counter() - started
    np.testing.assert_array_equal(points[vertex_ids], result[vertex_ids])
    np.testing.assert_array_equal(points[abs(points) == 0.5], result[abs(points) == 0.5])
    before = observe(output, name, target, points, cells, surface, True, periodic_auditor=checked)
    after = observe(output, name, target, result, cells, surface, True, periodic_auditor=checked)
    np.savez_compressed(output / "mesh.npz", points=result, tetrahedra=cells, surface=surface)
    assert all(w.digest(Path(f)) == h for f, h in hashes.items())
    w.save(
        output / "result.json",
        dict(
            case=name,
            controls=controls(),
            inputs=hashes,
            selected_wall_nodes=len(selected),
            free_periodic_classes=len(free),
            local_classes=len(local_classes),
            support_tetrahedra=int(support.sum()),
            geometry_solve_relative_residuals=residuals,
            trial_steps=trials,
            accepted_fraction=chosen,
            geometry_extension_seconds=seconds,
            before=dict(geometry=before[0], periodic=before[1], quality=before[2]),
            after=dict(geometry=after[0], periodic=after[1], quality=after[2]),
            admitted=after[0]["geometry_admitted"],
            elements=len(cells),
            nodes=len(points),
            corners_and_connectivity_unchanged=True,
            mechanical_solve_performed=False,
            comparative_workflow_timing=False,
            mesh_sha256=w.digest(output / "mesh.npz"),
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--fixtures", type=Path, required=True)
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args()
    w.imports()
    controls()
    experiment(a.source, a.fixtures, a.run)
