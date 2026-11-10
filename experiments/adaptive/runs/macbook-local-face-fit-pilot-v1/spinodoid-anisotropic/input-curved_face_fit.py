"""Experimental face-wise fitting of Tet10 wall nodes at fixed connectivity.

Only mid-edge nodes move, along fixed constrained normal directions. The fit
uses the implicit function, never benchmark samples, volume targets or tensors.
Off-surface interpolation nodes are intentional. This is not a production API.
"""

import argparse
import json
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

import comparison_fixtures as fixtures
from comparison_mesh_audit import observe
import microgen_workflow as w
from periodic_trace import periodic_pairs
import scorecard


def barycentric(degree):
    return (
        np.array(
            [(i, j, degree - i - j) for i in range(degree + 1) for j in range(degree + 1 - i)],
            dtype=float,
        )
        / degree
    )


def fit(target, points, surface, classes, *, degree=6, power=2, constrained=False, local=False):
    original = points.copy()
    corners = points[surface[:, :3]]
    caps = np.any(
        np.all(abs(corners + 0.5) < 1e-10, axis=1) | np.all(abs(corners - 0.5) < 1e-10, axis=1),
        axis=1,
    )
    wall = surface[~caps]
    mids = np.unique(wall[:, 3:])
    locked = abs(abs(points) - 0.5) < 1e-10
    directions = np.zeros_like(points)
    directions[mids] = scorecard.gradient(target.function, points[mids])
    directions[locked] = 0
    grouped = np.zeros_like(points)
    np.add.at(grouped, classes, directions)
    directions = grouped[classes]
    directions /= np.maximum(np.linalg.norm(directions, axis=1), 1e-30)[:, None]
    roots, inverse = np.unique(classes[mids], return_inverse=True)
    columns = np.full(len(points), -1, dtype=int)
    columns[mids] = inverse
    assert np.all(columns[wall[:, 3:]] >= 0)
    pairs = wall[:, np.array([[0, 1], [1, 2], [2, 0]])]
    lengths = np.linalg.norm(points[pairs[:, :, 1]] - points[pairs[:, :, 0]], axis=2)
    limits = np.full(len(roots), np.inf)
    np.minimum.at(limits, columns[wall[:, 3:]].ravel(), 0.02 * lengths.ravel())
    basis = w.tri_values(barycentric(degree))
    levels = np.array([v for v in target.levels if np.isfinite(v)])
    center = points[wall[:, :3]].mean(axis=1)
    selected = levels[np.argmin(abs(target.function(*center.T)[:, None] - levels), axis=1)]
    area = (
        np.linalg.norm(
            np.cross(
                points[wall[:, 1]] - points[wall[:, 0]], points[wall[:, 2]] - points[wall[:, 0]]
            ),
            axis=1,
        )
        / 2
    )
    weights = np.repeat(np.sqrt(area / len(basis)), len(basis))

    def residual(p, values):
        probes = np.einsum("sn,fnd->fsd", values, p[wall]).reshape(-1, 3)
        gradient = scorecard.gradient(target.function, probes)
        effective = gradient.copy()
        if constrained:
            effective[abs(abs(probes) - 0.5) < 1e-10] = 0
        raw = target.function(*probes.T) - np.repeat(selected, len(values))
        fixed = np.all(abs(abs(probes) - 0.5) < 1e-10, axis=1)
        assert abs(raw[fixed]).max(initial=0) < 1e-10, "Fixed corner is off the implicit wall"
        # All coordinates are fixed at cube corners. Remove only their roundoff
        # residual, which otherwise becomes enormous after dividing by zero.
        raw[fixed] = 0
        gradient[fixed] = 0
        norms = np.maximum(np.linalg.norm(effective, axis=1), 1e-12)
        distance = raw / norms
        return distance, gradient / norms[:, None]

    holdout_basis = w.tri_values(barycentric(11))
    held_before = residual(points, holdout_basis)[0]
    full_wall, full_levels = wall, selected
    active = np.ones(len(roots), dtype=bool)
    if local:
        initial = residual(points, basis)[0].reshape(len(wall), len(basis))
        marked = np.max(abs(initial), axis=1) >= 0.5 * abs(initial).max()
        active[:] = False
        active[np.unique(columns[wall[marked, 3:]])] = True
        support = np.any(active[columns[wall[:, 3:]]], axis=1)
        wall, selected, area = wall[support], selected[support], area[support]
        weights = np.repeat(np.sqrt(area / len(basis)), len(basis))
    records = []
    offsets = np.zeros(len(roots))
    started = time.perf_counter()
    for iteration in range(2 if power == 2 else 6):
        distance, normals = residual(points, basis)
        coeff = (
            np.einsum(
                "fsd,fed->fse", normals.reshape(len(wall), len(basis), 3), directions[wall[:, 3:]]
            )
            * basis[None, :, 3:]
        )
        row = np.repeat(np.arange(len(distance)), 3)
        col = np.broadcast_to(columns[wall[:, None, 3:]], coeff.shape).ravel()
        matrix = coo_matrix((coeff.ravel(), (row, col)), shape=(len(distance), len(roots))).tocsr()
        irls = np.maximum(abs(distance) / max(abs(distance).max(), 1e-30), 1e-3) ** (
            (power - 2) / 2
        )
        iteration_weights = weights * irls
        matrix = matrix.multiply(iteration_weights[:, None]).tocsr()
        rhs = -matrix.T @ (distance * iteration_weights)
        if local:
            diagonal = np.asarray(matrix.power(2).sum(axis=0)).ravel()
            step = rhs / np.maximum(diagonal, diagonal.max() * 1e-10)
            step[~active] = 0
        else:
            normal = (matrix.T @ matrix).tocsc()
            diagonal = normal.diagonal()
            normal += diags(np.maximum(diagonal, diagonal.max() * 1e-12) * 1e-10)
            step = spsolve(normal, rhs)
        step = np.clip(offsets + step, -limits, limits) - offsets
        objective = lambda d: float(np.sum(abs(d) ** power * weights**2) ** (1 / power))
        before = objective(distance)
        accepted = False
        for fraction in (1.0, 0.5, 0.25, 0.125):
            trial_offsets = offsets + fraction * step
            trial = original.copy()
            trial[mids] += directions[mids] * trial_offsets[inverse, None]
            trial_residual = residual(trial, basis)[0]
            after = objective(trial_residual)
            maximum_ok = not local or abs(trial_residual).max() <= abs(distance).max() + 1e-14
            if after < before - 1e-14 and maximum_ok:
                points, offsets, accepted = trial, trial_offsets, True
                break
        records.append(
            dict(
                iteration=iteration,
                accepted=accepted,
                objective_before=before,
                objective_after=after,
                step_fraction=fraction if accepted else 0.0,
            )
        )
    fitted_faces = len(wall)
    wall, selected = full_wall, full_levels
    held_after = residual(points, holdout_basis)[0]
    assert np.array_equal(points[locked], original[locked])
    assert np.array_equal(points[np.unique(surface[:, :3])], original[np.unique(surface[:, :3])])
    return points, dict(
        seconds=time.perf_counter() - started,
        wall_faces=len(wall),
        movable_nodes=len(mids),
        scalar_unknowns=len(roots),
        local_diagonal=local,
        active_scalar_unknowns=int(active.sum()),
        fitted_faces=fitted_faces,
        fit_degree=degree,
        holdout_degree=11,
        iterations=records,
        loss_power=power,
        cap_constrained_distance=constrained,
        max_displacement=float(np.linalg.norm(points - original, axis=1).max()),
        displacement_bound_fraction_of_shortest_edge=0.02,
        holdout_before=dict(
            max=float(abs(held_before).max()), rms=float(np.sqrt(np.mean(held_before**2)))
        ),
        holdout_after=dict(
            max=float(abs(held_after).max()), rms=float(np.sqrt(np.mean(held_after**2)))
        ),
        objective="Corner-area weighted p-norm of normalized implicit residual on symmetric barycentric samples",
        constraints="Fixed corners, fixed cap coordinates, common periodic-class scalar displacement and direction",
        benchmark_samples_used_for_fit=False,
        volume_or_stiffness_target_used=False,
    )


def controls():
    from scorecard import Case

    # Independent exact radial distances, including interpolation nodes.
    corners = np.array([[0.2, 0.0, 1.0], [0.0, 0.2, 1.0], [-0.2, -0.2, 1.0]])
    corners /= np.linalg.norm(corners, axis=1)[:, None]
    mids = corners[np.array([[0, 1], [1, 2], [2, 0]])].mean(axis=1)
    mids /= np.linalg.norm(mids, axis=1)[:, None]
    points = np.vstack((corners, mids)) * 0.3
    target = Case("sphere-control", lambda x, y, z: x * x + y * y + z * z - 0.09, (-np.inf, 0.0))
    surface = np.arange(6)[None]
    moved, _ = fit(target, points, surface, np.arange(6))
    basis = w.tri_values(barycentric(41))
    exact = lambda p: np.abs(np.linalg.norm(basis @ p, axis=1) - 0.3)
    assert exact(moved).max() < exact(points).max()
    assert np.array_equal(moved[:3], points[:3])
    # Cyclic face relabelling must produce the same physical result.
    permuted, _ = fit(target, points, surface[:, [1, 2, 0, 4, 5, 3]], np.arange(6))
    np.testing.assert_allclose(permuted, moved, atol=1e-12, rtol=0)
    plane = points.copy()
    plane[:, 2] = 0.2
    flat, _ = fit(
        Case("plane-control", lambda x, y, z: z - 0.2, (-np.inf, 0.0)), plane, surface, np.arange(6)
    )
    np.testing.assert_array_equal(flat, plane)
    return dict(
        sphere_exact_max_before=float(exact(points).max()),
        sphere_exact_max_after=float(exact(moved).max()),
        plane_unchanged=True,
        fixed_corners=True,
        face_permutation_invariant=True,
    )


def run(
    source, fixture_root, output, *, power=2, constrained=False, local=False, native_field=False
):
    w.imports()
    control = controls()
    saved = json.loads((source / "result.json").read_text())
    name = saved["case"]
    files = [
        Path(__file__),
        source / "result.json",
        source / "mesh.npz",
        Path(w.__file__),
        w.HERE / "comparison_mesh_audit.py",
        w.HERE / "periodic_trace.py",
        w.HERE / "scorecard.py",
        w.HERE / "comparison_fixtures.py",
        w.HERE / "quadratic_edge_removal_probe.py",
    ]
    for suffix in ("-samples.npz", "-specimen.json"):
        path = source / (name + suffix)
        if not path.exists():
            path = fixture_root / name / (name + suffix)
        files.append(path)
        shutil.copyfile(path, output / path.name)
    if name.startswith("spinodoid"):
        files.append(fixtures.coefficients_path(name, fixture_root))
    if name == "graded-gyroid-density-linear":
        files.append(fixture_root / "linear-density-calibration.json")
    hashes = {str(p): w.digest(p) for p in files}
    shutil.copyfile(__file__, output / "input-curved_face_fit.py")
    assert hashes[str(source / "mesh.npz")] == saved["mesh_sha256"]
    with np.load(source / "mesh.npz") as mesh:
        points, cells, surface = (mesh[k].copy() for k in ("points", "tetrahedra", "surface"))
    assert cells.shape[1] == 10 and max(len(points), len(cells)) <= w.CAP
    target = fixtures.case(name, fixture_root)
    classes = (
        periodic_pairs(points, surface)[1] if fixtures.periodic(name) else np.arange(len(points))
    )
    baseline = observe(
        output,
        name,
        target,
        points,
        cells,
        surface,
        fixtures.periodic(name),
        periodic_auditor=periodic_pairs,
    )
    np.testing.assert_allclose(
        baseline[0]["sampled_two_sided_max"],
        saved["geometry"]["sampled_two_sided_max"],
        atol=1e-12,
        rtol=0,
    )
    started = time.perf_counter()
    fit_target = target
    if native_field:
        import meshers
        from scorecard import Case

        compiled = meshers.compile_field(target.function, gradient=False)

        def native(x, y, z):
            return compiled.values(np.column_stack((x, y, z)))

        fit_target = Case(name, native, target.levels)
    moved, diagnostics = fit(
        fit_target, points, surface, classes, power=power, constrained=constrained, local=local
    )
    fit_seconds = time.perf_counter() - started
    diagnostics["native_field"] = native_field
    if native_field:
        probes = np.random.default_rng(20261110).uniform(-0.5, 0.5, (4096, 3))
        np.testing.assert_allclose(
            fit_target.function(*probes.T), target.function(*probes.T), rtol=1e-12, atol=1e-12
        )
    after = observe(
        output,
        name,
        target,
        moved,
        cells,
        surface,
        fixtures.periodic(name),
        periodic_auditor=periodic_pairs,
    )
    np.testing.assert_array_equal(list(baseline[2].values()), list(after[2].values()))
    np.testing.assert_array_equal(points[np.unique(cells[:, :4])], moved[np.unique(cells[:, :4])])
    np.savez_compressed(output / "mesh.npz", points=moved, tetrahedra=cells, surface=surface)
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    w.save(
        output / "result.json",
        dict(
            case=name,
            source=str(source),
            inputs=hashes,
            controls=control,
            nodes=len(points),
            elements=len(cells),
            fit_seconds=fit_seconds,
            diagnostics=diagnostics,
            before=dict(geometry=baseline[0], periodic=baseline[1], quality=baseline[2]),
            after=dict(geometry=after[0], periodic=after[1], quality=after[2]),
            admitted=after[0]["geometry_admitted"],
            mesh_sha256=w.digest(output / "mesh.npz"),
            source_unchanged=True,
            native_only=False,
            new_FE_solve=False,
            old_domain_preserved=False,
            reference_promoted=False,
            comparative_workflow_timing=False,
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path)
    parser.add_argument("--fixtures", type=Path)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--loss-power", type=int, choices=(2, 8), default=2)
    parser.add_argument("--constrained-distance", action="store_true")
    parser.add_argument("--selection", type=Path)
    parser.add_argument("--local", action="store_true")
    parser.add_argument("--native-field", action="store_true")
    args = parser.parse_args()
    assert not (args.run / "result.json").exists()
    if args.selection:
        selection = json.loads(args.selection.read_text())
        shutil.copyfile(args.selection, args.run / "selection.json")
        for entry in selection["cases"]:
            directory = args.run / entry["case"]
            directory.mkdir()
            run(
                Path(entry["source"]),
                Path(entry["fixtures"]),
                directory,
                power=8,
                constrained=True,
                local=args.local,
                native_field=args.native_field,
            )
    else:
        run(
            args.source,
            args.fixtures,
            args.run,
            power=args.loss_power,
            constrained=args.constrained_distance,
            local=args.local,
            native_field=args.native_field,
        )
