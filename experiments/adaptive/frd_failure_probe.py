"""Diagnose stranded wall corners and test a bounded cap-neighbor repair.

The candidate uses the fixed implicit field, linear corner quality and periodic
classes. Admission samples are observers only. It changes no physical fixture,
gate, topology or count. Re-elevation and complete admission are optional because
the short diagnostic mode is useful while inspecting a frozen failure.
"""

import argparse
import importlib
import importlib.metadata
import inspect
import json
from pathlib import Path
import shutil
import sys

import microgen_workflow as w
import numpy as np
from scipy.optimize import minimize

from nested_refinement import EDGES, positive_quality
from periodic_trace import periodic_pairs
import scorecard
import tpms_family_fixtures as fixtures


def frozen_inputs(source, fixture_root, name):
    """Check the frozen specimen and record code before any numerical work."""
    source_result = source / "result.json"
    source_mesh = source / "mesh.npz"
    saved = json.loads(source_result.read_text())
    assert saved["case"] == name
    assert w.digest(source_mesh) == saved["mesh_sha256"], "Source mesh differs from its frozen result"
    fixture_manifest = fixture_root / "result.json"
    manifest = json.loads(fixture_manifest.read_text())
    assert manifest["passed"]
    recorded = {(Path(path) if Path(path).is_absolute() else w.ROOT / path).resolve(): digest
                for path, digest in manifest["inputs"].items()}
    protocol_path = fixture_root / "protocol.json"
    assert w.digest(protocol_path) == recorded[protocol_path.resolve()], "Frozen fixture protocol changed"
    protocol = json.loads(protocol_path.read_text())
    inputs = {Path(__file__).resolve(), source_result.resolve(), source_mesh.resolve(),
              fixture_manifest.resolve(), protocol_path.resolve()}
    selected = {}
    for suffix in ("-definition.json", "-samples.npz", "-specimen.json"):
        fixture = fixture_root / name / (name + suffix)
        expected = recorded[fixture.resolve()]
        assert w.digest(fixture) == expected, f"Frozen fixture changed: {fixture}"
        candidate = source / (name + suffix)
        if candidate.is_file():
            assert w.digest(candidate) == expected, f"Source specimen differs from fixed fixture: {candidate}"
        else:
            candidate = fixture
        selected[suffix] = candidate.resolve()
        inputs.update((fixture.resolve(), candidate.resolve()))

    # The formula and its case wrapper must be the ones used to freeze the field.
    surface_functions = importlib.import_module("microgen.shape.surface_functions")
    for module in (surface_functions, fixtures, scorecard):
        path = Path(inspect.getsourcefile(module)).resolve()
        expected = protocol["source_hashes"][str(path)]
        assert w.digest(path) == expected, f"Frozen field/protocol source changed: {path}"
        inputs.add(path)

    helpers = ("microgen_workflow.py", "tpms_family_fixtures.py", "microgen_catalog.py",
               "scorecard.py", "nested_refinement.py", "periodic_trace.py",
               "comparison_mesh_audit.py", "quadratic_edge_removal_probe.py",
               "quadrature_audit.py", "facet_incidence_audit.py", "box_projection.py",
               "job_parts.py")
    inputs.update((w.HERE / helper).resolve() for helper in helpers)
    inputs.add((w.ROOT.parent / "periodic-microgen/experiments/periodic-fea/periodic.py").resolve())
    for module_name in ("fedoo.lib_elements.tetrahedron", "fedoo.lib_elements.element_base",
                        "scipy.optimize._minimize", "scipy.optimize._slsqp_py", "pueue_progress"):
        inputs.add(Path(inspect.getsourcefile(importlib.import_module(module_name))).resolve())
    inputs.update(path.resolve() for path in (w.ROOT / "crates/meshers-core/src").rglob("*.rs"))
    inputs.add((w.ROOT / "Cargo.lock").resolve())
    assert all(path.is_file() for path in inputs)
    hashes = {str(path): w.digest(path) for path in sorted(inputs)}
    provenance = dict(inputs=hashes, source_mesh_matches_frozen_result=True,
        fixture_files_match_frozen_manifest=True, field_sources_match_frozen_protocol=True,
        selected_specimen_files={suffix: str(path) for suffix, path in selected.items()},
        versions={package: importlib.metadata.version(package)
                  for package in ("numpy", "scipy", "fedoo", "microgen", "pyvista", "vtk", "autograd")},
        python=sys.version, executable=sys.executable)
    return hashes, selected, provenance


def project(case, points):
    """The observer's Newton projection, with exact cap-coordinate locks."""
    foot = points.copy()
    locked = np.isclose(abs(points), 0.5, atol=1e-12, rtol=0)
    levels = np.array([level for level in case.levels if np.isfinite(level)])
    values = case.function(*points.T)
    level = levels[np.argmin(abs(values[:, None] - levels), axis=1)]
    for _ in range(40):
        gradient = scorecard.gradient(case.function, foot)
        gradient[locked] = 0
        residual = case.function(*foot.T) - level
        step = residual[:, None] * gradient / np.maximum((gradient**2).sum(1), 1e-30)[:, None]
        foot -= step
        if np.linalg.norm(step, axis=1).max(initial=0) < 1e-12:
            break
    residual = abs(case.function(*foot.T) - level)
    good = (residual < 1e-9) & (abs(foot).max(axis=1) <= 0.5 + 1e-12)
    np.testing.assert_array_equal(foot[locked], points[locked])
    return foot, good


def cap_faces(points, surface):
    corners = points[surface[:, :3]]
    return np.any(np.all(np.isclose(corners, -0.5, atol=1e-12, rtol=0), axis=1)
                  | np.all(np.isclose(corners, 0.5, atol=1e-12, rtol=0), axis=1), axis=1)


def triangle_quality(xyz):
    """Also the maximum attainable MMG tet q over every possible fourth vertex.

    For a fixed triangle, the best apex lies above its centroid at squared height
    2/9 times the sum of its three squared edge lengths. Substitution gives this
    bound. Thus a boundary value below .1 cannot be fixed by moving its apex.
    """
    area2 = np.linalg.norm(np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]), axis=1)
    squared = sum(((xyz[:, a] - xyz[:, b])**2).sum(1) for a, b in ((0, 1), (1, 2), (2, 0)))
    return 2 * np.sqrt(3) * area2 / squared


def forward(case, points, wall):
    samples = np.einsum("sn,fnd->fsd", w.tri_values(scorecard.BARY, 2), points[wall])
    foot, good = project(case, samples.reshape(-1, 3))
    distance = np.linalg.norm(foot - samples.reshape(-1, 3), axis=1)
    index = int(distance.argmax())
    face, probe = divmod(index, len(scorecard.BARY))
    return dict(maximum=float(distance[index]), failed_faces=int(np.any(
        distance.reshape(len(wall), -1) > 0.002, axis=1).sum()),
        projection_failed=int((~good).sum()), wall_face=face, probe=probe,
        barycentric=scorecard.BARY[probe].tolist(), nodes=wall[face].tolist(),
        point=samples[face, probe].tolist(), foot=foot[index].tolist())


def corner_rows(case, points, cells, surface, wall, classes):
    nodes = np.unique(wall[:, :3])
    foot, good = project(case, points[nodes])
    assert good.all()
    distance = np.linalg.norm(foot - points[nodes], axis=1)
    rows = []
    for group in np.unique(classes[nodes[distance > 1e-8]]):
        selected = np.flatnonzero(classes[nodes] == group)
        copies = nodes[selected]
        incident = np.flatnonzero(np.any(np.isin(cells[:, :4], copies), axis=1))
        faces = np.flatnonzero(np.any(np.isin(surface[:, :3], copies), axis=1))
        candidate = points.copy()
        candidate[copies] += np.mean(foot[selected] - points[copies], axis=0)
        before = points[surface[faces, :3]]
        after = candidate[surface[faces, :3]]
        n0 = np.cross(before[:, 1] - before[:, 0], before[:, 2] - before[:, 0])
        n1 = np.cross(after[:, 1] - after[:, 0], after[:, 2] - after[:, 0])
        ratio = (n0 * n1).sum(1) / (n0 * n0).sum(1)
        bounds = triangle_quality(after)
        rows.append(dict(group=int(group), copies=copies.tolist(), incident=incident.tolist(),
            distance=float(distance[selected].max()),
            q_before=float(positive_quality(points[cells[incident, :4]]).min()),
            q_after_full_projection=float(positive_quality(candidate[cells[incident, :4]]).min()),
            minimum_oriented_face_ratio=float(ratio.min()),
            boundary_apex_bound=float(bounds.min()),
            intrinsically_blocking_boundary_faces=faces[bounds < 0.1].tolist()))
    return sorted(rows, key=lambda row: -row["distance"])


def settle(case, original, cells, rows, policy, progress):
    points = original.copy()
    trace = []
    for sweep in range(16):
        moved = 0
        for row in rows:
            copies = np.array(row["copies"])
            incident = np.array(row["incident"])
            foot, good = project(case, points[copies])
            assert good.all()
            delta = np.mean(foot - points[copies], axis=0)
            if np.linalg.norm(delta) < 1e-10:
                continue
            coordinates = points[cells[incident, :4]]
            before = positive_quality(coordinates).min()
            floor = 0.1 if policy == "full_absolute" else max(0.1, 0.5 * min(0.3, before))
            mask = np.isin(cells[incident, :4], copies)
            fractions = [1.0] if policy == "full_absolute" else 0.5 ** np.arange(9)
            for fraction in fractions:
                trial = coordinates.copy()
                trial[mask] += fraction * delta
                if positive_quality(trial).min() >= floor:
                    points[copies] += fraction * delta
                    moved += 1
                    break
        nodes = np.array([row["copies"][0] for row in rows])
        foot, good = project(case, points[nodes])
        distance = np.linalg.norm(foot - points[nodes], axis=1)
        trace.append(dict(sweep=sweep + 1, moved_classes=moved,
                          maximum_corner_distance=float(distance.max()),
                          classes_over_gate=int((distance > 0.002).sum())))
        progress.update((sweep + 1) / 16, "Settling wall corners", distance=float(distance.max()))
        if moved == 0:
            break
    return points, dict(trace=trace, minimum_quality=float(positive_quality(points[cells[:, :4]]).min()),
        remaining=[dict(group=row["group"], distance=float(d))
                   for row, d in zip(rows, distance, strict=True) if d > 0.002])


def cap_repair(case, points, cells, surface, wall, classes, rows):
    """At most four cap-only neighbor trials for each unresolved wall class.

    Each trial sets the wall class exactly on its level, then maximizes the worst
    local Tet4 q using the two free coordinates of one periodic cap-node class.
    Its radius is one quarter of the longest incident edge. Every incident tet
    of both moving classes participates, including the opposite periodic cap.
    """
    points = points.copy()
    wall_groups = set(classes[np.unique(wall[:, :3])])
    boundary_groups = set(classes[np.unique(surface[:, :3])])
    trials = []
    for row in rows:
        copies = np.array(row["copies"])
        foot, good = project(case, points[copies])
        assert good.all()
        distance = np.linalg.norm(foot - points[copies], axis=1).max()
        if distance <= 0.002:
            continue
        incident = np.array(row["incident"])
        base = points.copy()
        base[copies] += np.mean(foot - points[copies], axis=0)
        bad = incident[positive_quality(base[cells[incident, :4]]) < 0.1]
        neighbors = set(classes[np.unique(cells[bad, :4])]) & (boundary_groups - wall_groups)
        neighbors = sorted(neighbors, key=lambda group: np.linalg.norm(
            (points[group] - points[copies[0]] + 0.5) % 1 - 0.5))[:4]
        for group in neighbors:
            cap_nodes = np.flatnonzero(classes == group)
            locked = np.isclose(abs(points[cap_nodes[0]]), 0.5, atol=1e-12, rtol=0)
            axes = np.flatnonzero(~locked)
            if len(axes) != 2:
                continue
            support = np.flatnonzero(np.any(np.isin(cells[:, :4], np.r_[copies, cap_nodes]), axis=1))
            cap_support = np.flatnonzero(np.any(np.isin(cells[:, :4], cap_nodes), axis=1))
            local = base[cells[support, :4]]
            mask = np.isin(cells[support, :4], cap_nodes)
            edge_vectors = points[cells[cap_support][:, EDGES[:, 1]]] - points[cells[cap_support][:, EDGES[:, 0]]]
            radius = 0.25 * np.linalg.norm(edge_vectors, axis=2).max()
            value_scale = max(float(np.linalg.norm(scorecard.gradient(case.function, points[cap_nodes[:1]]))) * radius, 1e-8)

            def constraints(x):
                delta = np.zeros(3)
                delta[axes] = radius * x[:2]
                candidate = local.copy()
                candidate[mask] += delta
                values = case.function(*(points[cap_nodes] + delta).T)
                material = []
                if np.isfinite(case.levels[0]):
                    material.extend((values - case.levels[0]) / value_scale)
                if np.isfinite(case.levels[1]):
                    material.extend((case.levels[1] - values) / value_scale)
                return np.r_[positive_quality(candidate) - x[2], material]

            initial = np.r_[0.0, 0.0, positive_quality(local).min()]
            result = minimize(lambda x: -x[2], initial, method="SLSQP",
                              bounds=[(-1.0, 1.0), (-1.0, 1.0), (0.0, 1.0)],
                              constraints={"type": "ineq", "fun": constraints},
                              options={"maxiter": 100, "ftol": 1e-12})
            delta = np.zeros(3)
            delta[axes] = radius * result.x[:2]
            candidate = base.copy()
            candidate[cap_nodes] += delta
            q = positive_quality(candidate[cells[support, :4]])
            accepted = bool(q.min() >= 0.1 and constraints(result.x).min() >= -1e-10
                            and abs(candidate[np.r_[copies, cap_nodes]]).max() <= 0.5 + 1e-12)
            trials.append(dict(wall_group=row["group"], cap_group=int(group),
                cap_nodes=cap_nodes.tolist(), support_tets=len(support), radius=float(radius),
                delta=delta.tolist(), iterations=int(result.nit), optimizer_success=bool(result.success),
                minimum_quality=float(q.min()), accepted=accepted))
            if accepted:
                points = candidate
                break
    return points, trials


def elevate(case, points, cells, surface, wall, classes, progress):
    """Diagnostic Python reconstruction, not an established exact native replay.

    It uses the native edge order, 23 Jacobian sample locations, 1e-3 ratio guard,
    and eight halvings followed by flattening. Projection differs: finite
    differences and 40 unclipped Newton steps replace the native gradient,
    clipped steps and bracket fallback. The Python reach check also applies to
    converged Newton feet; native reach only bounds its fallback bracket.
    Periodic displacements are averaged through saved node classes rather than
    grouping current straight-midpoint coordinates. No curved fitting runs here.
    Full independent admission below judges the resulting mesh on its own.
    """
    from fedoo.lib_elements.tetrahedron import Tet10

    points = points.copy()
    pairs = np.sort(cells[:, EDGES].reshape(-1, 2), axis=1)
    edges, first, inverse = np.unique(pairs, axis=0, return_index=True, return_inverse=True)
    mids = cells[:, 4:].ravel()[first]
    assert np.array_equal(mids[inverse], cells[:, 4:].ravel())
    straight = points.copy()
    straight[mids] = points[edges].mean(axis=1)
    points[mids] = straight[mids]
    wall_mids = np.unique(wall[:, 3:])
    foot, good = project(case, straight[wall_mids])
    lengths = np.zeros(len(points))
    lengths[mids] = np.linalg.norm(points[edges[:, 1]] - points[edges[:, 0]], axis=1)
    good &= np.linalg.norm(foot - straight[wall_mids], axis=1) <= 0.5 * lengths[wall_mids]
    curved = wall_mids[good]
    points[curved] = foot[good]
    for group in np.unique(classes[curved]):
        selected = curved[classes[curved] == group]
        points[selected] = straight[selected] + (points[selected] - straight[selected]).mean(axis=0)
    bary = list(np.eye(4)) + [(np.eye(4)[a] + np.eye(4)[b]) / 2 for a, b in EDGES]
    bary += [(np.ones(4) - np.eye(4)[i]) / 3 for i in range(4)] + [np.full(4, 0.25)]
    bary += [np.full(4, b) + np.eye(4)[i] * (a - b)
             for a, b in ((0.7, 0.1), (0.1, 0.3)) for i in range(4)]
    derivatives = Tet10(15).shape_function_derivative(np.array(bary)[:, 1:])
    reductions = np.zeros(len(points), dtype=int)
    is_curved = np.zeros(len(points), dtype=bool)
    is_curved[curved] = True
    for iteration in range(128):
        low = np.full(len(cells), np.inf)
        high = np.full(len(cells), -np.inf)
        coordinates = points[cells]
        for derivative in derivatives:
            determinant = np.linalg.det(np.einsum("an,tnb->tab", derivative, coordinates))
            low = np.minimum(low, determinant)
            high = np.maximum(high, determinant)
        bad = low <= 1e-3 * high
        progress.update(min(0.95, (iteration + 1) / 16), "Rebuilding quadratic geometry", bad_tets=int(bad.sum()))
        if not bad.any():
            return points, dict(iterations=iteration, minimum_jacobian_ratio=float((low / high).min()),
                                partially_curved_edges=int(((reductions > 0) & is_curved).sum()))
        bad_nodes = np.unique(cells[bad, 4:])
        groups = np.unique(classes[bad_nodes[is_curved[bad_nodes]]])
        assert len(groups), "A linear corner star is invalid before curvature"
        selected = curved[np.isin(classes[curved], groups) & is_curved[curved]]
        flatten = selected[reductions[selected] >= 8]
        halve = selected[reductions[selected] < 8]
        points[halve] = 0.5 * (points[halve] + straight[halve])
        reductions[halve] += 1
        points[flatten] = straight[flatten]
        is_curved[flatten] = False
    raise RuntimeError("Bounded elevation did not converge")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--audit-candidate", action="store_true")
    args = parser.parse_args()
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    name = json.loads((args.source / "result.json").read_text())["case"]
    hashes, specimen_files, provenance = frozen_inputs(args.source, args.fixtures, name)
    w.save(args.run / "input-manifest.json", provenance)
    case = fixtures.case(name, args.fixtures)
    with np.load(args.source / "mesh.npz") as saved:
        points, cells, surface = (saved[key].copy() for key in ("points", "tetrahedra", "surface"))
    assert max(len(points), len(cells)) <= w.CAP
    _, classes, _ = periodic_pairs(points, surface, declared_bounds=([-0.5] * 3, [0.5] * 3))
    wall = surface[~cap_faces(points, surface)]
    parts = ["diagnose", "full_absolute", "settle_halves", "cap_repair"]
    if args.audit_candidate:
        parts += ["elevate", "admission"]
    with Job(parts=parts) as job:
        with job.part("diagnose") as part:
            part.update(0, "Finding corner and face witnesses", tets=len(cells), nodes=len(points))
            witness = forward(case, points, wall)
            rows = corner_rows(case, points, cells, surface, wall, classes)
            regular = np.array([[[0., 0., 0.], [1., 0., 0.], [0.5, np.sqrt(3) / 2, 0.]]])
            np.testing.assert_allclose(triangle_quality(regular), 1, atol=1e-14)
            part.update(1, "Classified frozen corner failure", distance=witness["maximum"], classes=len(rows))
        report = dict(case=name, inputs=hashes, input_integrity=provenance,
                      baseline_forward=witness, corners=rows,
                      elements=len(cells), nodes=len(points), policies={}, mechanical_solve_performed=False,
                      reconstruction_scope="Diagnostic Python reconstruction; bitwise/native equivalence unestablished",
                      exact_native_reproduction_verified=False)
        for policy in ("full_absolute", "settle_halves"):
            with job.part(policy) as part:
                candidate, result = settle(case, points, cells, rows, policy, part)
            report["policies"][policy] = result
            job.result(policy, maximum_corner_distance=result["trace"][-1]["maximum_corner_distance"],
                       remaining_classes=len(result["remaining"]))
        with job.part("cap_repair") as part:
            part.update(0, "Repairing blocked cap neighbors", unresolved=len(result["remaining"]))
            candidate, trials = cap_repair(case, candidate, cells, surface, wall, classes, rows)
            report["cap_trials"] = trials
            remaining = corner_rows(case, candidate, cells, surface, wall, classes)
            report["candidate_corner_maximum"] = max((row["distance"] for row in remaining), default=0.0)
            report["candidate_corner_minimum_quality"] = float(positive_quality(candidate[cells[:, :4]]).min())
            part.update(1, "Cap-neighbor trials complete", distance=report["candidate_corner_maximum"])
        np.savez_compressed(args.run / "candidate-corners.npz", points=candidate, tetrahedra=cells, surface=surface)
        report["candidate_corners_sha256"] = w.digest(args.run / "candidate-corners.npz")
        report["corner_file_scope"] = "Corner changes only; existing quadratic nodes are stale until re-elevation"
        if args.audit_candidate:
            with job.part("elevate") as part:
                candidate, report["elevation"] = elevate(case, candidate, cells, surface, wall, classes, part)
            with job.part("admission") as part:
                from functools import partial
                from comparison_mesh_audit import observe

                part.update(0, "Checking complete quadratic candidate", tets=len(cells), nodes=len(points))
                for suffix in ("-samples.npz", "-specimen.json"):
                    copied = args.run / (name + suffix)
                    shutil.copyfile(specimen_files[suffix], copied)
                    assert w.digest(copied) == hashes[str(specimen_files[suffix])]
                geometry, periodic, quality = observe(args.run, name, case, candidate, cells, surface, True,
                    periodic_auditor=partial(periodic_pairs, declared_bounds=([-0.5] * 3, [0.5] * 3)))
                report.update(geometry=geometry, periodic=periodic, quality=quality,
                              admitted=geometry["geometry_admitted"])
                np.savez_compressed(args.run / "mesh.npz", points=candidate, tetrahedra=cells, surface=surface)
                report["mesh_sha256"] = w.digest(args.run / "mesh.npz")
                part.update(1, "Complete admission checked", admitted=geometry["geometry_admitted"])
            job.result("candidate", admitted=geometry["geometry_admitted"],
                       distance=geometry.get("sampled_two_sided_max"), minimum_quality=quality["min"])
        assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
        w.save(args.run / "result.json", report)
    print(f"Saved FRD probe to {args.run.resolve()}")


if __name__ == "__main__":
    main()
