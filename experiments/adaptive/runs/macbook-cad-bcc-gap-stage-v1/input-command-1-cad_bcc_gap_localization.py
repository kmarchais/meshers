"""Read-only localization of an exactly reproduced BCC linear rejection.

Run inside macbook_stage.py with 240 s / 8 GiB bounds and --interrupt none.
The caller must inspect the capture first and supply its final node/tet counts.
This reads frozen snapshots, never remeshes, changes options or promotes a mesh.
"""

import argparse
from itertools import product
import json
import math
from pathlib import Path
import re
import shutil
import sys

import microgen_workflow as w

KINDS = ("corner-0", "corner-1", "corner-2", "centroid", "edge-01", "edge-12", "edge-20")


def finish_samples(np, points, surface, bounds):
    """Preserve finish's face/sample order, duplicates and exact cap predicates."""
    lo, hi = np.asarray(bounds)
    triangles = points[surface]
    caps = np.any(np.all(triangles == lo, axis=1) | np.all(triangles == hi, axis=1), axis=1)
    ids = np.flatnonzero(~caps)
    p = triangles[ids]
    samples = np.stack((p[:, 0], p[:, 1], p[:, 2], (p[:, 0]+p[:, 1]+p[:, 2])/3,
        (p[:, 0]+p[:, 1])/2, (p[:, 1]+p[:, 2])/2, (p[:, 2]+p[:, 0])/2), axis=1)
    locks = np.zeros(samples.shape, dtype=bool)
    locks[:, :3] = (p == lo) | (p == hi)
    for index, (a, b) in enumerate(((0, 1), (1, 2), (2, 0)), 4):
        locks[:, index] = ((p[:, a] == lo) & (p[:, b] == lo)) | ((p[:, a] == hi) & (p[:, b] == hi))
    return samples.reshape(-1, 3), locks.reshape(-1, 3), np.repeat(ids, 7), np.tile(np.arange(7), len(ids))


def distances(np, compiled, samples, locks, bounds, periodic, scale):
    canonical = samples.copy()
    lo, hi = np.asarray(bounds)
    for axis, enabled in enumerate(periodic):
        if enabled:
            canonical[canonical[:, axis] == hi[axis], axis] = lo[axis]
    values = compiled.values(canonical)/scale
    gradients = np.array([compiled.gradient(*point) for point in samples])/scale
    full = np.sqrt((gradients[:, 0]**2+gradients[:, 1]**2)+gradients[:, 2]**2)
    free = np.where(locks, 0., gradients)
    norms = np.sqrt((free[:, 0]**2+free[:, 1]**2)+free[:, 2]**2)
    errors = np.full(len(samples), np.inf)
    np.divide(abs(values), norms, out=errors, where=norms > 0)
    errors[values == 0] = 0
    return errors, values, gradients, full, norms, canonical


def controls(np, meshers, fixture):
    def cylinder(axis):
        return dict(center=[0., 0., 0.], axis=axis, radius=.1, height=2.)

    z = cylinder([0., 0., 1.])
    compiled = meshers.compile_field(fixture.field([z]))
    theta = .2
    points = np.array([[.1*math.cos(theta), -.1*math.sin(theta), 0.],
                       [.1*math.cos(theta), .1*math.sin(theta), 0.], [.1, 0., .1]])
    samples, locks, face_ids, kinds = finish_samples(np, points, np.array([[0, 1, 2]]), [[-.5]*3, [.5]*3])
    assert len(samples) == 7 and kinds.tolist() == list(range(7)) and face_ids.tolist() == [0]*7
    error, *_ = distances(np, compiled, samples, locks, [[-.5]*3, [.5]*3], [False]*3, 1.)
    expected = .1*(1-math.cos(theta))
    assert abs(error[4]-expected) < 1e-15 and max(error[:3]) < 1e-15
    # An actual box-plane lock removes part of a radial normal.
    y = cylinder([0., 1., 0.])
    compiled_y = meshers.compile_field(fixture.field([y]))
    cap_point = np.array([[.06, 0., .06]])
    constrained, *_ = distances(np, compiled_y, cap_point, np.array([[False, False, True]]),
        [[-.5]*3, [.5, .5, .06]], [False]*3, 1.)
    unconstrained, *_ = distances(np, compiled_y, cap_point, np.zeros((1, 3), bool),
        [[-.5]*3, [.5, .5, .06]], [False]*3, 1.)
    assert abs(constrained[0]-unconstrained[0]*math.sqrt(2)) < 1e-15
    x = cylinder([1., 0., 0.])
    junction = meshers.compile_field(fixture.field([x, y]))
    point = np.array([[.08, .08, 0.]])
    joint, *_ = distances(np, junction, point, np.zeros((1, 3), bool), [[-.5]*3, [.5]*3], [False]*3, 1.)
    np.testing.assert_array_equal(junction.gradient(*point[0]), [.5, .5, 0.])
    assert abs(joint[0]-.02*math.sqrt(2)) < 1e-15
    # Cap-only triangles are skipped, including triangles on an exact upper plane.
    cap = np.array([[0., 0., .5], [.1, 0., .5], [0., .1, .5]])
    assert len(finish_samples(np, cap, np.array([[0, 1, 2]]), [[-.5]*3, [.5]*3])[0]) == 0
    cap[2, 2] = .4
    _, cap_locks, _, _ = finish_samples(np, cap, np.array([[0, 1, 2]]), [[-.5]*3, [.5]*3])
    assert cap_locks[:, 2].tolist() == [True, True, False, False, True, False, False]
    trace = meshers.compile_field(lambda x, y, z: x*x+y-.2)
    _, _, trace_gradient, *_ = distances(np, trace, np.array([[.5, .1, 0.]]),
        np.array([[True, False, False]]), [[-.5]*3, [.5]*3], [True, False, False], .5)
    np.testing.assert_array_equal(trace_gradient, [[2., 2., 0.]])
    assert compiled(.02, 0., 1.) == 0 and compiled(.1, 0., 1.) == 0
    # Preserve the compiler's exact min/max tie rule in this SDF expression.
    # Residual zero still gives distance zero, irrespective of this half normal.
    np.testing.assert_array_equal(compiled.gradient(.02, 0., 1.), [0., 0., .5])
    return dict(regular_cylinder_chord=error[4], expected_chord=expected,
        cap_locked_distance=constrained[0], unconstrained_distance=unconstrained[0],
        two_branch_exact_tie_distance=joint[0], sample_order_and_cap_exclusion=True,
        finite_cylinder_cap_and_rim=True, periodic_gradient_keeps_physical_side=True)


def branch_rows(np, compiled_terms, records, point):
    values = np.array([term(*point) for term in compiled_terms])
    best = float(values.min())
    active = np.flatnonzero(values == best)
    order = np.argsort(values, kind="stable")
    rows = []
    for index in order[:min(4, len(order))]:
        record = records[index]
        delta = point-np.asarray(record["center"])
        axis = np.asarray(record["axis"])
        axial_position = float(delta@axis)
        radial = float(np.linalg.norm(delta-axial_position*axis)-record["radius"])
        axial = abs(axial_position)-record["height"]/2
        rows.append(dict(index=int(index), value=float(values[index]), radial_component=radial,
            axial_component=axial, component_gap=abs(radial-axial),
            dominant_component="radial" if radial > axial else "axial" if axial > radial else "rim-tie"))
    return dict(exact_minimum_indices=active.tolist(), exact_tie_count=len(active),
        second_branch_gap=float(values[order[1]]-best) if len(order) > 1 else None,
        near_minimum_count_1e_10=int((values-best <= 1e-10).sum()), nearest=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--expected-final-nodes", type=int, required=True)
    parser.add_argument("--expected-final-tets", type=int, required=True)
    args = parser.parse_args()
    for key in ("capture", "package", "run"):
        setattr(args, key, getattr(args, key).resolve())
    capture = json.loads((args.capture/"result.json").read_text())
    assert capture["complete"] and capture["input_hashes_verified"] and capture["status"] == "exact-original-failure-captured"
    assert not capture["native_admitted"]
    assert all(row["failure_reproduced"] and row["expected_stages"] for row in capture["controls"])
    observed = next(row for row in capture["controls"] if row["mode"] == "observed")
    stages = observed["snapshots"]
    assert stages[-1]["stage"] == "final-linear"
    assert (stages[-1]["nodes"], stages[-1]["elements"]) == (args.expected_final_nodes, args.expected_final_tets)
    assert all(row["cap_satisfied"] for row in stages), "Read actual resource cap failures before another diagnostic"
    sys.path.insert(0, str(args.package))
    import numpy as np
    import meshers
    import cad_lattice_fixture as fixture
    import scorecard
    from facet_incidence_audit import audit
    from periodic_trace import periodic_pairs
    from nested_refinement import EDGES

    assert Path(meshers.__file__).resolve().is_relative_to(args.package)
    assert w.digest(Path(meshers._meshers.__file__)) == observed["native_binary_sha256"]
    manifest = json.loads((args.capture/"manifest.json").read_text())
    inputs = dict(manifest["inputs"])
    for row in stages:
        inputs[row["path"]] = row["sha256"]
    for path in (Path(__file__).resolve(), Path(fixture.__file__), Path(scorecard.__file__),
                 args.capture/"result.json", args.capture/"manifest.json"):
        measured = w.digest(path)
        if str(path) in inputs:
            assert inputs[str(path)] == measured, f"Changed frozen helper: {path}"
        inputs[str(path)] = measured
    for path, expected in inputs.items():
        assert w.digest(Path(path)) == expected, f"Changed diagnostic input: {path}"
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"driver.py")
    w.save(args.run/"manifest.json", dict(inputs=inputs, expected_final_nodes=args.expected_final_nodes,
        expected_final_tets=args.expected_final_tets, capture=str(args.capture), numerical_thread_policy=w.NUMERICAL_THREAD_POLICY,
        scope="Read-only rejected-state localization; no mesh admission, remeshing or timing comparison"))
    check = controls(np, meshers, fixture)
    w.save(args.run/"controls.json", check)
    recipe = json.loads((args.capture/"primitive-receipt.json").read_text())["normalized_primitives"]
    compiled = meshers.compile_field(fixture.field(recipe))
    terms = [meshers.compile_field(fixture.field([record])) for record in recipe]
    options = observed["native_options"]
    bounds = np.asarray(options["bounds"]).reshape(3, 2).T
    scale = max(abs(compiled(*(bounds[0]+np.asarray(p)*(bounds[1]-bounds[0])))) for p in product((0., .23, .61, 1.), repeat=3))
    assert compiled.has_gradient and scale > 0
    result = dict(complete=False, native_admitted=False, controls=check, stages=[],
        sampled_error_kind="First-order field residual divided by unlocked gradient norm; not exact CAD distance",
        direct_CAD_distance_performed=False, finite_element_solve_performed=False, comparative_timing_claim=False)
    finite = lambda value: float(value) if np.isfinite(value) else None
    for row in stages:
        print("Localizing", row["stage"], flush=True)
        with np.load(row["path"]) as archive:
            points, cells, surface = [archive[key] for key in ("points", "tetrahedra", "surface")]
        assert len(points) == row["nodes"] and len(cells) == row["elements"]
        samples, locks, faces, kinds = finish_samples(np, points, surface, bounds)
        error, values, gradients, full, free, canonical = distances(np, compiled, samples, locks, bounds, options["periodic"], scale)
        counts = locks.sum(axis=1)
        classifications = []
        for kind in range(7):
            for count in range(4):
                mask = (kinds == kind) & (counts == count)
                if mask.any():
                    classifications.append(dict(kind=KINDS[kind], locked_axes=count, samples=int(mask.sum()),
                        maximum=finite(error[mask].max()), nonfinite=int((~np.isfinite(error[mask])).sum()),
                        above_tolerance=int((error[mask] > options["geometry_tolerance"]).sum())))
        # Preserve every sample; branch attribution is explicitly limited to the
        # worst 128, so it cannot support an all-violations regularity assertion.
        np.savez_compressed(args.run/(row["stage"]+"-samples.npz"), points=samples, locked=locks,
            face_ids=faces, sample_kind=kinds, distance=error, normalized_value=values,
            normalized_gradient=gradients, full_gradient_norm=full, free_gradient_norm=free)
        selected = np.argsort(-error, kind="stable")[:128]
        worst = []
        vertex_branches = {}
        for index in selected:
            nodes = surface[faces[index]]
            for node in nodes:
                if int(node) not in vertex_branches:
                    vertex_branches[int(node)] = branch_rows(np, terms, recipe, points[node])
            worst.append(dict(sample_index=int(index), face=int(faces[index]), face_nodes=surface[faces[index]].tolist(),
                kind=KINDS[kinds[index]], point=samples[index].tolist(), canonical_point=canonical[index].tolist(),
                locked=locks[index].tolist(), distance=finite(error[index]), normalized_value=float(values[index]),
                normalized_gradient=gradients[index].tolist(), full_gradient_norm=float(full[index]),
                free_gradient_norm=float(free[index]), physical_branches=branch_rows(np, terms, recipe, samples[index]),
                canonical_branches=branch_rows(np, terms, recipe, canonical[index]),
                face_vertex_branches=[vertex_branches[int(node)] for node in nodes]))
        w.save(args.run/(row["stage"]+"-worst.json"), worst)
        quality, _, volumes = scorecard.tet_metrics(points, cells)
        facets = audit(points, cells, surface)
        try:
            _, _, periodic = periodic_pairs(points, surface, declared_bounds=bounds)
            periodic["passed"] = True
        except ValueError as failure:
            periodic = dict(passed=False, error=str(failure))
        summary = dict(stage=row["stage"], maximum=finite(error.max()), samples=len(error),
            nonfinite=int((~np.isfinite(error)).sum()),
            above_tolerance=int((error > options["geometry_tolerance"]).sum()), categories=classifications,
            branch_attribution_samples=len(selected), branch_attribution_complete=len(selected) == len(samples),
            branch_attribution_covers_all_violations=int((error > options["geometry_tolerance"]).sum()) <= len(selected),
            nodes=len(points), elements=len(cells), q_min=float(quality.min()), q_passed=bool(quality.min() >= .1),
            minimum_signed_volume=float(volumes.min()), volume=float(volumes.sum()),
            independent_facets=facets, independent_periodicity=periodic, native_admitted=False)
        if row["stage"] == "final-linear":
            unique_edges = np.unique(np.sort(cells[:, EDGES].reshape(-1, 2), axis=1), axis=0)
            quadratic_nodes = len(points)+len(unique_edges)
            summary["quadratic_count_preflight"] = dict(unique_edges=len(unique_edges),
                prospective_nodes=quadratic_nodes, elements=len(cells),
                count_cap_satisfied=max(quadratic_nodes, len(cells)) <= w.CAP,
                scope="Connectivity count only; no elevation, curved Jacobian, geometry or FE claim")
        result["stages"].append(summary)
        w.save(args.run/"result.json", result)
    expected = float(re.search(r"sampled geometry error ([^ ]+)", observed["error"])[1])
    assert result["stages"][-1]["nonfinite"] == 0
    difference = abs(result["stages"][-1]["maximum"]-expected)
    assert difference <= 32*np.spacing(expected), f"finish reproduction differs by {difference}"
    for path, digest in inputs.items():
        assert w.digest(Path(path)) == digest, f"Diagnostic changed input: {path}"
    result.update(complete=True, input_hashes_verified=True, final_native_error=expected,
        error_reproduction_absolute_difference=difference, error_reproduction_tolerance=32*np.spacing(expected),
        status="linear-error-localized")
    w.save(args.run/"result.json", result)
    print("Original finish error reproduced; rejected-state categories archived", flush=True)


if __name__ == "__main__":
    main()
