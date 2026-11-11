"""One native curved Tet10 control on the frozen crossed-cylinder specimen.

The straight elevation is an exact geometry-preserving control. Native order2
uses the saved Tet4 recipe with only element order changed; curved_fit stays
off. Original corner/connectivity equality is measured, not assumed. No FEA
or matched-order comparison against the saved actual Tet4 mesh is claimed.
"""

import argparse
from functools import partial
import json
from pathlib import Path
import shutil
import sys
import time
import traceback
from types import SimpleNamespace

import cad_intersection_probe as probe
import cad_union_comparison as base
import microgen_workflow as w


def archive(run, name, arrays, part, *, native_pairs=None):
    import numpy as np
    from comparison_mesh_audit import observe
    from periodic_trace import periodic_pairs

    points, cells, surface = arrays
    output = run / name
    output.mkdir()
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(run / (base.NAME+suffix), output / (base.NAME+suffix))
    np.savez_compressed(output / "mesh.npz", points=points, tetrahedra=cells, surface=surface)
    result = dict(nodes=len(points), elements=len(cells), order=2, qualified=False,
        mesh_file=str(output / "mesh.npz"), mesh_sha256=w.digest(output / "mesh.npz"),
        status="archived-before-admission")
    w.save(output / "result.json", result)
    if max(len(points), len(cells)) > w.CAP:
        result.update(status="cap-rejected", error="Node or tetrahedron cap exceeded")
        w.save(output / "result.json", result)
        return result
    assert cells.shape[1] == 10 and surface.shape[1] == 6
    check = partial(periodic_pairs, declared_bounds=([-.5]*3, [.5]*3))
    try:
        part.update(.1, "Auditing full Tet10 geometry and Tri6 periodic caps", nodes=len(points), tets=len(cells))
        result["geometry"], result["periodic"], result["quality_quantiles"] = observe(
            output, base.NAME, base.target(), points, cells, surface, True, periodic_auditor=check)
        result["qualified"] = result["geometry"]["geometry_admitted"]
        result["jacobian_scope"] = ("Exact affine reduction to constant positive determinant" if
            result["geometry"]["exact_straight_reduction_verified"] else
            "Positive125 Duffy probes plus element nodes, barycenter and 15-point integration rule")
        if native_pairs is not None:
            checked, _, _ = check(points, surface)
            identical = all({tuple(map(int, row)) for row in supplied} ==
                            {tuple(map(int, row)) for row in independent}
                            for supplied, independent in zip(native_pairs, checked, strict=True))
            result["native_and_independent_pairs_agree"] = identical
            result["qualified"] &= identical
        result["status"] = "admitted" if result["qualified"] else "admission-rejected"
    except Exception as error:
        result.update(status="observer-failed", qualified=False, error=str(error), traceback=traceback.format_exc())
    part.update(1, "Finished full-order admission", admitted=result["qualified"])
    w.save(output / "result.json", result)
    return result


def localization(source, arrays, part):
    import numpy as np
    import scorecard

    points, cells, surface = arrays
    wall = probe.wall_mask(points, surface)
    groups = probe.face_groups(points, surface)
    weights = w.tri_values(scorecard.BARY, 2)
    positions = np.einsum("sn,fnd->fsd", weights, points[surface[wall]])
    feet, converged = scorecard.project(base.target(), positions.reshape(-1, 3))
    forward = np.linalg.norm(feet-positions.reshape(-1, 3), axis=1).reshape(positions.shape[:2])
    good = converged.reshape(positions.shape[:2])
    description = {}
    for name in ("smooth", "crossing", "ridge_touch"):
        selected = groups[name][wall]
        description[name] = dict(faces=int(selected.sum()),
            diagnostic_forward_max=float(forward[selected].max(initial=0)),
            converged_samples=int(good[selected].sum()), sample_count=int(good[selected].size))
    part.update(.25, "Localized quadratic wall patches", wall_faces=int(wall.sum()))
    with np.load(source / (base.NAME+"-samples.npz")) as samples:
        names = ("walls", "caps", "ridges", "rims")
        queries = np.concatenate([samples[name] for name in names])
        vertices, faces = w.tessellate(points, surface, 16)
        distances = scorecard.distance_to_mesh(vertices, faces, queries)
        reverse = {}
        offset = 0
        for name in names:
            count = len(samples[name])
            values = distances[offset:offset+count]
            worst = int(values.argmax())
            reverse[name] = dict(count=count, maximum=float(values.max()), p99=float(np.quantile(values, .99)),
                worst_point=samples[name][worst].tolist())
            offset += count
    part.update(1, "Observed frozen ridge, rim and wall samples", ridge_error=reverse["ridges"]["maximum"])
    return dict(face_groups=description, reverse_groups=reverse, reverse_subdivisions=16,
        projection_scope="Unconstrained Newton diagnostic by original corner branch labels; shared admission is separate")


def worker(args):
    np, _, meshers, *_ = w.imports()
    from microgen_tet10_control import straight_elevation
    from nested_refinement import EDGES
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    source = args.source
    saved = json.loads((source / "result.json").read_text())
    assert saved["qualified"] and w.digest(source / "mesh.npz") == saved["mesh_sha256"]
    with np.load(source / "mesh.npz") as original:
        old_points, old_cells, old_surface = (original[key] for key in ("points", "tetrahedra", "surface"))
    paths = {Path(__file__), Path(probe.__file__), Path(base.__file__), Path(base.fixture.__file__), Path(w.__file__),
             Path(meshers.__file__), Path(meshers._meshers.__file__)}
    paths.update(w.HERE/name for name in ("microgen_tet10_control.py", "scorecard.py", "periodic_trace.py",
        "comparison_mesh_audit.py", "quadratic_edge_removal_probe.py", "quadrature_audit.py",
        "facet_incidence_audit.py", "nested_refinement.py", "box_projection.py"))
    paths.update(source/name for name in ("result.json", "mesh.npz", base.NAME+"-specimen.json", base.NAME+"-samples.npz"))
    hashes = {str(path.resolve()): w.digest(path) for path in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(source/(base.NAME+suffix), args.run/(base.NAME+suffix))
    options = dict(base.NATIVE, order=2, curved_fit=False)
    w.save(args.run/"manifest.json", dict(inputs=hashes, native_options=options, worker_timeout_s=240,
        memory_bytes=w.MEMORY, numerical_thread_policy=w.NUMERICAL_THREAD_POLICY,
        node_cap=w.CAP, tetrahedron_cap=w.CAP, native_binary_sha256=w.digest(Path(meshers._meshers.__file__)),
        hypothesis="Quadratic wall geometry removes smooth-cylinder chord error at unchanged linear topology",
        scope="One fresh native order2 pilot and one exact affine elevation control. Not a matched-order actual comparison, repeated timing study or FEA result"))
    result = dict(inputs=hashes, source=str(source), complete=False, original_specimen_and_samples_unchanged=True,
        node_cap=w.CAP, tetrahedron_cap=w.CAP, baseline_nodes=len(old_points), baseline_elements=len(old_cells),
        native_options=options, solve_performed=False, matched_order_actual_comparison=False)
    with Job(parts={"straight-elevation": 1, "native-generation": 1, "native-admission": 3, "localization": 2}) as job:
        with job.part("straight-elevation") as part:
            part.update(0, "Counting edges before exact straight elevation", tets=len(old_cells))
            straight = straight_elevation(old_points, old_cells, old_surface, node_cap=w.CAP)
            result["straight_node_count"] = straight["nodes"]
            if not straight["within_node_cap"]:
                result.update(complete=True, status="straight-node-cap-rejected")
                w.save(args.run/"result.json", result)
                return
            arrays = tuple(straight[key] for key in ("points", "tetrahedra", "surface"))
            result["straight"] = archive(args.run, "straight-tet10", arrays, part)
            assert result["straight"]["qualified"]
            assert result["straight"]["geometry"]["exact_straight_reduction_verified"]
            for key in ("sampled_two_sided_max", "relative_volume_error", "q_min_corner"):
                np.testing.assert_allclose(result["straight"]["geometry"][key], saved["geometry"][key], atol=1e-12, rtol=1e-10)
        job.result("straight-elevation", nodes=straight["nodes"], exact_domain_preserved=True)
        w.save(args.run/"result.json", result)
        with job.part("native-generation") as part:
            part.update(0, "Generating one native curved Tet10 mesh", cells=options["cells"])
            started = time.perf_counter()
            try:
                mesh = meshers.generate(base.fixture.periodic_field, **options)
                result["native_generation_s"] = time.perf_counter()-started
                points, cells, surface = (np.asarray(getattr(mesh, name)) for name in ("points", "tetrahedra", "surface"))
                # Archive immediately, even if the subsequent checks reject it.
                np.savez_compressed(args.run/"native-returned.npz", points=points, tetrahedra=cells, surface=surface)
                result["native_returned_sha256"] = w.digest(args.run/"native-returned.npz")
                result["native_diagnostics"] = mesh.diagnostics
                assert cells.shape[1] == 10 and surface.shape[1] == 6
                equal = dict(corners_bitwise_equal=bool(np.array_equal(points[:len(old_points)], old_points)),
                    tetrahedron_corners_bitwise_equal=bool(np.array_equal(cells[:, :4], old_cells)),
                    surface_corners_bitwise_equal=bool(np.array_equal(surface[:, :3], old_surface)))
                result["linear_topology_reproduced"] = equal
                result["same_linear_topology"] = all(equal.values())
                displacement = points[cells[:, 4:]]-points[cells[:, EDGES]].mean(axis=2)
                result["maximum_quadratic_edge_displacement"] = float(np.linalg.norm(displacement, axis=2).max())
                part.update(1, "Native Tet10 returned and archived", nodes=len(points), tets=len(cells))
            except Exception as error:
                result.update(complete=True, status="native-failed", error_type=type(error).__name__,
                    error=str(error), traceback=traceback.format_exc(), native_generation_s=time.perf_counter()-started)
                part.update(1, "Native generation or order check failed", failed=True)
                job.result("native-generation", failed=True, error_type=type(error).__name__)
                w.save(args.run/"result.json", result)
                assert all(w.digest(Path(p)) == h for p, h in hashes.items())
                return
        job.result("native-generation", seconds=result["native_generation_s"], nodes=len(points),
                   same_linear_topology=result["same_linear_topology"])
        w.save(args.run/"result.json", result)
        with job.part("native-admission") as part:
            result["native"] = archive(args.run, "native-curved-tet10", (points, cells, surface), part,
                                       native_pairs=mesh.periodic_pairs)
        job.result("native-admission", admitted=result["native"]["qualified"],
                   minimum_jacobian=result["native"].get("geometry", {}).get("minimum_sampled_jacobian"))
        w.save(args.run/"result.json", result)
        with job.part("localization") as part:
            if result["native"]["qualified"]:
                part.update(0, "Localizing the native curved boundary", nodes=len(points), tets=len(cells))
                result["localization"] = localization(source, (points, cells, surface), part)
            else:
                part.update(1, "Rejected native mesh retained", admitted=False)
        job.result("localization", ridge_error=result.get("localization", {}).get("reverse_groups", {}).get("ridges", {}).get("maximum"))
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    result.update(complete=True, status="control-complete", input_hashes_verified=True)
    w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=probe.NATIVE)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "source"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=240)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    row = json.loads((args.run/"result.json").read_text())
    native = row.get("native", {})
    geometry = native.get("geometry", {})
    result("native-curved-tet10", status=row["status"], admitted=native.get("qualified", False),
        nodes=native.get("nodes"), same_linear_topology=row.get("same_linear_topology"),
        boundary=geometry.get("sampled_two_sided_max"), volume_error=geometry.get("relative_volume_error"),
        q_min=geometry.get("q_min_corner"), minimum_jacobian=geometry.get("minimum_sampled_jacobian"))


if __name__ == "__main__":
    main()
