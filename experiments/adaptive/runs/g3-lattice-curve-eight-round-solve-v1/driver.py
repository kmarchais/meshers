"""Bounded geometry-only curve-grading controls, or one charged live candidate."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

CURVE_LENGTH = .0025
SLOPE = .8


def curve_callback(mesh, function):
    import numpy as np
    from periodic_crease_audit import features, project_box_curve
    from scipy.spatial import cKDTree

    data = features({"points": mesh.points, "tetrahedra": mesh.tetrahedra, "surface": mesh.surface})
    clouds = []
    for edge, midpoint in zip(data["curve_edges"], data["curve_midpoints"], strict=True):
        first, second = mesh.points[list(edge)]
        middle = mesh.points[midpoint]
        speed = max(np.linalg.norm(-3 * first - second + 4 * middle),
                    np.linalg.norm(first + 3 * second - 4 * middle))
        divisions = max(2, int(np.ceil(speed / (CURVE_LENGTH / 8))))
        parameter = np.linspace(0, 1, divisions + 1)[:, None]
        cloud = ((1 - parameter) * (1 - 2 * parameter) * first +
                 parameter * (2 * parameter - 1) * second + 4 * parameter * (1 - parameter) * middle)
        clouds.append(cloud)
    counts = [len(cloud) for cloud in clouds]
    projected, projection = project_box_curve(function, np.concatenate(clouds))
    cursor, maximum_spacing = 0, 0.
    for count in counts:
        maximum_spacing = max(maximum_spacing,
            float(np.linalg.norm(np.diff(projected[cursor:cursor + count], axis=0), axis=1).max()))
        cursor += count
    assert maximum_spacing <= CURVE_LENGTH / 4
    canonical = np.unique(np.round((projected + .5) % 1, 12) % 1, axis=0)
    tree = cKDTree(canonical, boxsize=1)
    def evaluate(coordinates_x, coordinates_y, coordinates_z):
        query = (np.column_stack((coordinates_x, coordinates_y, coordinates_z)) + .5) % 1
        distance, _ = tree.query(query)
        return np.minimum(CURVE_LENGTH + SLOPE * distance, .1)
    queries = np.random.default_rng(20261019).uniform(-.5, .5, (4096, 3))
    values = evaluate(*queries.T)
    assert np.isfinite(values).all() and (values >= CURVE_LENGTH).all() and (values <= .1).all()
    for axis in range(3):
        low, high = queries.copy(), queries.copy()
        low[:, axis], high[:, axis] = -.5, .5
        np.testing.assert_array_equal(evaluate(*low.T), evaluate(*high.T))
    return evaluate, canonical, {"coarse_curve_edges": len(data["curve_edges"]),
        "projected_samples": len(projected), "periodic_unique_samples": len(canonical),
        "maximum_projected_sample_spacing": maximum_spacing, **projection}


def worker(args):
    started = time.perf_counter()
    np, _, meshers, *_ = workflow.imports()
    from periodic import periodic_pairs

    output = args.run / args.worker
    row = {"case": "smooth-lattice", "cells": 24, "order": 2, "solver_quadrature_points": 15,
        "adaptive": True, "adapt_iterations": args.adapt_iterations,
        "adapt_tolerance": .004, "native_tolerance": 1., "min_length": .002,
        "optimization_passes": 1, "gradation": SLOPE, "qualified": False, "status": "started",
        "route": args.worker, "solve_node_cap": 250000,
        "solver": "cached-elastic-block-amg-cg" if args.solve else None,
        "field_representation": "Raw smooth-union values with existing positive void clamp0.02; corrected compiler gradient"}
    try:
        case, _, _ = workflow.specimen("smooth-lattice")
        def field(coordinates_x, coordinates_y, coordinates_z):
            return np.minimum(case.function(coordinates_x, coordinates_y, coordinates_z), .02)
        queries = np.random.default_rng(20261019).uniform(-.5, .5, (4096, 3))
        expected, actual = case.function(*queries.T), field(*queries.T)
        np.testing.assert_array_equal(actual <= 0, expected <= 0)
        np.testing.assert_array_equal(actual[expected <= 0], expected[expected <= 0])
        common = dict(bounds=workflow.BOUNDS, cells=24, band=case.band, periodic=(True, True, True),
                      threads=1, optimize_passes=1, minimum_quality=.1, geometry_tolerance=1.,
                      max_tetrahedra=workflow.CAP, order=2)
        size_field = None
        if args.worker == "curve":
            phase = time.perf_counter()
            coarse = meshers.generate(field, adapt=False, **common)
            assert len(coarse.tetrahedra) <= workflow.CAP and len(coarse.points) <= 100000
            row["coarse_generation_s"] = time.perf_counter() - phase
            row.update(coarse_nodes=len(coarse.points), coarse_elements=len(coarse.tetrahedra))
            np.savez_compressed(output / "coarse-mesh.npz", points=coarse.points,
                                tetrahedra=coarse.tetrahedra, surface=coarse.surface)
            phase = time.perf_counter()
            row["coarse_geometry"] = workflow.geometry(output, "smooth-lattice",
                coarse.points, coarse.tetrahedra, coarse.surface)
            assert row["coarse_geometry"]["geometry_admitted"]
            row["coarse_audit_s"] = time.perf_counter() - phase
            phase = time.perf_counter()
            size_field, cloud, metadata = curve_callback(coarse, case.function)
            row["feature_preparation_s"] = time.perf_counter() - phase
            row["curve_feature"] = metadata
            np.savez_compressed(output / "curve-samples.npz", points=cloud)
        phase = time.perf_counter()
        mesh = meshers.generate(field, adapt=True, adapt_iterations=args.adapt_iterations, adapt_tolerance=.004,
            min_length=.002, max_length=.1, gradation=SLOPE, quality_guard=.3, size_field=size_field, **common)
        row.update(generation_s=time.perf_counter() - phase, nodes=len(mesh.points),
                   elements=len(mesh.tetrahedra), diagnostics=mesh.diagnostics)
        assert len(mesh.tetrahedra) <= workflow.CAP
        if size_field is not None:
            from periodic_crease_audit import features

            data = features({"points": mesh.points, "tetrahedra": mesh.tetrahedra, "surface": mesh.surface})
            locations = mesh.points[np.asarray(data["curve_edges"])]
            targets = size_field(*locations.reshape(-1, 3).T).reshape(-1, 2).mean(axis=1)
            lengths = np.linalg.norm(locations[:, 1] - locations[:, 0], axis=1)
            ratios = lengths / targets
            row["curve_size_realization"] = {
                "curve_edges": len(lengths),
                "edge_length_quantiles": np.quantile(lengths, [0, .1, .5, .9, 1]).tolist(),
                "length_target_ratio_quantiles": np.quantile(ratios, [0, .1, .5, .9, 1]).tolist(),
                "edges_above_split_ratio": int(np.count_nonzero(ratios > 4 / 3 + 1e-8)),
                "fraction_above_split_ratio": float(np.mean(ratios > 4 / 3 + 1e-8)),
                "scope": "Requested endpoint targets, not proof of native split admissibility"}
        np.savez_compressed(output / "mesh.npz", points=mesh.points, tetrahedra=mesh.tetrahedra, surface=mesh.surface)
        row["geometry"] = workflow.geometry(output, "smooth-lattice", mesh.points, mesh.tetrahedra, mesh.surface)
        _, _, row["periodic"] = periodic_pairs(mesh.points, mesh.surface)
        row["qualified"] = row["geometry"]["geometry_admitted"]
        row["status"] = "admitted-needs-reference" if row["qualified"] else "geometry-rejected"
        workflow.save(output / "result.json", row)
        if row["qualified"] and len(mesh.points) > row["solve_node_cap"]:
            row["status"] = "solve-node-cap"
        elif row["qualified"] and args.solve:
            from amg_reference import AuditedAMG
            from stress_recovery import solve, vector_first_reactions

            phase = time.perf_counter()
            def checkpoint(completed, columns, records):
                row.update(status="solving", completed_loads=completed,
                           partial_stiffness_columns=columns.tolist(), solver_records=records)
                workflow.save(output / "result.json", row)
            def save_field(completed, displacement):
                np.savez_compressed(output / f"load-{completed}.npz", displacement=displacement)
            fields, tensor, records = solve(mesh.points, mesh.tetrahedra, n_elm_gp=15,
                solver=AuditedAMG(elastic_blocks=True), checkpoint=checkpoint,
                field_checkpoint=save_field, reaction_reader=vector_first_reactions, release_reduced_cache=True)
            assert len(records) == 6 and max(record["relative_residual"] for record in records) < 2e-9
            row.pop("partial_stiffness_columns", None)
            row.update(solve_s=time.perf_counter() - phase, stiffness=tensor.tolist(),
                       solver_records=records, status="solved-needs-reference")
            np.savez_compressed(output / "fields.npz", displacement=fields)
    except Exception as error:
        row.update(qualified=False, status="failed", error=str(error), traceback=traceback.format_exc())
    row["total_audit_workflow_s"] = time.perf_counter() - started
    workflow.save(output / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--route", choices=("both", "baseline", "curve"), default="both")
    parser.add_argument("--solve", action="store_true")
    parser.add_argument("--adapt-iterations", type=int, choices=(4, 8), default=4)
    parser.add_argument("--worker", choices=("baseline", "curve"))
    args = parser.parse_args()
    if args.solve and args.route != "curve":
        parser.error("Select the single measured curve candidate before any live solve")
    if args.worker:
        worker(args)
        return
    screen = json.loads((args.screen / "result.json").read_text())
    assert screen["predeclared_curve_localization_screen_passed"]
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    import meshers._meshers
    import numpy as np
    from periodic_crease_audit import features

    parent = Path(screen["parent"])
    checked = features(np.load(parent / "mesh.npz"))
    assert len(checked["curve_edges"]) == screen["coarse_curve_edges"]
    assert int(checked["groups"]["curve_vertex_incident"].sum()) == screen["rules"][0]["groups"]["curve_vertex_incident"]["elements"]
    paths = [Path(__file__), args.screen / "manifest.json", args.screen / "result.json",
        parent / "mesh.npz", Path(meshers._meshers.__file__),
        workflow.ROOT / "crates/meshers-core/src/tet_adapt.rs",
        workflow.ROOT / "crates/meshers-python/python/meshers/expression.py"]
    paths.extend(workflow.HERE / name for name in ("periodic_crease_audit.py", "microgen_workflow.py",
                 "amg_reference.py", "stress_recovery.py", "tet10_integration.py"))
    source = workflow.HERE / "runs/g3-microgen-defaults-v4"
    paths.extend(source / ("smooth-lattice" + suffix) for suffix in ("-samples.npz", "-specimen.json"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 1200 if args.solve else 300, "tetrahedra_cap": workflow.CAP,
        "solve_node_cap": 250000, "adapt_iterations": args.adapt_iterations,
        "curve_length": CURVE_LENGTH, "grading_slope": SLOPE,
        "metric": "min(0.1, curve_length + grading_slope * torus_distance_to_projected_curve_samples)",
        "coarse_solve_count": 0, "live_feature_preparation_charged": True, "geometry_only": not args.solve,
        "policy_note": "One fixed curve-local metric with matching native gradation; no physical crease smoothing"})
    routes = ("baseline", "curve") if args.route == "both" else (args.route,)
    for route in routes:
        (args.run / route).mkdir()
        for suffix in ("-samples.npz", "-specimen.json"):
            shutil.copyfile(source / ("smooth-lattice" + suffix), args.run / route / ("smooth-lattice" + suffix))
    with Job(parts=list(routes)) as job:
        for route in routes:
            with job.part(route) as part:
                command = [sys.executable, __file__, *sys.argv[1:], "--worker", route]
                record = workflow.bounded(command, args.run / route / "worker.log", part,
                                          timeout=1200 if args.solve else 300)
                workflow.save(args.run / route / "worker.json", record)
                assert record["returncode"] == 0, f"Read {route}/worker.log"
                row = json.loads((args.run / route / "result.json").read_text())
                job.result(route, qualified=int(row["qualified"]), nodes=row.get("nodes", 0),
                           seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
