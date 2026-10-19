"""One charged live Tet10 stress-energy sizing pilot with fixed admission gates."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                  RAYON_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def energy_targets(mesh, fields, tensor):
    import fedoo as fd
    import numpy as np
    from periodic import periodic_pairs
    from stress_recovery import stresses, gram

    points, cells = mesh.points, mesh.tetrahedra
    element = fd.lib_elements.element_list.get_element("tet10")(15)
    inverse = np.linalg.inv(tensor)
    mass, volume = np.zeros(len(cells)), np.zeros(len(cells))
    for derivative, weight in zip(element.shape_function_derivative(element.xi_pg), element.w_pg, strict=True):
        determinant = np.linalg.det(np.einsum("an,tnb->tab", derivative, points[cells]))
        assert determinant.min() > 0 and weight > 0
        stress = stresses(points, cells, fields, derivative)
        mass += weight * determinant * np.einsum("ij,tji->t", inverse, gram(stress))
        volume += weight * determinant
    assert mass.min() > 0 and volume.min() > 0 and abs(mass.sum() - 6) < 1e-7
    selected = np.argsort(mass)[-int(np.ceil(.2 * len(cells))):]
    lengths = np.clip(.5 * np.cbrt(6 * volume[selected]), .01, .05)
    _, classes, _ = periodic_pairs(points, mesh.surface)
    targets = np.full(len(points), .1)
    np.minimum.at(targets, classes[cells[selected].ravel()], np.repeat(lengths, 10))
    targets = targets[classes]
    canonical = np.round((points + .5) % 1, 12) % 1
    unique, grouping = np.unique(canonical, axis=0, return_inverse=True)
    values = np.full(len(unique), .1)
    np.minimum.at(values, grouping, targets)
    assert (values >= .01).all() and (values <= .1).all()
    return unique, values, mass, volume, selected


def size_callback(points, values):
    import numpy as np
    from scipy.spatial import cKDTree

    tree = cKDTree(points, boxsize=1)
    def evaluate(coordinates_x, coordinates_y, coordinates_z):
        query = (np.column_stack((coordinates_x, coordinates_y, coordinates_z)) + .5) % 1
        distance, indices = tree.query(query, k=4)
        weights = 1 / np.maximum(distance, 1e-12)**2
        return np.sum(weights * values[indices], axis=1) / weights.sum(axis=1)
    queries = np.random.default_rng(20261019).uniform(-.5, .5, (4096, 3))
    for axis in range(3):
        low, high = queries.copy(), queries.copy()
        low[:, axis], high[:, axis] = -.5, .5
        np.testing.assert_array_equal(evaluate(*low.T), evaluate(*high.T))
    assert np.isfinite(evaluate(*queries.T)).all() and (evaluate(*queries.T) > 0).all()
    return evaluate


def worker(args):
    started = time.perf_counter()
    np, _, meshers, *_ = workflow.imports()
    from amg_reference import AuditedAMG
    from gradient_tie_probe import symmetric_lattice
    from periodic import periodic_pairs
    from stress_recovery import solve, vector_first_reactions

    output = args.run / args.worker
    row = {"case": "smooth-lattice", "cells": 24, "order": 2, "solver_quadrature_points": 15,
        "min_length": .002, "native_tolerance": 1., "optimization_passes": 1,
        "adaptive": args.worker == "final", "adapt_tolerance": .004,
        "qualified": False, "status": "started", "solver": "cached-elastic-block-amg-cg",
        "field_representation": "Existing domain-equivalent symmetric union with positive void clamp0.02",
        "pilot_policy": "Mark top20percent by all-load stress-energy mass; halve volume length; clip[.01,.05]; all other targets.1"}
    try:
        case, _, _ = workflow.specimen("smooth-lattice")
        def field(coordinates_x, coordinates_y, coordinates_z):
            return np.minimum(symmetric_lattice(coordinates_x, coordinates_y, coordinates_z), .02)
        queries = np.random.default_rng(20261019).uniform(-.5, .5, (4096, 3))
        expected, actual = case.function(*queries.T), field(*queries.T)
        np.testing.assert_array_equal(actual <= 0, expected <= 0)
        np.testing.assert_allclose(actual[expected <= 0], expected[expected <= 0], atol=1e-12, rtol=1e-12)
        size_field = None
        if args.worker == "final":
            coarse = json.loads((args.run / "coarse/result.json").read_text())
            assert coarse["qualified"] and coarse["status"] == "solved-needs-reference"
            targets = np.load(args.run / "coarse/targets.npz")
            size_field = size_callback(targets["points"], targets["values"])
            row["coarse_workflow_s"] = coarse["total_audit_workflow_s"]
        phase = time.perf_counter()
        mesh = meshers.generate(field, bounds=workflow.BOUNDS, cells=24, band=case.band,
            periodic=(True, True, True), threads=1, optimize_passes=1, minimum_quality=.1,
            max_tetrahedra=workflow.CAP, order=2, geometry_tolerance=1.,
            adapt=args.worker == "final", adapt_iterations=4, adapt_tolerance=.004,
            min_length=.002, max_length=.1, gradation=.3, quality_guard=.3, size_field=size_field)
        row.update(generation_s=time.perf_counter() - phase, nodes=len(mesh.points),
                   elements=len(mesh.tetrahedra), diagnostics=mesh.diagnostics)
        assert len(mesh.tetrahedra) <= workflow.CAP
        np.savez_compressed(output / "mesh.npz", points=mesh.points, tetrahedra=mesh.tetrahedra, surface=mesh.surface)
        row["geometry"] = workflow.geometry(output, "smooth-lattice", mesh.points, mesh.tetrahedra, mesh.surface)
        _, _, row["periodic"] = periodic_pairs(mesh.points, mesh.surface)
        row["qualified"] = row["geometry"]["geometry_admitted"]
        row["status"] = "admitted-needs-reference" if row["qualified"] else "geometry-rejected"
        node_cap = 100000 if args.worker == "coarse" else 200000
        row["solve_node_cap"] = node_cap
        workflow.save(output / "result.json", row)
        if row["qualified"] and len(mesh.points) <= node_cap:
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
            if args.worker == "coarse":
                phase = time.perf_counter()
                locations, values, mass, volume, selected = energy_targets(mesh, fields, tensor)
                size_callback(locations, values)
                np.savez_compressed(output / "targets.npz", points=locations, values=values,
                                    mass=mass, volume=volume, selected_cells=selected)
                row.update(indicator_and_target_s=time.perf_counter() - phase,
                    selected_elements=len(selected), target_quantiles=np.quantile(values, [0, .1, .5, .9, 1]).tolist())
        elif row["qualified"]:
            row["status"] = "solve-node-cap"
    except Exception as error:
        row.update(qualified=False, status="failed", error=str(error), traceback=traceback.format_exc())
    row["stage_workflow_s"] = time.perf_counter() - started
    row["total_audit_workflow_s"] = row["stage_workflow_s"] + row.get("coarse_workflow_s", 0.)
    workflow.save(output / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--worker", choices=("coarse", "final"))
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    screen = json.loads((args.screen / "result.json").read_text())
    assert screen["predeclared_energy_baseline_screen_passed"]
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    import meshers._meshers
    paths = [Path(__file__), args.screen / "manifest.json", args.screen / "result.json",
        Path(meshers._meshers.__file__), workflow.ROOT / "crates/meshers-core/src/tet_adapt.rs"]
    paths.extend(workflow.HERE / name for name in ("amg_reference.py", "gradient_tie_probe.py",
                 "stress_recovery.py", "tet10_integration.py", "microgen_workflow.py"))
    source = workflow.HERE / "runs/g3-microgen-defaults-v4"
    paths.extend(source / ("smooth-lattice" + suffix) for suffix in ("-samples.npz", "-specimen.json"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    for stage in ("coarse", "final"):
        (args.run / stage).mkdir()
        for suffix in ("-samples.npz", "-specimen.json"):
            shutil.copyfile(source / ("smooth-lattice" + suffix), args.run / stage / ("smooth-lattice" + suffix))
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeouts_s": {"coarse": 600, "final": 1200},
        "tetrahedra_cap": workflow.CAP, "coarse_node_cap": 100000, "final_node_cap": 200000,
        "pyamg": importlib.metadata.version("pyamg"), "selected_fraction": .2,
        "length_multiplier": .5, "target_clip": [.01, .05], "unmarked_target": .1,
        "oracle_input_to_live_indicator": False, "coarse_workflow_charged": True})
    started = time.perf_counter()
    with Job(parts={"coarse": 1, "final": 2}) as job:
        for stage, timeout in (("coarse", 600), ("final", 1200)):
            with job.part(stage) as part:
                record = workflow.bounded([sys.executable, __file__, *sys.argv[1:], "--worker", stage],
                    args.run / stage / "worker.log", part, timeout=timeout)
                workflow.save(args.run / stage / "worker.json", record)
                assert record["returncode"] == 0, f"Read {stage}/worker.log"
                row = json.loads((args.run / stage / "result.json").read_text())
                job.result(stage, nodes=row.get("nodes", 0), qualified=int(row["qualified"]),
                           seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
                if stage == "coarse":
                    assert row["status"] == "solved-needs-reference", "Coarse stage rejects; do not run final"
    workflow.save(args.run / "summary.json", {"controller_wall_s": time.perf_counter() - started,
        "final_status": row["status"], "qualified": row["qualified"],
        "charged_live_workflow_s": row["total_audit_workflow_s"],
        "scope": "One live non-recovery energy-sizing pilot; six coarse loads and field construction charged; reference analysis and certification excluded"})
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
