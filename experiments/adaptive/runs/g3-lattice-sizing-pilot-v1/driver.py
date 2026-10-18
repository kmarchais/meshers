"""Build and test spatial adaptation, then measure one two-arm sizing pilot."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

import microgen_workflow as workflow

ROOT = workflow.ROOT
os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
    RAYON_NUM_THREADS="1", CARGO_BUILD_JOBS="2", RUSTUP_TOOLCHAIN="1.95.0",
    VIRTUAL_ENV=str(ROOT / ".venv"), PATH=str(Path.home() / ".cargo/bin") + os.pathsep + os.environ["PATH"])
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def size_from_recovery(mesh, fields, tensor, nominal):
    import numpy as np
    import fedoo as fd
    from scipy.spatial import cKDTree
    from periodic import periodic_pairs
    from stress_recovery import stresses, recover, gram

    points, tetrahedra = mesh.points, mesh.tetrahedra
    _, classes, _ = periodic_pairs(points, mesh.surface)
    volumes = np.linalg.det(points[tetrahedra[:, 1:]] - points[tetrahedra[:, :1]]) / 6
    element = fd.lib_elements.element_list.get_element("tet4")(4)
    stress = stresses(points, tetrahedra, fields, element.shape_function_derivative([[0.25] * 3])[0])
    recovered = recover(stress, volumes, tetrahedra, classes)
    mass = np.zeros(len(tetrahedra))
    inverse = np.linalg.inv(tensor)
    for coordinate, weight in zip(element.xi_pg, element.w_pg, strict=True):
        interpolation = element.shape_function(coordinate[None])[0]
        difference = np.einsum("n,tnlij->tlij", interpolation, recovered[tetrahedra]) - stress
        mass += 6 * weight * volumes * np.einsum("ij,tji->t", inverse, gram(difference))
    length = np.cbrt(6 * volumes)
    targets = np.maximum(mass / (volumes * length**2), 1e-30)**(-1 / 5)
    nodal = np.full(len(points), np.inf)
    np.minimum.at(nodal, classes[tetrahedra.ravel()], np.repeat(targets, 4))
    nodal = nodal[classes]
    nodal = np.clip(nominal * nodal / np.median(nodal), 0.03, 0.1)
    canonical = np.round((points + 0.5) % 1, 12) % 1
    unique, grouping = np.unique(canonical, axis=0, return_inverse=True)
    values = np.full(len(unique), np.inf)
    np.minimum.at(values, grouping, nodal)
    tree = cKDTree(unique, boxsize=1)
    def evaluate(coordinates_x, coordinates_y, coordinates_z):
        query = (np.column_stack((coordinates_x, coordinates_y, coordinates_z)) + 0.5) % 1
        distance, indices = tree.query(query, k=4)
        weights = 1 / np.maximum(distance, 1e-12)**2
        return np.sum(weights * values[indices], axis=1) / weights.sum(axis=1)
    return evaluate, nodal, mass


def candidate(args):
    np, _, meshers, _, _, _ = workflow.imports()
    from periodic import periodic_pairs, tensor_difference
    from stress_recovery import solve

    row = {"case": args.case, "route": args.worker, "nominal_length": args.nominal,
           "cells": args.cells, "coarse_cells": args.coarse_cells, "min_length": args.min_length,
           "order": 1, "solver_quadrature_points": 4,
           "status": "started", "qualified": False}
    started = time.perf_counter()
    try:
        case, _, _ = workflow.specimen(args.case)
        field = case.function
        if args.case == "smooth-lattice":
            from gradient_tie_probe import symmetric_lattice
            queries = np.random.default_rng(20261018).uniform(-0.5, 0.5, (4096, 3))
            expected = case.function(*queries.T)
            symmetric = symmetric_lattice(*queries.T)
            np.testing.assert_allclose(symmetric, expected, rtol=1e-12, atol=1e-12)
            def field(coordinates_x, coordinates_y, coordinates_z):
                return np.minimum(symmetric_lattice(coordinates_x, coordinates_y, coordinates_z), 0.02)
            actual = field(*queries.T)
            np.testing.assert_array_equal(actual < 0, expected < 0)
            np.testing.assert_allclose(actual[expected <= 0], expected[expected <= 0], rtol=1e-12, atol=1e-12)
            row["lattice_field_representation"] = {"symmetric_minimum": True, "positive_void_ceiling": 0.02,
                "verification_points": len(queries), "material_mismatches": 0,
                "maximum_symmetric_field_difference": float(abs(symmetric - expected).max())}
        common = dict(bounds=workflow.BOUNDS, periodic=(True, True, True), band=case.band,
                      threads=1, optimize_passes=1, minimum_quality=0.1, max_tetrahedra=workflow.CAP)
        if args.worker == "recovery":
            phase = time.perf_counter()
            coarse = meshers.generate(field, cells=args.coarse_cells, geometry_tolerance=1.0, **common)
            row["coarse_generation_s"] = time.perf_counter() - phase
            phase = time.perf_counter()
            fields, tensor, records = solve(coarse.points, coarse.tetrahedra)
            row["coarse_solve_s"] = time.perf_counter() - phase
            row["coarse_solver_records"] = records
            phase = time.perf_counter()
            size_field, lengths, mass = size_from_recovery(coarse, fields, tensor, args.nominal)
            row["recovery_s"] = time.perf_counter() - phase
            row["target_quantiles"] = np.quantile(lengths, [0, 0.1, 0.5, 0.9, 1]).tolist()
            np.savez_compressed(args.run / "coarse-recovery.npz", points=coarse.points,
                tetrahedra=coarse.tetrahedra, displacement=fields, target_length=lengths, estimated_mass=mass)
        else:
            size_field = lambda coordinates_x, *_: np.full_like(coordinates_x, args.nominal)
        phase = time.perf_counter()
        mesh = meshers.generate(field, cells=args.cells, geometry_tolerance=0.002,
            adapt=True, min_length=args.min_length, max_length=0.1, gradation=0.3, adapt_iterations=4,
            quality_guard=0.3, size_field=size_field, **common)
        row.update(generation_s=time.perf_counter() - phase, nodes=len(mesh.points),
                   elements=len(mesh.tetrahedra), native_diagnostics=mesh.diagnostics)
        np.savez_compressed(args.run / (args.worker + "-mesh.npz"), points=mesh.points,
                            tetrahedra=mesh.tetrahedra, surface=mesh.surface)
        phase = time.perf_counter()
        row["geometry"] = workflow.geometry(args.run, args.case, mesh.points, mesh.tetrahedra, mesh.surface)
        _, _, row["periodic"] = periodic_pairs(mesh.points, mesh.surface)
        row["audit_s"] = time.perf_counter() - phase
        row["qualified"] = row["geometry"]["geometry_admitted"]
        if not row["qualified"]:
            row["status"] = "geometry-rejected"
        elif len(mesh.points) > 100_000:
            row["status"] = "solve-node-cap"
        else:
            phase = time.perf_counter()
            fields, tensor, records = solve(mesh.points, mesh.tetrahedra)
            row.update(solve_s=time.perf_counter() - phase, stiffness=tensor.tolist(), solver_records=records)
            reference = json.loads((args.source / "summary.json").read_text())["references"].get(args.case)
            row["status"] = "solved-needs-reference"
            if reference is not None:
                assert reference["qualified_inputs"]
                row["reference_disagreement"] = tensor_difference(tensor, np.array(reference["stiffness"]))
                row["status"] = "solved"
            np.savez_compressed(args.run / (args.worker + "-fields.npz"), displacement=fields)
    except Exception as error:
        row.update(status="failed", qualified=False, error=str(error), traceback=traceback.format_exc())
    row["total_workflow_s"] = time.perf_counter() - started
    workflow.save(args.run / (args.worker + ".json"), row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=workflow.HERE / "runs/g3-microgen-defaults-v4")
    parser.add_argument("--case", default="schwarz-solid")
    parser.add_argument("--nominal", type=float, default=0.065)
    parser.add_argument("--cells", type=int, default=16)
    parser.add_argument("--coarse-cells", type=int, default=12)
    parser.add_argument("--min-length", type=float, default=0.02)
    parser.add_argument("--worker-timeout", type=int, choices=(600, 1200), default=1200)
    parser.add_argument("--worker", choices=("uniform", "recovery"))
    parser.add_argument("--skip-build", action="store_true")
    args = parser.parse_args()
    if args.worker:
        candidate(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    inputs = [args.source / (args.case + suffix) for suffix in ("-samples.npz", "-specimen.json")]
    inputs.append(args.source / "summary.json")
    hashes = {str(path): workflow.digest(path) for path in inputs}
    for path in inputs[:2]:
        shutil.copyfile(path, args.run / path.name)
    source_paths = [ROOT / path for path in ("crates/meshers-core/src/geometry.rs",
        "crates/meshers-core/src/implicit.rs", "crates/meshers-core/src/tet_adapt.rs",
        "crates/meshers-python/src/lib.rs", "crates/meshers-python/python/meshers/__init__.py",
        "crates/meshers-python/tests/test_size_field.py", "experiments/adaptive/stress_recovery.py",
        "experiments/adaptive/gradient_tie_probe.py", "experiments/adaptive/microgen_workflow.py",
        "experiments/adaptive/sizing_pilot.py")]
    manifest = {"script_sha256": workflow.digest(Path(__file__)), "inputs": hashes,
        "source_sha256": {str(path): workflow.digest(path) for path in source_paths},
        "memory_bytes": workflow.MEMORY, "threads": 1, "rust": "1.95.0", "build_jobs": 2,
        "nominal": args.nominal, "case": args.case, "cells": args.cells,
        "coarse_cells": args.coarse_cells, "min_length": args.min_length,
        "worker_timeout_s": args.worker_timeout, "tetrahedra_cap": workflow.CAP,
        "solve_node_cap": 100_000}
    workflow.save(args.run / "manifest.json", manifest)
    tasks = [] if args.skip_build else [
        ("build", [str(ROOT / ".venv/bin/maturin"), "develop", "--release", "--locked", "--manifest-path",
                   str(ROOT / "crates/meshers-python/Cargo.toml"), "--features", "experimental-surfaces"], 1200),
        ("core-tests", ["cargo", "test", "--release", "--locked", "-p", "meshers-core", "--lib", "--", "--test-threads=1"], 900),
        ("python-tests", [sys.executable, "-m", "pytest", "-q", "crates/meshers-python/tests/test_size_field.py",
            "crates/meshers-python/tests/test_api.py", "crates/meshers-python/tests/test_adaptive.py",
            "crates/meshers-python/tests/test_fedoo.py"], 600)]
    for mode in ("uniform", "recovery"):
        tasks.append((mode, [sys.executable, __file__, "--worker", mode, "--run", str(args.run),
            "--source", str(args.source), "--case", args.case, "--nominal", str(args.nominal),
            "--cells", str(args.cells), "--coarse-cells", str(args.coarse_cells),
            "--min-length", str(args.min_length)], args.worker_timeout))
    workers = {}
    with Job(parts=[label for label, _, _ in tasks]) as job:
        for label, command, timeout in tasks:
            with job.part(label) as part:
                workers[label] = workflow.bounded(command, args.run / (label + ".log"), part, timeout=timeout)
                workflow.save(args.run / "workers.json", workers)
                assert workers[label]["returncode"] == 0, f"Read {label}.log"
                if label in ("uniform", "recovery"):
                    indexed = args.run / label
                    indexed.mkdir()
                    for destination, original in (("result.json", label + ".json"),
                            ("mesh.npz", label + "-mesh.npz"), ("fields.npz", label + "-fields.npz")):
                        source = args.run / original
                        if source.exists():
                            os.link(source, indexed / destination)
                job.result(label, peak_gib=workers[label]["peak_rss_bytes"] / 1024**3,
                           seconds=workers[label]["wall_s"])
    import meshers._meshers
    manifest["native_sha256"] = workflow.digest(Path(meshers._meshers.__file__))
    workflow.save(args.run / "manifest.json", manifest)
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())
    assert all(workflow.digest(Path(path)) == value for path, value in manifest["source_sha256"].items())
    workflow.save(args.run / "summary.json", {"inputs_unchanged": True,
        "candidates": [json.loads((args.run / (mode + ".json")).read_text()) for mode in ("uniform", "recovery")]})


if __name__ == "__main__":
    main()
