"""One bounded adaptive admission check on a fixed difficult specimen."""

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
os.environ.update(CARGO_BUILD_JOBS="2", RUSTUP_TOOLCHAIN="1.95.0",
    VIRTUAL_ENV=str(workflow.ROOT / ".venv"),
    PATH=str(Path.home() / ".cargo/bin") + os.pathsep + os.environ["PATH"])
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    np, _, meshers, _, _, _ = workflow.imports()
    from periodic import periodic_pairs

    row = {"case": args.case, "cells": args.cells, "min_length": args.min_length,
           "order": args.order, "native_tolerance": args.native_tolerance, "qualified": False}
    started = time.perf_counter()
    try:
        case, shape, _ = workflow.specimen(args.case)
        field = case.function
        if args.compiled_spinodoid:
            assert args.case == "spinodoid"
            phase = time.perf_counter()
            frep = shape._frep
            coefficients = [(*map(float, mode / frep.cell_size), float(amplitude), float(offset))
                for mode, amplitude, offset in zip(frep.modes, frep.amplitudes, frep.phases, strict=True)]
            def expanded(coordinates_x, coordinates_y, coordinates_z):
                value = float(frep.threshold)
                for mode_x, mode_y, mode_z, amplitude, offset in coefficients:
                    angle = 2 * np.pi * (mode_x * (coordinates_x + 0.5)
                        + mode_y * (coordinates_y + 0.5) + mode_z * (coordinates_z + 0.5)) + offset
                    value = value - amplitude * np.cos(angle)
                return value
            field = meshers.compile_field(expanded)
            queries = np.random.default_rng(20261018).uniform(-0.5, 0.5, (4096, 3))
            expected = case.function(*queries.T)
            actual = field.values(queries)
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
            assert np.array_equal(actual < 0, expected < 0)
            derivative_error = 0.0
            for point in queries[:32]:
                numeric = []
                for direction in np.eye(3):
                    plus, minus = point + 1e-6 * direction, point - 1e-6 * direction
                    numeric.append(float(case.function(*plus) - case.function(*minus)) / 2e-6)
                analytic = np.asarray(field.gradient(*point))
                np.testing.assert_allclose(analytic, numeric, rtol=2e-7, atol=2e-7)
                derivative_error = max(derivative_error, float(np.max(abs(analytic - numeric))))
            row["compiled_field"] = {"modes": len(coefficients), "expression_nodes": field.node_count,
                "maximum_value_difference": float(np.max(abs(actual - expected))),
                "maximum_gradient_difference": derivative_error, "verification_points": len(queries),
                "preparation_and_verification_s": time.perf_counter() - phase}
            workflow.save(args.run / "fourier-coefficients.json", {
                "coefficients": coefficients, "threshold": float(frep.threshold)})
        if args.clip_void is not None:
            assert args.case == "smooth-lattice" and args.clip_void > 0
            original = field
            def clipped(coordinates_x, coordinates_y, coordinates_z):
                return np.minimum(original(coordinates_x, coordinates_y, coordinates_z), args.clip_void)
            field = clipped
            queries = np.random.default_rng(20261018).uniform(-0.5, 0.5, (4096, 3))
            expected = original(*queries.T)
            actual = field(*queries.T)
            assert np.array_equal(actual < 0, expected < 0)
            np.testing.assert_array_equal(actual[expected <= 0], expected[expected <= 0])
            row["void_clamp"] = {"positive_ceiling": args.clip_void,
                "verification_points": len(queries), "material_mismatches": 0,
                "domain_identity": "min(phi, positive_constant) <= 0 iff phi <= 0"}
        phase = time.perf_counter()
        mesh = meshers.generate(field, bounds=workflow.BOUNDS, cells=args.cells,
            band=case.band, periodic=(True, True, True), threads=1, optimize_passes=1,
            minimum_quality=0.1, geometry_tolerance=args.native_tolerance, max_tetrahedra=workflow.CAP,
            adapt=True, min_length=args.min_length, max_length=0.1, gradation=0.3,
            adapt_iterations=4, quality_guard=0.3, adapt_tolerance=0.002, order=args.order)
        row.update(generation_s=time.perf_counter() - phase, nodes=len(mesh.points),
                   elements=len(mesh.tetrahedra), diagnostics=mesh.diagnostics)
        np.savez_compressed(args.run / "mesh.npz", points=mesh.points,
                            tetrahedra=mesh.tetrahedra, surface=mesh.surface)
        row["geometry"] = workflow.geometry(args.run, args.case, mesh.points, mesh.tetrahedra, mesh.surface)
        _, _, row["periodic"] = periodic_pairs(mesh.points, mesh.surface)
        row["qualified"] = row["geometry"]["geometry_admitted"]
        row["status"] = "admitted-needs-reference" if row["qualified"] else "geometry-rejected"
        if row["qualified"] and args.solve:
            if len(mesh.points) > 100_000:
                row["status"] = "solve-node-cap"
            else:
                from stress_recovery import solve
                phase = time.perf_counter()
                fields, tensor, records = solve(mesh.points, mesh.tetrahedra)
                row.update(solve_s=time.perf_counter() - phase, stiffness=tensor.tolist(),
                           solver_records=records, status="solved-needs-reference")
                np.savez_compressed(args.run / "fields.npz", displacement=fields)
    except Exception as error:
        row.update(status="failed", error=str(error), traceback=traceback.format_exc())
    row["total_audit_workflow_s"] = time.perf_counter() - started
    if args.solve:
        row["total_workflow_s"] = row["total_audit_workflow_s"]
    workflow.save(args.run / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--case", default="spinodoid")
    parser.add_argument("--cells", type=int, default=16)
    parser.add_argument("--min-length", type=float, default=0.01)
    parser.add_argument("--compiled-spinodoid", action="store_true")
    parser.add_argument("--order", type=int, default=1, choices=(1, 2))
    parser.add_argument("--native-tolerance", type=float, default=0.002)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument("--clip-void", type=float)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--python-tests", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    source = workflow.HERE / "runs/g3-microgen-defaults-v4"
    inputs = [source / (args.case + suffix) for suffix in ("-samples.npz", "-specimen.json")]
    hashes = {str(path): workflow.digest(path) for path in inputs}
    for path in inputs:
        shutil.copyfile(path, args.run / path.name)
    import meshers._meshers
    workflow.save(args.run / "manifest.json", {"script_sha256": workflow.digest(Path(__file__)),
        "native_sha256": workflow.digest(Path(meshers._meshers.__file__)), "inputs": hashes,
        "memory_bytes": workflow.MEMORY, "threads": 1, "timeout_s": 600})
    parts = ["build", "core-tests", "adaptive-admission"] if args.build else ["adaptive-admission"]
    if args.python_tests:
        parts.insert(-1, "python-tests")
    with Job(parts=parts) as job:
        if args.build:
            commands = [
                ("build", [str(workflow.ROOT / ".venv/bin/maturin"), "develop", "--release", "--locked",
                    "--manifest-path", str(workflow.ROOT / "crates/meshers-python/Cargo.toml"),
                    "--features", "experimental-surfaces"]),
                ("core-tests", ["cargo", "test", "--release", "--locked", "-p", "meshers-core", "--lib",
                    "--", "--test-threads=1"])]
            for label, command in commands:
                with job.part(label) as part:
                    result = workflow.bounded(command, args.run / (label + ".log"), part, timeout=1200)
                    workflow.save(args.run / (label + ".json"), result)
                    assert result["returncode"] == 0, f"Read {label}.log"
                    job.result(label, seconds=result["wall_s"], peak_gib=result["peak_rss_bytes"] / 1024**3)
        if args.python_tests:
            with job.part("python-tests") as part:
                command = [sys.executable, "-m", "pytest", "-q"] + [
                    str(workflow.ROOT / "crates/meshers-python/tests" / name) for name in
                    ("test_api.py", "test_adaptive.py", "test_size_field.py", "test_fedoo.py")]
                result = workflow.bounded(command, args.run / "python-tests.log", part, timeout=600)
                workflow.save(args.run / "python-tests.json", result)
                assert result["returncode"] == 0, "Read python-tests.log"
                job.result("python-tests", seconds=result["wall_s"], peak_gib=result["peak_rss_bytes"] / 1024**3)
        with job.part("adaptive-admission") as part:
            command = [sys.executable, __file__, "--run", str(args.run),
                "--case", args.case, "--cells", str(args.cells), "--min-length", str(args.min_length),
                "--order", str(args.order), "--native-tolerance", str(args.native_tolerance),
                "--worker"]
            if args.compiled_spinodoid:
                command.append("--compiled-spinodoid")
            if args.solve:
                command.append("--solve")
            if args.clip_void is not None:
                command.extend(["--clip-void", str(args.clip_void)])
            result = workflow.bounded(command, args.run / "worker.log", part, timeout=600)
            workflow.save(args.run / "worker.json", result)
            assert result["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("adaptive-admission", qualified=int(row["qualified"]),
                       nodes=row.get("nodes", 0), seconds=row["total_audit_workflow_s"])
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())
    manifest = json.loads((args.run / "manifest.json").read_text())
    manifest["native_sha256_after"] = workflow.digest(Path(meshers._meshers.__file__))
    manifest["tet_adapt_sha256"] = workflow.digest(workflow.ROOT / "crates/meshers-core/src/tet_adapt.rs")
    workflow.save(args.run / "manifest.json", manifest)


if __name__ == "__main__":
    main()
