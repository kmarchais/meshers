"""One bounded adaptive admission check on a fixed difficult specimen."""

import argparse
import importlib.metadata
import json
import os
import platform
from pathlib import Path
import shutil
import sys
import time
import traceback

import microgen_workflow as workflow
from job_parts import worker_part

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", RAYON_NUM_THREADS="1")
os.environ.update(CARGO_BUILD_JOBS="2", RUSTUP_TOOLCHAIN="1.95.0",
    VIRTUAL_ENV=str(Path(sys.executable).parent.parent),
    PATH=str(Path.home() / ".cargo/bin") + os.pathsep + os.environ["PATH"])
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    np, _, meshers, _, _, _ = workflow.imports()
    from periodic import periodic_pairs

    row = {"case": args.case, "cells": args.cells, "min_length": args.min_length,
           "max_length": args.max_length, "interior_refine_nodes": args.interior_refine_nodes,
           "order": args.order, "native_tolerance": args.native_tolerance,
           "optimization_passes": args.optimization_passes,
           "adaptive": not args.no_adapt, "adapt_tolerance": args.adapt_tolerance,
           "solver_quadrature_points": args.quadrature_points if args.order == 2 else 4,
           "solve_node_cap": args.solve_node_cap, "qualified": False}
    row["solver"] = "cached-elastic-block-amg-cg" if args.amg else "jacobi-cg"
    row["crease_grading_requested"] = getattr(args, "crease_grading", False)
    row["crease_bisection_requested"] = getattr(args, "crease_bisection", False)
    row["crease_longest_closure_requested"] = getattr(args, "crease_longest_closure", False)
    row["crease_quality_closure_requested"] = getattr(args, "crease_quality_closure", False)
    row["quadratic_quality_pass_requested"] = getattr(args, "quadratic_quality_pass", False)
    started = time.perf_counter()
    if args.order == 2 and args.quadrature_points == 4:
        row.update(timing_comparable=False,
            timing_scope="Reference study includes extra own-rule and fifteen-point fixed-field audits")
    try:
        case, shape, _ = workflow.specimen(args.case)
        field = case.function
        if args.symmetric_union:
            assert args.case == "smooth-lattice"
            from gradient_tie_probe import symmetric_lattice
            field = symmetric_lattice
            prepared = meshers.compile_field(field)
            queries = np.random.default_rng(20261018).uniform(-0.5, 0.5, (4096, 3))
            expected, actual = case.function(*queries.T), prepared.values(queries)
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
            assert np.array_equal(actual < 0, expected < 0)
            tied = np.array([[0., 0.17, -0.17], [0.17, 0., -0.17], [-0.17, 0.17, 0.]])
            import scorecard
            numeric = scorecard.gradient(case.function, tied)
            analytic = np.array([prepared.gradient(*point) for point in tied])
            np.testing.assert_allclose(analytic, numeric, rtol=1e-7, atol=1e-7)
            row["symmetric_union"] = {"identity": "min(a,b) = (a+b-abs(a-b))/2",
                "maximum_field_value_difference": float(abs(actual - expected).max()),
                "maximum_tie_gradient_difference": float(abs(analytic - numeric).max()),
                "value_verification_points": len(queries), "gradient_verification_points": len(tied)}
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
        size_field = None
        if args.cap_length is not None:
            def size_field(coordinates_x, coordinates_y, coordinates_z):
                distance = np.minimum.reduce([.5 - abs(coordinates_x),
                    .5 - abs(coordinates_y), .5 - abs(coordinates_z)])
                fraction = np.clip(distance / args.cap_band, 0, 1)
                smooth = fraction**2 * (3 - 2 * fraction)
                return args.cap_length + (.1 - args.cap_length) * smooth
            queries = np.random.default_rng(20261018).uniform(-.5, .5, (4096, 3))
            values = size_field(*queries.T)
            assert np.isfinite(values).all() and (values > 0).all()
            for axis in range(3):
                low, high = queries.copy(), queries.copy()
                low[:, axis], high[:, axis] = -.5, .5
                np.testing.assert_array_equal(size_field(*low.T), size_field(*high.T))
                np.testing.assert_allclose(size_field(*low.T), args.cap_length, atol=1e-15, rtol=0)
            np.testing.assert_allclose(size_field(*np.zeros((3, 1))), .1, atol=1e-15, rtol=0)
            row["cap_size_field"] = {"cap_target_length": args.cap_length,
                "transition_band": args.cap_band, "interior_maximum_length": .1,
                "formula": "cap + (0.1-cap)*smoothstep(min_box_plane_distance/band)",
                "periodic_verification_points_per_axis": len(queries)}
        generation_options = dict(bounds=workflow.BOUNDS, cells=args.cells,
            band=case.band, periodic=(True, True, True), threads=1, optimize_passes=args.optimization_passes,
            minimum_quality=0.1,
            geometry_tolerance=args.native_tolerance, max_tetrahedra=workflow.CAP,
            adapt=not args.no_adapt, min_length=args.min_length, max_length=args.max_length, gradation=0.3,
            adapt_iterations=4, quality_guard=0.3, adapt_tolerance=args.adapt_tolerance, order=args.order,
            size_field=size_field, interior_refine_nodes=args.interior_refine_nodes)
        mesh = None
        if getattr(args, "crease_grading", False) or getattr(args, "crease_bisection", False):
            seed = meshers.generate(field, **generation_options)
            assert max(len(seed.points), len(seed.tetrahedra)) <= workflow.CAP
            directory = args.run / "feature-seed"
            directory.mkdir(exist_ok=False)
            for suffix in ("-samples.npz", "-specimen.json"):
                shutil.copyfile(args.run / (args.case+suffix), directory / (args.case+suffix))
            np.savez_compressed(directory / "mesh.npz", points=seed.points,
                tetrahedra=seed.tetrahedra, surface=seed.surface)
            seed_geometry = workflow.geometry(directory, args.case, seed.points, seed.tetrahedra, seed.surface)
            _, _, seed_periodic = periodic_pairs(seed.points, seed.surface)
            assert seed_geometry["geometry_admitted"]
            if getattr(args, "crease_bisection", False):
                from crease_directional_refinement import generate
                mesh, metadata = generate(seed, case.function, args, seed_geometry["volume"])
                row.update(route="meshers-experimental-exact-map-directional-refinement",
                    native_only=False, saved_seed_replayed=False)
            else:
                from periodic_crease_grading import make_metric
                target, metadata, cloud = make_metric(seed, case.function, args.min_length, args.max_length, args.cells)
                np.savez_compressed(args.run / "crease-cloud.npz", canonical_points=cloud)
                generation_options["size_field"] = target
            metadata.update(seed_nodes=len(seed.points), seed_tets=len(seed.tetrahedra),
                seed_mesh_sha256=workflow.digest(directory / "mesh.npz"),
                seed_geometry=seed_geometry, seed_periodic=seed_periodic)
            row["crease_bisection" if getattr(args, "crease_bisection", False) else "crease_grading"] = metadata
            row.update(timing_comparable=False,
                timing_scope="Fresh native seed, geometry controls and declared grading/refinement stages are charged; no paired timing")
            del seed
        if mesh is None:
            mesh = meshers.generate(field, **generation_options)
        if args.quadratic_quality_pass:
            from types import SimpleNamespace
            from quadratic_quality_preflight import apply_patches,QUALITY_PASS_POLICY
            from quadratic_edge_removal_probe import feasible_patches,measure
            assert max(len(mesh.points),len(mesh.tetrahedra)) <= workflow.CAP
            directory = args.run / "quality-source"
            directory.mkdir(exist_ok=False)
            np.savez_compressed(directory / "mesh.npz",points=mesh.points,
                tetrahedra=mesh.tetrahedra,surface=mesh.surface)
            for suffix in ("-samples.npz","-specimen.json"):
                shutil.copyfile(args.run/(args.case+suffix),directory/(args.case+suffix))
            _,volume = measure(mesh.points[mesh.tetrahedra],positive_rule=False)
            before_volume = float(volume.sum())
            _,_,before_periodic = periodic_pairs(mesh.points,mesh.surface)
            patches,statistics = feasible_patches(mesh.points,mesh.tetrahedra,mesh.surface,before_volume)
            out,operations = apply_patches(mesh.points,mesh.tetrahedra,mesh.surface,patches)
            _,_,after_periodic = periodic_pairs(out['points'],out['surface'])
            assert before_periodic['periodic_pair_counts']==after_periodic['periodic_pair_counts']
            workflow.save(args.run / "quality-operations.json",operations)
            quality = dict(policy=QUALITY_PASS_POLICY,source_nodes=len(mesh.points),
                source_elements=len(mesh.tetrahedra),source_volume=before_volume,
                source_mesh_sha256=workflow.digest(directory / "mesh.npz"),
                operations_sha256=workflow.digest(args.run / "quality-operations.json"),
                accepted_disjoint_operations=len(operations),local_statistics=statistics,
                curved_boundary_coordinates_identical=True,old_approximation_space_preserved=not bool(operations),
                fe_fields_used=False,field_transfer_performed=False)
            mesh = SimpleNamespace(**out,diagnostics=dict(source="Curved disjoint quality edge removal",
                before_quality_pass=mesh.diagnostics,quadratic_quality_pass=quality))
            row.update(quadratic_quality_pass=quality,old_approximation_space_preserved=not bool(operations),
                route="meshers-experimental-curved-quality-edge-removal",native_only=False,
                saved_seed_replayed=False,timing_comparable=False,
                timing_scope="Fresh source generation and curved quality operation charged; no paired timing")
            if args.crease_bisection and operations:
                row['crease_bisection']['approximation_space_scope'] = 'Nested only before subsequent quality retopology'
        row.update(generation_s=time.perf_counter() - phase, nodes=len(mesh.points),
                   elements=len(mesh.tetrahedra), diagnostics=mesh.diagnostics)
        np.savez_compressed(args.run / "mesh.npz", points=mesh.points,
                            tetrahedra=mesh.tetrahedra, surface=mesh.surface)
        if len(mesh.points) > workflow.CAP or len(mesh.tetrahedra) > workflow.CAP:
            raise ValueError("Mesh node/tetrahedron admission cap exceeded")
        row["geometry"] = workflow.geometry(args.run, args.case, mesh.points, mesh.tetrahedra, mesh.surface)
        if args.quadratic_quality_pass:
            assert abs(row['geometry']['volume']/row['quadratic_quality_pass']['source_volume']-1)<1e-10
        _, _, row["periodic"] = periodic_pairs(mesh.points, mesh.surface)
        row["qualified"] = row["geometry"]["geometry_admitted"]
        row["status"] = "admitted-needs-reference" if row["qualified"] else "geometry-rejected"
        workflow.save(args.run / "result.json", row)
        if row["qualified"] and args.solve:
            if len(mesh.points) > args.solve_node_cap:
                row["status"] = "solve-node-cap"
            elif args.isolate_solve:
                import psutil
                row.update(solve_process_boundary="exec-replace", generation_process_id=os.getpid(),
                    generation_memory_before_exec_bytes=psutil.Process().memory_info().rss,
                    solve_exec_clock=time.perf_counter(), solve_mesh_sha256=workflow.digest(args.run / "mesh.npz"))
                workflow.save(args.run / "result.json", row)
                os.execv(sys.executable, [sys.executable, str(workflow.HERE / "fresh_mesh_solve.py"),
                    "--run", str(args.run), "--started", str(started)])
            else:
                from stress_recovery import solve, vector_first_reactions
                phase = time.perf_counter()
                def checkpoint(completed, partial_tensor, records):
                    row.update(status="solving", completed_loads=completed,
                        partial_stiffness_columns=partial_tensor.tolist(), solver_records=records)
                    workflow.save(args.run / "result.json", row)
                solver, save_field = None, None
                if args.amg:
                    from amg_reference import AuditedAMG
                    solver = AuditedAMG(elastic_blocks=True)
                    def save_field(completed, displacement):
                        np.savez_compressed(args.run / f"load-{completed}.npz", displacement=displacement)
                fields, tensor, records = solve(mesh.points, mesh.tetrahedra, checkpoint=checkpoint,
                    n_elm_gp=row["solver_quadrature_points"], solver=solver,
                    field_checkpoint=save_field, reaction_reader=vector_first_reactions if args.amg else None,
                    release_reduced_cache=args.amg)
                row.pop("partial_stiffness_columns", None)
                row.update(solve_s=time.perf_counter() - phase, stiffness=tensor.tolist(),
                           solver_records=records, status="solved-needs-reference")
                np.savez_compressed(args.run / "fields.npz", displacement=fields)
                if args.order == 2 and args.quadrature_points == 4:
                    from periodic import tensor_difference
                    from tet10_integration import audit_saved, integrated_tensors
                    workflow.save(args.run / "result.json", row)
                    audit_saved(args.run, args.run / "four-point-field-consistency.json")
                    _, higher_energy, _ = integrated_tensors(mesh.points, mesh.tetrahedra, fields)
                    row["higher_rule_fixed_field_difference"] = tensor_difference(higher_energy, tensor)
                    row["quadrature_sensitivity_limit"] = 1e-4
                    row["quadrature_sensitivity_passed"] = row["higher_rule_fixed_field_difference"][
                        "worst_relative_strain_energy_difference"] <= row["quadrature_sensitivity_limit"]
                    row["quadrature_scope"] = "Own-rule equilibrium and fifteen-point reintegration; not a fifteen-point re-solve"
                    if not row["quadrature_sensitivity_passed"]:
                        row.update(qualified=False, status="quadrature-sensitivity-rejected")
    except Exception as error:
        row.update(status="failed", qualified=False, error=str(error), traceback=traceback.format_exc())
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
    parser.add_argument("--max-length", type=float, default=0.1)
    parser.add_argument("--interior-refine-nodes", type=int, default=0)
    parser.add_argument("--crease-grading", action="store_true")
    parser.add_argument("--crease-bisection", action="store_true")
    parser.add_argument("--crease-longest-closure", action="store_true")
    parser.add_argument("--crease-quality-closure", action="store_true")
    parser.add_argument("--quadratic-quality-pass", action="store_true")
    parser.add_argument("--compiled-spinodoid", action="store_true")
    parser.add_argument("--order", type=int, default=1, choices=(1, 2))
    parser.add_argument("--native-tolerance", type=float, default=0.002)
    parser.add_argument("--optimization-passes", type=int, default=1)
    parser.add_argument("--adapt-tolerance", type=float, default=0.002)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument("--isolate-solve", action="store_true")
    parser.add_argument("--memory-gib", type=int, choices=(8, 10, 12), default=8)
    parser.add_argument("--clip-void", type=float)
    parser.add_argument("--no-adapt", action="store_true")
    parser.add_argument("--symmetric-union", action="store_true")
    parser.add_argument("--solve-node-cap", type=int, choices=(100_000, 150_000, 200_000, 250_000), default=100_000)
    parser.add_argument("--quadrature-points", type=int, choices=(4, 15), default=15)
    parser.add_argument("--worker-timeout", type=int, choices=(600, 1200, 1800), default=600)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--python-tests", action="store_true")
    parser.add_argument("--cap-length", type=float)
    parser.add_argument("--cap-band", type=float, default=.05)
    parser.add_argument("--amg", action="store_true")
    args = parser.parse_args()
    assert args.interior_refine_nodes >= 0
    assert not args.quadratic_quality_pass or args.order == 2
    assert not args.isolate_solve or (args.solve and args.amg and args.order == 2 and args.quadrature_points == 15)
    assert not (args.crease_grading and args.crease_bisection)
    assert not args.crease_longest_closure or args.crease_bisection
    assert not args.crease_quality_closure or (args.crease_bisection and not args.crease_longest_closure)
    assert not (args.crease_grading or args.crease_bisection) or (args.case == "smooth-lattice" and not args.no_adapt
        and args.order == 2 and args.cap_length is None and args.interior_refine_nodes == 0)
    assert args.cap_length is None or (not args.no_adapt and 0 < args.cap_length < .1 and 0 < args.cap_band < .5)
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
    python_sources = [Path(__file__), Path(workflow.__file__), workflow.HERE / "stress_recovery.py",
                      workflow.HERE / "tet10_integration.py",
                      workflow.HERE / "job_parts.py", workflow.HERE / "facet_incidence_audit.py",
                      workflow.HERE / "gradient_tie_probe.py", workflow.ROOT / "crates/meshers-python/python/meshers/expression.py"]
    python_sources.extend(workflow.ROOT.parent / "periodic-microgen/experiments/periodic-fea" / name
        for name in ("homogenization.py", "periodic.py"))
    source_hashes = {str(path): workflow.digest(path) for path in python_sources}
    if args.crease_grading or args.crease_bisection:
        source_hashes.update({str(workflow.HERE/name): workflow.digest(workflow.HERE/name)
            for name in ("periodic_crease_grading.py", "lattice_crease_mechanism.py", "periodic_crease_audit.py")})
    if args.crease_bisection:
        source_hashes.update({str(workflow.HERE/name): workflow.digest(workflow.HERE/name)
            for name in ("crease_directional_refinement.py", "test_crease_directional_refinement.py", "nested_refinement.py", "test_nested_refinement.py", "hp_shared_domain.py")})
    if args.crease_longest_closure:
        source_hashes.update({str(workflow.HERE/name): workflow.digest(workflow.HERE/name)
            for name in ("longest_edge_closure.py", "test_longest_edge_closure.py")})
    if args.crease_quality_closure:
        source_hashes.update({str(workflow.HERE/name): workflow.digest(workflow.HERE/name)
            for name in ("quality_split_closure.py", "test_quality_split_closure.py",
                "longest_edge_closure.py", "test_longest_edge_closure.py", "fixtures/quality-split-trap.json")})
    if args.isolate_solve:
        source_hashes[str(workflow.HERE / "fresh_mesh_solve.py")] = workflow.digest(workflow.HERE / "fresh_mesh_solve.py")
    if args.quadratic_quality_pass:
        source_hashes.update({str(workflow.HERE/name):workflow.digest(workflow.HERE/name)
            for name in ('quadratic_quality_preflight.py','quadratic_edge_removal_probe.py',
                'test_quadratic_quality_pass.py','nested_refinement.py','quadrature_audit.py')})
    if args.amg:
        source_hashes[str(workflow.HERE / "amg_reference.py")] = workflow.digest(workflow.HERE / "amg_reference.py")
    for index, path in enumerate(source_hashes):
        shutil.copyfile(path, args.run / f"input-{index}-{Path(path).name}")
    import meshers._meshers
    workflow.save(args.run / "manifest.json", {"script_sha256": workflow.digest(Path(__file__)),
        "native_sha256": workflow.digest(Path(meshers._meshers.__file__)), "inputs": hashes,
        "platform": platform.platform(), "architecture": platform.machine(),
        "python": sys.version, "executable": sys.executable,
        "virtual_env": os.environ["VIRTUAL_ENV"],
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY,
        "packages": {name: importlib.metadata.version(name) for name in
            ("meshers", "numpy", "scipy", "fedoo", "microgen", "mmgpy", "pyamg", "pyvista", "vtk")},
        "memory_bytes": args.memory_gib*1024**3, "threads": 1, "timeout_s": args.worker_timeout,
        "max_length": args.max_length,
        "interior_refine_nodes": args.interior_refine_nodes,
        "crease_grading_requested": args.crease_grading,
        "crease_bisection_requested": args.crease_bisection,
        "crease_longest_closure_requested": args.crease_longest_closure,
        "crease_quality_closure_requested": args.crease_quality_closure,
        "quadratic_quality_pass_requested": args.quadratic_quality_pass,
        "solve_process_boundary": "exec-replace" if args.isolate_solve else "same-process",
        "python_source_sha256": source_hashes,
        "solver": "cached-elastic-block-amg-cg" if args.amg else "jacobi-cg",
        "pyamg": importlib.metadata.version("pyamg") if args.amg else None})
    parts = ["build", "core-tests", "adaptive-admission"] if args.build else ["adaptive-admission"]
    if args.python_tests:
        parts.insert(-1, "python-tests")
    with Job(parts=parts) as job:
        if args.build:
            commands = [
                ("build", [str(Path(sys.executable).with_name("maturin")), "develop", "--release", "--locked",
                    "--manifest-path", str(workflow.ROOT / "crates/meshers-python/Cargo.toml"),
                    "--features", "experimental-surfaces"]),
                ("core-tests", ["cargo", "test", "--release", "--locked", "-p", "meshers-core", "--lib",
                    "--", "--test-threads=1"])]
            for label, command in commands:
                with worker_part(job, label) as part:
                    result = workflow.bounded(command, args.run / (label + ".log"), part, timeout=1200)
                    workflow.save(args.run / (label + ".json"), result)
                    assert result["returncode"] == 0, f"Read {label}.log"
                    job.result(label, seconds=result["wall_s"], peak_gib=result["peak_rss_bytes"] / 1024**3)
        if args.python_tests:
            with worker_part(job, "python-tests") as part:
                command = [sys.executable, "-m", "pytest", "-q"] + [
                    str(workflow.ROOT / "crates/meshers-python/tests" / name) for name in
                    ("test_api.py", "test_adaptive.py", "test_size_field.py", "test_fedoo.py")]
                result = workflow.bounded(command, args.run / "python-tests.log", part, timeout=600)
                workflow.save(args.run / "python-tests.json", result)
                assert result["returncode"] == 0, "Read python-tests.log"
                job.result("python-tests", seconds=result["wall_s"], peak_gib=result["peak_rss_bytes"] / 1024**3)
        with worker_part(job, "adaptive-admission", memory_bytes=args.memory_gib*1024**3) as part:
            command = [sys.executable, __file__, "--run", str(args.run),
                "--case", args.case, "--cells", str(args.cells), "--min-length", str(args.min_length),
                "--max-length", str(args.max_length),
                "--interior-refine-nodes", str(args.interior_refine_nodes),
                "--order", str(args.order), "--native-tolerance", str(args.native_tolerance),
                "--optimization-passes", str(args.optimization_passes),
                "--adapt-tolerance", str(args.adapt_tolerance),
                "--solve-node-cap", str(args.solve_node_cap),
                "--quadrature-points", str(args.quadrature_points),
                "--worker-timeout", str(args.worker_timeout),
                "--memory-gib", str(args.memory_gib),
                "--worker"]
            if args.compiled_spinodoid:
                command.append("--compiled-spinodoid")
            if args.crease_grading:
                command.append("--crease-grading")
            if args.crease_bisection:
                command.append("--crease-bisection")
            if args.crease_longest_closure:
                command.append("--crease-longest-closure")
            if args.crease_quality_closure:
                command.append("--crease-quality-closure")
            if args.quadratic_quality_pass:
                command.append("--quadratic-quality-pass")
            if args.amg:
                command.append("--amg")
            if args.cap_length is not None:
                command.extend(["--cap-length", str(args.cap_length), "--cap-band", str(args.cap_band)])
            if args.solve:
                command.append("--solve")
            if args.isolate_solve:
                command.append("--isolate-solve")
            if args.clip_void is not None:
                command.extend(["--clip-void", str(args.clip_void)])
            if args.no_adapt:
                command.append("--no-adapt")
            if args.symmetric_union:
                command.append("--symmetric-union")
            result = workflow.bounded(command, args.run / "worker.log", part, timeout=args.worker_timeout,
                memory_bytes=args.memory_gib*1024**3)
            workflow.save(args.run / "worker.json", result)
            assert result["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("adaptive-admission", qualified=int(row["qualified"]),
                       nodes=row.get("nodes", 0), seconds=row["total_audit_workflow_s"],
                       completed_loads=row.get("completed_loads", 0))
            assert row["qualified"], "Read result.json: geometry admission failed"
            if args.solve:
                assert (row["status"] == "solved-needs-reference"
                    and row.get("completed_loads") == 6 and (args.run / "fields.npz").exists()), \
                    "Read result.json: six fresh loads did not complete"
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())
    assert all(workflow.digest(Path(path)) == value for path, value in source_hashes.items())
    manifest = json.loads((args.run / "manifest.json").read_text())
    manifest["native_sha256_after"] = workflow.digest(Path(meshers._meshers.__file__))
    manifest["tet_adapt_sha256"] = workflow.digest(workflow.ROOT / "crates/meshers-core/src/tet_adapt.rs")
    workflow.save(args.run / "manifest.json", manifest)


if __name__ == "__main__":
    main()
