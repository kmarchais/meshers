"""Continue the frozen BCC fixture through CAD sampling and one native attempt."""

import argparse
from functools import partial
import json
from pathlib import Path
import shutil
import sys
import time
import traceback
from types import SimpleNamespace

import microgen_workflow as w
import cad_lattice_fixture as fixture
import cad_bcc_coverage_control as original
import cad_bcc_classifier_diagnosis as diagnosis

SOURCE = w.HERE/"runs/macbook-cad-bcc-coverage-v1"
DIAGNOSIS = w.HERE/"runs/macbook-cad-bcc-classifier-v1"
RESOLUTION = w.HERE/"runs/macbook-cad-bcc-membership-resolution-v1"


def compare_records(actual, expected):
    """Check identical entity structure and numerical support within 1e-12."""
    changes = []
    def visit(a, b, path):
        if isinstance(b, dict):
            assert isinstance(a, dict) and a.keys() == b.keys(), path
            for key in b:
                visit(a[key], b[key], path+"/"+key)
        elif isinstance(b, list):
            assert isinstance(a, list) and len(a) == len(b), path
            for index, (x, y) in enumerate(zip(a, b, strict=True)):
                visit(x, y, path+"/"+str(index))
        elif isinstance(b, float):
            difference = abs(float(a)-b)
            assert difference <= 1e-12, f"Changed CAD entity support: {path}: {difference}"
            if difference:
                changes.append(dict(path=path, absolute_difference=difference))
        else:
            assert a == b, f"Changed CAD entity structure: {path}"
    visit(actual, expected, "inventory")
    return dict(passed=True, numerical_absolute_tolerance=1e-12,
        maximum_numeric_difference=max((row["absolute_difference"] for row in changes), default=0.), differences=changes)


def worker(args):
    np, _, meshers, microgen, _, scorecard = w.imports()
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    old = json.loads((args.source/"result.json").read_text())
    proof = json.loads((args.resolution/"result.json").read_text())
    assert old["complete"] and old["input_hashes_verified"] and not old["native_attempted"]
    assert proof["complete"] and proof["input_hashes_verified"] and proof["resolution_conditions_met"]
    assert not proof["unresolved_ambiguity_ids"] and not proof["original_classifier_pass"]
    assert old["native_options"] == json.loads(json.dumps(original.OPTIONS))
    frozen_inputs = dict(proof["inputs"])
    for path, digest in old["inputs"].items():
        if "/python/meshers/" not in path:
            frozen_inputs[path] = digest
    for path, digest in frozen_inputs.items():
        assert w.digest(Path(path)) == digest, f"Changed frozen fixture or observer input: {path}"
    paths = {Path(path) for path in frozen_inputs}
    paths.update((Path(__file__).resolve(), Path(original.__file__).resolve(), Path(diagnosis.__file__).resolve(),
        Path(meshers.__file__), Path(meshers._meshers.__file__), args.resolution/"result.json", args.resolution/"queries.json"))
    paths.update(Path(meshers.__file__).parent/name for name in ("compiler.py", "expression.py", "intersection.py"))
    hashes = {str(path): w.digest(path) for path in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    for name in ("primitive-receipt.json", "raw-primitive-receipt.json", "membership-probes.npz", "membership-probes-input.npz", "actual-bcc.step"):
        shutil.copyfile(args.source/name, args.run/name)
        assert w.digest(args.run/name) == w.digest(args.source/name)
    options = dict(original.OPTIONS)
    result = dict(inputs=hashes, complete=False, fixture_verified=False, native_attempted=False, native_admitted=False,
        native_options=options, primitive_controls=old["controls"], controls_reused_by_verified_hash=True,
        original_classifier_pass=False, accepted_cap_resolution=str(args.resolution/"result.json"),
        membership_policy=proof["policy"], source_cad="Reconstructed unchanged public CAD; copied original STEP is provenance only",
        volume_reference_source=str(args.source/"result.json"), native_binary_sha256=w.digest(Path(meshers._meshers.__file__)),
        previous_native_binary_sha256=old["inputs"][str(Path(meshers._meshers.__file__))],
        comparative_timing_claim=False, actual_mmg_comparison=False, finite_element_solve_performed=False, production_change=False)
    w.save(args.run/"manifest.json", dict(inputs=hashes, native_options=options, worker_timeout_s=240, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, vtk_backend=vtkSMPTools.GetBackend(), maximum_nodes=w.CAP,
        maximum_tetrahedra=w.CAP, source_geometry="Original ordered finite cylinders, periodic split/fuse and cell clip",
        membership_policy=proof["policy"], geometry_gates=dict(distance=.002, absolute_relative_volume=.01, corner_quality=.1,
            positive125=True, full_Tri6_periodic_caps=True, coherent_oriented_facets=True)))
    w.save(args.run/"result.json", result)
    stage = "frozen-fixture"
    try:
        with Job(parts={"cad": 2, "membership": 2, "samples": 3, "native": 3, "admission": 4}) as job:
            with job.part("cad") as part:
                stage = "cad"
                part.update(0, "Reconstructing the unchanged public BCC CAD", resolved_probes=proof["resolved_ambiguities"])
                params = old["primitive_parameters"]
                shape = microgen.BodyCenteredCubic(strut_radius=params["radius"], cell_size=params["cell_size"], strut_joints=params["strut_joints"])
                np.testing.assert_array_equal(shape.vertices, params["vertices"])
                np.testing.assert_array_equal(shape.strut_vertex_pairs, params["strut_vertex_pairs"])
                recipe = json.loads((args.source/"primitive-receipt.json").read_text())
                assert fixture.primitives(shape) == recipe["original_ordered_primitives"]
                field = fixture.field(recipe["normalized_primitives"])
                started = time.perf_counter(); cad = shape.generate_cad()
                result["cad_generation_s"] = time.perf_counter()-started
                context, inventory = diagnosis.inventory(cad)
                w.save(args.run/"cad-inventory.json", inventory)
                frozen_inventory = json.loads((args.diagnosis/"fresh-inventory.json").read_text())
                result["cad_support_reproduction"] = compare_records(inventory, frozen_inventory)
                assert inventory["valid"] and inventory["solid_count"] == 1 and inventory["shells_closed"]
                assert not inventory["central_support_planes"] and cad._mesh_volume is None
                props = GProp_GProps()
                estimate = float(BRepGProp.VolumeProperties_s(cad.wrapped, props, 1e-10, True, True))
                regenerated_volume = float(abs(props.Mass()))
                volume = old["cad"]["volume"]
                assert np.isfinite(estimate) and estimate <= 1e-8 and abs(regenerated_volume-volume) <= 1e-12
                result["cad"] = dict(valid=True, volume_reference=volume, reconstructed_volume=regenerated_volume,
                    volume_difference=regenerated_volume-volume, quadrature_requested_relative_tolerance=1e-10,
                    quadrature_estimated_relative_error=estimate, reference_scope=old["cad"]["reference_scope"],
                    original_step_sha256=w.digest(args.run/"actual-bcc.step"), step_used_for_sampling=False)
                part.update(1, "Original CAD topology and numerical volume reproduced", faces=len(inventory["faces"]), volume=volume)
            with job.part("membership") as part:
                stage = "membership"
                with np.load(args.source/"membership-probes.npz") as saved:
                    queries, states, values = [saved[key] for key in ("points", "cad_states", "normalized_field")]
                np.testing.assert_array_equal(np.maximum(field(*queries.T), abs(queries).max(1)-.5), values)
                classifications = diagnosis.classify(cad, queries, part)
                assert all(row["reused"] == str(state) and row["fresh"] == str(state) for row, state in zip(classifications, states, strict=True))
                resolved_rows = json.loads((args.resolution/"queries.json").read_text())
                assert len(resolved_rows) == len(queries) == proof["original_queries"]
                np.testing.assert_array_equal([row["point"] for row in resolved_rows], queries)
                assert all(row["original_cad_state"] == str(state) for row, state in zip(resolved_rows, states, strict=True))
                assert all(row["resolution"] != "unsupported-ambiguity" for row in resolved_rows)
                w.save(args.run/"membership-reproduction.json", dict(queries=len(queries), original_states_reproduced=True,
                    original_classifier_pass=False, accepted_resolution_sha256=w.digest(args.resolution/"queries.json"),
                    resolved_ambiguities=proof["resolved_ambiguities"], proposed_state_counts=proof["proposed_state_counts"], classifications=classifications))
                result["membership"] = dict(original_queries=len(queries), unchanged_points=True, original_states_reproduced=True,
                    original_classifier_pass=False, accepted_resolution_conditions_met=True, resolution_sha256=w.digest(args.resolution/"result.json"))
                job.result("membership", queries=len(queries), retained_original_on_failures=proof["original_mismatch_count"], resolved=proof["resolved_ambiguities"])
            with job.part("samples") as part:
                stage = "samples"
                samples, coverage = fixture.cad_samples(cad, part)
                sample_path = args.run/(original.NAME+"-samples.npz")
                np.savez_compressed(sample_path, **samples)
                wall_values, cap_values = field(*samples["walls"].T), field(*samples["caps"].T)
                wall_error, cap_maximum = float(abs(wall_values).max()), float(cap_values.max())
                worst = int(abs(wall_values).argmax())
                result.update(exact_sample_coverage=coverage, maximum_field_at_cad_wall=wall_error,
                    maximum_field_at_cad_cap=cap_maximum, worst_wall_point=samples["walls"][worst].tolist(),
                    samples_sha256=w.digest(sample_path), samples_archived_before_checks=True)
                w.save(args.run/"cad-sample-check.json", dict(coverage=coverage, wall_residual_max=wall_error,
                    cap_field_max=cap_maximum, required_field_tolerance=1e-7, samples_sha256=w.digest(sample_path)))
                w.save(args.run/"result.json", result)
                assert wall_error <= 1e-7, "Original exact CAD wall sample check failed; no native attempt"
                assert cap_maximum <= 1e-7, "Original exact CAD cap sample check failed; no native attempt"
                specimen_path = args.run/(original.NAME+"-specimen.json")
                w.save(specimen_path, dict(volume=volume, volume_reference=old["cad"]["reference_scope"],
                    volume_relative_error_estimate=old["cad"]["quadrature_estimated_relative_error"], wall_samples=len(samples["walls"]),
                    cap_samples=len(samples["caps"]), periodic_axes=[True]*3, declared_bounds=[[-.5]*3, [.5]*3],
                    primitive_receipt_sha256=w.digest(args.run/"primitive-receipt.json")))
                result.update(fixture_verified=True, specimen_sha256=w.digest(specimen_path),
                    fixture_check_scope="Accepted narrow cap ambiguity resolution plus original exact CAD surface/edge samples; finite checks, not a formal equivalence proof")
                w.save(args.run/"fixture-before-meshing.json", dict(samples_sha256=result["samples_sha256"],
                    specimen_sha256=result["specimen_sha256"], native_attempted=False, native_options=options))
                job.result("fixture", verified=True, wall_residual=wall_error, cap_field_max=cap_maximum)
            with job.part("native") as part:
                stage = "native"
                result["native_attempted"] = True
                w.save(args.run/"result.json", result)
                part.update(0, "Generating the original single native Tet10 recipe", cells=options["cells"])
                started = time.perf_counter(); mesh = meshers.generate(field, **options)
                result["native_generation_s"] = time.perf_counter()-started
                points, cells, surface = [np.asarray(getattr(mesh, key)) for key in ("points", "tetrahedra", "surface")]
                np.savez_compressed(args.run/"native-returned.npz", points=points, tetrahedra=cells, surface=surface)
                result["native"] = dict(nodes=len(points), elements=len(cells), order=2, qualified=False,
                    diagnostics=mesh.diagnostics, mesh_sha256=w.digest(args.run/"native-returned.npz"))
                w.save(args.run/"result.json", result)
                assert cells.shape[1] == 10 and surface.shape[1] == 6
                assert max(len(points), len(cells)) <= w.CAP, "Original node/tet cap rejected"
                part.update(1, "Native output archived before admission", nodes=len(points), tets=len(cells))
            with job.part("admission") as part:
                stage = "admission"
                from comparison_mesh_audit import observe
                from periodic_trace import periodic_pairs
                part.update(0, "Auditing unchanged geometry, quality and full periodic gates", nodes=len(points), tets=len(cells))
                case = scorecard.Case(original.NAME, field, (-np.inf, 0), volume)
                check = partial(periodic_pairs, declared_bounds=([-.5]*3, [.5]*3))
                geometry, periodic, quality = observe(args.run, original.NAME, case, points, cells, surface, True, periodic_auditor=check)
                result["native"].update(geometry=geometry, periodic=periodic, quality_quantiles=quality)
                w.save(args.run/"result.json", result)
                independent, _, _ = check(points, surface)
                match = all({tuple(map(int, pair)) for pair in a} == {tuple(map(int, pair)) for pair in b}
                    for a, b in zip(mesh.periodic_pairs, independent, strict=True))
                result["native"].update(native_and_independent_pairs_agree=match,
                    qualified=bool(geometry["geometry_admitted"] and match))
                result["native_admitted"] = result["native"]["qualified"]
                part.update(1, "Single native admission completed", admitted=result["native_admitted"])
            job.result("native", admitted=result["native_admitted"], nodes=len(points), tets=len(cells),
                boundary=geometry.get("sampled_two_sided_max"), relative_volume_error=geometry["relative_volume_error"], q_min=quality["min"])
        result["status"] = "native-admitted" if result["native_admitted"] else "native-admission-rejected"
    except Exception as error:
        result.update(status="control-stopped", failed_stage=stage, error_type=type(error).__name__, error=str(error), traceback=traceback.format_exc())
        if stage == "native" and "native_generation_s" not in result:
            result["failed_native_elapsed_s"] = time.perf_counter()-started
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    if result["fixture_verified"]:
        assert w.digest(args.run/(original.NAME+"-samples.npz")) == result["samples_sha256"]
        assert w.digest(args.run/(original.NAME+"-specimen.json")) == result["specimen_sha256"]
    result.update(complete=True, input_hashes_verified=True)
    w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in (("source", SOURCE), ("diagnosis", DIAGNOSIS), ("resolution", RESOLUTION)):
        parser.add_argument("--"+name, type=Path, default=default)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("source", "diagnosis", "resolution", "run"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import report, result
    report(0, "Continuing frozen BCC CAD surface and native control", eta_s=240, memory_ceiling_gib=8)
    receipt = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=240)
    w.save(args.run/"worker.json", receipt)
    assert receipt["returncode"] == 0, "Read worker.log and worker.json"
    row = json.loads((args.run/"result.json").read_text())
    result("bcc-verified-continuation", status=row["status"], fixture_verified=row["fixture_verified"],
        native_attempted=row["native_attempted"], admitted=row["native_admitted"], failed_stage=row.get("failed_stage"))


if __name__ == "__main__":
    main()
