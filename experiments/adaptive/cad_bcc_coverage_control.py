"""One public BCC CAD fixture with field checks and one native Tet10 attempt."""

import argparse
from functools import partial
import importlib.metadata
import inspect
import json
from pathlib import Path
import shutil
import sys
import time
import traceback
from types import SimpleNamespace

import microgen_workflow as w
import cad_lattice_fixture as fixture
import cad_union_comparison as previous

NAME = "cad-bcc-flat-struts-radius01"
OPTIONS = dict(previous.NATIVE, order=2, curved_fit=False)


def controls():
    import numpy as np
    from scipy.spatial.transform import Rotation
    import cad_coverage_preflight as crossed

    shape = crossed.make_shape()
    originals = fixture.primitives(shape)
    images, _ = fixture.periodic_images(originals, ([-.5]*3, [.5]*3))
    merged, merges, _ = fixture.coalesce(images)
    f = fixture.field(merged)
    rng = np.random.default_rng(20261111)
    points = rng.uniform(-.5, .5, (2048, 3))
    difference = float(abs(f(*points.T)-crossed.periodic_field(*points.T)).max())
    assert difference < 1e-13
    rotation = Rotation.from_rotvec([.2, -.3, .1]).as_matrix()
    rotated = [dict(row, center=(rotation@np.asarray(row["center"])).tolist(),
        axis=(rotation@np.asarray(row["axis"])).tolist(), rotation=(rotation@np.asarray(row["rotation"])).tolist(),
        endpoints=(np.asarray(row["endpoints"])@rotation.T).tolist()) for row in originals]
    raw = fixture.field(originals, preserve_rotations=True)
    transformed = fixture.field(rotated, preserve_rotations=True)
    rigid_error = float(abs(raw(*points.T)-transformed(*(points@rotation.T).T)).max())
    assert rigid_error < 1e-13
    halves = [dict(id=index, center=[center, 0., 0.], axis=[1., 0., 0.], rotation=np.eye(3).tolist(),
        height=.5, radius=.1, endpoints=[[center-.25, 0., 0.], [center+.25, 0., 0.]])
        for index, center in enumerate((-.25, .25))]
    chain, _, _ = fixture.coalesce(halves)
    assert len(chain) == 1 and chain[0]["height"] == 1.
    assert float(fixture.field(halves, preserve_rotations=True)(0., 0., 0.)) == 0.
    assert float(fixture.field(chain)(0., 0., 0.)) == -.1
    separated = [dict(row) for row in halves]
    separated[1] = dict(separated[1], center=[.25+1e-16, 0., 0.],
        endpoints=[[1e-16, 0., 0.], [.5+1e-16, 0., 0.]])
    try:
        fixture.coalesce(separated)
    except ValueError:
        gap_rejected = True
    else:
        raise AssertionError("Ambiguous positive cap gap was merged")
    return dict(crossed_fixture_maximum_field_difference=difference, rigid_rotation_maximum_difference=rigid_error,
        crossed_primitive_images=len(images), crossed_normalized_terms=len(merged),
        hidden_internal_cap_zero_removed=True, positive_gap_rejected=gap_rejected,
        scope="Field representation controls, not CAD volume or universal equivalence proofs")


def membership(cad, normalized, raw, former_caps, part, run):
    import numpy as np
    import microgen.cad
    from OCP.BRepClass3d import BRepClass3d_SolidClassifier
    from OCP.TopAbs import TopAbs_IN, TopAbs_ON
    from OCP.gp import gp_Pnt
    from itertools import product

    rng = np.random.default_rng(20261111)
    random = rng.uniform(-.5, .5, (512, 3))
    probes = [np.zeros(3)]
    epsilon = 1e-5
    for point in product((-.5, .5), repeat=3):
        probes.extend(np.asarray(point)+epsilon*np.asarray(delta) for delta in product((-1, 0, 1), repeat=3))
    seen = set()
    for cap in former_caps:
        center, axis = np.asarray(cap["center"]), np.asarray(cap["axis"])
        if abs(center).max() > .5+1e-12 or tuple(center) in seen:
            continue
        seen.add(tuple(center))
        transverse = np.cross(axis, np.eye(3)[int(abs(axis).argmin())])
        transverse /= np.linalg.norm(transverse)
        other = np.cross(axis, transverse)
        for lateral in (np.zeros(3), .5*cap["radius"]*transverse, -.5*cap["radius"]*transverse,
                        .5*cap["radius"]*other, -.5*cap["radius"]*other):
            probes.extend(center+lateral+sign*epsilon*axis for sign in (-1, 0, 1))
    grid = np.linspace(-.5, .5, 9)
    cap_trace_error = []
    for axis in range(3):
        low = np.zeros((81, 3)); low[:, axis] = -.5
        low[:, [i for i in range(3) if i != axis]] = np.array(list(product(grid, repeat=2)))
        high = low.copy(); high[:, axis] = .5
        cap_trace_error.append(float(abs(normalized(*low.T)-normalized(*high.T)).max()))
        probes.extend(low); probes.extend(high)
    probes = np.unique(np.vstack((random, probes)), axis=0)
    assert len(probes) <= 8000
    classifiers = [BRepClass3d_SolidClassifier(solid) for solid in microgen.cad.enumerate_solids(cad)]
    values = normalized(*probes.T)
    raw_values = raw(*probes.T)
    bounded = np.maximum(values, abs(probes).max(axis=1)-.5)
    raw_bounded = np.maximum(raw_values, abs(probes).max(axis=1)-.5)
    np.savez_compressed(run/"membership-probes-input.npz", points=probes, normalized_field=bounded,
        raw_field=raw_bounded, original_random_points=random)
    w.save(run/"membership-preflight.json", dict(periodic_scalar_trace_errors=cap_trace_error,
        probes=len(probes), input_sha256=w.digest(run/"membership-probes-input.npz")))
    assert max(cap_trace_error) < 1e-12, "Adapter scalar traces are not periodic"
    states, mismatches = [], []
    for index, point in enumerate(probes):
        states_at_point = []
        for classifier in classifiers:
            classifier.Perform(gp_Pnt(*point.tolist()), 1e-9)
            states_at_point.append(classifier.State())
        state = "inside" if TopAbs_IN in states_at_point else "boundary" if TopAbs_ON in states_at_point else "outside"
        states.append(state)
        good = bounded[index] < 0 if state == "inside" else abs(bounded[index]) <= 1e-8 if state == "boundary" else bounded[index] > 0
        if not good:
            mismatches.append(dict(index=index, point=point.tolist(), cad_state=state,
                normalized_field=float(bounded[index]), raw_field=float(raw_bounded[index])))
        if index % 128 == 0:
            part.update(index/len(probes), "Comparing exact CAD membership and cap junctions", probes=index, mismatches=len(mismatches))
    np.savez_compressed(run/"membership-probes.npz", points=probes, cad_states=np.asarray(states),
        normalized_field=bounded, raw_field=raw_bounded, original_random_points=random)
    inside = np.asarray(states) == "inside"
    row = dict(passed=not mismatches, random_points=512, total_points=len(probes), former_cap_centers=len(seen),
        solid_count=len(classifiers),
        mismatch_count=len(mismatches), mismatches=mismatches, periodic_scalar_trace_errors=cap_trace_error,
        raw_nonnegative_at_cad_interior=int((inside & (raw_bounded >= 0)).sum()), classifier_tolerance=1e-9,
        boundary_residual_tolerance=1e-8, scope="Finite independent checks; raw internal cap zeros are field representation defects, not CAD volume mismatches")
    part.update(1, "CAD membership comparison completed", probes=len(probes), mismatches=len(mismatches))
    return row


def worker(args):
    np, _, meshers, microgen, _, scorecard = w.imports()
    import microgen.cad
    import OCP
    import OCP.OCP
    from OCP.BRepCheck import BRepCheck_Analyzer
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    shape = microgen.BodyCenteredCubic(strut_radius=.1, cell_size=1., strut_joints=False)
    sources = {Path(__file__).resolve(), Path(fixture.__file__).resolve(), Path(previous.__file__).resolve(),
        Path(previous.fixture.__file__).resolve(), Path(w.__file__).resolve(), Path(microgen.cad.__file__),
        Path(inspect.getfile(type(shape))), Path(inspect.getfile(type(shape).__mro__[1])),
        Path(inspect.getfile(microgen.Cylinder)), Path(inspect.getfile(microgen.periodic_split_and_translate)),
        Path(inspect.getfile(microgen.Box)), Path(inspect.getfile(microgen.Phase)), Path(inspect.getfile(microgen.Rve)),
        Path(inspect.getfile(microgen.operations.fuse_shapes)), Path(meshers.__file__), Path(meshers._meshers.__file__),
        Path(OCP.OCP.__file__)}
    sources.update(Path(meshers.__file__).parent/name for name in ("compiler.py", "expression.py", "intersection.py"))
    for package in ("microgen", "cadquery-ocp-novtk", "numpy", "scipy"):
        distribution = importlib.metadata.distribution(package)
        sources.update(Path(distribution.locate_file(path)) for path in distribution.files or []
            if Path(path).name in ("METADATA", "RECORD"))
    sources.update(w.HERE/name for name in ("comparison_mesh_audit.py", "scorecard.py", "quadrature_audit.py",
        "quadratic_edge_removal_probe.py", "nested_refinement.py", "facet_incidence_audit.py", "periodic_trace.py", "box_projection.py"))
    hashes = {str(path): w.digest(path) for path in sorted(sources)}
    for index, path in enumerate(sorted(sources)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    result = dict(inputs=hashes, complete=False, fixture_verified=False, native_attempted=False, native_admitted=False,
        native_options=OPTIONS, primitive_parameters=dict(api="microgen.BodyCenteredCubic", radius=.1, cell_size=1.,
            strut_joints=False, strut_vertex_pairs=shape.strut_vertex_pairs.tolist(), vertices=shape.vertices.tolist()),
        numerical_threads=1, vtk_backend=vtkSMPTools.GetBackend(), production_change=False,
        fixture_check_scope="Endpoint-anchored constructive identity plus finite independent CAD checks; not a formal equivalence proof",
        actual_mmg_comparison=False, finite_element_solve_performed=False, comparative_timing_claim=False)
    versions = {name: importlib.metadata.version(name) for name in ("microgen", "cadquery-ocp-novtk", "numpy", "scipy")}
    assert versions["cadquery-ocp-novtk"] == previous.CAD_PIN
    w.save(args.run/"manifest.json", dict(inputs=hashes, versions=versions, worker_timeout_s=240, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, native_options=OPTIONS, maximum_nodes=w.CAP,
        maximum_tetrahedra=w.CAP, source_geometry="Public finite cylinders, periodic split/fuse and cell clip; no spheres or smoothing"))
    w.save(args.run/"result.json", result)
    stage = "controls"
    try:
        with Job(parts={"controls": 1, "cad": 3, "membership": 2, "samples": 3, "native": 3, "admission": 4}) as job:
            with job.part("controls") as part:
                part.update(0, "Checking exact primitive field controls", struts=shape.strut_number)
                result["controls"] = controls()
                original = fixture.primitives(shape)
                images, proof = fixture.periodic_images(original, ([-.5]*3, [.5]*3))
                w.save(args.run/"raw-primitive-receipt.json", dict(original_ordered_primitives=original,
                    ordered_images=images, image_range=proof))
                normalized, merges, former_caps = fixture.coalesce(images)
                f, raw = fixture.field(normalized), fixture.field(images, preserve_rotations=True)
                w.save(args.run/"primitive-receipt.json", dict(original_ordered_primitives=original, ordered_images=images,
                    image_range=proof, normalized_primitives=normalized, merges=merges, former_caps=former_caps,
                    scope="Constructive coalescing identity within recorded arithmetic bounds; finite CAD validation follows"))
                part.update(1, "Metadata and image ranges recorded", images=len(images), terms=len(normalized))
            job.result("controls", struts=shape.strut_number, images=len(images), normalized_terms=len(normalized))
            stage = "cad"
            with job.part("cad") as part:
                part.update(0, "Generating the actual public BCC CAD solid", struts=shape.strut_number)
                started = time.perf_counter(); cad = shape.generate_cad()
                result["cad_generation_s"] = time.perf_counter()-started
                assert BRepCheck_Analyzer(cad.wrapped).IsValid(), "Public BCC CAD is invalid"
                assert cad._mesh_volume is None
                props = GProp_GProps()
                estimate = float(BRepGProp.VolumeProperties_s(cad.wrapped, props, 1e-10, True, True))
                volume = float(abs(props.Mass()))
                default_volume = cad.volume()
                assert volume > 0 and np.isfinite(estimate) and estimate <= 1e-8
                np.testing.assert_allclose(volume, default_volume, rtol=1e-8, atol=1e-12)
                cad.export_step(str(args.run/"actual-bcc.step"))
                result["cad"] = dict(valid=True, volume=volume, default_volume=default_volume,
                    quadrature_requested_relative_tolerance=1e-10, quadrature_estimated_relative_error=estimate,
                    step_sha256=w.digest(args.run/"actual-bcc.step"),
                    reference_scope="OCCT adaptive volume integration of valid exact CAD surfaces; numerical estimate, not certified analytic volume")
                part.update(1, "Actual BCC CAD and volume recorded", volume=volume)
            stage = "membership"
            with job.part("membership") as part:
                result["membership"] = membership(cad, f, raw, former_caps, part, args.run)
                w.save(args.run/"result.json", result)
                assert result["membership"]["passed"], "Field/CAD membership mismatch; no native attempt"
            job.result("membership", passed=True, probes=result["membership"]["total_points"])
            stage = "samples"
            with job.part("samples") as part:
                samples, coverage = fixture.cad_samples(cad, part)
                wall_error = float(abs(f(*samples["walls"].T)).max())
                assert wall_error <= 1e-7, "Exact CAD wall samples disagree with normalized field"
                assert (f(*samples["caps"].T) <= 1e-7).all()
                np.savez_compressed(args.run/(NAME+"-samples.npz"), **samples)
                w.save(args.run/(NAME+"-specimen.json"), dict(volume=volume,
                    volume_reference=result["cad"]["reference_scope"], volume_relative_error_estimate=estimate,
                    wall_samples=len(samples["walls"]), cap_samples=len(samples["caps"]), periodic_axes=[True]*3,
                    declared_bounds=[[-.5]*3, [.5]*3], primitive_receipt_sha256=w.digest(args.run/"primitive-receipt.json")))
                result.update(fixture_verified=True, exact_sample_coverage=coverage, maximum_field_at_cad_wall=wall_error,
                    samples_sha256=w.digest(args.run/(NAME+"-samples.npz")), specimen_sha256=w.digest(args.run/(NAME+"-specimen.json")))
            stage = "native"
            with job.part("native") as part:
                result["native_attempted"] = True
                w.save(args.run/"result.json", result)
                part.update(0, "Generating the single native Tet10 attempt", cells=16)
                started = time.perf_counter(); mesh = meshers.generate(f, **OPTIONS)
                result["native_generation_s"] = time.perf_counter()-started
                points, cells, surface = (np.asarray(getattr(mesh, key)) for key in ("points", "tetrahedra", "surface"))
                np.savez_compressed(args.run/"native-returned.npz", points=points, tetrahedra=cells, surface=surface)
                result["native"] = dict(nodes=len(points), elements=len(cells), order=2, qualified=False,
                    diagnostics=mesh.diagnostics, mesh_sha256=w.digest(args.run/"native-returned.npz"))
                part.update(1, "Native output archived before admission", nodes=len(points), tets=len(cells))
                assert cells.shape[1] == 10 and surface.shape[1] == 6
                assert max(len(points), len(cells)) <= w.CAP, "Native node/tet cap rejected"
            stage = "admission"
            with job.part("admission") as part:
                from comparison_mesh_audit import observe
                from periodic_trace import periodic_pairs
                part.update(0, "Auditing all unchanged geometry and periodic gates", nodes=len(points), tets=len(cells))
                case = scorecard.Case(NAME, f, (-np.inf, 0), volume)
                check = partial(periodic_pairs, declared_bounds=([-.5]*3, [.5]*3))
                geometry, periodic, quality = observe(args.run, NAME, case, points, cells, surface, True, periodic_auditor=check)
                independent, _, _ = check(points, surface)
                match = all({tuple(map(int, pair)) for pair in a} == {tuple(map(int, pair)) for pair in b}
                    for a, b in zip(mesh.periodic_pairs, independent, strict=True))
                result["native"].update(geometry=geometry, periodic=periodic, quality_quantiles=quality,
                    native_and_independent_pairs_agree=match, qualified=bool(geometry["geometry_admitted"] and match))
                result["native_admitted"] = result["native"]["qualified"]
                part.update(1, "Native admission complete", admitted=result["native_admitted"])
            job.result("native-admission", admitted=result["native_admitted"], boundary=geometry.get("sampled_two_sided_max"),
                relative_volume_error=geometry["relative_volume_error"], q_min=quality["min"])
        result["status"] = "native-admitted" if result["native_admitted"] else "native-admission-rejected"
    except Exception as error:
        result.update(status="control-stopped", failed_stage=stage, error_type=type(error).__name__,
            error=str(error), traceback=traceback.format_exc())
        if stage == "native" and "native_generation_s" not in result:
            result["failed_native_elapsed_s"] = time.perf_counter()-started
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    result.update(complete=True, input_hashes_verified=True)
    w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(); args.run = args.run.resolve()
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import report, result
    report(0, "Starting bounded exact BCC fixture and native control", eta_s=240, memory_ceiling_gib=8)
    receipt = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=240)
    w.save(args.run/"worker.json", receipt)
    assert receipt["returncode"] == 0, "Read worker.log and worker.json"
    row = json.loads((args.run/"result.json").read_text())
    result("bcc-coverage", status=row["status"], fixture_verified=row["fixture_verified"],
        native_attempted=row["native_attempted"], admitted=row["native_admitted"], failed_stage=row.get("failed_stage"))


if __name__ == "__main__":
    main()
