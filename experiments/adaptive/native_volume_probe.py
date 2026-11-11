"""Private native volume-correction payoff experiment on nine frozen Tet10 meshes.

Build and study are separate parent-owned queued stages. Direct Rust fixture
evaluators use analytic gradients; measured stage costs are not production JIT
workflow timings. All geometry and FEA admission policies remain unchanged.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import microgen_workflow as w


RUNS = w.HERE/"runs"
SOURCE = RUNS/"macbook-native-face-fit-nine-v1"
SELECTION = w.ROOT/"handoff/macbook/face-fit-selection-20261110.json"
TEST = "quadratic::fit::volume_probe::volume_bias_payoff"
CONTROL = "quadratic::fit::volume_probe::volume_probe_controls"
FLOAT_ENCODING = "ieee754-binary64-u64"
TRANSPORT_EVIDENCE = RUNS/"macbook-pmy-saved-trace-replay-v1/result.json"
CASES = ("gyroid", "schwarz-solid", "spinodoid", "smooth-lattice", "schwarz-d-solid",
         "graded-gyroid-linear", "graded-gyroid-periodic", "spinodoid-anisotropic",
         "graded-gyroid-density-linear")


def encode_floats(values):
    import numpy as np
    values = np.asarray(values, dtype=np.float64)
    assert np.isfinite(values).all(), "Nonfinite scientific input"
    return values.view(np.uint64).tolist()


def decode_floats(bits):
    import numpy as np
    # Mixed small integers and values above i64::MAX must never pass through
    # NumPy's inferred dtype, which can choose float64 and lose integer bits.
    raw = np.asarray(bits, dtype=object)
    assert all(isinstance(value, (int, np.integer)) and 0 <= value <= 0xffffffffffffffff
               for value in raw.flat), "Expected unsigned integer bit patterns"
    values = np.asarray(bits, dtype=np.uint64).view(np.float64)
    assert np.isfinite(values).all(), "Nonfinite scientific output"
    return values


def transport_controls():
    import numpy as np
    evidence = json.loads(TRANSPORT_EVIDENCE.read_text())["native_input_parse"]["example_differences"][0]
    assert evidence["node"] == 9009 and evidence["axis"] == 0
    assert evidence["python_bits"] == 13821171956431104683
    assert evidence["rust_bits"] == 13821171956431104684
    bits = np.array([0, 0x8000000000000000, 1, 0x000fffffffffffff,
                     0x0010000000000000, 0x7fefffffffffffff, 0xffefffffffffffff,
                     evidence["python_bits"]], dtype=np.uint64)
    payload = json.loads(json.dumps(encode_floats(bits.view(np.float64))))
    np.testing.assert_array_equal(decode_floats(payload).view(np.uint64), bits)
    return dict(encoding=FLOAT_ENCODING, python_roundtrip_bits_exact=True,
        evidence_source=str(TRANSPORT_EVIDENCE), evidence_sha256=w.digest(TRANSPORT_EVIDENCE),
        known_parser_difference=evidence, controls="positive/negative zero, subnormal endpoints, minimum normal and signed finite extremes")


def run_child(command, log, *, env=None, timeout=300):
    """Children stay in the outer bounded worker's process group.

    The outer RSS watchdog measures the whole family and kills the whole group
    on its900s/8GiB limit. A child has its own300s wall timeout.
    """
    start = time.perf_counter()
    with log.open("w") as stream:
        try:
            done = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                  env=env, timeout=timeout, check=False)
            code, failure = done.returncode, None
        except subprocess.TimeoutExpired:
            code, failure = -1, "child-timeout"
    record = dict(command=command, returncode=code, failure=failure,
                  elapsed_s=time.perf_counter()-start, timeout_s=timeout,
                  memory_enforcement="ancestor process-tree RSS watchdog, 8GiB")
    w.save(log.with_suffix(".json"), record)
    assert code == 0, f"Read {log}"
    return record


def build(directory):
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    os.environ.update(CARGO_BUILD_JOBS="2", RUSTUP_TOOLCHAIN="1.95.0")
    transport = transport_controls()
    sources = [w.ROOT/"Cargo.toml", w.ROOT/"Cargo.lock", w.ROOT/"crates/meshers-core/Cargo.toml",
               *sorted((w.ROOT/"crates/meshers-core/src").rglob("*.rs")), Path(__file__), TRANSPORT_EVIDENCE]
    hashes = {str(path.resolve()): w.digest(path) for path in sources}
    for path in sources:
        destination = directory/"sources"/path.relative_to(w.ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
    command = ["cargo", "+1.95.0", "test", "--locked", "--release", "-p", "meshers-core",
               "--lib", "--features", "experimental-quads", "--no-run", "--message-format=json"]
    run_child(command, directory/"cargo.log")
    artifacts = []
    for line in (directory/"cargo.log").read_text().splitlines():
        try:
            entry = json.loads(line)
        except json.JSONDecodeError:
            continue
        if entry.get("reason") == "compiler-artifact" and entry.get("executable"):
            if entry["target"]["name"] == "meshers_core" and entry["profile"]["test"]:
                artifacts.append(Path(entry["executable"]))
    assert len(artifacts) == 1, artifacts
    executable = directory/"volume-probe-tests"
    shutil.copy2(artifacts[0], executable)
    info = dict(executable=str(executable), executable_sha256=w.digest(executable),
        inputs=hashes, rust=subprocess.check_output(["rustc", "+1.95.0", "--version"], text=True).strip(),
        mode="release private core unit-test executable", build_jobs=2,
        production_defaults_changed=False, compiled_python_wheel_changed=False, float_transport=transport)
    w.save(directory/"binary.json", info)
    run_child([str(executable), "--exact", CONTROL, "--nocapture", "--test-threads=1"], directory/"controls.log")
    run_child([str(executable), "curved_fit", "--nocapture", "--test-threads=1"], directory/"existing-fit.log")
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    print(f"Private native probe built and controls passed: {directory/'binary.json'}")


def prepare(name, directory, entry):
    import numpy as np
    import comparison_fixtures as fixtures
    from gradient_tie_probe import symmetric_lattice
    import meshers

    source = SOURCE/(name+"-repeat-0-fit")
    original = SOURCE/(name+"-repeat-0-baseline")
    files = [source/"mesh.npz", original/"mesh.npz", source/"result.json", original/"result.json",
             source/(name+"-samples.npz"), source/(name+"-specimen.json")]
    saved = json.loads((source/"result.json").read_text())
    baseline = json.loads((original/"result.json").read_text())
    assert w.digest(source/"mesh.npz") == saved["mesh_sha256"]
    assert w.digest(original/"mesh.npz") == baseline["mesh_sha256"]
    with np.load(source/"mesh.npz") as archive, np.load(original/"mesh.npz") as old:
        points, cells, surface = (archive[key].copy() for key in ("points", "tetrahedra", "surface"))
        originals = old["points"].copy()
        np.testing.assert_array_equal(cells, old["tetrahedra"])
        np.testing.assert_array_equal(surface, old["surface"])
    assert points.shape == originals.shape and cells.shape[1] == 10
    assert max(len(points), len(cells)) <= w.CAP
    np.testing.assert_array_equal(points[np.unique(cells[:, :4])], originals[np.unique(cells[:, :4])])
    fixture_root = Path(entry["fixtures"]).resolve()
    target = fixtures.case(name, fixture_root)
    field = target.function
    if name == "smooth-lattice":
        def field(x, y, z):
            return np.minimum(symmetric_lattice(x, y, z), .02)
    compiled = meshers.compile_field(field)
    assert compiled.has_gradient, "Analytic fixture identity requires the current compiled gradient"
    with np.load(source/(name+"-samples.npz")) as samples:
        arrays = [samples[key][::max(1, len(samples[key])//128)] for key in ("walls", "caps")]
    arrays.append(points[::max(1, len(points)//128)])
    if name == "smooth-lattice":
        # Same six tie controls as gradient_tie_probe.py, not specimen samples.
        arrays.append(np.array([[0., .17, -.17], [0., -.17448645899797485, .1744864589979749],
            [.17, 0., -.17], [-.17, .17, 0.], [0., -.17, -.17], [.17, 0., .17]]))
    probes = np.concatenate(arrays)
    values = compiled.values(probes)
    gradients = np.array([compiled.gradient(*point) for point in probes])
    assert len(probes) > 0 and np.isfinite(values).all() and np.isfinite(gradients).all()
    np.testing.assert_allclose(values, field(*probes.T), atol=1e-11, rtol=1e-11)
    fixture = dict(name=name, periodic=[fixtures.periodic(name)]*3)
    if name.startswith("spinodoid"):
        path = fixtures.coefficients_path(name, fixture_root)
        frozen = json.loads(path.read_text())
        fixture.update(coefficients_bits=encode_floats(frozen["coefficients"]), threshold_bits=encode_floats(frozen["threshold"]))
        files.append(path)
    if name == "graded-gyroid-density-linear":
        path = fixture_root/"linear-density-calibration.json"
        frozen = json.loads(path.read_text())
        fixture.update(knots_bits=encode_floats(frozen["z"]), widths_bits=encode_floats(frozen["half_width"]))
        files.append(path)
    directory.mkdir()
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(source/(name+suffix), directory/(name+suffix))
    w.save(directory/"native-input.json", dict(float_encoding=FLOAT_ENCODING, fixture=fixture,
        points_bits=encode_floats(points), original_bits=encode_floats(originals),
        tetrahedra=cells.tolist(), surface=surface.tolist(), probes_bits=encode_floats(probes),
        values_bits=encode_floats(values), gradients_bits=encode_floats(gradients)))
    record = dict(case=name, source=str(source), original=str(original), fixtures=str(fixture_root),
        inputs={str(path.resolve()): w.digest(path) for path in files},
        native_input_sha256=w.digest(directory/"native-input.json"),
        field_identity_probes=len(probes), original_specimen_samples_unchanged=True,
        float_encoding=FLOAT_ENCODING, transport_scope="Every scientific input float is an exact u64 bit pattern; metadata and integer connectivity retain ordinary JSON types")
    w.save(directory/"preparation.json", record)
    return record


def audit(directory):
    np, *_ = w.imports()
    import comparison_fixtures as fixtures
    from comparison_mesh_audit import observe
    from periodic_trace import periodic_pairs
    from volume_fidelity_control import barycentric
    from gyroid_volume_rankings import comparison
    import scorecard
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    prepared = json.loads((directory/"preparation.json").read_text())
    name = prepared["case"]
    source = Path(prepared["source"])
    native = json.loads((directory/"native-output.json").read_text())
    assert native["float_encoding"] == FLOAT_ENCODING
    after = decode_floats(native.pop("points_bits"))
    with np.load(source/"mesh.npz") as old:
        before, cells, surface = (old[key].copy() for key in ("points", "tetrahedra", "surface"))
    assert before.shape == after.shape and np.isfinite(after).all()
    with np.load(Path(prepared["original"])/"mesh.npz") as old:
        original = old["points"].copy()
    corners = np.unique(cells[:, :4])
    np.testing.assert_array_equal(before[corners].view(np.uint64), after[corners].view(np.uint64))
    locked = abs(before) == .5
    np.testing.assert_array_equal(before[locked].view(np.uint64), after[locked].view(np.uint64))
    walls = surface[~np.any(np.all(abs(before[surface[:, :3]]+.5)<1e-10, axis=1)
                            |np.all(abs(before[surface[:, :3]]-.5)<1e-10, axis=1), axis=1)]
    limits = np.full(len(before), np.inf)
    lengths = np.linalg.norm(original[walls[:, [1, 2, 0]]]-original[walls[:, :3]], axis=2)
    np.minimum.at(limits, walls[:, 3:].ravel(), .02*lengths.ravel())
    assert np.all(np.linalg.norm(after-original, axis=1) <= limits+1e-12)
    fixed = np.ones(len(before), dtype=bool)
    fixed[np.unique(walls[:, 3:])] = False
    np.testing.assert_array_equal(before[fixed].view(np.uint64), after[fixed].view(np.uint64))
    np.savez_compressed(directory/"mesh.npz", points=after, tetrahedra=cells, surface=surface)
    target = fixtures.case(name, Path(prepared["fixtures"]))
    required = fixtures.periodic(name)
    part = Job().part("audit-"+name)
    standard, dense = [], []
    for i, points in enumerate((before, after)):
        part.update(i/4, "independent standard geometry audit", force=True, points=len(before), tets=len(cells))
        standard.append(observe(directory, name, target, points, cells, surface, required, periodic_auditor=periodic_pairs))
    original_bary = scorecard.BARY
    scorecard.BARY = np.unique(np.vstack((original_bary, barycentric(13))), axis=0)
    try:
        for i, points in enumerate((before, after)):
            part.update((i+2)/4, "independent dense geometry audit", force=True, points=len(before), tets=len(cells))
            dense.append(observe(directory, name, target, points, cells, surface, required, periodic_auditor=periodic_pairs))
    finally:
        scorecard.BARY = original_bary
    saved = json.loads((source/"result.json").read_text())
    for key in ("volume", "sampled_two_sided_max", "relative_volume_error"):
        np.testing.assert_allclose(standard[0][0][key], saved["geometry"][key], atol=1e-12, rtol=1e-12)
    np.testing.assert_array_equal(list(standard[0][2].values()), list(standard[1][2].values()))
    reference = json.loads((directory/(name+"-specimen.json")).read_text())
    evidence = dict(frozen_reference=reference, frozen_scores_unchanged=True)
    rankings = json.loads((RUNS/"macbook-gyroid-volume-rankings-v1/result.json").read_text())
    if name in [row["case"] for row in rankings["cases"]]:
        row = next(row for row in rankings["cases"] if row["case"] == name)
        lo, hi = row["empirical_reference_interval"]
        estimate = row["estimates"]["gauss128"]
        evidence.update(gauss128_estimate=estimate, empirical_reference_interval=[lo, hi],
            before_absolute_error=abs(standard[0][0]["volume"]-estimate),
            after_absolute_error=abs(standard[1][0]["volume"]-estimate),
            empirical_ranking=comparison(standard[0][0]["volume"], standard[1][0]["volume"], lo, hi),
            interval_scope="Empirical sensitivity interval, not a certified bound or promoted reference")
    elif name in ("schwarz-solid", "schwarz-d-solid"):
        evidence.update(exact_reference=.5, before_absolute_error=abs(standard[0][0]["volume"]-.5),
                        after_absolute_error=abs(standard[1][0]["volume"]-.5))
    repeats = native["repeats"]
    timing = {key: dict(median=float(np.median([r["timing"][key] for r in repeats])),
                       values=[r["timing"][key] for r in repeats]) for key in repeats[0]["timing"]}
    report = dict(case=name, native=native, timing=timing, nodes=len(before), elements=len(cells),
        before=dict(geometry=standard[0][0], periodic=standard[0][1]),
        after=dict(geometry=standard[1][0], periodic=standard[1][1]),
        dense_before=dense[0][0], dense_after=dense[1][0], volume_evidence=evidence,
        admitted=standard[1][0]["geometry_admitted"], dense_admitted=dense[1][0]["geometry_admitted"],
        standard_boundary_nonworse=standard[1][0]["sampled_two_sided_max"] <= standard[0][0]["sampled_two_sided_max"]+1e-12,
        dense_boundary_nonworse=dense[1][0]["sampled_two_sided_max"] <= dense[0][0]["sampled_two_sided_max"]+1e-12,
        counts_corners_caps_corner_quality_unchanged=True, mesh_sha256=w.digest(directory/"mesh.npz"),
        scope="Native correction stage with direct Rust analytic evaluators; not production JIT generation timing. No FEA solve or inherited stiffness score. Added native cap-distance and Jacobian-ratio guards may change the candidate versus Python.")
    assert all(w.digest(Path(path)) == digest for path, digest in prepared["inputs"].items())
    w.save(directory/"result.json", report)
    part.update(1., "independent audits complete", force=True, admitted=report["admitted"])


def study(args):
    w.imports()
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    info = json.loads(args.build_info.read_text())
    assert w.digest(Path(info["executable"])) == info["executable_sha256"]
    selection = json.loads(SELECTION.read_text())
    entries = {row["case"]: row for row in selection["cases"]}
    paths = [Path(__file__), args.build_info, SELECTION, Path(info["executable"]),
             RUNS/"macbook-gyroid-volume-rankings-v1/result.json", TRANSPORT_EVIDENCE]
    paths += [w.HERE/name for name in ("comparison_fixtures.py", "comparison_mesh_audit.py", "volume_fidelity_control.py",
        "gyroid_volume_rankings.py", "gradient_tie_probe.py", "scorecard.py", "microgen_workflow.py", "periodic_trace.py",
        "box_projection.py", "facet_incidence_audit.py", "quadratic_edge_removal_probe.py", "linear_density_fixture.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in paths}
    plan = {"prepare": 2}
    for name in CASES:
        plan["native-"+name] = 1
        plan["audit-"+name] = 5
    rows = []
    with Job(parts=plan) as job:
        with job.part("prepare") as part:
            transport = transport_controls()
            prepared = []
            for i, name in enumerate(CASES):
                part.update(i/len(CASES), "preparing frozen arrays and compiled field controls", force=True, cases=i)
                prepared.append(prepare(name, args.run/name, entries[name]))
            part.update(1., "all nine frozen inputs prepared", force=True, cases=9)
        w.save(args.run/"manifest.json", dict(inputs=hashes, build=info, preparation=prepared,
            numerical_threads=w.NUMERICAL_THREAD_POLICY, whole_timeout_s=900, child_timeout_s=300,
            memory_bytes=w.MEMORY, warmups=1, timed_repetitions=3, reference_promoted=False, float_transport=transport,
            method="Signed-envelope field-only bias step; no target volume. Native cap-constrained distance maximum and min/max Jacobian ratio guards supplement the earlier Python prototype."))
        for name in CASES:
            directory = args.run/name
            with job.part("native-"+name) as part:
                part.update(0., "starting private native payoff test", force=True, repeats=4, threads=1)
                env = dict(os.environ, MESHERS_VOLUME_INPUT=str(directory/"native-input.json"),
                    MESHERS_VOLUME_OUTPUT=str(directory/"native-output.json"), MESHERS_VOLUME_PROGRESS_PART="native-"+name)
                native_worker = run_child([info["executable"], "--ignored", "--exact", TEST,
                    "--nocapture", "--test-threads=1"], directory/"native.log", env=env)
            with job.part("audit-"+name) as part:
                part.update(0., "starting independent geometry audits", force=True, threads=1)
                audit_worker = run_child([sys.executable, str(Path(__file__).resolve()), "--audit", str(directory)], directory/"audit.log")
            row = json.loads((directory/"result.json").read_text())
            row.update(native_worker=native_worker, audit_worker=audit_worker)
            rows.append(row)
            w.save(args.run/"result.json", dict(complete=len(rows)==9, rows=rows, inputs=hashes,
                production_defaults_changed=False, production_workflow_speed_claim=False, reference_promoted=False))
            job.result(name, native_stage_s=row["timing"]["total_s"]["median"], admitted=row["admitted"],
                       dense_admitted=row["dense_admitted"], boundary_nonworse=row["standard_boundary_nonworse"] and row["dense_boundary_nonworse"])
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path)
    parser.add_argument("--build-info", type=Path)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--study-worker", action="store_true")
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    if args.build:
        return build(args.build)
    if args.audit:
        return audit(args.audit.resolve())
    assert args.run is not None and args.build_info is not None
    args.run, args.build_info = args.run.resolve(), args.build_info.resolve()
    if args.study_worker:
        return study(args)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run/"input-native_volume_probe.py")
    command = [sys.executable, str(Path(__file__).resolve()), "--study-worker", "--run", str(args.run), "--build-info", str(args.build_info)]
    record = w.bounded(command, args.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=900)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, f"Read {args.run/'worker.log'}"
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import result
    saved = json.loads((args.run/"result.json").read_text())
    for row in saved["rows"]:
        result(row["case"]+" private native volume correction", native_stage_s=row["timing"]["total_s"]["median"],
            admitted=row["admitted"], dense_admitted=row["dense_admitted"],
            boundary_nonworse=row["standard_boundary_nonworse"] and row["dense_boundary_nonworse"],
            jsonl=str(args.run/"reported-results.jsonl"))


if __name__ == "__main__":
    main()
