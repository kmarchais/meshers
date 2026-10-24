"""One bounded repeated preferred-star-only lattice preflight, with no solves."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", RAYON_NUM_THREADS="1")

import numpy as np
import microgen_workflow as workflow
from nested_preflight import fixture_inputs, curve_vertex_star_edges
from job_parts import worker_part

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job

MAX_ROUNDS = 8


def checkpoint_rounds(record):
    return record["total_rounds"] if "total_rounds" in record else record["completed_rounds"]


def worker(run, source, continued=None, box_projection_fallback=False):
    workflow.imports()
    from nested_refinement import refine, interpolation_check
    from periodic_crease_audit import features
    from hp_shared_domain import compose, periodic_jumps

    parent = json.loads((source / "result.json").read_text())
    assert parent["case"] == "smooth-lattice" and parent["qualified"] and parent["order"] == 2
    assert parent["geometry"]["geometry_admitted"] and len(parent["solver_records"]) == 6
    original = dict(np.load(source / "mesh.npz"))
    original_fields = np.load(source / "fields.npz")["displacement"]
    current, fields = original, original_fields
    root_parents = np.arange(len(original["tetrahedra"]))
    root_maps = np.repeat(np.eye(4)[None], len(root_parents), axis=0)
    offset = 0
    if continued:
        previous = json.loads((continued / "result.json").read_text())
        assert previous["qualified"] and (previous.get("preferred_only_repeated") or previous.get("projection_reaudit"))
        assert not previous["equilibrium_solve_performed"]
        assert Path(previous["source"]).resolve() == source.resolve()
        assert previous["stop_reason"] == "round-limit" or (previous.get("projection_reaudit") and previous["stop_reason"] == "auditor-updated-geometry-passed")
        assert bool(previous.get("box_projection_fallback", False)) == box_projection_fallback, "Change auditor only through a separate same-mesh re-audit"
        current = dict(np.load(continued / "mesh.npz"))
        transfer = np.load(continued / "prolongation.npz")
        fields, root_parents, root_maps = (transfer[key] for key in ("displacement", "parents", "parent_barycentric_maps"))
        interpolation_check(original["points"], original["tetrahedra"], original_fields,
            {**current, "displacement": fields, "parents": root_parents, "parent_barycentric_maps": root_maps})
        offset = checkpoint_rounds(previous)
    pin = int(np.linalg.norm(original["points"], axis=1).argmin())
    rounds, stop = [], "round-limit"
    for level in range(offset+1, offset+MAX_ROUNDS+1):
        feature = features(current)
        preferred = curve_vertex_star_edges(current["tetrahedra"], feature["curve_endpoints"])
        child = refine(current["points"], current["tetrahedra"], current["surface"], fields,
            preferred_edges=preferred, preferred_only=True)
        immediate = interpolation_check(current["points"], current["tetrahedra"], fields, child)
        selected = np.bincount(child["parents"], minlength=len(current["tetrahedra"])) > 1
        outside = int((selected & ~feature["groups"]["curve_vertex_incident"]).sum())
        assert outside == 0
        root_parents, root_maps = compose(root_parents, root_maps, child["parents"], child["parent_barycentric_maps"])
        composed = dict(child, parents=root_parents, parent_barycentric_maps=root_maps)
        reproduction = interpolation_check(original["points"], original["tetrahedra"], original_fields, composed)
        geometry = workflow.geometry(run, "smooth-lattice", child["points"], child["tetrahedra"], child["surface"],
            box_projection_fallback=box_projection_fallback)
        volume_change = geometry["volume"]/parent["geometry"]["volume"]-1
        assert abs(volume_change) < 1e-10
        periodic = periodic_jumps(child["points"], child["surface"], child["displacement"])
        np.testing.assert_array_equal(child["points"][pin], original["points"][pin])
        assert int(np.linalg.norm(child["points"], axis=1).argmin()) == pin
        directory = run / f"level-{level}"
        directory.mkdir()
        np.savez_compressed(directory / "mesh.npz", **{key: child[key] for key in ("points", "tetrahedra", "surface")})
        np.savez_compressed(directory / "prolongation.npz", **{key: child[key]
            for key in ("displacement", "parents", "parent_barycentric_maps")})
        row = dict(level=level, nodes=len(child["points"]), elements=len(child["tetrahedra"]),
            preferred_edges=len(preferred), source_curve_vertex_elements=int(feature["groups"]["curve_vertex_incident"].sum()),
            bisected_curve_vertex_elements=int(selected.sum()), noncurve_bisected_elements=outside,
            added_nodes=len(child["points"])-len(current["points"]),
            immediate_reproduction=immediate, original_reproduction=reproduction,
            geometry=geometry, relative_volume_change=volume_change, periodic=periodic,
            selection=child["diagnostics"], mesh_sha256=workflow.digest(directory / "mesh.npz"),
            prolongation_sha256=workflow.digest(directory / "prolongation.npz"))
        workflow.save(directory / "result.json", row)
        rounds.append(row)
        workflow.save(run / "checkpoint.json", dict(completed_rounds=len(rounds), total_rounds=level, rounds=rounds))
        print(f"Preferred-only round {level}: {len(child['points'])} nodes, {int(selected.sum())} curve parents split, geometry {geometry['geometry_admitted']}", flush=True)
        current = {key: child[key] for key in ("points", "tetrahedra", "surface")}
        fields = child["displacement"]
        if not geometry["geometry_admitted"]:
            stop = "geometry-rejected"
            break
        if not selected.any():
            stop = "no-admissible-preferred-star"
            break
        if len(child["points"]) == 250000 or len(child["tetrahedra"]) == workflow.CAP:
            stop = "count-cap"
            break
    np.savez_compressed(run / "mesh.npz", **current)
    np.savez_compressed(run / "prolongation.npz", displacement=fields, parents=root_parents, parent_barycentric_maps=root_maps)
    last = rounds[-1]
    workflow.save(run / "result.json", dict(case="smooth-lattice", source=str(source), order=2,
        preferred_only_repeated=True, qualified=last["geometry"]["geometry_admitted"],
        status="nested-preflight-passed" if last["geometry"]["geometry_admitted"] else "geometry-rejected",
        nodes=last["nodes"], elements=last["elements"], geometry=last["geometry"],
        relative_volume_change=last["relative_volume_change"], interpolation=last["original_reproduction"],
        rounds=rounds, completed_rounds=len(rounds), total_rounds=offset+len(rounds), maximum_rounds=MAX_ROUNDS, stop_reason=stop,
        continued_from=str(continued) if continued else None,
        box_projection_fallback=box_projection_fallback,
        original_pin=pin, periodic=last["periodic"], equilibrium_solve_performed=False,
        reference_promoted=False, native_candidate_promoted=False,
        scope="Repeated preferred curve-vertex-incident edge orbits only, no bulk fallback; unchanged original curved domain and included old displacement space. Saved-mesh geometry/resource preflight, no sizing sweep, accuracy or live workflow claim."))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--continue-from", type=Path, help="Passed round-limited geometric checkpoint; original fields remain unchanged")
    parser.add_argument("--box-projection-fallback", action="store_true", help="Explicit constrained target method, with historical default retained")
    args = parser.parse_args()
    if args.worker:
        worker(args.run, args.source, args.continue_from, args.box_projection_fallback)
        return
    fixtures = fixture_inputs(args.source)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    for path in fixtures:
        shutil.copyfile(path, args.run / path.name)
    paths = [Path(__file__), *fixtures]
    paths.extend(args.source / name for name in ("manifest.json", "result.json", "mesh.npz", "fields.npz"))
    paths.extend(workflow.HERE / name for name in ("nested_refinement.py", "nested_preflight.py", "hp_shared_domain.py",
        "periodic_crease_audit.py", "facet_incidence_audit.py", "scorecard.py", "microgen_workflow.py",
        "test_nested_refinement.py", "job_parts.py", "test_job_parts.py", "test_bounded.py"))
    paths.append(workflow.HERE / "test_nested_star_budget.py")
    if args.continue_from:
        paths.extend(args.continue_from / name for name in ("manifest.json", "result.json", "mesh.npz", "prolongation.npz"))
        prior = json.loads((args.continue_from / "manifest.json").read_text())["inputs"]
        if (args.continue_from / "analysis-inputs.json").exists():
            prior.update(json.loads((args.continue_from / "analysis-inputs.json").read_text()))
        assert all(workflow.digest(args.source / name) == prior[str(args.source / name)] for name in ("result.json", "mesh.npz", "fields.npz"))
    if args.box_projection_fallback:
        paths.extend(workflow.HERE / name for name in ("box_projection.py", "test_box_projection.py", "test_projection_failure_probe.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    (args.run / "requirements.txt").write_text(subprocess.check_output(["uv", "pip", "freeze", "--python", sys.executable], text=True))
    workflow.save(args.run / "manifest.json", dict(inputs=hashes, platform=platform.platform(), architecture=platform.machine(),
        numerical_thread_policy=workflow.NUMERICAL_THREAD_POLICY, maximum_rounds=MAX_ROUNDS,
        memory_bytes=workflow.MEMORY, maximum_nodes=250000, maximum_tets=workflow.CAP,
        stage_timeouts_s=dict(controls=120, preflight=600), no_live_workflow_timing=True,
        continue_from=str(args.continue_from) if args.continue_from else None,
        box_projection_fallback=args.box_projection_fallback,
        packages={name: importlib.metadata.version(name) for name in ("numpy", "scipy", "fedoo", "meshers")}))
    with Job(parts=["controls", "preferred-only-preflight"]) as job:
        for stage in ("controls", "preferred-only-preflight"):
            with worker_part(job, stage) as part:
                command = ([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_nested_refinement.py",
                    "experiments/adaptive/test_job_parts.py", "experiments/adaptive/test_bounded.py", "experiments/adaptive/test_nested_star_budget.py"] if stage == "controls"
                    else [sys.executable, __file__, "--run", str(args.run), "--source", str(args.source), "--worker"])
                if stage != "controls" and args.continue_from:
                    command.extend(["--continue-from", str(args.continue_from)])
                if stage != "controls" and args.box_projection_fallback:
                    command.append("--box-projection-fallback")
                record = workflow.bounded(command, args.run / (stage+".log"), part, timeout=120 if stage == "controls" else 600)
                workflow.save(args.run / (stage+"-worker.json"), record)
                assert record["returncode"] == 0, f"Read {stage}.log"
                job.result(stage, seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"]/1024**3)
    assert all(workflow.digest(Path(path)) == expected for path, expected in hashes.items())


if __name__ == "__main__":
    main()
