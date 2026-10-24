"""Re-audit an unchanged rejected mesh with an explicit constrained projector."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import microgen_workflow as workflow


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    from nested_refinement import interpolation_check
    from hp_shared_domain import periodic_jumps

    previous = json.loads((args.source / "result.json").read_text())
    assert not previous["qualified"] and previous["stop_reason"] == "geometry-rejected"
    assert previous["geometry"]["projection_failed"] > 0
    assert workflow.geometry_is_admitted(dict(previous["geometry"], projection_failed=0))
    original = Path(previous["source"])
    paths = [Path(__file__), Path(workflow.__file__)]
    paths.extend(args.source / name for name in ("manifest.json", "result.json", "mesh.npz", "prolongation.npz"))
    paths.extend(original / name for name in ("result.json", "mesh.npz", "fields.npz"))
    paths.extend(workflow.HERE / name for name in ("box_projection.py", "test_box_projection.py",
        "test_projection_failure_probe.py", "nested_refinement.py", "hp_shared_domain.py", "scorecard.py", "facet_incidence_audit.py"))
    for suffix in ("-samples.npz", "-specimen.json"):
        path = args.source / (previous["case"]+suffix)
        assert workflow.digest(path) == workflow.digest(original / path.name)
        paths.append(path)
        shutil.copyfile(path, args.output.parent / path.name)
    inputs = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_box_projection.py",
        "experiments/adaptive/test_projection_failure_probe.py", "experiments/adaptive/test_facet_admission.py"], check=True, timeout=60)
    for name in ("mesh.npz", "prolongation.npz"):
        shutil.copyfile(args.source / name, args.output.parent / name)
        assert workflow.digest(args.source / name) == workflow.digest(args.output.parent / name)
    mesh = dict(np.load(args.output.parent / "mesh.npz"))
    transfer = dict(np.load(args.output.parent / "prolongation.npz"))
    root = np.load(original / "mesh.npz")
    fields = np.load(original / "fields.npz")["displacement"]
    reproduction = interpolation_check(root["points"], root["tetrahedra"], fields, {**mesh, **transfer})
    legacy = workflow.geometry(args.output.parent, previous["case"], mesh["points"], mesh["tetrahedra"], mesh["surface"])
    assert legacy["projection_failed"] == previous["geometry"]["projection_failed"] and not legacy["geometry_admitted"]
    geometry = workflow.geometry(args.output.parent, previous["case"], mesh["points"], mesh["tetrahedra"], mesh["surface"], box_projection_fallback=True)
    for key in ("sampled_two_sided_max", "relative_volume_error", "q_min_corner", "minimum_sampled_jacobian", "forward_samples", "reverse_samples"):
        np.testing.assert_allclose(geometry[key], legacy[key], atol=1e-13)
    periodic = periodic_jumps(mesh["points"], mesh["surface"], transfer["displacement"])
    workflow.save(args.output, dict(projection_reaudit=True, audited_source=str(args.source), source=str(original),
        source_rejection_preserved=True, box_projection_fallback=True, case=previous["case"], order=2,
        qualified=geometry["geometry_admitted"], status="nested-preflight-passed" if geometry["geometry_admitted"] else "geometry-rejected",
        stop_reason="auditor-updated-geometry-passed" if geometry["geometry_admitted"] else "geometry-rejected",
        total_rounds=previous["total_rounds"], nodes=len(mesh["points"]), elements=len(mesh["tetrahedra"]),
        geometry=geometry, legacy_geometry=legacy, periodic=periodic, interpolation=reproduction,
        original_pin=previous["original_pin"], relative_volume_change=previous["relative_volume_change"],
        equilibrium_solve_performed=False, reference_promoted=False, native_candidate_promoted=False,
        mesh_sha256=workflow.digest(args.output.parent / "mesh.npz"),
        prolongation_sha256=workflow.digest(args.output.parent / "prolongation.npz"),
        scope="Hash-identical rejected geometry, fields and samples re-audited with opt-in box-constrained feet. Historical default and rejection retained; thresholds unchanged. No refinement, equilibrium solve, reference or workflow promotion."))
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


if __name__ == "__main__":
    main()
