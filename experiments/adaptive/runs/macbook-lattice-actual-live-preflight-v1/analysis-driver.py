"""Fresh actual lattice route with explicitly charged repair and star refinement."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

import microgen_workflow as workflow


def prepare(source, output):
    from microgen_tet10_control import straight_elevation
    from nested_preflight import curve_vertex_star_edges
    from nested_refinement import refine, interpolation_check
    from periodic_crease_audit import features
    from homogenization import strain_tensor
    from hp_shared_domain import periodic_jumps
    from interior_relocation_trial import validate_repair_source

    repair = json.loads((source / "result.json").read_text())
    generation = Path(repair["source"])
    parent = json.loads((generation / "result.json").read_text())
    assert repair["case"] == "smooth-lattice" and repair["qualified"]
    assert repair["repair"]["remaining_bad_elements"] == 0 and repair["repair"]["boundary_coordinates_unchanged"]
    assert parent["stage_audit_enabled"] and not parent.get("saved_seed_replay", False)
    validate_repair_source(parent)
    assert parent["mmg_calls"] == repair["parent_mmg_calls"]
    assert workflow.digest(workflow.HERE / "microgen_tuned.py") == json.loads((generation / "manifest.json").read_text())["script_sha256"]
    paths = [Path(__file__)]
    paths.extend(workflow.HERE / name for name in ("microgen_tet10_control.py", "nested_preflight.py",
        "nested_refinement.py", "periodic_crease_audit.py", "hp_shared_domain.py", "test_straight_elevation.py",
        "test_nested_refinement.py", "microgen_tuned.py", "interior_relocation_trial.py"))
    paths.extend(directory / name for directory in (source, generation) for name in ("manifest.json", "result.json", "mesh.npz"))
    inputs = {}
    for directory in (source, generation):
        for name, expected in json.loads((directory / "manifest.json").read_text())["inputs"].items():
            path = Path(name)
            assert workflow.digest(path) == expected, str(path)
            inputs[str(path)] = expected
    for name in ("vtk-seed.json", "mmg-stage-1.json", "mmg-stage-2.json"):
        path = generation / name
        stage = json.loads(path.read_text())
        assert 0 < stage["elements"] <= workflow.CAP and stage["nodes"] <= 250000
        paths.append(path)
    for suffix in ("-specimen.json", "-samples.npz"):
        path = source / ("smooth-lattice" + suffix)
        assert workflow.digest(path) == inputs[str(generation / path.name)]
        paths.append(path)
        shutil.copyfile(path, output.parent / path.name)
    import microgen.remesh
    paths.extend((Path(microgen.remesh.__file__),
        Path(importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin/mmg3d_O3"))))
    inputs.update({str(path): workflow.digest(path) for path in paths})
    workflow.save(output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, output.parent / "analysis-driver.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_straight_elevation.py",
        "experiments/adaptive/test_nested_refinement.py"], check=True, timeout=60)
    start = time.perf_counter()
    original = np.load(source / "mesh.npz")
    elevation = straight_elevation(original["points"], original["tetrahedra"], original["surface"])
    if not elevation["within_node_cap"]:
        workflow.save(output, dict(qualified=False, status="straight-lift-node-cap", nodes=elevation["nodes"], source=str(source)))
        return
    points, cells, surface = (elevation[name] for name in ("points", "tetrahedra", "surface"))
    from periodic import periodic_pairs
    from meshers import Mesh
    matching, _, _ = periodic_pairs(points, surface)
    converted = Mesh(points, cells, surface, np.zeros(len(surface), dtype=np.uint8), matching,
        {"source": "Fresh actual VTK51/two-pass MMG, explicit fixed-boundary repair and straight elevation"}).to_fedoo()
    np.testing.assert_array_equal(converted.nodes, points)
    np.testing.assert_array_equal(converted.elements, cells)
    assert converted.elm_type == "tet10"
    fields = np.array([points @ strain_tensor(load).T for load in np.eye(6)])
    feature = features(dict(points=points, tetrahedra=cells, surface=surface))
    lift_s = time.perf_counter() - start
    phase = time.perf_counter()
    refined = refine(points, cells, surface, fields,
        preferred_edges=curve_vertex_star_edges(cells, feature["curve_endpoints"]))
    reproduction = interpolation_check(points, cells, fields, refined)
    refinement_s = time.perf_counter() - phase
    np.savez_compressed(output.parent / "elevated-mesh.npz", points=points, tetrahedra=cells, surface=surface)
    np.savez_compressed(output.parent / "mesh.npz", **{key: refined[key] for key in ("points", "tetrahedra", "surface")})
    np.savez_compressed(output.parent / "prolongation.npz", **{key: refined[key]
        for key in ("displacement", "parents", "parent_barycentric_maps")})
    geometry = workflow.geometry(output.parent, "smooth-lattice", refined["points"], refined["tetrahedra"], refined["surface"])
    periodic = periodic_jumps(refined["points"], refined["surface"], refined["displacement"])
    volume_change = geometry["volume"] / repair["geometry"]["volume"] - 1
    assert abs(volume_change) < 1e-10
    prepared_s = time.perf_counter() - start
    workflow.save(output, dict(case="smooth-lattice", route=repair["route"]+"-straight-tet10-curve-star",
        source=str(source), generation_source=str(generation), nodes=len(refined["points"]), elements=len(refined["tetrahedra"]),
        order=2, status="fresh-actual-preflight-passed" if geometry["geometry_admitted"] else "geometry-rejected",
        qualified=geometry["geometry_admitted"], geometry=geometry, periodic=periodic, reproduction=reproduction,
        refinement=refined["diagnostics"], relative_volume_change=volume_change,
        lifted_nodes=len(points), lift_and_conversion_s=lift_s, refinement_s=refinement_s, total_preparation_s=prepared_s,
        generation_and_audits_s=parent["total_workflow_s"], repair_and_audits_s=repair["total_saved_mesh_study_s"],
        charged_preparation_workflow_s=parent["total_workflow_s"]+repair["total_saved_mesh_study_s"]+prepared_s,
        fresh_actual_generation=True, saved_seed_replay=False, timing_comparable=False,
        extra_steps=["SLSQP interior repair", "manual straight Tet10/direct fedoo", "curve-vertex-incident star refinement"],
        control_fields="Six affine macrostrain fields for interpolation and periodic controls; no parent FE solve or equilibrium-field replay",
        field_solve_performed=False, reference_promoted=False,
        scope="Fresh Mac actual route through geometry preflight. Diagnostic intermediate audits are charged, so no production timing ratio; no stiffness or continuum claim."))
    assert all(workflow.digest(Path(path)) == digest for path, digest in inputs.items())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    workflow.imports()
    prepare(args.source, args.output)


if __name__ == "__main__":
    main()
