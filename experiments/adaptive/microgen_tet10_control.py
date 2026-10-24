"""Straight Tet10 and direct-fedoo control on an admitted genuine Microgen mesh."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import time

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def straight_elevation(points, tetrahedra, surface, node_cap=250000):
    """Count before allocating, then preserve the exact affine parent domain."""
    import numpy as np

    assert tetrahedra.shape[1] == 4 and surface.shape[1] == 3
    edges = np.array(((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)))
    unique, inverse = np.unique(np.sort(tetrahedra[:, edges].reshape(-1, 2), axis=1), axis=0, return_inverse=True)
    nodes = len(points) + len(unique)
    result = {"nodes": nodes, "within_node_cap": nodes <= node_cap}
    if nodes > node_cap:
        return result
    elevated = np.column_stack((tetrahedra, inverse.reshape(-1, 6) + len(points)))
    keys = unique[:, 0] * len(points) + unique[:, 1]
    surface_edges = np.sort(surface[:, edges[:3]].reshape(-1, 2), axis=1)
    surface_keys = surface_edges[:, 0] * len(points) + surface_edges[:, 1]
    mapped = np.searchsorted(keys, surface_keys)
    np.testing.assert_array_equal(keys[mapped], surface_keys)
    elevated_surface = np.column_stack((surface, mapped.reshape(-1, 3) + len(points)))
    elevated_points = np.concatenate((points, points[unique].mean(axis=1)))
    np.testing.assert_array_equal(elevated_points[elevated[:, 4:]], points[tetrahedra[:, edges]].mean(axis=2))
    result.update(points=elevated_points, tetrahedra=elevated, surface=elevated_surface)
    return result


def worker(source, run, use_amg=False):
    np, _, meshers, *_ = workflow.imports()
    from periodic import periodic_pairs, tensor_difference
    from stress_recovery import solve, vector_first_reactions
    from tet10_integration import audit_saved, integrated_tensors

    parent = json.loads((source / "result.json").read_text())
    assert parent["qualified"] and parent["geometry"]["geometry_admitted"]
    parent_calls = parent["mmg_calls"] if "mmg_calls" in parent else parent["parent_mmg_calls"]
    assert parent["route"].startswith("microgen-vtk") and parent_calls
    original = np.load(source / "mesh.npz")
    points, tetrahedra, surface = (original[name] for name in ("points", "tetrahedra", "surface"))
    assert tetrahedra.shape[1] == 4 and len(tetrahedra) <= workflow.CAP
    started = time.perf_counter()
    elevation = straight_elevation(points, tetrahedra, surface)
    nodes = elevation["nodes"]
    row = {"case": parent["case"], "route": "microgen-vtk-mmg-straight-tet10-direct-fedoo",
        "parent_route": parent["route"],
        "parent_mmg_calls": parent_calls, "source": str(source),
        "nodes": nodes, "elements": len(tetrahedra), "order": 2, "solver_quadrature_points": 4,
        "qualified": False, "status": "lift-node-cap", "live_workflow_timing": False,
        "solver": "cached-elastic-block-amg-cg" if use_amg else "jacobi-cg",
        "release_completed_reduced_cache": use_amg,
        "scope": "Saved genuine Microgen mesh; manual straight elevation and direct fedoo; no native Meshers generation or public Microgen Tet10 conversion"}
    workflow.save(run / "result.json", row)
    if not elevation["within_node_cap"]:
        return
    elevated_points, elevated, elevated_surface = (elevation[key] for key in ("points", "tetrahedra", "surface"))
    matching, _, periodic = periodic_pairs(elevated_points, elevated_surface)
    tags = np.zeros(len(surface), dtype=np.uint8)
    for axis in range(3):
        for side, value in enumerate((-0.5, 0.5)):
            tags[np.all(points[surface, axis] == value, axis=1)] = 1 + 2 * axis + side
    converted = meshers.Mesh(elevated_points, elevated, elevated_surface, tags, matching,
        {"source": "Actual Microgen VTK/MMG; straight midpoint elevation only"}).to_fedoo()
    np.testing.assert_array_equal(converted.nodes, elevated_points)
    np.testing.assert_array_equal(converted.elements, elevated)
    assert converted.elm_type == "tet10"
    row["lift_and_conversion_s"] = time.perf_counter() - started
    np.savez_compressed(run / "mesh.npz", points=elevated_points, tetrahedra=elevated, surface=elevated_surface)
    row["geometry"] = workflow.geometry(run, parent["case"], elevated_points, elevated, elevated_surface,
        reverse_subdivisions=(1,))
    row.update(periodic=periodic, qualified=row["geometry"]["geometry_admitted"],
        status="admitted-needs-reference" if row["geometry"]["geometry_admitted"] else "geometry-rejected")
    workflow.save(run / "result.json", row)
    if not row["qualified"]:
        return
    phase = time.perf_counter()
    def checkpoint(completed, partial_tensor, records):
        row.update(status="solving", completed_loads=completed,
            partial_stiffness_columns=partial_tensor.tolist(), solver_records=records)
        workflow.save(run / "result.json", row)
    def save_field(completed, displacement):
        np.savez_compressed(run / f"load-{completed}.npz", displacement=displacement)
    solver = None
    if use_amg:
        from amg_reference import AuditedAMG
        solver = AuditedAMG(elastic_blocks=True)
    fields, tensor, records = solve(converted.nodes, converted.elements, checkpoint=checkpoint,
        n_elm_gp=4, field_checkpoint=save_field, solver=solver,
        reaction_reader=vector_first_reactions if use_amg else None,
        release_reduced_cache=use_amg)
    row.pop("partial_stiffness_columns", None)
    row.update(solve_s=time.perf_counter() - phase, stiffness=tensor.tolist(), solver_records=records,
        status="solved-needs-reference", quadrature_exact_for_affine_geometry=True)
    np.savez_compressed(run / "fields.npz", displacement=fields)
    _, higher_energy, _ = integrated_tensors(elevated_points, elevated, fields)
    row["higher_rule_fixed_field_difference"] = tensor_difference(higher_energy, tensor)
    assert row["higher_rule_fixed_field_difference"]["worst_relative_strain_energy_difference"] < 1e-7
    row.update(quadrature_sensitivity_limit=1e-4, quadrature_sensitivity_passed=True,
        higher_rule_consistency_limit=1e-7)
    workflow.save(run / "result.json", row)
    audit_saved(run, run / "field-consistency.json")
    row["total_saved_mesh_study_s"] = time.perf_counter() - started
    workflow.save(run / "result.json", row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--amg", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.source, args.run, args.amg)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    case = json.loads((args.source / "result.json").read_text())["case"]
    files = [args.source / name for name in ("result.json", "mesh.npz", case + "-samples.npz", case + "-specimen.json")]
    for path in files[2:]:
        shutil.copyfile(path, args.run / path.name)
    files += [Path(__file__), Path(workflow.__file__), workflow.HERE / "stress_recovery.py",
        workflow.HERE / "tet10_integration.py", workflow.ROOT / "crates/meshers-python/python/meshers/__init__.py"]
    if args.amg:
        files.extend(workflow.HERE / name for name in ("amg_reference.py", "test_cubic_reference.py"))
    hashes = {str(path): workflow.digest(path) for path in files}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 1800, "solve_node_cap": 250_000, "live_workflow_timing": False,
        "solver": "cached-elastic-block-amg-cg" if args.amg else "jacobi-cg",
        "pyamg": importlib.metadata.version("pyamg") if args.amg else None,
        "cg_rtol": 1e-12, "cg_atol": 1e-14})
    parts = ["four-point-controls", "straight-tet10"] if args.amg else ["straight-tet10"]
    with Job(parts=parts) as job:
        if args.amg:
            with job.part("four-point-controls") as part:
                record = workflow.bounded([sys.executable, "-m", "pytest", "-q",
                    str(workflow.HERE / "test_cubic_reference.py")],
                    args.run / "four-point-controls.log", part, timeout=120)
                workflow.save(args.run / "four-point-controls.json", record)
                assert record["returncode"] == 0, "Read four-point-controls.log"
                job.result("four-point-controls", seconds=record["wall_s"])
        with job.part("straight-tet10") as part:
            command = [sys.executable, __file__, "--worker", "--source", str(args.source),
                "--run", str(args.run)]
            if args.amg:
                command.append("--amg")
            record = workflow.bounded(command, args.run / "worker.log", part, timeout=1800)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("straight-tet10", seconds=record["wall_s"], nodes=row["nodes"],
                qualified=int(row["qualified"]), peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
