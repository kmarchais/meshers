"""Bound a four-point Tet10 re-solve against a saved fifteen-point control."""

import argparse
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

LIMIT = 1e-4


def worker(source, run, node_cap):
    np, fd, *_ = workflow.imports()
    from homogenization import MATERIAL, strain_tensor
    from periodic import periodic_pairs, tensor_difference
    from stress_recovery import solve
    from tet10_integration import audit_saved, integrated_tensors

    audit_saved(source, run / "control-field-consistency.json")
    control = json.loads((source / "result.json").read_text())
    assert control["qualified"] and control["geometry"]["geometry_admitted"]
    mesh = np.load(source / "mesh.npz")
    points, tetrahedra = mesh["points"], mesh["tetrahedra"]
    assert tetrahedra.shape[1] == 10 and len(points) <= node_cap
    element = fd.lib_elements.element_list.get_element("tet10")(4)
    quadrature = (element.xi_pg, element.w_pg)
    affine_fields = np.array([element.xi_nd @ strain_tensor(engineering).T for engineering in np.eye(6)])
    patch_stress, patch_energy, patch_volume = integrated_tensors(
        element.xi_nd, np.arange(10)[None], affine_fields, quadrature)
    np.testing.assert_allclose(patch_stress, MATERIAL / 6, atol=1e-13, rtol=1e-13)
    np.testing.assert_allclose(patch_energy, MATERIAL / 6, atol=1e-13, rtol=1e-13)
    assert abs(patch_volume - 1 / 6) < 1e-14
    expected = np.asarray(control["stiffness"])
    saved = np.load(source / "fields.npz")["displacement"]
    _, reintegrated, volume = integrated_tensors(points, tetrahedra, saved, quadrature)
    precheck = tensor_difference(reintegrated, expected)
    row = {"source": str(source), "quadrature_points": 4, "control_quadrature_points": 15,
        "predeclared_relative_stiffness_limit": LIMIT, "fixed_field_precheck": precheck,
        "relative_volume_shift": float(volume / control["geometry"]["volume"] - 1),
        "scope": "Saved-mesh quadrature sensitivity, not a live workflow or continuum accuracy certificate",
        "status": "precheck-rejected", "passed": False}
    workflow.save(run / "result.json", row)
    if precheck["worst_relative_strain_energy_difference"] > LIMIT or abs(row["relative_volume_shift"]) > LIMIT:
        return
    started = time.perf_counter()
    fields, tensor, records = solve(points, tetrahedra, n_elm_gp=4)
    row.update(solve_s=time.perf_counter() - started, stiffness=tensor.tolist(), solver_records=records)
    stress, energy, _ = integrated_tensors(points, tetrahedra, fields, quadrature)
    consistency = max(float(np.linalg.norm(stress - tensor) / np.linalg.norm(tensor)),
        float(np.linalg.norm(energy - tensor) / np.linalg.norm(tensor)))
    assert consistency < 1e-7 and len(records) == 6
    pairs, _, _ = periodic_pairs(points, mesh["surface"])
    jump = 0.0
    for load, engineering in enumerate(np.eye(6)):
        for axis, matching in enumerate(pairs):
            difference = fields[load, matching[:, 1]] - fields[load, matching[:, 0]]
            jump = max(jump, float(abs(difference - strain_tensor(engineering)[:, axis]).max(initial=0)))
    assert jump < 1e-9
    difference = tensor_difference(tensor, expected)
    row.update(status="solved", own_rule_consistency_relative=consistency, periodic_jump_max=jump,
        control_difference=difference, passed=difference["worst_relative_strain_energy_difference"] <= LIMIT)
    workflow.save(run / "result.json", row)
    np.savez_compressed(run / "fields.npz", displacement=fields)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--node-cap", type=int, choices=(100_000, 200_000), default=100_000)
    parser.add_argument("--worker-timeout", type=int, choices=(600, 1200), default=600)
    args = parser.parse_args()
    if args.worker:
        worker(args.source, args.run, args.node_cap)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    files = [args.source / name for name in ("result.json", "mesh.npz", "fields.npz")]
    files += [Path(__file__), Path(workflow.__file__), workflow.HERE / "stress_recovery.py",
        workflow.HERE / "tet10_integration.py"]
    hashes = {str(path): workflow.digest(path) for path in files}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "node_cap": args.node_cap, "worker_timeout_s": args.worker_timeout,
        "predeclared_relative_stiffness_limit": LIMIT})
    with Job(parts=["four-point-control"]) as job:
        with job.part("four-point-control") as part:
            record = workflow.bounded([sys.executable, __file__, "--worker", "--source", str(args.source),
                "--run", str(args.run), "--node-cap", str(args.node_cap)],
                args.run / "worker.log", part, timeout=args.worker_timeout)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("four-point-control", seconds=record["wall_s"], passed=int(row["passed"]),
                peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
