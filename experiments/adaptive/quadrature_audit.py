"""Audit immutable fields with positive Duffy rules, without claiming a re-solve."""

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow
from job_parts import worker_part

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def duffy_rule(order):
    import numpy as np

    abscissae, weights = np.polynomial.legendre.leggauss(order)
    abscissae, weights = (abscissae + 1) / 2, weights / 2
    first, second, third = np.meshgrid(abscissae, abscissae, abscissae, indexing="ij")
    weight_first, weight_second, weight_third = np.meshgrid(weights, weights, weights, indexing="ij")
    locations = np.column_stack((first.ravel(), ((1 - first) * second).ravel(),
        ((1 - first) * (1 - second) * third).ravel()))
    weights = (weight_first * weight_second * weight_third * (1 - first)**2 * (1 - second)).ravel()
    assert (weights > 0).all() and (locations >= 0).all() and (locations.sum(axis=1) <= 1).all()
    for degree_first in range(2 * order - 2):
        for degree_second in range(2 * order - 2 - degree_first):
            for degree_third in range(2 * order - 2 - degree_first - degree_second):
                degrees = np.array([degree_first, degree_second, degree_third])
                exact = math.prod(math.factorial(int(degree)) for degree in degrees) / math.factorial(int(degrees.sum()) + 3)
                assert abs(weights @ np.prod(locations**degrees, axis=1) - exact) < 1e-14
    return locations, weights


def worker(source, run):
    np, fd, *_ = workflow.imports()
    from homogenization import MATERIAL, strain_tensor
    from periodic import tensor_difference
    from tet10_integration import audit_saved, integrated_tensors

    audit_saved(source, run / (source.name + "-consistency.json"))
    row = json.loads((source / "result.json").read_text())
    mesh = np.load(source / "mesh.npz")
    fields = np.load(source / "fields.npz")["displacement"]
    expected = np.asarray(row["stiffness"])
    affine_points = fd.lib_elements.element_list.get_element("tet10")(15).xi_nd
    affine_fields = np.array([affine_points @ strain_tensor(engineering).T for engineering in np.eye(6)])
    rules, previous_energy = [], None
    for order in (4, 5):
        quadrature = duffy_rule(order)
        control_stress, control_energy, control_volume = integrated_tensors(
            affine_points, np.arange(10)[None], affine_fields, quadrature)
        np.testing.assert_allclose(control_stress, MATERIAL / 6, atol=1e-13, rtol=1e-13)
        np.testing.assert_allclose(control_energy, MATERIAL / 6, atol=1e-13, rtol=1e-13)
        assert abs(control_volume - 1 / 6) < 1e-14
        stress, energy, volume = integrated_tensors(mesh["points"], mesh["tetrahedra"], fields, quadrature)
        symmetrized_stress = (stress + stress.T) / 2
        rules.append({"points": order**3, "stress_integrated_stiffness": stress.tolist(),
            "fixed_field_energy_tensor": energy.tolist(), "volume": volume,
            "stress_vs_solved_tensor": tensor_difference(symmetrized_stress, expected),
            "energy_vs_solved_tensor": tensor_difference(energy, expected),
            "energy_vs_previous_rule": None if previous_energy is None else tensor_difference(energy, previous_energy),
            "stress_energy_relative_difference": float(np.linalg.norm(stress - energy) / np.linalg.norm(energy)),
            "stress_skew_relative": float(np.linalg.norm(stress - stress.T) / np.linalg.norm(stress)),
            "relative_volume_difference": float(volume / row["geometry"]["volume"] - 1)})
        previous_energy = energy
    workflow.save(run / (source.name + ".json"), {"source": str(source), "rules": rules,
        "controls": "Positive rules; exact tetrahedral monomials through degree 2*n-3; affine Tet10 elasticity patch",
        "scope": "Reintegration of fixed displacements, not higher-quadrature equilibrium re-solves or continuum bounds"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--sources", nargs="+", type=Path)
    parser.add_argument("--worker", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.run)
        return
    assert args.sources and len({source.name for source in args.sources}) == len(args.sources)
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    from tet10_integration import __file__ as integration_source
    sources = [Path(__file__), Path(integration_source), Path(workflow.__file__),
               workflow.HERE / "job_parts.py"]
    inputs = {str(source / name): workflow.digest(source / name) for source in args.sources
        for name in ("result.json", "mesh.npz", "fields.npz")}
    source_hashes = {str(source): workflow.digest(source) for source in sources}
    workflow.save(args.run / "manifest.json", {"inputs": inputs, "source_sha256": source_hashes,
        "numerical_thread_policy": workflow.NUMERICAL_THREAD_POLICY, "platform": sys.platform,
        "memory_bytes": workflow.MEMORY, "threads": 1, "worker_timeout_s": 300,
        "quadrature_points": [15, 64, 125], "fixed_displacement_fields": True})
    with Job(parts=[source.name for source in args.sources]) as job:
        for source in args.sources:
            with worker_part(job, source.name) as part:
                record = workflow.bounded([sys.executable, __file__, "--run", str(args.run),
                    "--worker", str(source)], args.run / (source.name + ".log"), part, timeout=300)
                workflow.save(args.run / (source.name + "-worker.json"), record)
                assert record["returncode"] == 0, f"Read {source.name}.log"
                row = json.loads((args.run / (source.name + ".json")).read_text())
                job.result(source.name, seconds=record["wall_s"],
                    energy_shift=row["rules"][-1]["energy_vs_solved_tensor"]["worst_relative_strain_energy_difference"],
                    peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == value for path, value in {**inputs, **source_hashes}.items())


if __name__ == "__main__":
    main()
