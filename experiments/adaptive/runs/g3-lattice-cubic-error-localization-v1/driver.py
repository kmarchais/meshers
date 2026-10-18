"""Localize same-domain enrichment energy for a retained-model worst strain."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import microgen_workflow as workflow

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import Job


def worker(args):
    np, fd, *_ = workflow.imports()
    from cubic_reference import register_cubic
    from quadrature_audit import duffy_rule

    parent = json.loads((args.parent / "result.json").read_text())
    cubic = json.loads((args.cubic / "result.json").read_text())
    assert parent["qualified"] and cubic["qualified"]
    assert parent["geometry"]["geometry_admitted"] and cubic["quadrature_sensitivity_passed"]
    assert cubic["geometry_source"] == str(args.parent)
    comparison = json.loads(args.comparison.read_text())
    candidate = next(row for row in comparison["candidates"] if row["source"].endswith("g3-lattice-adapt24-tol0028-v1"))
    strain = np.asarray(candidate["worst_model"]["engineering_strain_in_saved_tensor_basis"])
    original = np.load(args.parent / "mesh.npz")
    mesh = np.load(args.prepared / "mesh.npz")
    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    assert len(points) <= 250000 and len(cells) <= workflow.CAP
    element = register_cubic()()
    parent_values = fd.lib_elements.element_list.get_element("tet10")().shape_function(element.xi_nd)
    original_field = np.einsum("l,lnc->nc", strain,
        np.load(args.parent / "fields.npz")["displacement"])
    embedded = np.zeros((len(points), 3))
    assigned = np.zeros(len(points), dtype=bool)
    for index, weights in enumerate(parent_values):
        destinations = cells[:, index]
        values = np.einsum("n,tnc->tc", weights, original_field[original["tetrahedra"]])
        known = assigned[destinations]
        np.testing.assert_allclose(embedded[destinations[known]], values[known], atol=1e-12, rtol=0)
        embedded[destinations] = values
        assigned[destinations] = True
    assert assigned.all()
    correction = embedded - np.einsum("l,lnc->nc", strain,
        np.load(args.cubic / "fields.npz")["displacement"])
    expected = float(strain @ (np.asarray(parent["stiffness"]) - np.asarray(cubic["stiffness"])) @ strain)
    assert expected > 0
    coordinates, displacement = points[cells], correction[cells]
    results, previous = [], None
    triangle_points = points[surface[:, :3]]
    cap = np.zeros(len(surface), dtype=bool)
    for axis in range(3):
        for side in (-.5, .5):
            cap |= np.all(abs(triangle_points[:, :, axis] - side) < 1e-10, axis=1)
    wall_incident = np.isin(cells[:, :4], np.unique(surface[~cap, :3])).any(axis=1)
    cap_incident = np.isin(cells[:, :4], np.unique(surface[cap, :3])).any(axis=1)
    for label, quadrature in (("own15", (element.xi_pg, element.w_pg)),
        ("positive64", duffy_rule(4)), ("positive125", duffy_rule(5))):
        locations, weights = quadrature
        energy, volume = np.zeros(len(cells)), np.zeros(len(cells))
        for derivative, weight in zip(element.shape_function_derivative(locations), weights, strict=True):
            jacobian = np.einsum("an,tnb->tab", derivative, coordinates)
            determinants = np.linalg.det(jacobian)
            assert determinants.min() > 0 and weight > 0
            gradients = np.einsum("tba,an->tbn", np.linalg.inv(jacobian), derivative)
            gradient = np.einsum("tni,tjn->tij", displacement, gradients)
            deformation = (gradient + gradient.swapaxes(-1, -2)) / 2
            stress = 2 * deformation + np.trace(deformation, axis1=-2, axis2=-1)[:, None, None] * np.eye(3)
            energy += weight * determinants * np.einsum("tij,tij->t", deformation, stress)
            volume += weight * determinants
        assert energy.min() >= 0 and volume.min() > 0
        total = float(energy.sum())
        identity_error = abs(total / expected - 1)
        if label == "own15":
            assert identity_error < 1e-6
        count = int(np.ceil(.2 * len(cells)))
        top = np.argsort(energy)[-count:]
        volume_top = np.argsort(volume)[-count:]
        results.append({"quadrature": label, "enrichment_energy": total,
            "galerkin_identity_relative_error": identity_error,
            "relative_change_from_previous_rule": None if previous is None else abs(total / previous - 1),
            "highest_energy_twenty_percent_capture": float(energy[top].sum() / total),
            "largest_volume_twenty_percent_capture": float(energy[volume_top].sum() / total),
            "wall_vertex_incident_element_fraction": float(wall_incident.mean()),
            "wall_vertex_incident_energy_fraction": float(energy[wall_incident].sum() / total),
            "cap_vertex_incident_element_fraction": float(cap_incident.mean()),
            "cap_vertex_incident_energy_fraction": float(energy[cap_incident].sum() / total)})
        previous = total
        np.savez_compressed(args.run / (label + ".npz"), energy=energy, volume=volume)
    workflow.save(args.run / "result.json", {"source_parent": str(args.parent),
        "source_cubic": str(args.cubic), "engineering_strain": strain.tolist(),
        "worst_model_witness": candidate["worst_model"], "expected_galerkin_energy": expected,
        "rules": results, "scope": "Same-domain discrete enrichment energy localization; no continuum or local stress certificate; no meshing improvement inferred"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--cubic", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), args.comparison, args.prepared / "mesh.npz"]
    paths.extend(source / name for source in (args.parent, args.cubic) for name in ("result.json", "fields.npz"))
    paths.append(args.parent / "mesh.npz")
    paths.extend(workflow.HERE / name for name in ("cubic_reference.py", "quadrature_audit.py", "microgen_workflow.py"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    prior = json.loads((args.cubic / "manifest.json").read_text())["sources"]
    for path in (args.prepared / "mesh.npz", args.parent / "mesh.npz", args.parent / "fields.npz", args.parent / "result.json"):
        assert prior[str(path)] == hashes[str(path)]
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": 300, "threads": 1})
    with Job(parts=["same-domain-error"]) as job:
        with job.part("same-domain-error") as part:
            record = workflow.bounded([sys.executable, __file__, *sys.argv[1:], "--worker"],
                args.run / "worker.log", part, timeout=300)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            job.result("same-domain-error", seconds=record["wall_s"],
                peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == value for path, value in hashes.items())


if __name__ == "__main__":
    main()
