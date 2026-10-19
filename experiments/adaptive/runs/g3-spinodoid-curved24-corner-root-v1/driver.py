"""Probe selective restoration of native-damped quadratic edge proposals."""

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
    from periodic import periodic_pairs
    from scorecard import BARY, tet_metrics

    saved = json.loads((args.source / "result.json").read_text())
    assert saved["case"] == "spinodoid" and saved["order"] == 2
    mesh = np.load(args.source / "mesh.npz")
    points, cells, surface = (mesh[name].copy() for name in ("points", "tetrahedra", "surface"))
    assert len(cells) <= workflow.CAP and len(points) <= 150000
    assert cells.shape[1] == 10 and surface.shape[1] == 6
    case, shape, _ = workflow.specimen("spinodoid")
    frep = shape._frep
    frozen = json.loads((args.source / "fourier-coefficients.json").read_text())
    coefficients = [(*map(float, mode / frep.cell_size), float(amplitude), float(phase))
        for mode, amplitude, phase in zip(frep.modes, frep.amplitudes, frep.phases, strict=True)]
    np.testing.assert_array_equal(coefficients, frozen["coefficients"])
    assert frep.threshold == frozen["threshold"] and shape.resolution == 16
    baseline = workflow.geometry(args.run, "spinodoid", points, cells, surface, detail_count=20)
    assert abs(baseline["sampled_two_sided_max"] - saved["geometry"]["sampled_two_sided_max"]) < 1e-11
    assert baseline["projection_failed"] == saved["geometry"]["projection_failed"]
    _, classes, _ = periodic_pairs(points, surface)
    cap = np.zeros(len(surface), dtype=bool)
    for axis in range(3):
        for side in (-.5, .5):
            cap |= np.all(abs(points[surface[:, :3], axis] - side) < 1e-12, axis=1)
    wall = surface[~cap]
    endpoints = {}
    for edge, (first, second) in enumerate(((0, 1), (1, 2), (2, 0))):
        for face in wall:
            node = int(face[3 + edge])
            pair = tuple(sorted((int(face[first]), int(face[second]))))
            assert node not in endpoints or endpoints[node] == pair
            endpoints[node] = pair
    nodes = np.array(sorted(endpoints))
    residuals = abs(case.function(*points[nodes].T))
    candidates = nodes[residuals > 1e-8]
    groups = {}
    unresolved = []
    recovered = []
    for node in candidates:
        straight = points[list(endpoints[int(node)])].mean(axis=0)
        displacement = points[node] - straight
        if np.linalg.norm(displacement) <= 1e-14:
            unresolved.append(int(node))
            continue
        targets = straight + (2. ** np.arange(1, 9))[:, None] * displacement
        inside = (abs(targets) <= .5 + 1e-14).all(axis=1)
        roots = inside & (abs(case.function(*targets.T)) < 1e-10)
        if not roots.any():
            unresolved.append(int(node))
            continue
        reductions = int(np.flatnonzero(roots)[0]) + 1
        target = targets[reductions - 1]
        canonical = np.where(straight == .5, -.5, straight)
        key = tuple(np.rint((canonical + .5) * 1e11).astype(np.int64))
        entry = {"node": int(node), "straight": straight, "target": target, "reductions": reductions}
        groups.setdefault(key, []).append(entry)
        recovered.append(int(node))
    assert len(recovered) <= saved["diagnostics"]["quadratic"]["partially_curved_edges"]
    element = fd.lib_elements.element_list.get_element("tet10")(15)
    vertices = fd.lib_elements.element_list.get_element("tet4")().xi_nd
    barycentric = []
    for first, second in ((.7, .1), (.1, .3)):
        for index in range(4):
            values = np.full(4, second)
            values[index] = first
            barycentric.append(values)
    locations = np.concatenate((element.xi_nd, element.xi_pg, [[.25] * 3],
        np.array([[0 if index == missing else 1 / 3 for index in range(4)]
            for missing in range(4)]) @ vertices, np.array(barycentric) @ vertices))
    derivatives = element.shape_function_derivative(locations)

    def minimum_ratio(tetrahedra):
        determinants = np.array([np.linalg.det(np.einsum("an,tnb->tab", derivative, points[tetrahedra]))
            for derivative in derivatives])
        return float((determinants.min(axis=0) / determinants.max(axis=0)).min())

    before_ratio = minimum_ratio(cells)
    assert before_ratio > 1e-3
    trials = []
    changed = []
    for key, members in sorted(groups.items()):
        member_nodes = np.array([member["node"] for member in members])
        previous = points[member_nodes].copy()
        incident = cells[np.isin(cells, member_nodes).any(axis=1)]
        maximum_reductions = max(member["reductions"] for member in members)
        selected = None
        for reductions in range(maximum_reductions):
            for member in members:
                points[member["node"]] = member["straight"] + 2.**(-reductions) * (member["target"] - member["straight"])
            ratio = minimum_ratio(incident)
            trials.append({"nodes": member_nodes.tolist(), "remaining_halvings": reductions,
                "minimum_incident_jacobian_ratio": ratio, "accepted": ratio > 1e-3})
            if ratio > 1e-3:
                selected = reductions
                changed.extend(member_nodes.tolist())
                break
        if selected is None:
            points[member_nodes] = previous
    corner_trials = []
    if args.project_corners:
        seen = set()
        edge_order = np.array(((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)))
        corner_count = int(cells[:, :4].max()) + 1
        interpolation = workflow.tri_values(BARY)
        for probe in baseline["worst_boundary_samples"]:
            support = np.flatnonzero(abs(interpolation[probe["probe_index"]]) > 1e-14)
            if probe["failed"] or probe["distance"] <= .002 or len(support) != 1:
                continue
            node = probe["node_ids"][int(support[0])]
            if node >= corner_count or int(classes[node]) in seen:
                continue
            seen.add(int(classes[node]))
            members = np.flatnonzero(classes == classes[node])
            assert (members < corner_count).all()
            original = points.copy()
            direction = np.asarray(probe["foot"]) - points[node]
            incident = cells[np.isin(cells[:, :4], members).any(axis=1)]
            for fraction in (1., .5, .25):
                displacement = np.zeros_like(points)
                displacement[members] = fraction * direction
                displacement[cells[:, 4:].ravel()] = displacement[cells[:, :4][:, edge_order]].mean(axis=2).reshape(-1, 3)
                points[:] = original + displacement
                quality = float(tet_metrics(points, incident[:, :4])[0].min())
                ratio = minimum_ratio(incident)
                determinants = np.linalg.det(points[incident[:, 1:4]] - points[incident[:, :1]])
                accepted = quality >= .1 and ratio > 1e-3 and determinants.min() > 0
                corner_trials.append({"nodes": members.tolist(), "origin": original[node].tolist(),
                    "root": probe["foot"], "fraction": fraction,
                    "minimum_corner_quality": quality, "minimum_jacobian_ratio": ratio,
                    "accepted": bool(accepted)})
                if accepted:
                    break
                points[:] = original
    after_ratio = minimum_ratio(cells)
    assert after_ratio > 1e-3
    _, _, periodic = periodic_pairs(points, surface)
    after = workflow.geometry(args.run, "spinodoid", points, cells, surface, detail_count=20)
    weights = workflow.tri_values(BARY)
    supported = sum(bool(set(probe["node_ids"][index] for index in np.flatnonzero(
        abs(weights[probe["probe_index"]]) > 1e-14)) & set(recovered))
        for probe in baseline["worst_boundary_samples"])
    np.savez_compressed(args.run / "mesh.npz", points=points, tetrahedra=cells, surface=surface)
    workflow.save(args.run / "result.json", {"source": str(args.source),
        "case": "spinodoid", "nodes": len(points), "elements": len(cells),
        "native_partially_curved_edges": saved["diagnostics"]["quadratic"]["partially_curved_edges"],
        "wall_mid_nodes_with_residual": len(candidates), "recovered_damped_proposals": len(recovered),
        "unresolved_mid_nodes": unresolved, "restored_nodes": changed,
        "worst_twenty_probes_supported_on_damped_nodes": supported,
        "baseline_minimum_jacobian_ratio": before_ratio, "after_minimum_jacobian_ratio": after_ratio,
        "corner_projection_trials": corner_trials,
        "baseline_geometry": baseline, "geometry": after, "periodic": periodic, "trials": trials,
        "geometry_admitted": after["geometry_admitted"],
        "scope": "Saved-mesh selective curvature diagnostic only; unchanged corner topology and no live native implementation or stiffness result"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--project-corners", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), Path(workflow.__file__), workflow.HERE / "scorecard.py"]
    paths.extend(args.source / name for name in ("mesh.npz", "result.json", "fourier-coefficients.json",
        "spinodoid-samples.npz", "spinodoid-specimen.json"))
    for name in ("spinodoid-samples.npz", "spinodoid-specimen.json"):
        shutil.copyfile(args.source / name, args.run / name)
    hashes = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "worker_timeout_s": 300, "threads": 1, "jacobian_ratio_floor": 1e-3})
    with Job(parts=["selective-curve-restoration"]) as job:
        with job.part("selective-curve-restoration") as part:
            record = workflow.bounded([sys.executable, __file__, *sys.argv[1:], "--worker"],
                args.run / "worker.log", part, timeout=300)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            result = json.loads((args.run / "result.json").read_text())
            job.result("selective-curve-restoration", seconds=record["wall_s"],
                restored_nodes=len(result["restored_nodes"]), admitted=int(result["geometry_admitted"]))
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
