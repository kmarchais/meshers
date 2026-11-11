"""One fixed-count ridge projection control with the original quality floor.

Only endpoints of branch-crossing wall edges are considered. A full move to
the analytic ridge must stay within one quarter of the shortest incident wall
edge and preserve the original global minimum tetrahedron quality. Cell caps
and connectivity stay bitwise fixed. This is a saved-mesh feasibility test.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import cad_intersection_probe as probe
import cad_union_comparison as base
import microgen_workflow as w


def local_metrics(points, cells):
    import numpy as np

    p = points[cells]
    determinant = np.einsum("ij,ij->i", p[:, 1]-p[:, 0],
                           np.cross(p[:, 2]-p[:, 0], p[:, 3]-p[:, 0]))
    squared = sum(((p[:, a]-p[:, b])**2).sum(1) for a in range(4) for b in range(a))
    quality = (12 * (abs(determinant)/2)**(2/3) / squared)**1.5
    return quality, determinant


def repair(points, cells, surface, part):
    import numpy as np

    moved = points.copy()
    original_quality, determinant = local_metrics(points, cells)
    assert (determinant > 0).all()
    floor = float(original_quality.min())
    wall = surface[probe.wall_mask(points, surface)]
    edges = np.unique(np.sort(wall[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1), axis=0)
    branch = probe.branch_value(points)
    crossing = edges[((branch[edges[:, 0]] < -1e-12) & (branch[edges[:, 1]] > 1e-12)) |
                     ((branch[edges[:, 1]] < -1e-12) & (branch[edges[:, 0]] > 1e-12))]
    candidates = np.unique(crossing)
    cap_nodes = np.flatnonzero(np.any(abs(abs(points)-.5) < 1e-10, axis=1))
    assert not np.intersect1d(candidates, cap_nodes).size, "This control does not move cap classes"
    incident = [[] for _ in points]
    neighbors = [[] for _ in points]
    for i, tet in enumerate(cells):
        for node in tet:
            incident[node].append(i)
    for a, b in edges:
        neighbors[a].append(int(b))
        neighbors[b].append(int(a))
    target = probe.ridge_target(points[candidates])
    length = np.array([np.linalg.norm(points[neighbors[node]]-points[node], axis=1).min() for node in candidates])
    ratio = np.linalg.norm(target-points[candidates], axis=1)/length
    order = np.lexsort((candidates, ratio))
    attempts = []
    accepted = []
    for step, index in enumerate(order):
        node = int(candidates[index])
        selected = np.asarray(incident[node], dtype=int)
        current_branch = probe.branch_value(moved[np.r_[node, neighbors[node]]])
        still_crosses = ((current_branch[0] < -1e-12) & (current_branch[1:] > 1e-12)).any() or \
                        ((current_branch[0] > 1e-12) & (current_branch[1:] < -1e-12)).any()
        row = dict(node=node, requested_displacement=float(np.linalg.norm(target[index]-points[node])),
                   shortest_original_wall_edge=float(length[index]), displacement_ratio=float(ratio[index]))
        if not still_crosses:
            row["status"] = "crossings-already-covered"
        elif ratio[index] > .25:
            row["status"] = "move-bound-rejected"
        else:
            old = moved[node].copy()
            moved[node] = target[index]
            quality, determinant = local_metrics(moved, cells[selected])
            row.update(candidate_incident_q_min=float(quality.min()),
                       candidate_minimum_determinant=float(determinant.min()))
            if determinant.min() <= 0 or quality.min() < floor-1e-13:
                row["status"] = "quality-rejected"
                moved[node] = old
            else:
                row["status"] = "accepted"
                row["point"] = moved[node].tolist()
                accepted.append(node)
        attempts.append(row)
        part.update((step+1)/max(len(order), 1), "Testing complete ridge projections",
                    accepted=len(accepted), candidates=len(order), quality_floor=floor)
    np.testing.assert_array_equal(moved[cap_nodes], points[cap_nodes])
    quality, determinant = local_metrics(moved, cells)
    assert determinant.min() > 0 and quality.min() >= floor-1e-13
    if accepted:
        moved_ridge = moved[accepted]
        np.testing.assert_allclose(moved_ridge[:, 0]**2+moved_ridge[:, 2]**2, base.fixture.RADIUS**2, atol=1e-15)
        np.testing.assert_allclose(moved_ridge[:, 1]**2+moved_ridge[:, 2]**2, base.fixture.RADIUS**2, atol=1e-15)
    quantiles = [0, .01, .1, .5, .9, .99, 1]
    keys = ("min", "p01", "p10", "median", "p90", "p99", "max")
    statuses = {name: sum(row["status"] == name for row in attempts) for name in sorted({row["status"] for row in attempts})}
    return moved, dict(original_global_quality_floor=floor, displacement_limit_shortest_wall_edge=.25,
        original_crossing_edges=len(crossing), candidates=len(candidates), accepted_nodes=len(accepted),
        status_counts=statuses, attempts=attempts, cap_nodes_bitwise_unchanged=True,
        original_quality=dict(zip(keys, map(float, np.quantile(original_quality, quantiles)))),
        final_quality=dict(zip(keys, map(float, np.quantile(quality, quantiles)))),
        maximum_accepted_displacement=max((row["requested_displacement"] for row in attempts if row["status"] == "accepted"), default=0),
        maximum_accepted_displacement_ratio=max((row["displacement_ratio"] for row in attempts if row["status"] == "accepted"), default=0))


def worker(args):
    import numpy as np
    w.imports()
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    diagnosis = json.loads((args.diagnosis / "result.json").read_text())
    assert diagnosis["complete"] and diagnosis["input_hashes_verified"]
    source = Path(diagnosis["rows"]["native"]["source"])
    assert w.digest(source / "mesh.npz") == diagnosis["rows"]["native"]["input_mesh_sha256"]
    paths = {Path(__file__), Path(probe.__file__), Path(base.__file__), Path(base.fixture.__file__),
        Path(w.__file__), args.diagnosis / "result.json"}
    paths.update(w.HERE / name for name in ("scorecard.py", "periodic_trace.py", "comparison_mesh_audit.py",
        "facet_incidence_audit.py", "quadratic_edge_removal_probe.py", "nested_refinement.py", "box_projection.py"))
    paths.update(source / name for name in ("mesh.npz", "result.json", base.NAME+"-samples.npz", base.NAME+"-specimen.json"))
    hashes = {str(p.resolve()): w.digest(p) for p in sorted(paths)}
    for p in paths:
        if p.suffix == ".py":
            shutil.copyfile(p, args.run / ("input-"+p.name))
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(source / (base.NAME+suffix), args.run / (base.NAME+suffix))
    w.save(args.run / "manifest.json", dict(inputs=hashes, worker_timeout_s=180, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        hypothesis="Ridge-only projection improves ridge fidelity but does not repair the dominant smooth-cylinder chord error",
        scope="Saved-mesh fixed-connectivity feasibility control. No generation time, FEA, or overall superiority claim"))
    result = dict(inputs=hashes, source=str(source), complete=False, original_samples_unchanged=True,
        geometry_definition_unchanged=True, connectivity_unchanged=True, gate_thresholds_unchanged=True)
    with np.load(source / "mesh.npz") as mesh:
        points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    with Job(parts={"ridge-projection": 2, "admission": 2, "localization": 2}) as job:
        with job.part("ridge-projection") as part:
            part.update(0, "Preparing fixed-count ridge control", nodes=len(points), tets=len(cells))
            moved, result["repair"] = repair(points, cells, surface, part)
        job.result("ridge-projection", accepted=result["repair"]["accepted_nodes"],
                   q_min=result["repair"]["final_quality"]["min"])
        with job.part("admission") as part:
            part.update(0, "Archiving then checking unchanged full gates", nodes=len(points), tets=len(cells))
            result["candidate"] = base.archive_and_audit(args.run, "ridge-control", (moved, cells, surface))
            part.update(1, "Completed fixed-sample audit", admitted=result["candidate"]["qualified"])
        job.result("admission", admitted=result["candidate"]["qualified"],
                   boundary=result["candidate"].get("geometry", {}).get("sampled_two_sided_max"))
        w.save(args.run / "result.json", result)
        with job.part("localization") as part:
            if result["candidate"]["qualified"]:
                result["localization"] = probe.describe(args.run / "ridge-control", part)
            else:
                part.update(1, "Rejected candidate retained", admitted=False)
        job.result("localization", ridge_error=result.get("localization", {}).get("reverse_groups", {}).get("ridges", {}).get("maximum"))
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    result.update(complete=True, input_hashes_verified=True, nodes=len(points), elements=len(cells))
    w.save(args.run / "result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--diagnosis", type=Path, default=probe.RUNS / "macbook-cad-intersection-probe-v1")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "diagnosis"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run / "worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=180)
    w.save(args.run / "worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import result
    row = json.loads((args.run / "result.json").read_text())
    geometry = row["candidate"].get("geometry", {})
    result("ridge-control", admitted=row["candidate"]["qualified"], accepted_nodes=row["repair"]["accepted_nodes"],
           boundary=geometry.get("sampled_two_sided_max"), absolute_volume_error=abs(geometry.get("relative_volume_error", 0)),
           q_min=row["repair"]["final_quality"]["min"],
           ridge_error=row.get("localization", {}).get("reverse_groups", {}).get("ridges", {}).get("maximum"))


if __name__ == "__main__":
    main()
