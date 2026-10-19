"""Audit periodic wall-normal jumps and same-domain FE error near box curves."""

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


def features(mesh):
    import numpy as np

    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    assert cells.shape[1] == 10 and surface.shape[1] == 6
    triangle = points[surface[:, :3]]
    cap = np.any(np.all(workflow.on_box_plane(abs(triangle), .5), axis=1), axis=1)
    wall_edges, cap_edges, mids = set(), set(), {}
    for face, is_cap in zip(surface, cap, strict=True):
        target = cap_edges if is_cap else wall_edges
        for offset, (first, second) in enumerate(((0, 1), (1, 2), (2, 0))):
            edge = tuple(sorted((int(face[first]), int(face[second]))))
            target.add(edge)
            midpoint = int(face[offset + 3])
            assert edge not in mids or mids[edge] == midpoint
            mids[edge] = midpoint
    curve_edges = sorted(wall_edges & cap_edges)
    assert curve_edges
    endpoints = np.unique(np.asarray(curve_edges).ravel())
    curve_nodes = np.unique(np.concatenate((endpoints, [mids[edge] for edge in curve_edges])))
    curve_vertex = np.isin(cells[:, :4], endpoints).any(axis=1)
    curve_set = set(curve_edges)
    curve_edge = np.array([any(tuple(sorted((int(cell[first]), int(cell[second])))) in curve_set
        for first, second in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))) for cell in cells[:, :4]])
    cap_vertex = np.isin(cells[:, :4], np.unique(surface[cap, :3])).any(axis=1)
    wall_vertex = np.isin(cells[:, :4], np.unique(surface[~cap, :3])).any(axis=1)
    chord = np.linalg.norm(points[np.asarray(curve_edges)[:, 1]] - points[np.asarray(curve_edges)[:, 0]], axis=1)
    return {"curve_edges": curve_edges, "curve_nodes": curve_nodes, "curve_endpoints": endpoints,
        "groups": {"curve_vertex_incident": curve_vertex, "curve_edge_incident": curve_edge,
            "cap_without_curve_vertex": cap_vertex & ~curve_vertex,
            "wall_without_curve_vertex": wall_vertex & ~curve_vertex,
            "neither_wall_nor_cap_vertex": ~wall_vertex & ~cap_vertex},
        "curve_edge_chord_quantiles": np.quantile(chord, [0, .1, .5, .9, 1]).tolist()}


def gradient(function, points, step):
    import numpy as np

    return np.column_stack([(function(*(points + step * direction).T) -
        function(*(points - step * direction).T)) / (2 * step) for direction in np.eye(3)])


def normal_jumps(function, points):
    import numpy as np

    origin = points.copy()
    locked = workflow.on_box_plane(abs(points), .5)
    for iteration in range(40):
        derivative = gradient(function, points, 1e-6)
        derivative[locked] = 0
        residual = function(*points.T)
        step = residual[:, None] * derivative / np.maximum((derivative**2).sum(axis=1), 1e-30)[:, None]
        points -= step
        if np.linalg.norm(step, axis=1).max() < 1e-12:
            break
    residual = abs(function(*points.T)) / np.maximum(np.linalg.norm(gradient(function, points, 1e-6), axis=1), 1e-30)
    assert residual.max() < 1e-9 and abs(points).max() <= .5 + 1e-9
    rows, all_angles, all_changes = [], [], []
    maximum_trace_mismatch = 0.
    for axis in range(3):
        selected = workflow.on_box_plane(points[:, axis], -.5)
        low = points[selected].copy()
        high = low.copy()
        high[:, axis] = .5
        assert len(low)
        mismatch = float(abs(function(*low.T) - function(*high.T)).max())
        maximum_trace_mismatch = max(maximum_trace_mismatch, mismatch)
        assert mismatch < 1e-9
        angles = []
        for epsilon in (2e-6, 1e-6):
            branches = []
            for branch in (low, high):
                probe = branch.copy()
                on_plane = workflow.on_box_plane(abs(probe), .5)
                probe[on_plane] -= np.sign(probe[on_plane]) * epsilon
                derivative = gradient(function, probe, epsilon / 8)
                norm = np.linalg.norm(derivative, axis=1)
                assert norm.min() > 1e-8
                branches.append(derivative / norm[:, None])
            cosine = np.einsum("ni,ni->n", *branches)
            angles.append(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
        change = abs(angles[1] - angles[0])
        rows.append({"axis": axis, "paired_curve_nodes": len(low), "scalar_trace_mismatch": mismatch,
            "normal_jump_angle_deg_quantiles": np.quantile(angles[1], [0, .1, .5, .9, 1]).tolist(),
            "fraction_above_ten_degrees": float(np.mean(angles[1] >= 10)),
            "maximum_inward_limit_angle_change_deg": float(change.max())})
        all_angles.extend(angles[1])
        all_changes.extend(change)
    return {"axes": rows, "maximum_scalar_trace_mismatch": maximum_trace_mismatch,
        "maximum_projection_distance": float(np.linalg.norm(points - origin, axis=1).max()),
        "maximum_projected_field_distance": float(residual.max()),
        "median_normal_jump_deg": float(np.median(all_angles)),
        "fraction_above_ten_degrees": float(np.mean(np.asarray(all_angles) >= 10)),
        "maximum_inward_limit_angle_change_deg": float(max(all_changes))}


def worker(args):
    np, *_ = workflow.imports()
    from scipy.linalg import eigvalsh

    parent = json.loads((args.parent / "result.json").read_text())
    fine = json.loads((args.fine / "result.json").read_text())
    screen = json.loads((args.gram / "result.json").read_text())
    assert Path(screen["parent"]).resolve() == args.parent.resolve()
    assert parent["qualified"] and fine["qualified"]
    mesh = np.load(args.parent / "mesh.npz")
    refined = np.load(args.fine / "mesh.npz")
    coarse_features, fine_features = features(mesh), features(refined)
    case, _, _ = workflow.specimen("smooth-lattice")
    normals = normal_jumps(case.function, mesh["points"][coarse_features["curve_nodes"]].copy())
    inverse = np.linalg.inv(np.asarray(parent["stiffness"]))
    rules = []
    for label in ("own15", "positive125"):
        archive = np.load(args.gram / (label + ".npz"))
        errors = archive["enrichment_gram"]
        assert errors.shape == (len(mesh["tetrahedra"]), 6, 6)
        total = errors.sum(axis=0)
        assert np.linalg.eigvalsh(total).min() > 0
        trace = float(np.trace(inverse @ total))
        groups = {}
        for name, mask in coarse_features["groups"].items():
            selected = errors[mask].sum(axis=0)
            groups[name] = {"elements": int(mask.sum()), "element_fraction": float(mask.mean()),
                "whitened_trace_fraction": float(np.trace(inverse @ selected) / trace),
                "minimum_mixed_load_capture": float(eigvalsh(selected, total).min()),
                "maximum_mixed_load_capture": float(eigvalsh(selected, total).max())}
        rules.append({"quadrature": label, "total_whitened_correction": trace, "groups": groups})
    primary = rules[-1]["groups"]["curve_vertex_incident"]
    passed = (primary["element_fraction"] <= .1 and primary["whitened_trace_fraction"] >= .6
        and primary["minimum_mixed_load_capture"] >= .5
        and normals["fraction_above_ten_degrees"] >= .2 and normals["median_normal_jump_deg"] >= 5
        and normals["maximum_inward_limit_angle_change_deg"] <= .1)
    workflow.save(args.run / "result.json", {"parent": str(args.parent), "fine": str(args.fine),
        "normals": normals, "rules": rules, "predeclared_curve_localization_screen_passed": bool(passed),
        "coarse_curve_edges": len(coarse_features["curve_edges"]),
        "fine_curve_edges": len(fine_features["curve_edges"]),
        "coarse_curve_edge_chord_quantiles": coarse_features["curve_edge_chord_quantiles"],
        "fine_curve_edge_chord_quantiles": fine_features["curve_edge_chord_quantiles"],
        "scope": "One-sided periodic-extension wall normals and fixed-domain full FE enrichment localization; no stress-singularity proof, total geometry-error bound, live sizing benefit or stiffness certification"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--gram", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, args.run / "driver.py")
    paths = [Path(__file__), Path(workflow.__file__), args.gram / "manifest.json", args.gram / "result.json",
             args.gram / "own15.npz", args.gram / "positive125.npz"]
    paths.extend(source / name for source in (args.parent, args.fine) for name in ("result.json", "mesh.npz"))
    hashes = {str(path): workflow.digest(path) for path in paths}
    frozen = json.loads((args.gram / "manifest.json").read_text())["inputs"]
    assert all(workflow.digest(Path(path)) == digest for path, digest in frozen.items())
    workflow.save(args.run / "manifest.json", {"inputs": hashes, "memory_bytes": workflow.MEMORY,
        "threads": 1, "worker_timeout_s": 300, "maximum_curve_element_fraction": .1,
        "minimum_trace_capture": .6, "minimum_every_mixed_load_capture": .5,
        "minimum_normal_jump_fraction_above_ten_degrees": .2, "minimum_median_normal_jump_deg": 5,
        "maximum_inward_limit_angle_change_deg": .1,
        "primary_group": "Tets incident on a corner vertex of a wall/cap intersection edge"})
    with Job(parts=["periodic-crease-audit"]) as job:
        with job.part("periodic-crease-audit") as part:
            record = workflow.bounded([sys.executable, __file__, *sys.argv[1:], "--worker"],
                args.run / "worker.log", part, timeout=300)
            workflow.save(args.run / "worker.json", record)
            assert record["returncode"] == 0, "Read worker.log"
            row = json.loads((args.run / "result.json").read_text())
            job.result("periodic-crease-audit", screen_passed=int(row["predeclared_curve_localization_screen_passed"]),
                seconds=record["wall_s"], peak_gib=record["peak_rss_bytes"] / 1024**3)
    assert all(workflow.digest(Path(path)) == digest for path, digest in hashes.items())


if __name__ == "__main__":
    main()
