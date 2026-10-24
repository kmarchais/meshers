"""Locate fixed-sample reverse-distance witnesses without modifying the seed."""

import argparse
import json
from pathlib import Path
import shutil

import microgen_workflow as workflow

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.source, args.output = args.source.resolve(), args.output.resolve()
    assert not args.output.exists()
    _, _, _, _, pv, scorecard = workflow.imports()
    source_result = json.loads((args.source / "result.json").read_text())
    name = source_result["case"]
    paths = [Path(__file__), Path(workflow.__file__), Path(scorecard.__file__),
        args.source / "mesh.npz", args.source / "result.json",
        args.source / f"{name}-samples.npz", args.source / f"{name}-specimen.json",
        workflow.HERE / "runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json"]
    inputs = {str(path): workflow.digest(path) for path in paths}
    workflow.save(args.output.parent / "analysis-inputs.json", inputs)
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    # A known plane checks that returned cells, feet and distances agree.
    plane = pv.PolyData(np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]), [3, 0, 1, 2])
    cell, feet = plane.find_closest_cell([[.2, .2, .3], [.2, .2, -.4]], return_closest_point=True)
    np.testing.assert_array_equal(cell, [0, 0])
    np.testing.assert_allclose(feet, [[.2, .2, 0.], [.2, .2, 0.]], atol=1e-15)
    mesh = np.load(args.source / "mesh.npz")
    points, surface = mesh["points"], mesh["surface"]
    samples = np.load(args.source / f"{name}-samples.npz")
    queries = np.concatenate([samples["walls"], samples["caps"]])
    packed = np.column_stack((np.full(len(surface), 3), surface)).ravel()
    poly = pv.PolyData(points, packed)
    indices, closest = poly.find_closest_cell(queries, return_closest_point=True)
    distances = np.linalg.norm(closest - queries, axis=1)
    saved = source_result["geometry"]
    np.testing.assert_allclose(distances.max(), saved["sampled_reverse_max"], rtol=1e-12)
    np.testing.assert_allclose(distances, scorecard.distance_to_mesh(points, surface, queries), atol=1e-14, rtol=1e-12)
    case, _, _ = workflow.specimen(name)
    finite_levels = np.array([level for level in case.levels if np.isfinite(level)])
    rows = []
    for index in np.argsort(distances)[-32:][::-1]:
        face = surface[indices[index]]
        coordinates = points[face]
        gradient = scorecard.gradient(case.function, queries[index:index + 1])[0]
        edges = coordinates[[1, 2, 0]] - coordinates
        normal = np.cross(edges[0], -edges[2])
        normal /= np.linalg.norm(normal)
        level = finite_levels[np.argmin(abs(case.function(*queries[index]) - finite_levels))]
        rows.append({"sample_index": int(index), "sample_kind": "wall" if index < len(samples["walls"]) else "cap",
            "distance": float(distances[index]), "query": queries[index].tolist(),
            "query_field_residual": float(case.function(*queries[index]) - level),
            "closest_point": closest[index].tolist(), "surface_face": int(indices[index]),
            "node_ids": face.tolist(), "coordinates": coordinates.tolist(),
            "node_field_residuals": (case.function(*coordinates.T) - level).tolist(),
            "maximum_edge_length": float(np.linalg.norm(edges, axis=1).max()),
            "query_distance_to_box": float((.5 - abs(queries[index])).min()),
            "closest_face_is_cap": bool(np.any(np.all(workflow.on_box_plane(coordinates, -.5), axis=0)
                | np.all(workflow.on_box_plane(coordinates, .5), axis=0))),
            "absolute_normal_alignment": float(abs(normal @ gradient) / np.linalg.norm(gradient))})
    geometry = workflow.geometry(args.source, name, points, mesh["tetrahedra"], surface, detail_count=32)
    for key in ("sampled_forward_max", "sampled_reverse_max", "projection_failed", "relative_volume_error"):
        np.testing.assert_allclose(geometry[key], saved[key], rtol=1e-12)
    workflow.save(args.output, {"source": str(args.source), "case": name,
        "sample_count": len(queries), "wall_samples": len(samples["walls"]),
        "reverse_maximum": float(distances.max()), "reverse_samples_above_gate": int((distances > .002).sum()),
        "wall_samples_above_gate": int((distances[:len(samples["walls"])] > .002).sum()),
        "cap_samples_above_gate": int((distances[len(samples["walls"]):] > .002).sum()),
        "reverse_distance_quantiles": dict(zip(("p50", "p90", "p99", "p999", "maximum"),
            np.quantile(distances, [.5, .9, .99, .999, 1]).tolist(), strict=True)),
        "worst_reverse_samples": rows, "forward_diagnostics": geometry["worst_boundary_samples"],
        "projection_failed": geometry["projection_failed"], "known_plane_control_passed": True,
        "mesh_modified": False, "samples_regenerated": False, "qualified": False,
        "mmg_calls": [], "stiffness_solve_performed": False})
    assert all(workflow.digest(Path(path)) == expected for path, expected in inputs.items())
    print(f"Reverse distance {distances.max():.9g}; {(distances > .002).sum()} fixed samples exceed gate")


if __name__ == "__main__":
    main()
