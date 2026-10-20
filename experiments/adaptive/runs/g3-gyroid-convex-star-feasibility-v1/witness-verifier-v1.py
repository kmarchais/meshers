"""Extended-precision replay of a saved single-vertex tangent witness."""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

import microgen_workflow as workflow


def tangent_upper(points, cells, node, target, position, weights, origin, scale):
    coordinates = (np.asarray(points[cells], dtype=np.longdouble) - np.asarray(origin, dtype=np.longdouble)) / np.longdouble(scale)
    weights = np.asarray(weights, dtype=np.longdouble)
    position = np.asarray(position, dtype=np.longdouble)
    assert weights.min() >= 0 and abs(weights.sum() - 1) < 1e-12
    values, derivatives = [], []
    for vertices, cell in zip(coordinates, cells, strict=True):
        moving = int(np.flatnonzero(cell == node).item())
        squared_initial = sum(np.dot(vertices[first] - vertices[second], vertices[first] - vertices[second])
            for first, second in itertools.combinations(range(4), 2))
        def determinant(location):
            triangle = vertices.copy()
            triangle[moving] = location
            edges = triangle[1:] - triangle[0]
            return np.dot(edges[0], np.cross(edges[1], edges[2]))
        baseline = determinant(np.zeros(3, dtype=np.longdouble))
        determinant_gradient = np.array([determinant(direction) - baseline for direction in np.eye(3, dtype=np.longdouble)])
        moved = vertices.copy()
        moved[moving] = position
        squared = sum(np.dot(moved[first] - moved[second], moved[first] - moved[second])
            for first, second in itertools.combinations(range(4), 2))
        squared_gradient = sum(2 * (position - moved[other]) for other in range(4) if other != moving)
        coefficient = np.longdouble(str(target)) / (12 * np.sqrt(np.longdouble(3)))
        values.append(determinant(position) / squared_initial**np.longdouble("1.5")
            - coefficient * (squared / squared_initial)**np.longdouble("1.5"))
        derivatives.append(determinant_gradient / squared_initial**np.longdouble("1.5")
            - np.longdouble("1.5") * coefficient * np.sqrt(squared / squared_initial) * squared_gradient / squared_initial)
    value = weights @ np.asarray(values)
    gradient = weights @ np.asarray(derivatives)
    lower = (np.longdouble("-.5") - np.asarray(origin, dtype=np.longdouble)) / np.longdouble(scale)
    upper = (np.longdouble(".5") - np.asarray(origin, dtype=np.longdouble)) / np.longdouble(scale)
    return float(max(value + gradient @ (np.asarray(corner) - position)
        for corner in itertools.product(*zip(lower, upper, strict=True))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    result = json.loads((args.run / "result.json").read_text())
    manifest = json.loads((args.run / "manifest.json").read_text())
    assert all(workflow.digest(Path(path)) == expected for path, expected in manifest["inputs"].items())
    row = result["diagnostic"]
    source = Path(result["source"])
    mesh = np.load(source / "mesh.npz")
    points, cells, surface = (mesh[name] for name in ("points", "tetrahedra", "surface"))
    node = row["node"]
    selected = np.flatnonzero((cells == node).any(axis=1))
    assert len(selected) == row["incident_tets"]
    np.testing.assert_array_equal(points[node], row["physical_origin"])
    assert np.finfo(np.longdouble).nmant > np.finfo(float).nmant
    upper = tangent_upper(points, cells[selected], node, row["quality_target"], row["normalized_expansion_point"],
        row["weights"], row["physical_origin"], row["scale"])
    assert abs(upper - row["checked_tangent_upper_bound"]) < 1e-10
    obstructed = upper < -row["required_negative_margin"]
    assert obstructed == row["fixed_neighbor_obstruction"]
    boundary = np.zeros(len(points), dtype=bool)
    boundary[np.unique(surface)] = True
    active, neighbors = [], set()
    for local, weight in enumerate(row["weights"]):
        if weight <= 1e-12:
            continue
        cell = cells[selected[local]]
        opposite = cell[cell != node]
        free = [int(vertex) for vertex in opposite if not boundary[vertex]]
        neighbors.update(free)
        active.append({"element": int(selected[local]), "weight": weight,
            "fixed_opposite_vertices": opposite.tolist(), "other_interior_vertices": free})
    workflow.save(args.output, {"source": str(args.run), "source_result_sha256": workflow.digest(args.run / "result.json"),
        "verifier_sha256": workflow.digest(Path(__file__)), "passed": True,
        "longdouble_mantissa_bits": np.finfo(np.longdouble).nmant + 1, "extended_precision_upper_bound": upper,
        "checked_fixed_neighbor_obstruction": obstructed, "active_constraints": active,
        "other_interior_vertices_in_active_constraints": sorted(neighbors),
        "scope": "Independent scalar-triple-product and eight-box-corner replay; fixed-neighbor floating-point witness, not interval certification or coupled-relocation impossibility"})


if __name__ == "__main__":
    main()
