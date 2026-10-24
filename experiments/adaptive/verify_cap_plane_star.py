"""Independent Decimal replay of cap-only supporting-plane certificates."""

import argparse
from decimal import Decimal, localcontext
import json
from pathlib import Path
import shutil
import subprocess
import sys

import microgen_workflow as workflow

import numpy as np


def replay(points, diagnostic, target):
    triangles = diagnostic["incident_cap_triangles"]
    assert all(weight == 0 for weight in diagnostic["weights"][len(triangles):]), "Cap-only certificate required"
    with localcontext() as context:
        context.prec = 80
        number = lambda value: Decimal.from_float(float(value))
        axes = [axis for axis in range(3) if axis != diagnostic["axis"]]
        scale = number(diagnostic["scale"])
        position = list(map(number, diagnostic["normalized_expansion_point"]))
        weights = list(map(number, diagnostic["weights"][:len(triangles)]))
        assert min(weights) >= 0 and sum(weights) > 0
        total = sum(weights)
        weights = [weight / total for weight in weights]
        value, gradient = Decimal(0), [Decimal(0), Decimal(0)]
        for triangle, weight in zip(triangles, weights, strict=True):
            if not weight:
                continue
            origin = [number(points[triangle[0], axis]) for axis in axes]
            first, second = [[(number(points[node, axis]) - center) / scale
                for axis, center in zip(axes, origin, strict=True)] for node in triangle[1:]]
            a = [x - y for x, y in zip(first, position, strict=True)]
            b = [x - y for x, y in zip(second, position, strict=True)]
            determinant = a[0] * b[1] - a[1] * b[0]
            squared = sum(x * x for x in [*a, *b]) + sum((x - y)**2 for x, y in zip(first, second, strict=True))
            value += weight * (2 * Decimal(3).sqrt() * determinant - target * squared)
            derivative = [first[1] - second[1], second[0] - first[0]]
            for axis in range(2):
                gradient[axis] += weight * (2 * Decimal(3).sqrt() * derivative[axis]
                    - target * (4 * position[axis] - 2 * (first[axis] + second[axis])))
        lower = [max((Decimal("-.5") - number(points[node, axis])) / scale for node in diagnostic["nodes"])
            for axis in axes]
        upper = [min((Decimal(".5") - number(points[node, axis])) / scale for node in diagnostic["nodes"])
            for axis in axes]
        return value - sum(g * x for g, x in zip(gradient, position, strict=True)) + sum(
            max(g * low, g * high) for g, low, high in zip(gradient, lower, upper, strict=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    shutil.copyfile(__file__, args.output.parent / "analysis-driver.py")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "experiments/adaptive/test_verify_cap_plane_star.py"],
        check=True, timeout=30)
    inputs = json.loads((args.run / "analysis-inputs.json").read_text())
    assert all(workflow.digest(Path(path)) == expected for path, expected in inputs.items())
    result = json.loads((args.run / "result.json").read_text())
    source = Path(result["source"])
    mesh = np.load(source / "mesh.npz")
    points = mesh["points"]
    from cap_corner_probe import polygon_corners
    from cap_quality_bound import caps

    faces, _ = caps(mesh)
    rows = []
    for diagnostic in result["diagnostics"]:
        checked = replay(points, diagnostic, Decimal.from_float(diagnostic["target_cap_quality"]))
        assert abs(float(checked) - diagnostic["checked_tangent_upper_bound"]) < 1e-10
        gate = replay(points, diagnostic, Decimal(".1"))
        axis = diagnostic["axis"]
        axes = [other for other in range(3) if other != axis]
        boundary = {}
        for side in (-.5, .5):
            selected = faces[np.all(workflow.on_box_plane(points[faces, axis], side), axis=1)]
            corners, ambiguous = polygon_corners(points[:, axes], selected)
            boundary[side] = {row["node"] for row in corners} | set(ambiguous)
        active, free_neighbors = [], set()
        for triangle, weight in zip(diagnostic["incident_cap_triangles"], diagnostic["weights"], strict=False):
            if weight <= 1e-12:
                continue
            side = float(points[triangle[0], axis])
            side = -.5 if side < 0 else .5
            free = [node for node in triangle[1:] if node not in boundary[side]]
            free_neighbors.update(free)
            active.append({"triangle": triangle, "weight": weight, "other_free_cap_vertices": free})
        rows.append({"nodes": diagnostic["nodes"], "target_replay_upper_bound_decimal80": str(checked),
            "gate_target": .1, "gate_replay_upper_bound_decimal80": str(gate),
            "checked_gate_obstruction": gate < Decimal("-1e-8"), "active_cap_constraints": active,
            "other_free_cap_vertices_in_support": sorted(free_neighbors)})
    workflow.save(args.output, {"source": str(args.run), "source_result_sha256": workflow.digest(args.run / "result.json"),
        "source_mesh_sha256": workflow.digest(source / "mesh.npz"), "verifier_sha256": workflow.digest(Path(__file__)),
        "decimal_precision": 80, "weights_normalized_in_decimal": True, "passed": True, "diagnostics": rows,
        "scope": "Independent cap-area and squared-edge derivatives and box support at the original 0.1 gate; fixed-neighbor periodic translation only, not coupled cap motion or interval certification"})
    print("Verified cap-only tangent certificate at the unchanged quality gate with 80-digit arithmetic")


if __name__ == "__main__":
    main()
