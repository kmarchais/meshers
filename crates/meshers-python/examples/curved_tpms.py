"""Generate a cylinder or torus TPMS with retained curves and matching ends."""

import argparse

import meshers
import numpy as np


def gyroid(x, y, z):
    a, b, c = (
        2 * np.pi * (x / 0.75 + 0.125),
        2 * np.pi * (y / 0.75 + 0.125),
        2 * np.pi * (z + 0.125),
    )
    return np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=["cylinder", "torus"])
    parser.add_argument("output", help="Output .vtkhdf filename")
    parser.add_argument("--cells", type=int, default=24)
    args = parser.parse_args()
    constraints = {
        "sheet_upper": lambda x, y, z: gyroid(x, y, z) - 0.55,
        "sheet_lower": lambda x, y, z: -gyroid(x, y, z) - 0.55,
        "wall": lambda x, y, _z: np.sqrt(x * x + y * y) - 0.75,
    }
    mapping, transforms = None, None
    if args.kind == "torus":
        radius = 9 / np.pi
        angle = np.pi / 3

        def mapping(x, y, z):
            return (radius + x) * np.cos(z / radius), (radius + x) * np.sin(z / radius), -y

        rotation = np.eye(4)
        rotation[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
        transforms = {2: rotation}
    mesh = meshers.generate_intersection(
        constraints,
        bounds=(-0.8, 0.8, -0.8, 0.8, 0, 3),
        cells=(args.cells, args.cells, round(args.cells * 5 / 3)),
        coordinate_map=mapping,
        periodic=(False, False, True),
        periodic_transforms=transforms,
        optimize_passes=12,
        geometry_tolerance=0.02,
    )
    mesh.write_vtkhdf(args.output)
    print(f"{len(mesh.tetrahedra):,} tetrahedra; {len(mesh.feature_edges):,} feature edges")
    print(mesh.diagnostics)


if __name__ == "__main__":
    main()
