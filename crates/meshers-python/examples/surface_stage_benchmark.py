"""Measure the experimental surface stages on periodic TPMS cell bands."""

import argparse
import time

import meshers
import numpy as np


def split_p(x, y, z):
    x, y, z = (2 * np.pi * v for v in (x, y, z))
    return (
        1.1
        * (
            np.sin(2 * x) * np.cos(y) * np.sin(z)
            + np.sin(2 * y) * np.cos(z) * np.sin(x)
            + np.sin(2 * z) * np.cos(x) * np.sin(y)
        )
        - 0.2
        * (
            np.cos(2 * x) * np.cos(2 * y)
            + np.cos(2 * y) * np.cos(2 * z)
            + np.cos(2 * z) * np.cos(2 * x)
        )
        - 0.4 * (np.cos(2 * x) + np.cos(2 * y) + np.cos(2 * z))
    )


def gyroid(x, y, z):
    x, y, z = (2 * np.pi * v for v in (x, y, z))
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def matching_caps(surface):
    pairs = []
    for axis in range(3):
        faces = []
        for side in range(2):
            triangles = surface.triangles[surface.labels == 2 + 2 * axis + side]
            points = surface.points[triangles].copy()
            points[:, :, axis] = 0
            vertices = np.rint(points * 1e9).astype(np.int64)
            faces.append({tuple(sorted(map(tuple, face))) for face in vertices})
        pairs.append(len(faces[0]) == len(faces[1]) and faces[0] == faces[1])
    return pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("geometry", choices=("gyroid", "split_p"))
    parser.add_argument("--cells", type=int, default=24)
    parser.add_argument("--smooth", type=int, default=10)
    parser.add_argument("--improve", type=int, default=12)
    parser.add_argument("--polish", type=int, default=40)
    parser.add_argument("--periodic", action="store_true")
    parser.add_argument("--weak", action="store_true")
    args = parser.parse_args()
    field = meshers.compile_field(globals()[args.geometry])
    t = time.perf_counter()
    surface = meshers.generate_surface(
        field,
        bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
        cells=args.cells,
        band=(-0.25, 0.25),
        periodic=(args.periodic,) * 3,
        smoothing_iterations=args.smooth,
        improvement_rounds=args.improve,
        polish_passes=args.polish,
    )
    print(
        args.geometry,
        f"{args.cells}->{surface.diagnostics['background_cells']}",
        args.smooth,
        args.improve,
        args.polish,
        args.periodic,
        len(surface.points),
        len(surface.triangles),
        f"native={surface.diagnostics['seconds']:.3f}s",
        f"total={time.perf_counter() - t:.3f}s",
        f"min_angle={surface.diagnostics['minimum_angle_degrees']:.3f}",
        f"matching_caps={matching_caps(surface)}",
        flush=True,
    )
    vertices = surface.points[surface.triangles]
    angles = []
    for axis in range(3):
        a = vertices[:, (axis + 1) % 3] - vertices[:, axis]
        b = vertices[:, (axis + 2) % 3] - vertices[:, axis]
        angles.append(
            np.degrees(
                np.arctan2(
                    np.linalg.norm(np.cross(a, b), axis=1), np.einsum("ij,ij->i", a, b)
                )
            )
        )
    weakest = np.min(angles, axis=0)
    print(
        "labels:",
        {
            int(label): (
                round(float(np.min(weakest[surface.labels == label])), 2),
                int(np.count_nonzero(weakest[surface.labels == label] < 10)),
            )
            for label in np.unique(surface.labels)
        },
        flush=True,
    )
    if args.weak:
        order = np.argsort(weakest)
        for fi in order[:12]:
            face = surface.triangles[fi]
            lengths = np.linalg.norm(
                surface.points[face] - surface.points[np.roll(face, 1)], axis=1
            )
            masks = [
                {
                    int(v)
                    for v in surface.labels[np.any(surface.triangles == vertex, axis=1)]
                }
                for vertex in face
            ]
            print(
                int(fi),
                int(surface.labels[fi]),
                round(float(weakest[fi]), 3),
                np.round(lengths, 5),
                np.round(surface.points[face], 5),
                masks,
                flush=True,
            )


if __name__ == "__main__":
    main()
