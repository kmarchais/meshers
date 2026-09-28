"""Compare graded TPMS paths in fresh processes at matched grid-point counts."""

import argparse
import json
import time

import meshers
import numpy as np

from tpms_scaling import peak_megabytes


GRADE_AXES = "x"


def gyroid(x, y, z):
    x, y, z = (2 * np.pi * v for v in (x, y, z))
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def thickness(x, y, z):
    if GRADE_AXES == "x":
        return 0.6 + 0.1 * x
    return 0.6 + (0.1 / 3) * (x + y + z)


def normalized_field(x, y, z):
    return gyroid(x, y, z) / (0.5 * thickness(x, y, z))


def volume_constraints():
    return {
        "upper": lambda x, y, z: gyroid(x, y, z) - 0.5 * thickness(x, y, z),
        "lower": lambda x, y, z: -gyroid(x, y, z) - 0.5 * thickness(x, y, z),
    }


def minimum_angle(points, triangles):
    vertices = points[triangles]
    angles = []
    for axis in range(3):
        a = vertices[:, (axis + 1) % 3] - vertices[:, axis]
        b = vertices[:, (axis + 2) % 3] - vertices[:, axis]
        angles.append(
            np.degrees(
                np.arctan2(
                    np.linalg.norm(np.cross(a, b), axis=1),
                    np.einsum("ij,ij->i", a, b),
                )
            )
        )
    return float(np.min(angles))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mode",
        choices=(
            "microgen_surface",
            "microgen_volume",
            "microgen_meshers_volume",
            "meshers_surface",
            "meshers_volume",
            "meshers_volume_slab_tile",
        ),
    )
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--grid-points-per-cell", type=int, default=16)
    parser.add_argument("--polish-passes", type=int, default=10)
    parser.add_argument("--grade-axes", choices=("x", "xyz"), default="x")
    args = parser.parse_args()
    global GRADE_AXES
    GRADE_AXES = args.grade_axes
    periodic = tuple(axis not in GRADE_AXES for axis in "xyz")
    if args.mode == "meshers_volume_slab_tile" and GRADE_AXES != "x":
        parser.error("slab tiling requires a grade that is periodic in y and z")
    u = args.repeats
    count = args.grid_points_per_cell * u
    bounds = (-u / 2, u / 2) * 3
    start = time.perf_counter()
    if args.mode.startswith("microgen"):
        from microgen import Tpms
        from microgen.shape.surface_functions import gyroid as microgen_gyroid

        shape = Tpms(
            microgen_gyroid,
            offset=thickness,
            repeat_cell=u,
            resolution=args.grid_points_per_cell,
        )
        setup_seconds = time.perf_counter() - start
        if args.mode == "microgen_surface":
            result = shape.generate_surface_mesh()
            triangles = result.faces.reshape(-1, 4)[:, 1:]
            output = {
                "points": result.n_points,
                "elements": result.n_cells,
                "open_edges": result.n_open_edges,
                "minimum_angle": minimum_angle(result.points, triangles),
            }
        elif args.mode == "microgen_volume":
            result = shape._generate_legacy_volume_mesh()
            output = {
                "points": result.n_points,
                "elements": result.n_cells,
                "cell_types": np.unique(result.celltypes, return_counts=True)[0].tolist(),
            }
        else:
            result = shape.generate_meshers(
                periodic=periodic,
                minimum_quality=0,
                geometry_tolerance=0.01,
                optimize_passes=4,
                max_tetrahedra=10_000_000,
            )
            output = {
                "points": len(result.points),
                "elements": len(result.tetrahedra),
                "quality": result.diagnostics["minimum_mmg_quality"],
                "sampled_error": result.diagnostics["sampled_surface_error"],
            }
    else:
        setup_seconds = 0.0
        if args.mode == "meshers_surface":
            result = meshers.generate_surface(
                normalized_field,
                bounds=bounds,
                cells=count - 1,
                band=(-1, 1),
                periodic=periodic,
                polish_passes=args.polish_passes,
            )
            output = {
                "points": len(result.points),
                "elements": len(result.triangles),
                "minimum_angle": result.diagnostics["minimum_angle_degrees"],
                "actual_cells": result.diagnostics["background_cells"],
            }
        elif args.mode == "meshers_volume":
            result = meshers.generate_intersection(
                volume_constraints(),
                bounds=bounds,
                cells=count - 1,
                periodic=periodic,
                geometry_tolerance=0.01,
                minimum_quality=0,
                max_tetrahedra=10_000_000,
                optimize_passes=4,
            )
            output = {
                "points": len(result.points),
                "elements": len(result.tetrahedra),
                "quality": result.diagnostics["minimum_mmg_quality"],
                "sampled_error": result.diagnostics["sampled_surface_error"],
            }
        else:
            slab_bounds = (-u / 2, u / 2, -0.5, 0.5, -0.5, 0.5)
            slab = meshers.generate(
                normalized_field,
                bounds=slab_bounds,
                cells=(
                    count - 1,
                    args.grid_points_per_cell - 1,
                    args.grid_points_per_cell - 1,
                ),
                band=(-1, 1),
                periodic=periodic,
                geometry_tolerance=0.01,
                minimum_quality=0,
                max_tetrahedra=10_000_000,
                optimize_passes=4,
                threads=1,
            )
            generation_seconds = time.perf_counter() - start
            result = meshers.tile_periodic(
                slab,
                bounds=slab_bounds,
                repeats=(1, u, u),
            )
            output = {
                "points": len(result.points),
                "elements": len(result.tetrahedra),
                "slab_elements": len(slab.tetrahedra),
                "generation_seconds": round(generation_seconds, 3),
                "slab_quality": slab.diagnostics["minimum_mmg_quality"],
                "slab_error": slab.diagnostics["sampled_surface_error"],
            }
    output.update(
        mode=args.mode,
        repeats=u,
        grade_axes=GRADE_AXES,
        grid_points_per_cell=args.grid_points_per_cell,
        setup_seconds=round(setup_seconds, 3),
        total_seconds=round(time.perf_counter() - start, 3),
        peak_megabytes=peak_megabytes(),
    )
    print(json.dumps(output), flush=True)


if __name__ == "__main__":
    main()
