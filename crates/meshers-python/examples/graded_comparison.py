"""Compare TPMS paths in fresh processes at matched grid-point counts."""

import argparse
import json
import time

import meshers
import numpy as np
from tpms_scaling import peak_megabytes, split_p

GRADE_AXES = "x"
GEOMETRY = "gyroid"
UNIFORM_THICKNESS = 0.5


def gyroid(x, y, z):
    x, y, z = (2 * np.pi * v for v in (x, y, z))
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def thickness(x, y, z):
    if GRADE_AXES == "none":
        return UNIFORM_THICKNESS
    if GRADE_AXES == "x":
        return 0.6 + 0.1 * x
    return 0.6 + (0.1 / 3) * (x + y + z)


def field_value(x, y, z):
    return (gyroid if GEOMETRY == "gyroid" else split_p)(x, y, z)


def normalized_field(x, y, z):
    return field_value(x, y, z) / (0.5 * thickness(x, y, z))


def volume_constraints():
    return {
        "upper": lambda x, y, z: field_value(x, y, z) - 0.5 * thickness(x, y, z),
        "lower": lambda x, y, z: -field_value(x, y, z) - 0.5 * thickness(x, y, z),
    }


def triangle_quality(points, triangles):
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
    values = np.concatenate(angles)
    return {
        "minimum_angle": float(np.min(values)),
        "angle_p01": float(np.quantile(values, 0.01)),
    }


def tetra_quality(points, tetrahedra):
    vertices = points[tetrahedra]
    determinant = np.linalg.det(vertices[:, 1:] - vertices[:, :1])
    edge_sum = sum(
        np.sum((vertices[:, i] - vertices[:, j]) ** 2, axis=1)
        for i in range(4)
        for j in range(i + 1, 4)
    )
    quality = np.sqrt(432 * determinant**2 / edge_sum**3)
    return {
        "minimum_mmg_quality": float(np.min(quality)),
        "quality_p01": float(np.quantile(quality, 0.01)),
    }


def face_keys(points, triangles):
    coordinates = np.rint(points * 1e9).astype(np.int64)
    return {tuple(sorted(map(tuple, coordinates[face]))) for face in triangles}


def periodic_mismatch(points, triangles, periodic):
    """Count unmatched boundary nodes and cap triangles after box translation."""
    node_mismatch = []
    triangle_mismatch = []
    for axis, enabled in enumerate(periodic):
        if not enabled:
            continue
        transverse = [i for i in range(3) if i != axis]
        lo, hi = np.min(points[:, axis]), np.max(points[:, axis])

        def cap(bound):
            nodes = points[np.abs(points[:, axis] - bound) < 1e-8][:, transverse]
            node_keys = set(map(tuple, np.rint(nodes * 1e8).astype(np.int64)))
            faces = triangles[
                np.all(np.abs(points[triangles, axis] - bound) < 1e-8, axis=1)
            ]
            face_keys_2d = {
                tuple(
                    sorted(
                        map(
                            tuple,
                            np.rint(points[face][:, transverse] * 1e8).astype(np.int64),
                        )
                    )
                )
                for face in faces
            }
            return node_keys, face_keys_2d

        low_nodes, low_faces = cap(lo)
        high_nodes, high_faces = cap(hi)
        node_mismatch.append(len(low_nodes ^ high_nodes))
        triangle_mismatch.append(len(low_faces ^ high_faces))
    return {
        "periodic_node_mismatch": node_mismatch,
        "periodic_triangle_mismatch": triangle_mismatch,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mode",
        choices=(
            "microgen_surface",
            "microgen_volume",
            "microgen_mmg_volume",
            "microgen_mmgpy_volume",
            "microgen_meshers_volume",
            "meshers_surface",
            "meshers_volume",
            "meshers_volume_slab_tile",
        ),
    )
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--grid-points-per-cell", type=int, default=16)
    parser.add_argument("--polish-passes", type=int, default=10)
    parser.add_argument("--grade-axes", choices=("none", "x", "xyz"), default="x")
    parser.add_argument("--geometry", choices=("gyroid", "split_p"), default="gyroid")
    parser.add_argument("--uniform-thickness", type=float, default=0.5)
    parser.add_argument("--optimize-passes", type=int, default=4)
    parser.add_argument("--snap", type=float, default=0.2)
    parser.add_argument("--minimum-quality", type=float, default=0)
    args = parser.parse_args()
    global GRADE_AXES, GEOMETRY, UNIFORM_THICKNESS
    GRADE_AXES = args.grade_axes
    GEOMETRY = args.geometry
    UNIFORM_THICKNESS = args.uniform_thickness
    if not 0 < UNIFORM_THICKNESS < 2:
        parser.error("uniform thickness must lie between 0 and 2")
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
        from microgen.shape.surface_functions import split_p as microgen_split_p

        if args.mode == "microgen_mmg_volume":
            from microgen.remesh import remesh_keeping_boundaries_for_fem
        if args.mode == "microgen_mmgpy_volume":
            import mmgpy
            from microgen import BoxMesh

    import_seconds = time.perf_counter() - start
    runtime_start = time.perf_counter()
    if args.mode.startswith("microgen"):
        shape = Tpms(
            microgen_gyroid if GEOMETRY == "gyroid" else microgen_split_p,
            offset=UNIFORM_THICKNESS if GRADE_AXES == "none" else thickness,
            repeat_cell=u,
            resolution=args.grid_points_per_cell,
        )
        setup_seconds = time.perf_counter() - runtime_start
        if args.mode == "microgen_surface":
            result = shape.generate_surface_mesh()
            triangles = result.faces.reshape(-1, 4)[:, 1:]
            output = {
                "points": result.n_points,
                "elements": result.n_cells,
                "open_edges": result.n_open_edges,
                **triangle_quality(result.points, triangles),
            }
        elif args.mode in (
            "microgen_volume",
            "microgen_mmg_volume",
            "microgen_mmgpy_volume",
        ):
            result = shape._generate_legacy_volume_mesh()
            legacy_seconds = time.perf_counter() - runtime_start
            if args.mode == "microgen_mmg_volume":
                result = remesh_keeping_boundaries_for_fem(
                    result, periodic=GRADE_AXES == "none"
                )
                output = {
                    "points": result.n_points,
                    "elements": result.n_cells,
                    "legacy_seconds": round(legacy_seconds, 3),
                    "mmg_seconds": round(
                        time.perf_counter() - runtime_start - legacy_seconds, 3
                    ),
                    **tetra_quality(result.points, result.cells_dict[10]),
                }
            elif args.mode == "microgen_mmgpy_volume":
                box = BoxMesh.from_pyvista(result.triangulate())
                boundary, _ = box.boundary_elements(box.rve)
                merged = box.to_pyvista().merge(boundary)
                mmg_mesh = mmgpy.from_pyvista(merged)
                required = face_keys(mmg_mesh.get_vertices(), mmg_mesh.get_triangles())
                mmg_mesh.set_required_triangles(
                    np.arange(boundary.n_cells, dtype=np.int32)
                )
                mmg_result = mmg_mesh.remesh(verbose=-1)
                if mmg_result["return_code"] != 0:
                    raise RuntimeError(f"mmgpy remeshing failed: {mmg_result}")
                output = {
                    "points": len(mmg_mesh.get_vertices()),
                    "elements": len(mmg_mesh.get_tetrahedra()),
                    "required_faces_preserved": required.issubset(
                        face_keys(mmg_mesh.get_vertices(), mmg_mesh.get_triangles())
                    ),
                    "legacy_seconds": round(legacy_seconds, 3),
                    "mmgpy_seconds": round(
                        time.perf_counter() - runtime_start - legacy_seconds, 3
                    ),
                    **tetra_quality(mmg_mesh.get_vertices(), mmg_mesh.get_tetrahedra()),
                }
            else:
                output = {
                    "points": result.n_points,
                    "elements": result.n_cells,
                    "cell_types": np.unique(result.celltypes, return_counts=True)[
                        0
                    ].tolist(),
                    "legacy_seconds": round(legacy_seconds, 3),
                }
        else:
            result = shape.generate_meshers(
                periodic=periodic,
                minimum_quality=args.minimum_quality,
                geometry_tolerance=0.01,
                optimize_passes=args.optimize_passes,
                snap=args.snap,
                max_tetrahedra=10_000_000,
            )
            output = {
                "points": len(result.points),
                "elements": len(result.tetrahedra),
                "quality": result.diagnostics["minimum_mmg_quality"],
                "sampled_error": result.diagnostics["sampled_surface_error"],
                **tetra_quality(result.points, result.tetrahedra),
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
                "actual_cells": result.diagnostics["background_cells"],
                **triangle_quality(result.points, result.triangles),
            }
        elif args.mode == "meshers_volume":
            options = dict(
                bounds=bounds,
                cells=count - 1,
                periodic=periodic,
                geometry_tolerance=0.01,
                minimum_quality=args.minimum_quality,
                max_tetrahedra=10_000_000,
                optimize_passes=args.optimize_passes,
                snap=args.snap,
            )
            if GRADE_AXES == "none":
                result = meshers.generate(
                    gyroid if GEOMETRY == "gyroid" else split_p,
                    band=(-0.5 * UNIFORM_THICKNESS, 0.5 * UNIFORM_THICKNESS),
                    **options,
                )
            else:
                result = meshers.generate_intersection(volume_constraints(), **options)
            output = {
                "points": len(result.points),
                "elements": len(result.tetrahedra),
                "quality": result.diagnostics["minimum_mmg_quality"],
                "sampled_error": result.diagnostics["sampled_surface_error"],
                "initial_quality": result.diagnostics.get(
                    "initial_minimum_mmg_quality"
                ),
                "accepted_vertex_moves": result.diagnostics.get(
                    "accepted_vertex_moves"
                ),
                "accepted_reconnections": result.diagnostics.get(
                    "accepted_reconnections"
                ),
                "quality_optimization_passes": result.diagnostics.get(
                    "quality_optimization_passes"
                ),
                **tetra_quality(result.points, result.tetrahedra),
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
                minimum_quality=args.minimum_quality,
                max_tetrahedra=10_000_000,
                optimize_passes=args.optimize_passes,
                snap=args.snap,
                threads=1,
            )
            generation_seconds = time.perf_counter() - runtime_start
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
    runtime_seconds = round(time.perf_counter() - runtime_start, 3)
    if any(periodic):
        if args.mode == "microgen_mmgpy_volume":
            check_points = mmg_mesh.get_vertices()
            check_triangles = mmg_mesh.get_triangles()
        elif args.mode in ("microgen_volume", "microgen_mmg_volume"):
            check_surface = result.extract_surface(
                algorithm="dataset_surface"
            ).triangulate()
            check_points = check_surface.points
            check_triangles = check_surface.faces.reshape(-1, 4)[:, 1:]
        elif args.mode == "microgen_surface":
            check_points = result.points
            check_triangles = triangles
        else:
            check_points = result.points
            check_triangles = (
                result.triangles if args.mode == "meshers_surface" else result.surface
            )
        output.update(periodic_mismatch(check_points, check_triangles, periodic))
        if args.mode in ("microgen_meshers_volume", "meshers_volume"):
            output["periodic_pair_counts"] = [
                len(pairs) for pairs in result.periodic_pairs
            ]
    output.update(
        mode=args.mode,
        geometry=GEOMETRY,
        uniform_thickness=UNIFORM_THICKNESS if GRADE_AXES == "none" else None,
        repeats=u,
        grade_axes=GRADE_AXES,
        grid_points_per_cell=args.grid_points_per_cell,
        setup_seconds=round(setup_seconds, 3),
        import_seconds=round(import_seconds, 3),
        runtime_seconds=runtime_seconds,
        total_seconds=round(time.perf_counter() - start, 3),
        peak_megabytes=peak_megabytes(),
    )
    print(json.dumps(output), flush=True)


if __name__ == "__main__":
    main()
