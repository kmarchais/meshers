# Copyright (c) 2026 Kevin Marchais

"""Mesh periodic implicit solids and measure the result independently of the mesher.

Every solid fills the centred unit cell and is periodic by construction. Each record holds
tetrahedron shape, wall distance in both directions, volume against a grid reference and the
exact matching of periodic cap nodes, next to the diagnostics the mesher reports about itself.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import importlib.metadata
import itertools
import json
import operator
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import meshers
import numpy as np
import pyvista as pv
from numpy.typing import NDArray

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from benchmark_hardware import hardware_snapshot

Coordinates = NDArray[np.float64]
Indices = NDArray[np.intp]
Field = Callable[[Coordinates], Coordinates]

SPECTRUM = Path(__file__).with_name("implicit-coverage") / "spinodoid-spectrum.csv"
PAIRS = tuple(itertools.combinations(range(4), 2))
DIFFERENCE_STEP = 1e-6
ROOT_BISECTIONS = 30
REFERENCE_GRID = 64
VOLUME_GRID = 160
CHUNK = 1 << 18
IMAGE_REACH = 1e-9
GOOD_QUALITY = 0.1
GOOD_DIHEDRALS = (10.0, 165.0)
GOOD_VOLUME_ERROR = 0.005
# Allowed drift from the baseline before a record counts as a regression.
QUALITY_DRIFT = 0.03
DISTANCE_DRIFT = 0.15
VOLUME_DRIFT = 0.002
SIZE_DRIFT = 0.15


def _wrapped(points: Coordinates, centre: tuple[float, float, float]) -> Coordinates:
    """Return the offsets from `centre` to its nearest periodic image."""
    delta = points - np.asarray(centre)
    return delta - np.round(delta)


def spheres(centres: Coordinates, radius: float) -> Field:
    """Return the signed distance to a sharp union of periodic spheres."""
    return lambda p: np.min(
        [np.linalg.norm(_wrapped(p, tuple(c)), axis=1) - radius for c in centres],
        axis=0,
    )


def fibre(axis: int, centre: tuple[float, float, float], radius: float) -> Field:
    """Return the signed distance to a periodic cylinder along `axis`."""
    others = [a for a in range(3) if a != axis]
    return lambda p: np.linalg.norm(_wrapped(p, centre)[:, others], axis=1) - radius


def union(*fields: Field) -> Field:
    """Return the sharp union of solids."""
    return lambda p: np.min([f(p) for f in fields], axis=0)


def complement(field: Field) -> Field:
    """Return the solid filling the void of `field`."""
    return lambda p: -field(p)


def gyroid(p: Coordinates) -> Coordinates:
    """Return the Gyroid field of the unit cell."""
    s, c = np.sin(2 * np.pi * p), np.cos(2 * np.pi * p)
    return (s * np.roll(c, -1, axis=1)).sum(axis=1)


def spinodoid(fraction: float) -> Field:
    """Return the archived cosine spectrum shifted to the requested solid fraction."""
    with SPECTRUM.open(encoding="utf8") as stream:
        rows = np.array([[float(v) for v in row.values()] for row in csv.DictReader(stream)])

    def raw(p: Coordinates) -> Coordinates:
        return (rows[:, 3] * np.cos(2 * np.pi * p @ rows[:, :3].T + rows[:, 4])).sum(axis=1)

    level = np.quantile(raw(np.random.default_rng(7).uniform(-0.5, 0.5, (200_000, 3))), fraction)
    return lambda p: raw(p) - level


def _cases() -> dict[str, Field]:
    random = np.random.default_rng(3).uniform(-0.5, 0.5, (8, 3))
    graze = 0.5 - 0.2 / np.sqrt(2)
    return {
        "sphere": spheres(np.zeros((1, 3)), 0.35),
        "sphere-corner": spheres(np.full((1, 3), 0.5), 0.35),
        "sphere-tangent-face": spheres(np.zeros((1, 3)), 0.49),
        "spheres-random": spheres(random, 0.17),
        "spheres-random-matrix": complement(spheres(random, 0.17)),
        "fibre": fibre(0, (0, 0.1, -0.2), 0.2),
        "fibre-through-edge": fibre(0, (0, 0.45, 0.42), 0.2),
        "fibre-grazing-edge": fibre(0, (0, graze + 0.004, graze + 0.004), 0.2),
        "fibre-inside-face": fibre(0, (0, 0, 0.5 - 0.2 - 0.003), 0.2),
        "fibre-across-face": fibre(0, (0, 0, 0.5 - 0.2 + 0.003), 0.2),
        "fibres-xyz": union(
            fibre(0, (0, 0.2, 0.2), 0.12),
            fibre(1, (-0.2, 0, -0.2), 0.12),
            fibre(2, (0.2, -0.2, 0), 0.12),
        ),
        "gyroid-skeleton": lambda p: gyroid(p) - 0.3,
        "gyroid-sheet": lambda p: abs(gyroid(p)) - 0.4,
        "gyroid-thin-sheet": lambda p: abs(gyroid(p)) - 0.15,
        "spinodoid-50": spinodoid(0.5),
        "spinodoid-30": spinodoid(0.3),
        "spinodoid-70": spinodoid(0.7),
    }


CASES = _cases()


def _chunked(field: Field, points: Coordinates) -> Coordinates:
    return np.concatenate([field(points[i : i + CHUNK]) for i in range(0, len(points), CHUNK)])


def gradient(field: Field, points: Coordinates) -> Coordinates:
    """Return the central-difference gradient of `field`."""
    step = DIFFERENCE_STEP * np.eye(3)
    return np.stack([(field(points + e) - field(points - e)) for e in step], axis=1) / (
        2 * DIFFERENCE_STEP
    )


def tet_shapes(points: Coordinates, tets: Indices) -> tuple[Coordinates, Coordinates, Coordinates]:
    """Return signed volume, MMG quality and the six dihedral angles of each tetrahedron."""
    p = points[tets]
    volume = np.linalg.det(p[:, 1:] - p[:, :1]) / 6
    lengths = sum(((p[:, i] - p[:, j]) ** 2).sum(axis=1) for i, j in PAIRS)
    quality = 72 * np.sqrt(3) * volume / np.maximum(lengths, 1e-300) ** 1.5
    normals = []
    for i in range(4):
        a, b, c = (p[:, j] for j in range(4) if j != i)
        n = np.cross(b - a, c - a)
        n *= np.sign(np.einsum("ij,ij->i", n, a - p[:, i]))[:, None]
        normals.append(n / np.maximum(np.linalg.norm(n, axis=1), 1e-300)[:, None])
    cosines = np.stack(
        [-np.einsum("ij,ij->i", normals[i], normals[j]) for i, j in PAIRS],
        axis=1,
    )
    return volume, quality, np.degrees(np.arccos(np.clip(cosines, -1, 1)))


def wall_facets(points: Coordinates, surface: Indices) -> Indices:
    """Return the boundary triangles that do not lie in a face of the cell."""
    p = points[surface]
    # Cap nodes sit on the cell faces exactly.
    cap = np.isin(abs(p), 0.5).all(axis=1).any(axis=1)
    return surface[~cap]


def mesh_to_wall(field: Field, points: Coordinates, wall: Indices) -> Coordinates:
    """Return, per wall facet, an upper bound on the distance from its samples to the wall.

    Each sample walks down the field gradient to a root. The length of that walk bounds the
    distance from above whether or not the field is a distance function.

    Returns:
        One distance per wall facet, infinite where a sample found no root.
    """
    p = points[wall]
    weights = np.array(
        [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [1, 1, 0], [0, 1, 1], [1, 0, 1]],
        dtype=float,
    )
    weights /= weights.sum(axis=1, keepdims=True)
    start = np.einsum("sk,nkd->nsd", weights, p).reshape(-1, 3)
    at = start.copy()
    for _ in range(12):
        value, slope = field(at), gradient(field, at)
        at -= (value / np.maximum((slope**2).sum(axis=1), 1e-300))[:, None] * slope
    distance = np.linalg.norm(at - start, axis=1)
    scale = np.maximum(np.linalg.norm(gradient(field, at), axis=1), 1e-300)
    distance[abs(field(at)) / scale > DIFFERENCE_STEP] = np.inf
    return distance.reshape(len(wall), -1).max(axis=1)


def wall_reference(field: Field, grid: int = REFERENCE_GRID) -> Coordinates:
    """Return wall points where the field changes sign along the lines of a regular grid."""
    axis = (np.arange(grid) + 0.5) / grid - 0.5
    nodes = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    values = _chunked(field, nodes)
    found = []
    for a in range(3):
        last = nodes[:, a] > axis[-1] - 0.5 / grid
        lo = nodes[~last]
        hi = lo.copy()
        hi[:, a] += 1 / grid
        lo_value = values[~last]
        hi_value = np.roll(values.reshape(grid, grid, grid), -1, axis=a).reshape(-1)[~last]
        crossing = (lo_value < 0) != (hi_value < 0)
        lo, hi, negative = lo[crossing], hi[crossing], lo_value[crossing] < 0
        for _ in range(ROOT_BISECTIONS):
            middle = (lo + hi) / 2
            same = (field(middle) < 0) == negative
            lo[same], hi[~same] = middle[same], middle[~same]
        found.append((lo + hi) / 2)
    return np.concatenate(found)


def wall_to_mesh(reference: Coordinates, points: Coordinates, wall: Indices) -> Coordinates:
    """Return the distance from each reference wall point to the meshed wall."""
    faces = np.hstack([np.full((len(wall), 1), 3), wall]).ravel()
    _, closest = pv.PolyData(points, faces).find_closest_cell(reference, return_closest_point=True)
    return np.linalg.norm(np.asarray(closest) - reference, axis=1)


def reference_volume(field: Field, grid: int = VOLUME_GRID) -> float:
    """Return the solid fraction of the unit cell sampled at the centres of a regular grid."""
    axis = (np.arange(grid) + 0.5) / grid - 0.5
    nodes = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    return float((_chunked(field, nodes) < 0).mean())


def cap_mismatch(points: Coordinates, surface: Indices) -> tuple[int, int, float]:
    """Return how the nodes of cap triangles miss their images on the opposite face.

    A node is unmatched when no node of the opposite cap lies within `IMAGE_REACH` of its
    translate. It is inexact when its image is there but not bitwise identical in the plane
    of the face.

    Returns:
        The unmatched count, the inexact count and the largest in-plane offset of an
        inexact node.
    """
    unmatched, inexact, offset = 0, 0, 0.0
    for a in range(3):
        others = [b for b in range(3) if b != a]
        faces = [
            points[np.unique(surface[np.isin(points[surface][:, :, a], side).all(axis=1)])][
                :,
                others,
            ]
            for side in (-0.5, 0.5)
        ]
        shared = {row.tobytes() for row in faces[0]} & {row.tobytes() for row in faces[1]}
        for face, other in zip(faces, faces[::-1], strict=True):
            loose = np.array([row for row in face if row.tobytes() not in shared])
            if not len(loose):
                continue
            gaps = (
                np.linalg.norm(loose[:, None] - other[None], axis=2).min(axis=1)
                if len(other)
                else np.full(len(loose), np.inf)
            )
            unmatched += int((gaps > IMAGE_REACH).sum())
            inexact += int((gaps <= IMAGE_REACH).sum())
            offset = max(offset, float(gaps[gaps <= IMAGE_REACH].max(initial=0.0)))
    return unmatched, inexact, offset


@dataclass
class Counter:
    """Count the points a field is asked to evaluate."""

    field: Field
    points: int = 0

    def __call__(self, x: Coordinates, y: Coordinates, z: Coordinates) -> Coordinates:
        """Return the field at separate coordinate arrays."""
        self.points += len(x)
        return self.field(np.stack([x, y, z], axis=1))


def _percentiles(values: Coordinates) -> dict[str, float]:
    low, first, median = np.percentile(values, [0, 1, 50])
    return {"minimum": float(low), "first_percentile": float(first), "median": float(median)}


def audit(field: Field, mesh: meshers.Mesh, reference: Coordinates, volume: float) -> dict:
    """Return the independent measurements of a meshed solid."""
    points = np.asarray(mesh.points)
    tets = np.asarray(mesh.tetrahedra)[:, :4]
    surface = np.asarray(mesh.surface)[:, :3]
    signed, quality, dihedral = tet_shapes(points, tets)
    wall = wall_facets(points, surface)
    forward = mesh_to_wall(field, points, wall)
    spread = np.percentile(forward, [10, 90])
    mismatch = cap_mismatch(points, surface)
    return {
        "tetrahedra": len(tets),
        "wall_facets": len(wall),
        "inverted": int((signed <= 0).sum()),
        "quality": _percentiles(quality),
        "dihedral_degrees": [float(dihedral.min()), float(dihedral.max())],
        "mesh_to_wall": float(forward.max()),
        "mesh_to_wall_spread": float(spread[1] / max(spread[0], 1e-300)),
        "wall_to_mesh": float(wall_to_mesh(reference, points, wall).max()),
        "relative_volume_error": float(signed.sum() / volume - 1),
        "unmatched_cap_nodes": mismatch[0],
        "inexact_cap_nodes": mismatch[1],
        "cap_node_offset": mismatch[2],
    }


def verdict(record: dict, tolerance: float) -> list[str]:
    """Return the names of the acceptance limits a measured mesh misses."""
    checks = {
        "inverted": record["inverted"] == 0,
        "quality": record["quality"]["minimum"] >= GOOD_QUALITY,
        "dihedral": GOOD_DIHEDRALS[0] <= record["dihedral_degrees"][0]
        and record["dihedral_degrees"][1] <= GOOD_DIHEDRALS[1],
        "mesh_to_wall": record["mesh_to_wall"] <= tolerance,
        "wall_to_mesh": record["wall_to_mesh"] <= tolerance,
        "volume": abs(record["relative_volume_error"]) <= GOOD_VOLUME_ERROR,
        "periodicity": record["unmatched_cap_nodes"] == 0,
        "exact_periodicity": record["inexact_cap_nodes"] == 0,
    }
    return [name for name, passed in checks.items() if not passed]


def suite() -> list[tuple[str, int, float]]:
    """Return the fixed (case, cells, tolerance) configurations of the regression suite.

    Every case runs at sixteen cells. The cases that depend most on where the lattice falls
    also run coarser and finer, and two smooth ones at the tolerance Microgen asks for.
    """
    sensitive = (
        "sphere-corner",
        "spheres-random",
        "fibre-grazing-edge",
        "fibre-across-face",
        "gyroid-thin-sheet",
        "spinodoid-50",
        "spinodoid-30",
    )
    return [
        *((name, 16, 0.02) for name in CASES),
        *((name, cells, 0.02) for name in sensitive for cells in (12, 24)),
        *((name, 24, 0.01) for name in ("sphere", "gyroid-sheet")),
    ]


def regressions(record: dict, reference: dict) -> list[str]:
    """Return what a record does worse than the baseline record of the same configuration."""
    if "rejected" in reference:
        return []
    if "rejected" in record:
        return ["now rejected: " + record["rejected"][:80]]
    found = [
        f"now misses {name}"
        for name in record["missed_limits"]
        if name not in reference["missed_limits"]
    ]
    tolerance = record["tolerance"]
    if record["quality"]["minimum"] < reference["quality"]["minimum"] - QUALITY_DRIFT:
        found.append("minimum quality")
    found.extend(
        key
        for key in ("mesh_to_wall", "wall_to_mesh")
        if record[key] > reference[key] + DISTANCE_DRIFT * tolerance
    )
    if (
        abs(record["relative_volume_error"])
        > abs(reference["relative_volume_error"]) + VOLUME_DRIFT
    ):
        found.append("volume")
    if record["tetrahedra"] > (1 + SIZE_DRIFT) * reference["tetrahedra"]:
        found.append("tetrahedra")
    found.extend(
        key
        for key in ("adaptation_abandoned", "defeaturing_moves")
        if record["mesher"][key] and not reference["mesher"][key]
    )
    return found


_key = operator.itemgetter("case", "cells", "adapt", "tolerance")


def measure(name: str, cells: int, *, adapt: bool, tolerance: float, quality: float) -> dict:
    """Return the record of one meshed case, or of its rejection by the mesher."""
    field = CASES[name]
    counter = Counter(field)
    record: dict = {"case": name, "cells": cells, "adapt": adapt, "tolerance": tolerance}
    start = time.perf_counter()
    try:
        mesh = meshers.generate(
            counter,
            bounds=(-0.5, 0.5) * 3,
            cells=cells,
            periodic=(True,) * 3,
            minimum_quality=quality,
            geometry_tolerance=tolerance,
            threads=1,
            adapt=adapt,
            compile_expression=False,
        )
    except meshers.MeshingError as error:
        record.update(rejected=str(error))
    # Generation only: the independent measures below are not the mesher's cost.
    record.update(seconds=time.perf_counter() - start, field_points=counter.points)
    if "rejected" not in record:
        record.update(
            audit(field, mesh, wall_reference(field), reference_volume(field)),
            mesher={
                key: mesh.diagnostics.get(key)
                for key in (
                    "minimum_mmg_quality",
                    "sampled_surface_error",
                    "sampled_coverage_error",
                    "adaptation_abandoned",
                    "balanced_wall_vertices",
                    "improved_elements",
                    "defeaturing_moves",
                    "refined_lattice_edges",
                    "removed_bridging_elements",
                    "reverted_snap_vertices",
                )
            },
        )
        record["missed_limits"] = verdict(record, tolerance)
    return record


def _line(record: dict) -> str:
    head = (
        f"{record['case']:22s} c={record['cells']:<3d} {'adapt' if record['adapt'] else 'grid ':5s}"
    )
    if "rejected" in record:
        return f"{head} REJECTED {record['rejected'][:90]}"
    q = record["quality"]
    low, high = record["dihedral_degrees"]
    return (
        f"{head} tets={record['tetrahedra']:<7d} q={q['minimum']:.3f}/{q['first_percentile']:.3f}"
        f"/{q['median']:.3f} dih={low:5.1f}-{high:5.1f} out={record['mesh_to_wall']:.4f}"
        f" in={record['wall_to_mesh']:.4f} vol={100 * record['relative_volume_error']:+.2f}%"
        f" moved={record['mesher']['defeaturing_moves']}"
        f"{' FALLBACK' if record['mesher']['adaptation_abandoned'] else ''}"
        f" {record['seconds']:.0f}s"
        f" {'ok' if not record['missed_limits'] else 'MISSES ' + ','.join(record['missed_limits'])}"
    )


def main() -> int:
    """Return a failing code when a mesh is rejected, misses a limit or regresses."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", choices=sorted(CASES), help="default: all")
    parser.add_argument("--cells", type=int, nargs="+", default=[16])
    parser.add_argument("--mode", choices=("grid", "adapt", "both"), default="grid")
    parser.add_argument("--tolerance", type=float, default=0.01)
    parser.add_argument("--quality", type=float, default=GOOD_QUALITY)
    parser.add_argument("--suite", action="store_true", help="run the regression configurations")
    parser.add_argument("--baseline", type=Path, help="fail on records worse than this file")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    modes = {"grid": [False], "adapt": [True], "both": [False, True]}[args.mode]
    configurations = (
        [c for c in suite() if not args.case or c[0] in args.case]
        if args.suite
        else [(name, cells, args.tolerance) for name in args.case or CASES for cells in args.cells]
    )
    reference = (
        {_key(r): r for r in json.loads(args.baseline.read_text(encoding="utf8"))}
        if args.baseline
        else {}
    )
    records, regressed = [], 0
    for (name, cells, tolerance), adapt in itertools.product(configurations, modes):
        record = measure(name, cells, adapt=adapt, tolerance=tolerance, quality=args.quality)
        records.append(record)
        if not args.baseline:
            worse = []
        elif _key(record) in reference:
            worse = regressions(record, reference[_key(record)])
        else:
            # A record with nothing to be compared against proves nothing.
            worse = ["no baseline record"]
        regressed += bool(worse)
        sys.stdout.write(
            _line(record) + (" REGRESSION " + "; ".join(worse) if worse else "") + "\n",
        )
        sys.stdout.flush()
    if args.output:
        with args.output.open("x", encoding="utf8") as stream:
            json.dump(records, stream, indent=2)
            stream.write("\n")
        # The machine and build the records were measured on, next to them.
        snapshot = {
            "date": dt.datetime.now(dt.UTC).date().isoformat(),
            "meshers": importlib.metadata.version("meshers"),
            "hardware": hardware_snapshot(),
        }
        machine = args.output.with_suffix(".machine.json")
        machine.write_text(json.dumps(snapshot, indent=2) + "\n", encoding="utf8")
    if args.baseline:
        return int(regressed > 0)
    return int(any("rejected" in r or r["missed_limits"] for r in records))


if __name__ == "__main__":
    raise SystemExit(main())
