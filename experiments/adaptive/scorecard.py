"""Deterministic surface scorecard: counts, triangle quality, two-sided distance, volume.

No timing: the GPU queue is busy, so only resolution-independent metrics are measured.
Run: uv run --no-sync python experiments/adaptive/scorecard.py [case ...]
"""

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.stats import qmc

BENCHMARK = Path(__file__).resolve().parents[3] / "microgen-benchmark/benchmarks/microgen"
sys.path.insert(0, str(BENCHMARK))
TAU = 2 * np.pi


def gyroid(x, y, z):
    a, b, c = TAU * x, TAU * y, TAU * z
    return (np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)) / 0.55


def schwarz_p(x, y, z):
    return (np.cos(TAU * x) + np.cos(TAU * y) + np.cos(TAU * z)) / 0.6


def sphere(x, y, z):
    return np.sqrt(x * x + y * y + z * z) - 0.35


AXES = (0.42, 0.25, 0.12)


def ellipsoid(x, y, z):
    # Scaled to an approximate distance near the surface; curvature varies 60x.
    a, b, c = AXES
    return (np.sqrt((x / a) ** 2 + (y / b) ** 2 + (z / c) ** 2) - 1) * c


CORNERS = np.array([[i, j, k] for i in (-0.5, 0.5) for j in (-0.5, 0.5) for k in (-0.5, 0.5)])


def lattice(x, y, z):
    # Eight capsule struts from the center to the corners, radius 0.12 (issue #4).
    d = None
    for corner in CORNERS:
        t = np.clip((x * corner[0] + y * corner[1] + z * corner[2]) / (corner @ corner), 0, 1)
        r = np.sqrt((x - t * corner[0]) ** 2 + (y - t * corner[1]) ** 2 + (z - t * corner[2]) ** 2)
        d = r if d is None else np.minimum(d, r)
    return d - 0.12


@dataclass(frozen=True)
class Case:
    name: str
    function: object
    levels: tuple  # wall levels; (-inf, 0) is the negative sublevel set
    volume: float | None = None  # analytic volume when known

    @property
    def band(self):
        return None if np.isinf(self.levels[0]) else self.levels


CASES = {
    "sphere": Case("sphere", sphere, (-np.inf, 0.0), 4 / 3 * np.pi * 0.35**3),
    "ellipsoid": Case("ellipsoid", ellipsoid, (-np.inf, 0.0), 4 / 3 * np.pi * np.prod(AXES)),
    "gyroid": Case("gyroid", gyroid, (-1.0, 1.0)),
    "schwarz-p": Case("schwarz-p", schwarz_p, (-1.0, 1.0)),
    "lattice": Case("lattice", lattice, (-np.inf, 0.0)),
}
LO, HI = -0.5, 0.5


def gradient(f, p, h=1e-6):
    g = np.empty_like(p)
    for a in range(3):
        e = np.zeros(3)
        e[a] = h
        g[:, a] = (f(*(p + e).T) - f(*(p - e).T)) / (2 * h)
    return g


def project(case, p, iterations=40):
    """Newton foot point on the nearest wall level; returns points and convergence."""
    f = case.function
    value = f(*p.T)
    levels = np.array([v for v in case.levels if np.isfinite(v)])
    level = levels[np.argmin(np.abs(value[:, None] - levels[None]), axis=1)]
    q = p.copy()
    active = np.arange(len(q))
    for _ in range(iterations):
        x = q[active]
        r = f(*x.T) - level[active]
        g = gradient(f, x)
        g2 = np.maximum((g * g).sum(1), 1e-30)
        step = (r / g2)[:, None] * g
        q[active] = x - step
        active = active[np.linalg.norm(step, axis=1) > 1e-13]
        if not len(active):
            break
    residual = np.abs(f(*q.T) - level) / np.maximum(np.linalg.norm(gradient(f, q), axis=1), 1e-30)
    return q, residual < 1e-10


def on_cap(points):
    return np.isclose(points, LO, atol=1e-9) | np.isclose(points, HI, atol=1e-9)


def triangle_quality(points, faces):
    a, b, c = (points[faces[:, i]] for i in range(3))
    area2 = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    lengths = np.stack([((b - a) ** 2).sum(1), ((c - b) ** 2).sum(1), ((a - c) ** 2).sum(1)], 1)
    shape = 2 * np.sqrt(3) * area2 / lengths.sum(1)
    angles = []
    for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
        e1, e2 = v - u, w - u
        cosine = (e1 * e2).sum(1) / np.sqrt((e1 * e1).sum(1) * (e2 * e2).sum(1))
        angles.append(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
    return shape, np.min(angles, axis=0)


def mesh_volume(points, faces):
    a, b, c = (points[faces[:, i]] for i in range(3))
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6)


def reference_volume(case, power=23):
    if case.volume is not None:
        return case.volume, 0.0
    cache = Path(__file__).with_name("reference-volumes.json")
    known = json.loads(cache.read_text()) if cache.exists() else {}
    if case.name in known:
        return tuple(known[case.name])
    # Four scrambled Sobol sets; their spread estimates the sampling error.
    estimates = []
    for seed in range(4):
        sobol = qmc.Sobol(3, scramble=True, seed=seed)
        inside = 0
        for _ in range(2 ** (power - 20)):
            x = sobol.random(2**20) - 0.5
            v = case.function(*x.T)
            inside += int(((v >= case.levels[0]) & (v <= case.levels[1])).sum())
        estimates.append(inside / 2**power)
    known[case.name] = (float(np.mean(estimates)), float(np.std(estimates) / 2))
    cache.write_text(json.dumps(known, indent=1))
    return tuple(known[case.name])


def wall_samples(case, count=2**17):
    """Points on the true walls inside the box, for the surface-to-mesh direction."""
    x = qmc.Sobol(3, scramble=True, seed=11).random(count) - 0.5
    q, converged = project(case, x)
    keep = converged & np.all((q >= LO - 1e-12) & (q <= HI + 1e-12), axis=1)
    return q[keep]


def distance_to_mesh(points, faces, queries):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

    poly = vtk.vtkPolyData()
    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(points), deep=True))
    poly.SetPoints(pts)
    cells = vtk.vtkCellArray()
    flat = np.hstack([np.full((len(faces), 1), 3), faces]).astype(np.int64).ravel()
    cells.SetCells(len(faces), numpy_to_vtkIdTypeArray(flat, deep=True))
    poly.SetPolys(cells)
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(poly)
    locator.BuildLocator()
    closest = [0.0, 0.0, 0.0]
    cell, sub, d2 = vtk.reference(0), vtk.reference(0), vtk.reference(0.0)
    out = np.empty(len(queries))
    for i, q in enumerate(queries):
        locator.FindClosestPoint(q, closest, cell, sub, d2)
        out[i] = np.sqrt(float(d2))
    return out


BARY = np.array(
    [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1 / 3, 1 / 3, 1 / 3]]
    + [[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]]
    + [[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]]
    + [[0.25, 0.75, 0], [0, 0.25, 0.75], [0.75, 0, 0.25]]
)


def score(case, points, faces, walls=None):
    """All metrics for one closed triangle surface."""
    points = np.asarray(points, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    corners = points[faces]
    cap = np.zeros(len(faces), dtype=bool)
    for a in range(3):
        for side in (LO, HI):
            cap |= np.all(np.isclose(corners[:, :, a], side, atol=1e-9), axis=1)
    wall = faces[~cap]
    samples = np.einsum("sk,fkd->fsd", BARY, points[wall]).reshape(-1, 3)
    foot, converged = project(case, samples)
    forward = np.linalg.norm(foot - samples, axis=1)
    walls = wall_samples(case) if walls is None else walls
    reverse = distance_to_mesh(points, faces, walls)
    shape, angle = triangle_quality(points, faces)
    reference, uncertainty = reference_volume(case)
    volume = mesh_volume(points, faces)
    return {
        "vertices": len(points),
        "triangles": len(faces),
        "wall_triangles": len(wall),
        "forward_max": float(forward.max()),
        "forward_p99": float(np.quantile(forward, 0.99)),
        "forward_unconverged": int((~converged).sum()),
        "reverse_max": float(reverse.max()),
        "reverse_p99": float(np.quantile(reverse, 0.99)),
        "hausdorff": float(max(forward.max(), reverse.max())),
        "volume_error": (volume - reference) / reference,
        "volume_reference_uncertainty": uncertainty / reference,
        "shape_min": float(shape.min()),
        "shape_p05": float(np.quantile(shape, 0.05)),
        "shape_mean": float(shape.mean()),
        "angle_min_deg": float(angle.min()),
        "below_30deg": float((angle < 30).mean()),
    }


# Mesh producers -------------------------------------------------------------


def flying_edges(case, n):
    from types import SimpleNamespace

    from microgen_benchmark import flying_edges as fe

    implicit = case.band is None
    proxy = SimpleNamespace(
        scale=1.0,
        field=case.function,
        implicit=implicit,
        band=(-1, 0) if implicit else case.levels,
    )
    points, faces, _, _ = fe.generate(proxy, n)
    return points, faces


def meshers_surface(case, n, quality=False):
    import meshers

    mesh = meshers.generate_surface(
        case.function,
        bounds=(LO, HI) * 3,
        cells=n,
        band=case.band,
        quality=quality,
        geometry_tolerance=0,
    )
    return mesh.points, mesh.triangles


def meshers_volume(case, n, passes=1):
    import meshers

    mesh = meshers.generate(
        case.function,
        bounds=(LO, HI) * 3,
        cells=n,
        band=case.band if case.band is not None else None,
        optimize_passes=passes,
        geometry_tolerance=1.0,
        threads=1,
    )
    used, faces = np.unique(mesh.surface, return_inverse=True)
    return mesh.points[used], faces.reshape(-1, 3)


PRODUCERS = {
    "flying-edges": flying_edges,
    "meshers-surface": meshers_surface,
    "meshers-surface-quality": lambda c, n: meshers_surface(c, n, quality=True),
    "meshers-volume": meshers_volume,
}


def main():
    names = sys.argv[1:] or list(CASES)
    grids = (16, 24, 32, 48)
    rows = []
    for name in names:
        case = CASES[name]
        walls = wall_samples(case)
        for producer, make in PRODUCERS.items():
            for n in grids:
                try:
                    points, faces = make(case, n)
                    row = {"case": name, "method": producer, "cells": n} | score(
                        case, points, faces, walls
                    )
                except Exception as error:  # noqa: BLE001 - failures are results
                    row = {"case": name, "method": producer, "cells": n, "error": str(error)}
                rows.append(row)
                with Path(__file__).with_name("baseline.jsonl").open("a") as log:
                    log.write(json.dumps(row) + "\n")
    out = Path(__file__).with_name("baseline.json")
    out.write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()


def tet_metrics(points, tets):
    """MMG quality and dihedral angles of every tetrahedron."""
    p = points[tets]
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    l2 = sum(((p[:, i] - p[:, j]) ** 2).sum(1) for i, j in edges)
    volume = np.einsum("ij,ij->i", p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])) / 6
    q = (12 * (3 * np.abs(volume)) ** (2 / 3) / l2) ** 1.5
    normals = []
    for i in range(4):
        f = [j for j in range(4) if j != i]
        n = np.cross(p[:, f[1]] - p[:, f[0]], p[:, f[2]] - p[:, f[0]])
        side = np.einsum("ij,ij->i", n, p[:, i] - p[:, f[0]])
        n *= -np.sign(side)[:, None]
        normals.append(n / np.linalg.norm(n, axis=1)[:, None])
    dihedral = np.full(len(tets), 180.0)
    for i in range(4):
        for j in range(i + 1, 4):
            c = np.einsum("ij,ij->i", normals[i], normals[j])
            dihedral = np.minimum(dihedral, 180 - np.degrees(np.arccos(np.clip(c, -1, 1))))
    return q, dihedral, volume


def score_volume(case, mesh, walls=None):
    used, faces = np.unique(mesh.surface, return_inverse=True)
    row = score(case, mesh.points[used], faces.reshape(-1, 3), walls)
    q, dihedral, volume = tet_metrics(mesh.points, mesh.tetrahedra)
    reference, _ = reference_volume(case)
    return row | {
        "tetrahedra": len(mesh.tetrahedra),
        "nodes": len(mesh.points),
        "tet_volume_error": float(volume.sum() / reference - 1),
        "q_min": float(q.min()),
        "q_p01": float(np.quantile(q, 0.01)),
        "q_mean": float(q.mean()),
        "q_below_0.3": float((q < 0.3).mean()),
        "dihedral_min_deg": float(dihedral.min()),
    }


@dataclass
class TetMesh:
    points: np.ndarray
    tetrahedra: np.ndarray
    surface: np.ndarray


def boundary_faces(tets):
    faces = np.concatenate([tets[:, [1, 2, 3]], tets[:, [0, 3, 2]], tets[:, [0, 1, 3]], tets[:, [0, 2, 1]]])
    keys = np.sort(faces, axis=1)
    _, index, counts = np.unique(keys, axis=0, return_index=True, return_counts=True)
    return faces[index[counts == 1]]


def mmg_adapt(case, seed, tolerance, max_length, hgrad=1.3):
    """MMG3D adaptation of a seed mesh with tagged wall and box patches."""
    import mmgpy
    import pyvista as pv

    points = np.asarray(seed.points, dtype=float)
    tets = np.asarray(seed.tetrahedra, dtype=np.int64)
    triangles = boundary_faces(tets)
    grid = pv.UnstructuredGrid(
        np.concatenate([np.c_[np.full(len(tets), 4), tets].ravel(), np.c_[np.full(len(triangles), 3), triangles].ravel()]),
        np.concatenate([np.full(len(tets), pv.CellType.TETRA), np.full(len(triangles), pv.CellType.TRIANGLE)]).astype(np.uint8),
        points,
    )
    native = mmgpy.from_pyvista(grid)
    corners = points[triangles]
    refs = np.ones(len(triangles), dtype=np.int64)
    if case.band is not None:
        value = case.function(*corners.mean(axis=1).T)
        refs = np.where(np.abs(value - case.levels[0]) < np.abs(value - case.levels[1]), 1, 2)
    for axis in range(3):
        for side, reference in ((LO, 3 + 2 * axis), (HI, 4 + 2 * axis)):
            refs[np.all(np.isclose(corners[:, :, axis], side, atol=1e-12), axis=1)] = reference
    native.set_triangles(triangles.astype(np.int32), refs)
    native.remesh(hausd=tolerance, hgrad=hgrad, hmax=max_length, verbose=-1)
    out = mmgpy.to_pyvista(native)
    tet = out.cells_dict[pv.CellType.TETRA].astype(np.int64)
    used, inverse = np.unique(tet, return_inverse=True)
    tet = inverse.reshape(-1, 4)
    pts = np.asarray(out.points)[used]
    return TetMesh(pts, tet, boundary_faces(tet))
