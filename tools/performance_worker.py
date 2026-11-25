# Copyright (c) 2026 Kevin Marchais

"""Measure an installed wheel after the benchmark driver configures its numerical threads."""

from __future__ import annotations

import hashlib
import importlib
import importlib.metadata
import inspect
import math
import platform
import sys
import time
from itertools import combinations_with_replacement
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict, cast

import meshers
import numpy as np

if TYPE_CHECKING:
    from meshers import Mesh
    from numpy.typing import NDArray

    FloatArray = NDArray[np.float64]
    Indices = NDArray[np.int64]

if sys.platform == "win32":
    try:
        import psutil
    except ModuleNotFoundError as error:
        if error.name != "psutil":
            raise
        psutil = None
else:
    import resource

CAP = 250_000
_MAX_JACOBIAN_DEPTH = 12
_TET4_NODES = 4
_TRI3_NODES = 3
# Face opposite each tetrahedron corner, in Tri6 order (ab, bc, ca).
_TET_FACES = np.array([
    [1, 2, 3, 5, 9, 8],
    [0, 3, 2, 7, 9, 6],
    [0, 1, 3, 4, 8, 7],
    [0, 2, 1, 6, 5, 4],
])
_TRIANGLES_PER_BOUNDARY_EDGE = 2
_MIN_CELL_QUALITY = 0.01
_MAX_SURFACE_ERROR = 0.1
_INNER_CYLINDER, _BOTTOM_PLANE, _START_PLANE, _END_PLANE = 1, 3, 5, 6


class _FixtureOptions(TypedDict):
    cells: int
    optimize_passes: int
    geometry_tolerance: float
    max_tetrahedra: int


class _VolumeOptions(_FixtureOptions, total=False):
    order: int
    adapt: bool
    min_length: float
    max_length: float
    adapt_tolerance: float
    gradation: float
    compile_expression: bool
    compile: bool


def _check_curved_jacobians(corners: FloatArray) -> None:
    """Bound the cubic Jacobian determinant over each reference tetrahedron.

    The four physical Jacobian matrices interpolate affinely. Cubic Bernstein
    coefficients bound their determinant over the whole tetrahedron. Subdivide
    inconclusive cells along the edge with the largest Jacobian variation.

    A fold or an unresolved bound raises RuntimeError.
    """
    triples = list(combinations_with_replacement(range(4), 3))
    powers = np.array([[triple.count(i) for i in range(4)] for triple in triples])
    bary = powers / 3
    factors = np.array([6 / math.prod(math.factorial(int(n)) for n in row) for row in powers])
    basis = factors * np.prod(bary[:, None, :] ** powers, axis=2)
    to_coefficients = np.linalg.inv(basis).T

    def check(batch: FloatArray, depth: int) -> None:
        matrices = np.einsum("pi,tijk->tpjk", bary, batch)
        values = np.linalg.det(matrices)
        if not np.isfinite(values).all() or np.any(values <= 0):
            msg = "Nonpositive sampled curved Jacobian"
            raise RuntimeError(msg)
        coefficients = values @ to_coefficients
        scale = np.linalg.norm(batch, axis=2).max(axis=1).prod(axis=1)
        ambiguous = np.any(coefficients <= 128 * np.finfo(float).eps * scale[:, None], axis=1)
        for cell in batch[ambiguous]:
            if depth >= _MAX_JACOBIAN_DEPTH:
                msg = "Unresolved curved Jacobian positivity"
                raise RuntimeError(msg)
            i, j = np.unravel_index(
                np.argmax(np.linalg.norm(cell[:, None] - cell[None, :], axis=(2, 3))),
                (4, 4),
            )
            children = np.stack([cell, cell])
            midpoint = (cell[i] + cell[j]) / 2
            children[0, i] = midpoint
            children[1, j] = midpoint
            check(children, depth + 1)

    for start in range(0, len(corners), 1024):
        check(corners[start : start + 1024], 0)


def cell_geometry(vertices: FloatArray) -> FloatArray:
    """Integrate Tet4/Tet10 volume and sample Jacobians independently of Meshers.

    Returns:
        Integrated physical volumes per element.

    Raises:
        RuntimeError: A sampled Jacobian is nonpositive or its lower bound is unresolved.
    """
    corner = np.linalg.det(vertices[:, 1:4] - vertices[:, :1])
    if vertices.shape[1] == _TET4_NODES:
        return corner / 6
    edges = np.array([(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)])
    bary = np.vstack([
        np.full((1, 4), 0.25),
        np.full((4, 4), 1 / 6) + np.eye(4) / 3,
        np.eye(4),
        (np.eye(4)[edges[:, 0]] + np.eye(4)[edges[:, 1]]) / 2,
    ])
    dl = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    corner_range = range(5, 5 + _TET4_NODES)
    determinants = []
    corner_maps = []
    vertices = np.subtract(vertices, vertices[:, :1])
    for index, location in enumerate(bary):
        derivative = np.vstack([
            (4 * location - 1)[:, None] * dl,
            4
            * (
                location[edges[:, 0], None] * dl[edges[:, 1]]
                + location[edges[:, 1], None] * dl[edges[:, 0]]
            ),
        ])
        jacobians = np.einsum("tni,nj->tij", vertices, derivative)
        if index in corner_range:  # The four reference vertices follow the five quadrature points.
            corner_maps.append(jacobians)
        determinant = np.linalg.det(jacobians)
        if not np.isfinite(determinant).all() or np.any(determinant <= 0):
            msg = "Nonpositive sampled curved Jacobian"
            raise RuntimeError(msg)
        determinants.append(determinant)
    _check_curved_jacobians(np.stack(corner_maps, axis=1))
    return np.array([-2 / 15, *([3 / 40] * 4)]) @ np.array(determinants[:5])


def peak_rss_bytes() -> int:
    """Read peak resident memory.

    Returns:
        The peak working set or resident size, in bytes.

    Raises:
        ModuleNotFoundError: The Windows peak working set requires psutil.
    """
    if sys.platform == "win32":
        if psutil is None:
            msg = "Windows benchmark memory sampling requires psutil"
            raise ModuleNotFoundError(msg, name="psutil")
        return psutil.Process().memory_info().peak_wset
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (
        1 if sys.platform == "darwin" else 1024
    )


def mapped_cap_error(mesh: Mesh) -> float:
    """Measure chord and plane errors on the four caps of the cylindrical fixture.

    Returns:
        Largest sampled distance to the corresponding cap surface.

    Raises:
        RuntimeError: The fixture is missing a cap or has unexpected cap labels.
    """
    caps = (_INNER_CYLINDER, _BOTTOM_PLANE, _START_PLANE, _END_PLANE)
    if set(np.unique(mesh.boundary_tags)) != {0, *caps}:
        msg = "Unexpected mapped-intersection cap labels"
        raise RuntimeError(msg)
    errors = []
    for tag in caps:
        tri = mesh.points[mesh.surface[mesh.boundary_tags == tag, :3]]
        samples = np.concatenate([
            tri.reshape(-1, 3),
            ((tri + np.roll(tri, 1, axis=1)) / 2).reshape(-1, 3),
            tri.mean(axis=1),
        ])
        x, y, z = samples.T
        if tag == _INNER_CYLINDER:
            distance = np.hypot(x, y) - 3
        elif tag == _BOTTOM_PLANE:
            distance = z
        else:
            angle = 0 if tag == _START_PLANE else 1 / 3
            distance = y * np.cos(angle) - x * np.sin(angle)
        errors.append(float(np.max(abs(distance))))
    return max(errors)


def surface_measurements(case: str, mesh: Mesh) -> dict[str, float]:
    """Return comparable wall and cap measurements for the installed wheel."""
    if case == "mapped_intersection":
        return {
            "surface_error": max(mesh.diagnostics["surface_errors"].values()),
            "mapped_cap_error": mapped_cap_error(mesh),
        }
    return {"surface_error": mesh.diagnostics["sampled_surface_error"]}


def _volume_options(case: str, options: _FixtureOptions) -> _VolumeOptions:
    volume_options: _VolumeOptions = {**options}
    # Use the keyword supported by each wheel, outside the timed call.
    if "compile_expression" in inspect.signature(meshers.generate).parameters:
        volume_options["compile_expression"] = case != "callback_gyroid"
    else:
        volume_options["compile"] = case != "callback_gyroid"
    if case == "gyroid_tet10":
        volume_options["order"] = 2
    if case.startswith("adaptive_"):
        volume_options.update(
            adapt=True,
            min_length=0.005,
            max_length=0.1,
            adapt_tolerance=0.002,
            geometry_tolerance=0.003,
            gradation=0.6,
        )
    return volume_options


def _measure_generation(case: str) -> tuple[Mesh, float, int]:
    def gyroid(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
        a, b, c = 2 * np.pi * x, 2 * np.pi * y, 2 * np.pi * z
        return np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)

    def spinodoid(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
        # Fixed periodic Fourier fixture; not a stochastic field generator.
        return cast(
            "FloatArray",
            sum(
                a * np.sin(2 * np.pi * (i * x + j * y + k * z) + phase)
                for i, j, k, a, phase in [
                    (1, 1, 0, 0.6, 0.2),
                    (0, 1, 1, 0.5, 0.7),
                    (1, 0, -1, 0.45, 1.1),
                    (2, 1, 0, 0.25, 0.4),
                    (0, 2, -1, 0.2, 0.9),
                    (1, -1, 2, 0.15, 1.4),
                ]
            ),
        )

    def lattice(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
        a, b, c = (np.sin(np.pi * v) ** 2 for v in (x, y, z))

        def union(f: FloatArray, g: FloatArray) -> FloatArray:
            return np.minimum(f, g) - np.maximum(0.08 - np.abs(f - g), 0) ** 2 / 0.32

        return union(union(a + b - 0.24, b + c - 0.24), a + c - 0.24)

    options: _FixtureOptions = {
        "cells": 24,
        "optimize_passes": 2,
        "geometry_tolerance": 0.1,
        "max_tetrahedra": CAP,
    }
    volume_options = _volume_options(case, options)
    start = time.perf_counter()
    if case == "mapped_intersection":

        def mapping(
            x: FloatArray,
            y: FloatArray,
            z: FloatArray,
        ) -> tuple[FloatArray, FloatArray, FloatArray]:
            return (3 + x) * np.cos(z / 3), (3 + x) * np.sin(z / 3), -y

        mesh = meshers.generate_intersection(
            {"right": lambda x, _y, _z: x - 0.501, "back": lambda _x, y, _z: y - 0.751},
            coordinate_map=mapping,
            **options,
        )
    else:
        field = (
            spinodoid
            if "spinodoid" in case
            else lattice
            if "lattice" in case
            else (
                lambda x, y, z: (
                    np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y) + np.cos(2 * np.pi * z)
                )
            )
            if case == "schwarz_p"
            else gyroid
        )
        if case == "graded":

            def field(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:
                return gyroid(x, y, z) / (0.35 + 0.15 * x)

        mesh = meshers.generate(
            field,
            band=(
                (-1, 1)
                if case == "graded"
                else None
                if case
                in {
                    "schwarz_p",
                    "spinodoid",
                    "smooth_lattice",
                    "adaptive_spinodoid",
                    "adaptive_lattice",
                }
                else (-0.5, 0.5)
            ),
            periodic=(False,) * 3 if case == "graded" else (True,) * 3,
            threads=1,
            batch_size=4096,
            **volume_options,
        )
    elapsed = time.perf_counter() - start
    # Capture OS high-water RSS before allocating arrays for independent checks.
    return mesh, elapsed, peak_rss_bytes()


def _sorted_trace(faces: Indices) -> Indices:
    # Attach each Tri6 midpoint to its opposite corner before reordering.
    order = np.argsort(faces[:, :3], axis=1)
    keys = np.take_along_axis(faces[:, :3], order, axis=1)
    if faces.shape[1] == _TET_FACES.shape[1]:
        mids = np.take_along_axis(faces[:, [4, 5, 3]], order, axis=1)
        keys = np.column_stack((keys, mids))
    return keys[np.lexsort(keys.T[::-1])]


def _volume_cap(mesh: Mesh, on_plane: NDArray[np.bool_]) -> Indices:
    # A nondegenerate tet in the box meets its supporting plane in a face
    # precisely when three corners lie there. Only gather those cells.
    corners = on_plane[mesh.tetrahedra[:, :4]]
    cells = np.flatnonzero(np.count_nonzero(corners, axis=1) == _TRI3_NODES)
    opposite = np.argmin(corners[cells], axis=1)
    slots = _TET_FACES[opposite]
    if mesh.tetrahedra.shape[1] == _TET4_NODES:
        slots = slots[:, :3]
    return mesh.tetrahedra[cells[:, None], slots]


def check_periodic_axis(mesh: Mesh, axis: int, transform: FloatArray) -> None:
    """Check the solver's cap connectivity and pairs on a periodic unit-box axis.

    Classify volume faces in parameter space for mapped meshes, then verify
    every supplied pair against the fixture's rigid map in physical space.
    Wall-only contacts may have partners without belonging to a material face.

    Raises:
        RuntimeError: The pair table, tagged boundary or opposite trace is inconsistent.
    """
    coordinates = mesh.points if mesh.parameters is None else mesh.parameters
    pairs = mesh.periodic_pairs[axis]
    if (
        pairs.shape[1:] != (2,)
        or not len(pairs)
        or np.any(pairs < 0)
        or np.any(pairs >= len(mesh.points))
    ):
        msg = f"Invalid periodic pair table on axis {axis}"
        raise RuntimeError(msg)
    source, target = pairs.T
    if (
        len(np.unique(source)) != len(source)
        or len(np.unique(target)) != len(target)
        or not np.all(coordinates[source, axis] == 0)
        or not np.all(coordinates[target, axis] == 1)
        or not np.allclose(
            mesh.points[source] @ transform[:3, :3].T + transform[:3, 3],
            mesh.points[target],
            atol=1e-10,
            rtol=0,
        )
    ):
        msg = f"Invalid periodic correspondence on axis {axis}"
        raise RuntimeError(msg)
    correspondence = np.full(len(mesh.points), -1, dtype=np.int64)
    correspondence[source] = target
    traces = []
    for side in (0, 1):
        cap = _volume_cap(mesh, coordinates[:, axis] == side)
        tagged = mesh.surface[mesh.boundary_tags == 2 * axis + 1 + side]
        volume_trace = _sorted_trace(cap)
        if (
            not len(cap)
            or cap.shape != tagged.shape
            or not np.array_equal(volume_trace, _sorted_trace(tagged))
        ):
            msg = f"Tagged cap differs from volume trace on axis {axis}, side {side}"
            raise RuntimeError(msg)
        if side == 0:
            cap = correspondence[cap]
            if np.any(cap < 0):
                msg = f"Missing periodic cap pair on axis {axis}"
                raise RuntimeError(msg)
            volume_trace = _sorted_trace(cap)
        traces.append(volume_trace)
    if not np.array_equal(traces[0], traces[1]):
        msg = f"Opposite cap traces differ on axis {axis}"
        raise RuntimeError(msg)


def sample(case: str) -> dict[str, float]:
    """Measure a fresh process, then independently validate its mesh.

    Returns:
        Timing, peak resident memory and independently computed mesh metrics.

    Raises:
        RuntimeError: The mesh fails the benchmark's validity or quality checks.
    """
    mesh, elapsed, peak = _measure_generation(case)
    if max(len(mesh.points), len(mesh.tetrahedra)) > CAP:
        msg = "Benchmark node/tetrahedron cap exceeded"
        raise RuntimeError(msg)
    points = mesh.points[mesh.tetrahedra]
    determinants = np.linalg.det(points[:, 1:4] - points[:, :1])
    volumes = cell_geometry(points)
    if not len(points) or not np.all(np.isfinite(mesh.points)) or not np.all(determinants > 0):
        msg = "Nonempty, finite, positively oriented mesh required"
        raise RuntimeError(msg)
    edges = sum(
        np.sum((points[:, i] - points[:, j]) ** 2, axis=1)
        for i in range(4)
        for j in range(i + 1, 4)
    )
    quality = np.sqrt(432 * determinants**2 / edges**3)
    surface = mesh.surface[:, :3]
    boundary_edges = np.sort(
        np.concatenate([surface[:, [0, 1]], surface[:, [1, 2]], surface[:, [2, 0]]]),
        axis=1,
    )
    _, incidence = np.unique(boundary_edges, axis=0, return_counts=True)
    if not np.all(incidence == _TRIANGLES_PER_BOUNDARY_EDGE):
        msg = "Closed manifold boundary required"
        raise RuntimeError(msg)
    if case not in {"graded", "mapped_intersection"}:
        for axis in range(3):
            transform = np.eye(4)
            transform[axis, 3] = 1
            check_periodic_axis(mesh, axis, transform)
    record = {
        "seconds": elapsed,
        "peak_rss_bytes": peak,
        "tetrahedra": len(points),
        "nodes": len(mesh.points),
        "volume": float(volumes.sum()),
        "quality_min": float(quality.min()),
        "quality_p05": float(np.quantile(quality, 0.05)),
        "quality_p10": float(np.quantile(quality, 0.1)),
        "quality_median": float(np.median(quality)),
        "quality_mean": float(quality.mean()),
        "quality_p90": float(np.quantile(quality, 0.9)),
        "quality_max": float(quality.max()),
        **surface_measurements(case, mesh),
    }
    if (
        record["quality_min"] < _MIN_CELL_QUALITY
        or max(record["surface_error"], record.get("mapped_cap_error", 0)) > _MAX_SURFACE_ERROR
    ):
        msg = "Benchmark geometry or quality floor failed"
        raise RuntimeError(msg)
    return record


def runtime_identity() -> dict[str, str]:
    """Return versions and the actual loaded native binary fingerprint."""
    native_module = importlib.import_module("meshers._meshers")
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "meshers": importlib.metadata.version("meshers"),
        "native_sha256": hashlib.sha256(Path(native_module.__file__).read_bytes()).hexdigest(),
    }
