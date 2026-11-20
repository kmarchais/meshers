# Copyright (c) 2026 Kevin Marchais

"""Experimental surface-only meshing for printing and previewing implicit solids."""

import time
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Unpack

import numpy as np

from . import _meshers
from ._options import MeshingOptions, check_keywords
from .expression import TraceError, compile_field

_generate_surface = getattr(_meshers, "generate_surface", None)
SURFACE_MESHING_AVAILABLE = _generate_surface is not None


@dataclass(frozen=True)
class SurfaceMesh:
    """Triangle arrays with closed-manifold checks and generation diagnostics."""

    points: np.ndarray
    triangles: np.ndarray
    boundary_tags: np.ndarray
    diagnostics: dict


class _SurfaceOptions(MeshingOptions, total=False):
    """Options for triangle extraction and adaptive remeshing."""

    gradient: Callable | None
    quality: bool
    adaptive: bool
    iterations: int
    balanced: bool
    fast: bool


def _remeshing(
    options: _SurfaceOptions,
    tolerance: float,
    extent: float,
) -> tuple[float, float, float, float, int, bool, bool] | None:
    if not options.get("adaptive"):
        return None
    if tolerance <= 0:
        msg = "adaptive remeshing needs a positive geometry_tolerance"
        raise ValueError(msg)
    minimum, maximum = options.get("min_length"), options.get("max_length")
    return (
        tolerance,
        tolerance / 2 if minimum is None else float(minimum),
        extent / 4 if maximum is None else float(maximum),
        float(options.get("gradation", 0.3)),
        int(options.get("iterations", 8)),
        bool(options.get("balanced", False)),
        bool(options.get("fast", False)),
    )


def generate_surface(
    field: Callable | _meshers.CompiledField,
    *,
    bounds: Sequence[float] = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
    cells: int | Sequence[int] = 24,
    band: Sequence[float] | None = None,
    geometry_tolerance: float | None = None,
    **options: Unpack[_SurfaceOptions],
) -> SurfaceMesh:
    """Generate a closed boundary without tetrahedra.

    ``quality=True`` smooths and improves triangles while retaining cap features.
    This experimental API checks manifold topology and finite positive triangle
    shapes. It does not certify intersections or geometric approximation errors.
    Sampling can miss features below the background spacing.
    ``adaptive=True`` remeshes towards curvature-adaptive edge lengths so that the
    sampled wall distance meets ``geometry_tolerance`` with as few triangles as
    possible; ``cells`` then only sets the extraction grid that resolves topology.
    ``balanced=True`` also moves wall vertices outward by half the local chord sag,
    so triangles deviate to both sides of the surface: about half the triangles at
    the same tolerance and almost no volume bias, but vertices are then off the
    surface. Use it for printing and display, not for curved quadratic FEA meshes.
    ``fast=True`` (with ``adaptive=True``) only coarsens the extraction and splits
    where the tolerance is exceeded, without flips or smoothing: faster, a few more
    triangles, and lower but valid triangle shapes. ``iterations`` then counts
    coarsening passes; two are enough.
    Set ``geometry_tolerance=0`` for a quick preview without distance sampling;
    topology and finite triangle checks still run. Diagnostics mark that omission.

    Other Parameters:
        gradient (Callable | None): Optional analytic field gradient. Defaults to None.
        quality (bool): Smooth and improve triangles. Defaults to False.
        batch_size (int): Maximum callback batch size. Defaults to 4096.
        cancel (CancellationToken | None): Cooperative cancellation token. Defaults to None.
        adaptive (bool): Remesh toward the sampled distance tolerance. Defaults to False.
        min_length (float | None): Minimum adaptive edge length. Defaults to half the tolerance.
        max_length (float | None): Maximum adaptive edge length. Defaults to a quarter of the
            shortest box side.
        gradation (float): Maximum target-length growth per unit distance. Defaults to 0.3.
        iterations (int): Adaptive remeshing rounds. Defaults to 8.
        balanced (bool): Allow wall vertices to move outward to balance chord error.
            Defaults to False.
        fast (bool): Use the coarsen-and-split adaptive path. Defaults to False.

    Returns:
        Triangle mesh and measured generation diagnostics.

    Raises:
        RuntimeError: Surface support was omitted from the native build.
        ValueError: Bounds or geometry tolerance are invalid.
        MeshingError: Generation fails or sampled error exceeds tolerance.
    """
    check_keywords(options, _SurfaceOptions)
    if _generate_surface is None:
        msg = "Surface-only meshing is unavailable in this meshers build"
        raise RuntimeError(msg)
    domain = np.asarray(bounds, dtype=float)
    if domain.shape != (6,):
        msg = "bounds must contain six coordinates"
        raise ValueError(msg)
    started = time.perf_counter()
    gradient = options.get("gradient")
    if callable(field) and not isinstance(field, _meshers.CompiledField):
        try:
            field = compile_field(field, gradient=gradient)
            gradient = None
        except TraceError as error:
            warnings.warn(f"Using NumPy callbacks: {error}", RuntimeWarning, stacklevel=2)
    preparation = time.perf_counter() - started
    extent = float(np.min(domain[1::2] - domain[::2]))
    tolerance = 0.01 * extent if geometry_tolerance is None else float(geometry_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0:
        msg = "geometry_tolerance must be finite and nonnegative"
        raise ValueError(msg)
    quality_passes = (4, 4, 8) if options.get("quality") else (0, 0, 0)
    result = _generate_surface(
        field,
        gradient,
        [domain[::2].tolist(), domain[1::2].tolist()],
        cells,
        list(band) if band is not None else [-np.inf, 0.0],
        *quality_passes,
        options.get("batch_size", 4096),
        options.get("cancel"),
        tolerance > 0,
        _remeshing(options, tolerance, extent),
    )
    labels = np.asarray(result["labels"])
    if tolerance > 0 and result["maximum_sampled_surface_error"] > tolerance:
        msg = (
            f"Sampled surface error {result['maximum_sampled_surface_error']} "
            f"exceeds tolerance {tolerance}; increase resolution"
        )
        raise _meshers.MeshingError(msg)
    field_labels = 2
    tags = np.where(labels < field_labels, 0, labels - 1).astype(np.uint8)
    return SurfaceMesh(
        np.asarray(result["points"]),
        np.asarray(result["triangles"], dtype=np.int64),
        tags,
        {
            key: result[key]
            for key in (
                "extraction_seconds",
                "compute_seconds",
                "minimum_triangle_shape",
                "volume",
                "maximum_sampled_surface_error",
                "geometry_checked",
            )
        }
        | {"preparation_seconds": preparation, "experimental": True}
        | ({"adaptive": result["adaptive"]} if "adaptive" in result else {}),
    )
