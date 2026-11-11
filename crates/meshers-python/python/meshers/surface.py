"""Experimental surface-only meshing for printing and previewing implicit solids."""

import time
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from . import _meshers
from .expression import TraceError, compile_field


@dataclass(frozen=True)
class SurfaceMesh:
    """Triangle arrays with closed-manifold checks and generation diagnostics."""

    points: np.ndarray
    triangles: np.ndarray
    boundary_tags: np.ndarray
    diagnostics: dict


def generate_surface(  # ruff: ignore[too-many-arguments] - match the volume API options.
    field: Callable | _meshers.CompiledField,
    *,
    bounds: Sequence[float] = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
    cells: int | Sequence[int] = 24,
    band: Sequence[float] | None = None,
    gradient: Callable | None = None,
    quality: bool = False,
    geometry_tolerance: float | None = None,
    batch_size: int = 4096,
    cancel: _meshers.CancellationToken | None = None,
    adaptive: bool = False,
    min_length: float | None = None,
    max_length: float | None = None,
    gradation: float = 0.3,
    iterations: int = 8,
    balanced: bool = False,
    fast: bool = False,
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

    Returns:
        Triangle mesh and measured generation diagnostics.

    Raises:
        RuntimeError: Surface support was omitted from the native build.
        ValueError: Bounds or geometry tolerance are invalid.
        MeshingError: Generation fails or sampled error exceeds tolerance.
    """
    if not hasattr(_meshers, "generate_surface"):
        raise RuntimeError("Surface-only meshing is unavailable in this meshers build")
    domain = np.asarray(bounds, dtype=float)
    if domain.shape != (6,):
        raise ValueError("bounds must contain six coordinates")
    started = time.perf_counter()
    if callable(field) and not isinstance(field, _meshers.CompiledField):
        try:
            field = compile_field(field, gradient=gradient)
            gradient = None
        except TraceError as error:
            warnings.warn(f"Using NumPy callbacks: {error}", RuntimeWarning, stacklevel=2)
    preparation = time.perf_counter() - started
    tolerance = (
        0.01 * float(np.min(domain[1::2] - domain[::2]))
        if geometry_tolerance is None
        else float(geometry_tolerance)
    )
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("geometry_tolerance must be finite and nonnegative")
    options = (4, 4, 8) if quality else (0, 0, 0)
    extent = float(np.min(domain[1::2] - domain[::2]))
    remeshing = None
    if adaptive:
        if tolerance <= 0:
            raise ValueError("adaptive remeshing needs a positive geometry_tolerance")
        remeshing = (
            tolerance,
            tolerance / 2 if min_length is None else float(min_length),
            extent / 4 if max_length is None else float(max_length),
            float(gradation),
            int(iterations),
            bool(balanced),
            bool(fast),
        )
    result = _meshers.generate_surface(
        field,
        gradient,
        [domain[::2].tolist(), domain[1::2].tolist()],
        cells,
        list(band) if band is not None else [-np.inf, 0.0],
        *options,
        batch_size,
        cancel,
        tolerance > 0,
        remeshing,
    )
    labels = np.asarray(result["labels"])
    if tolerance > 0 and result["maximum_sampled_surface_error"] > tolerance:
        raise _meshers.MeshingError(
            f"Sampled surface error {result['maximum_sampled_surface_error']} "
            f"exceeds tolerance {tolerance}; increase resolution"
        )
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
