# Copyright (c) 2026 Kevin Marchais

"""Bound the distance from one sample to the seeded spinodoid wall inside its unit cell."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from implicit_coverage import load_reference
from numpy.typing import NDArray
from reporting import write_json_line
from scipy.optimize import brentq
from scipy.stats import qmc

if TYPE_CHECKING:
    from microgen import Spinodoid

CELL_HALF_WIDTH = 0.5
SPECTRUM = Path(__file__).with_name("implicit-coverage") / "spinodoid-spectrum.csv"


Coordinates = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class Spectrum:
    """Fourier coefficients of the pinned periodic field."""

    modes: NDArray[np.int64]
    cell_size: Coordinates | float
    amplitudes: Coordinates
    phases: Coordinates
    threshold: float


def reference_spectrum(shape: Spinodoid) -> Spectrum:
    """Read the frozen unit-cell spectrum and check it against the public signed field.

    Returns:
        The coefficients used to bound derivatives of the historical specimen.

    Raises:
        ValueError: The archived coefficients do not describe the reference field.
    """
    rows = np.loadtxt(SPECTRUM, delimiter=",", skiprows=1)
    modes = rows[:, :3].astype(np.int64)
    amplitudes, phases = rows[:, 3], rows[:, 4]
    # The 16-cell FFT used by this specimen has no frequency above 8 per axis.
    # Sampling on 17 points per axis distinguishes all signed frequencies.
    n = shape.resolution + 1
    if np.any(2 * np.abs(modes) >= n):
        msg = "Spectrum exceeds the reference sampling bandwidth"
        raise ValueError(msg)
    axis = np.arange(n) / n
    points = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    values = -np.cos(2 * np.pi * (points @ modes.T) + phases) @ amplitudes
    if not np.allclose(values, shape.evaluate(*points.T), rtol=0, atol=1e-13):
        msg = "Archived spectrum differs from the public Microgen field"
        raise ValueError(msg)
    return Spectrum(modes=modes, amplitudes=amplitudes, phases=phases, cell_size=1.0, threshold=0.0)


def distance_bounds(
    source: Spectrum,
    point: Coordinates,
    free: NDArray[np.bool_],
) -> tuple[float | None, float | None]:
    """Return a Taylor lower bound and the shortest root found along sampled rays.

    Fixed coordinates restrict both bounds to the specified cap or edge. A missing upper
    bound means the ray search found no crossing within its search radius.
    """
    wave = np.asarray(2 * np.pi * source.modes / source.cell_size, dtype=np.float64)
    amplitudes = source.amplitudes
    angles = wave @ (point + CELL_HALF_WIDTH) + source.phases
    residual = source.threshold - np.cos(angles) @ amplitudes
    if residual == 0:
        return 0.0, 0.0
    if not free.any():
        return None, None

    gradient = (amplitudes * np.sin(angles)) @ wave
    lower = taylor_radius(wave[:, free], amplitudes, angles, gradient[free], residual)
    directions = qmc.Sobol(3, seed=71).random_base2(8) - 0.5
    directions = np.vstack([directions, np.eye(3), -np.eye(3), gradient, -gradient])
    directions[:, ~free] = 0
    upper = ray_distance(source, point, wave, angles, directions)
    return float(lower), upper


def taylor_radius(
    k: Coordinates,
    amplitudes: Coordinates,
    angles: Coordinates,
    gradient: Coordinates,
    residual: float,
) -> float:
    """Return a radius within which the Fourier sum cannot reach zero."""
    hessian = (k.T * (amplitudes * np.cos(angles))) @ k
    linear = np.linalg.norm(gradient)
    quadratic = np.linalg.norm(hessian, ord=2) / 2
    cubic = abs(amplitudes) @ (np.linalg.norm(k, axis=1) ** 3) / 6

    def bound(d: float) -> float:
        return ((cubic * d + quadratic) * d + linear) * d - abs(residual)

    diameter = np.sqrt(len(gradient))
    return brentq(bound, 0, diameter) if bound(diameter) >= 0 else diameter


def ray_distance(
    source: Spectrum,
    point: Coordinates,
    wave: Coordinates,
    angles: Coordinates,
    directions: Coordinates,
) -> float | None:
    """Return the shortest zero found on an in-box ray, or None if none crossed."""
    amplitudes = source.amplitudes
    residual = source.threshold - np.cos(angles) @ amplitudes
    lengths = np.linalg.norm(directions, axis=1)
    upper = None
    for direction in directions[lengths > 0] / lengths[lengths > 0, None]:
        moving = direction != 0
        end = np.copysign(CELL_HALF_WIDTH, direction[moving])
        limit = np.min((end - point[moving]) / direction[moving])
        slope = wave @ direction
        for d in np.geomspace(0.00001, 0.25, 30):
            if d > limit:
                break
            if residual * (source.threshold - np.cos(angles + d * slope) @ amplitudes) <= 0:
                root = brentq(
                    lambda t, s=slope: source.threshold - np.cos(angles + t * s) @ amplitudes,
                    0,
                    d,
                )
                upper = root if upper is None else min(upper, root)
                break
    return upper


def main() -> None:
    """Report distance bounds for a point reported by the native mesher."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--point", required=True, nargs=3, type=float)
    parser.add_argument("--lock-caps", action="store_true", help="Restrict motion to incident caps")
    args = parser.parse_args()
    point = np.array(args.point)
    if not np.all(np.isfinite(point) & (abs(point) <= CELL_HALF_WIDTH)):
        parser.error("Supply finite coordinates in the centered unit cell")
    if version("microgen") != "2.0.0b1":
        parser.error("Use the frozen spinodoid from microgen==2.0.0b1")
    cases = load_reference(args.reference).cases
    source = reference_spectrum(cases.spinodoid_source(1))
    free = abs(point) < CELL_HALF_WIDTH if args.lock_caps else np.ones(3, dtype=bool)
    lower, upper = distance_bounds(source, point, free)
    write_json_line({
        "point": args.point,
        "free_axes": free.tolist(),
        "lower": lower,
        "upper": upper,
    })


if __name__ == "__main__":
    main()
