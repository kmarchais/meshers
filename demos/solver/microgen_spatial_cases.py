# Copyright (c) 2026 Kevin Marchais

"""Raw-level recipes for the rotated graded cube and octagonally bored infill."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pyvista as pv
from microgen.shape import GradedInfill, Infill, from_field
from microgen.shape.surface_functions import gyroid
from scipy.spatial.transform import Rotation

if TYPE_CHECKING:
    from collections.abc import Callable

    from microgen import Tpms
    from numpy.typing import NDArray

HALF = 0.5
BOUNDS = (-HALF, HALF) * 3
ROT = Rotation.from_euler("zy", [28, 17], degrees=True).as_matrix()
CENTER = np.array([0.37, -0.21, 0.19])
SHIFT = (0.07, -0.02, 0.05)
SKIN, CORE, TRANSITION, SMOOTHNESS = 1.45, 0.48, 0.42, 0.19
OFFSET = 0.70
BORE_APOTHEM = 0.185
BORE_CENTER = (0.065, -0.04)
VTK_RESOLUTION = 28


def gyroid_value(x: NDArray, y: NDArray, z: NDArray) -> NDArray:
    """Return the unnormalized trigonometric gyroid."""
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def bore_envelope() -> pv.PolyData:
    """Return an outward-oriented cube with a shifted regular octagonal through-hole."""
    angles = np.arange(8) * np.pi / 4
    radial = np.column_stack((np.cos(angles), np.sin(angles)))
    outer = HALF * radial / np.max(abs(radial), axis=1)[:, None]
    inner = BORE_CENTER + radial * (BORE_APOTHEM / np.cos(np.pi / 8))
    plane = np.vstack([outer, inner])
    points = np.vstack([np.column_stack([plane, np.full(16, z)]) for z in (-HALF, HALF)])
    faces = []
    for i in range(8):
        j = (i + 1) % 8
        faces.extend([
            [i, j, j + 16, i + 16],
            [i + 8, i + 24, j + 24, j + 8],
            [i, i + 8, j + 8, j],
            [i + 16, j + 16, j + 24, i + 24],
        ])
    return pv.PolyData.from_regular_faces(points, np.asarray(faces)).triangulate()


def bore_level(x: NDArray, y: NDArray) -> NDArray:
    """Return a positive level inside the through-hole and zero on its eight faces."""
    a = np.arange(8) * np.pi / 4 + np.pi / 8
    values = [
        (x - BORE_CENTER[0]) * cx + (y - BORE_CENTER[1]) * cy
        for cx, cy in zip(np.cos(a), np.sin(a), strict=True)
    ]
    maximum = values[0]
    for value in values[1:]:
        maximum = np.maximum(maximum, value)
    return BORE_APOTHEM - maximum


@dataclass
class Specimen:
    """A pinned upstream object, its pointwise recipe and a rigid local meshing frame."""

    name: str
    upstream: Tpms
    raw: Callable
    part: str
    depth_max: float = 0.0

    def world(self, points: NDArray) -> NDArray:
        """Return physical coordinates for cube-local points."""
        return points @ ROT.T + CENTER if self.name == "distance-infill" else points

    def local(self, points: NDArray) -> NDArray:
        """Return cube-local coordinates for physical points."""
        return (points - CENTER) @ ROT if self.name == "distance-infill" else points

    def offset(self, x: NDArray, y: NDArray, z: NDArray) -> NDArray | float:
        """Return the raw level width inside the cube in the local frame."""
        if self.name != "distance-infill":
            return OFFSET
        depth = HALF - np.maximum(np.abs(x), np.maximum(np.abs(y), np.abs(z)))
        ratio = np.clip(depth / self.depth_max, 0, 1)
        return SKIN + (CORE - SKIN) * (1 + np.tanh((ratio - TRANSITION) / SMOOTHNESS)) / 2

    def level(self, x: NDArray, y: NDArray, z: NDArray) -> NDArray:
        """Return the signed raw level that Meshers must represent."""
        if self.name == "bored-infill":
            xw, yw, zw = (
                (q + s) * 2 * np.pi / c
                for q, s, c in zip(
                    (x, y, z),
                    self.upstream.phase_shift,
                    self.upstream.cell_size,
                    strict=True,
                )
            )
            return np.maximum(OFFSET / 2 - self.raw(xw, yw, zw), bore_level(x, y))
        physical = [
            row[0] * x + row[1] * y + row[2] * z + c for row, c in zip(ROT, CENTER, strict=True)
        ]
        angles = [
            (q + shift) * (2 * np.pi / size)
            for q, shift, size in zip(
                physical,
                self.upstream.phase_shift,
                self.upstream.cell_size,
                strict=True,
            )
        ]
        value = self.raw(*angles)
        half = self.offset(x, y, z) / 2
        return np.abs(value) - half

    def grid(self, options: dict) -> pv.UnstructuredGrid:
        """Return a physical Tet4 grid using only public Microgen volume entry points."""
        shape = from_field(self.level, bounds=BOUNDS)
        grid = shape.generate_volume_mesh(backend="meshers", periodic=(False,) * 3, **options)
        grid.points = self.world(grid.points)
        return grid


def specimen(name: str, phase_shift: tuple[float, ...] = SHIFT) -> Specimen:
    """Return the rotated graded cube or the bored cube constructed with upstream objects.

    Raises:
        ValueError: Unknown specimen.
    """
    if name == "distance-infill":
        obj = pv.Cube(points_dtype="float64").triangulate()
        obj.points = obj.points @ ROT.T + CENTER
        origin, end = np.asarray(obj.bounds).reshape(3, 2).T
        world = np.stack(np.meshgrid(*np.linspace(origin, end, 20).T, indexing="ij"), axis=-1)
        depth_max = float(HALF - np.max(np.abs((world - CENTER) @ ROT), axis=-1).min())
        obj_tpms = GradedInfill(
            obj=obj,
            surface_function=gyroid,
            repeat_cell=(2, 2, 2),
            phase_shift=phase_shift,
            resolution=VTK_RESOLUTION,
            offset_skin=SKIN,
            offset_core=CORE,
            transition=TRANSITION,
            smoothness=SMOOTHNESS,
        )
        return Specimen(name, obj_tpms, gyroid_value, "sheet", depth_max)
    if name == "bored-infill":
        shape = Infill(
            obj=bore_envelope(),
            surface_function=gyroid,
            offset=OFFSET,
            repeat_cell=(1, 1, 1),
            resolution=VTK_RESOLUTION,
            phase_shift=phase_shift,
        )
        return Specimen(name, shape, gyroid_value, "upper skeletal")
    msg = f"Unknown infill specimen {name}"
    raise ValueError(msg)
