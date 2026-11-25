# Copyright (c) 2026 Kevin Marchais

"""Upstream Microgen objects in a common, centred unit box."""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from functools import partial
from typing import TYPE_CHECKING

import microgen as mg
import numpy as np
from microgen.shape import implicit_ops as ops
from microgen.shape import surface_functions as sf
from scipy.spatial.transform import Rotation

if TYPE_CHECKING:
    from periodic_bench import Coordinates, Field

UNIT = (-0.5, 0.5) * 3
SHAPES = (
    "sphere",
    "ellipsoid",
    "capsule",
    "cylinder",
    "box",
    "polyhedron",
    "two-balls",
    "lens",
    "difference",
    "complement",
    "spherical-shell",
    "smooth-union",
    "smooth-intersection",
    "smooth-difference",
    "blend",
    "batch-union",
    "repeat",
    "cubic-lattice",
    "octet-truss",
    "spinodoid",
)
SURFACES = (
    "gyroid",
    "schwarz_p",
    "schwarz_d",
    "neovius",
    "schoen_iwp",
    "schoen_frd",
    "fischer_koch_s",
    "pmy",
    "honeycomb",
    "lidinoid",
    "split_p",
    "honeycomb_gyroid",
    "honeycomb_schwarz_p",
    "honeycomb_schwarz_d",
    "honeycomb_schoen_iwp",
    "honeycomb_lidinoid",
)
PARTS = ("sheet", "upper", "lower")
CASES = (
    *SHAPES,
    *(f"{fn}-{part}" for fn, part in itertools.product(SURFACES, PARTS)),
    "graded-iwp",
)


@dataclass
class Case:
    """A public meshing call and an independent occupancy query in the same frame."""

    shape: mg.Shape
    options: dict
    field: Field
    periodic: bool = False
    expected_volume: float | None = None
    expected_components: int | None = None
    hull_vertices: Coordinates | None = None
    two_sphere_operation: str | None = None
    component_targets: Coordinates | None = None


def _tpms(name: str) -> Case:
    graded = name == "graded-iwp"
    family, part = ("schoen_iwp", "sheet") if graded else name.rsplit("-", 1)
    surface = getattr(sf, family)
    offset = (lambda _x, _y, z: 1.1 + 0.2 * z) if graded else (1.1 if part == "sheet" else 0.0)
    shape = mg.Tpms(surface_function=surface, offset=offset, resolution=12)

    def field(points: Coordinates) -> Coordinates:
        raw = surface(*(2 * np.pi * points.T))
        if part == "sheet":
            return np.abs(raw) - (0.55 + 0.1 * points[:, 2] if graded else 0.55)
        return raw if part == "lower" else -raw

    options = {"type_part": "sheet" if part == "sheet" else f"{part} skeletal"}
    return Case(shape, options, field, periodic=not graded)


def _primitive(name: str) -> tuple[mg.Shape, float | None]:
    rotation = Rotation.from_euler("zy", [21, 32], degrees=True)
    pose = {"orientation": rotation, "center": (0.02, 0.01, -0.015)}
    ball = 4 * np.pi / 3
    cases = {
        "sphere": (partial(mg.Sphere, radius=0.30, **pose), ball * 0.30**3),
        "ellipsoid": (
            partial(mg.Ellipsoid, radii=(0.35, 0.24, 0.13), **pose),
            ball * 0.35 * 0.24 * 0.13,
        ),
        "capsule": (
            partial(mg.Capsule, height=0.40, radius=0.12, **pose),
            np.pi * 0.12**2 * 0.40 + ball * 0.12**3,
        ),
        "cylinder": (
            partial(mg.Cylinder, height=0.55, radius=0.19, **pose),
            np.pi * 0.19**2 * 0.55,
        ),
        "box": (partial(mg.Box, dim=(0.45, 0.32, 0.23), **pose), 0.45 * 0.32 * 0.23),
        "polyhedron": (lambda: mg.Polyhedron().scale(0.23), 8 / 3 * 0.23**3),
    }
    factory, volume = cases[name]
    return factory(), volume


def _composition(name: str) -> tuple[mg.Shape, float | None, int | None]:
    ball = 4 * np.pi / 3
    left = mg.Sphere(radius=0.23, center=(-0.12, 0, 0))
    right = mg.Sphere(radius=0.23, center=(0.12, 0, 0))
    isolated = mg.Sphere(radius=0.12, center=(0.24, 0.03, 0))
    lens = np.pi * (4 * 0.23 + 0.24) * (2 * 0.23 - 0.24) ** 2 / 12
    cases = {
        "two-balls": (lambda: left | isolated, ball * (0.23**3 + 0.12**3), 2),
        "lens": (lambda: left & right, lens, 1),
        "difference": (lambda: left - right, ball * 0.23**3 - lens, 1),
        "complement": (lambda: ~left, 1 - ball * 0.23**3, 1),
        "spherical-shell": (
            lambda: ops.shell(mg.Sphere(radius=0.29), 0.12),
            ball * (0.35**3 - 0.23**3),
            1,
        ),
        "smooth-union": (lambda: left.smooth_union(right, 0.11), None, 1),
        "smooth-intersection": (lambda: left.smooth_intersection(right, 0.11), None, 1),
        "smooth-difference": (lambda: left.smooth_difference(right, 0.11), None, 1),
        "blend": (lambda: ops.blend(left, right, 0.37), None, 1),
        "batch-union": (lambda: ops.batch_smooth_union([left, right, isolated]), None, 1),
        "repeat": (
            lambda: ops.repeat(mg.Sphere(radius=0.135), (0.5, 0.5, 0.5)),
            8 * ball * 0.135**3,
            27,
        ),
    }
    factory, volume, components = cases[name]
    return factory(), volume, components


def _hull_vertices(name: str) -> Coordinates | None:
    if name == "polyhedron":
        return 0.23 * np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
    if name == "box":
        half = 0.5 * np.array([0.45, 0.32, 0.23])
        rotation = Rotation.from_euler("zy", [21, 32], degrees=True)
        corners = half * np.array(list(itertools.product((-1, 1), repeat=3)))
        return rotation.apply(corners) + np.array([0.02, 0.01, -0.015])
    return None


def _component_targets(name: str) -> Coordinates | None:
    ball = 4 * np.pi / 3
    if name == "two-balls":
        return ball * np.array([0.23**3, 0.12**3])
    if name == "repeat":
        return ball * 0.135**3 * np.repeat([1, 0.5, 0.25, 0.125], [1, 6, 12, 8])
    return None


def specimen(name: str) -> Case:
    """Construct one specimen using upstream methods.

    Returns:
        The upstream object, audit field and any analytic expectations.
    """
    if name not in SHAPES:
        return _tpms(name)
    if name == "spinodoid":
        shape = mg.Spinodoid(
            k0=4 * np.pi,
            bandwidth=3,
            offset=0,
            resolution=16,
            anisotropy=(1, 1, 2),
            mode_amplitude_threshold=0.05,
            seed=20261124,
            center=(-0.5,) * 3,
        )
        return Case(shape, {}, lambda p: shape.evaluate(*p.T), periodic=True)
    if name in {"cubic-lattice", "octet-truss"}:
        lattice = mg.Cubic if name == "cubic-lattice" else mg.OctetTruss
        shape = lattice(strut_radius=0.085, strut_joints=True).generate_implicit()
        volume, components = None, 1
    elif name in SHAPES[:6]:
        shape, volume = _primitive(name)
        components = 1
    else:
        shape, volume, components = _composition(name)
    return Case(
        shape,
        {"bounds": UNIT},
        lambda p: shape.evaluate(*p.T),
        name == "repeat",
        volume,
        components,
        _hull_vertices(name),
        name if name in {"lens", "difference"} else None,
        component_targets=_component_targets(name),
    )
