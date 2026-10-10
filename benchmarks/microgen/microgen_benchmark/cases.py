"""Fixed physical geometry definitions shared by all backends."""

from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class Case:
    id: str
    name: str
    function: str
    offset: float
    description: str
    formula: str
    scale: float = 1.0
    grading: float = 0.0

    def record(self):
        return asdict(self)

    def width(self, x, y, z):
        return self.offset * (1 + self.grading * x / self.scale)

    def field(self, x, y, z):
        a, b, c = 2 * np.pi * x / self.scale, 2 * np.pi * y / self.scale, 2 * np.pi * z / self.scale
        if self.function == "gyroid":
            raw = np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)
        elif self.function == "schwarz_p":
            raw = np.cos(a) + np.cos(b) + np.cos(c)
        elif self.function == "schwarz_d":
            raw = (
                np.sin(a) * np.sin(b) * np.sin(c)
                + np.sin(a) * np.cos(b) * np.cos(c)
                + np.cos(a) * np.sin(b) * np.cos(c)
                + np.cos(a) * np.cos(b) * np.sin(c)
            )
        elif self.function == "neovius":
            raw = 3 * np.cos(a) + np.cos(b) + np.cos(c) + 4 * np.cos(a) * np.cos(b) * np.cos(c)
        elif self.function == "solid":
            raw = x * 0
        else:
            raise ValueError(f"Unknown field {self.function}")
        return raw / (0.5 * self.width(x, y, z))

    def shape(self, cells):
        from microgen import Tpms
        from microgen.shape import surface_functions

        raw = (
            (lambda x, y, z: x * 0)
            if self.function == "solid"
            else getattr(surface_functions, self.function)
        )
        offset = self.width if self.grading else self.offset
        return Tpms(raw, offset=offset, cell_size=self.scale, resolution=cells + 1)


CASES = [
    Case(
        "solid-cell",
        "Solid cell",
        "solid",
        1,
        "Full cube; affine patch and scaling control.",
        "f = 0",
    ),
    Case(
        "gyroid",
        "Gyroid",
        "gyroid",
        0.5,
        "Uniform sheet, one unit cell.",
        "sin X cos Y + sin Y cos Z + sin Z cos X",
    ),
    Case(
        "schwarz-p",
        "Schwarz P",
        "schwarz_p",
        0.6,
        "Uniform sheet with narrow junctions.",
        "cos X + cos Y + cos Z",
    ),
    Case(
        "diamond",
        "Schwarz D / Diamond",
        "schwarz_d",
        0.8,
        "Sheet with tetrahedral junction structure.",
        "sin X sin Y sin Z + sin X cos Y cos Z + cos X sin Y cos Z + cos X cos Y sin Z",
    ),
    Case(
        "neovius",
        "Microgen Neovius variant",
        "neovius",
        1,
        "Released Microgen formula; unequal cosine coefficients, not the symmetric textbook form.",
        "3 cos X + cos Y + cos Z + 4 cos X cos Y cos Z",
    ),
    Case(
        "graded-gyroid",
        "Graded Gyroid",
        "gyroid",
        0.5,
        "Thickness varies by 40% from left to right.",
        "sin X cos Y + sin Y cos Z + sin Z cos X",
        grading=0.4,
    ),
    Case(
        "thin-gyroid",
        "Thin Gyroid",
        "gyroid",
        0.25,
        "Thin-wall resolution stress test.",
        "sin X cos Y + sin Y cos Z + sin Z cos X",
    ),
    Case(
        "scaled-gyroid",
        "Gyroid at scale 0.001",
        "gyroid",
        0.5,
        "Same dimensionless shape at millimetre scale.",
        "sin X cos Y + sin Y cos Z + sin Z cos X",
        scale=0.001,
    ),
]

PROFILES = {
    "meshers-fast": {"optimize_passes": 0},
    "meshers-balanced": {"optimize_passes": 4},
    "meshers-quality": {"optimize_passes": 12},
}
WORKFLOWS = ["vtk", "mmg", *PROFILES]
