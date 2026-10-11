"""Fixed physical geometry definitions shared by all backends."""

from dataclasses import asdict, dataclass
from functools import lru_cache
from itertools import product

import numpy as np


@lru_cache(maxsize=8)
def spinodoid_source(scale):
    """Freeze spectral resolution independently of meshing resolution."""
    from microgen import Spinodoid

    return Spinodoid(
        k0=4 * np.pi / scale,
        bandwidth=3 / scale,
        offset=0,
        cell_size=scale,
        resolution=16,
        seed=20261010,
        mode_amplitude_threshold=0.05,
    )


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

    @property
    def implicit(self):
        return self.function in {"bcc_lattice", "spinodoid"}

    @property
    def band(self):
        # tanh keeps the implicit solid above the inactive lower isovalue.
        return (-1, 0) if self.implicit else (-1, 1)

    def material(self, x, y, z):
        value = self.field(x, y, z)
        return (value >= self.band[0]) & (value <= self.band[1])

    def record(self):
        return asdict(self)

    def width(self, x, y, z):
        return self.offset * (1 + self.grading * x / self.scale)

    def field(self, x, y, z):
        if self.function == "spinodoid":
            # Microgen's spectral cell is [0,L]^3, our specimen is centered.
            source = spinodoid_source(self.scale)
            arrays = np.broadcast_arrays(x, y, z)
            # Bound memory during independent Sobol and boundary sampling.
            flat = [a.ravel() + self.scale / 2 for a in arrays]
            values = np.empty(flat[0].size)
            for start in range(0, len(values), 4096):
                sl = slice(start, start + 4096)
                values[sl] = source.evaluate(*(a[sl] for a in flat))
            return np.tanh(values.reshape(arrays[0].shape))
        if self.function == "bcc_lattice":
            arrays = np.broadcast_arrays(x / self.scale, y / self.scale, z / self.scale)
            p = np.stack(arrays, axis=-1)
            distance = np.full(arrays[0].shape, np.inf)
            for corner in product((-0.5, 0.5), repeat=3):
                end = np.asarray(corner)
                t = np.clip(np.sum(p * end, axis=-1) / np.dot(end, end), 0, 1)
                distance = np.minimum(distance, np.linalg.norm(p - t[..., None] * end, axis=-1))
            return np.tanh(distance / self.offset - 1)
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

        if self.implicit:
            from microgen.shape import from_field

            return from_field(self.field, bounds=(-self.scale / 2, self.scale / 2) * 3)

        raw = (
            (lambda x, y, z: x * 0)
            if self.function == "solid"
            else getattr(surface_functions, self.function)
        )
        offset = self.width if self.grading else self.offset
        return Tpms(raw, offset=offset, cell_size=self.scale, resolution=cells + 1)

    def volume_grid(self, shape, cells):
        if self.implicit:
            return shape.generate_volume_mesh(resolution=cells + 1)
        return shape.grid_sheet

    def verify_field(self, shape, cells):
        if self.implicit:
            axis = np.linspace(-self.scale / 2, self.scale / 2, cells + 1)
            sample = np.array(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1)
            expected, actual = self.field(*sample), shape.evaluate(*sample)
        else:
            sample = np.asarray(shape.grid.points).T
            expected = self.field(*sample) * 0.5 * self.width(*sample)
            actual = np.asarray(shape.grid.point_data["surface"])
            offset = np.broadcast_to(np.asarray(shape.offset).ravel(order="F"), (sample.shape[1],))
            if not np.allclose(offset, self.width(*sample), rtol=1e-12, atol=1e-12):
                raise ValueError("Actual Microgen offset differs from the declared graded field")
        if not np.allclose(expected, actual, rtol=1e-12, atol=1e-12):
            raise ValueError("Case formula and Microgen's sampled surface field disagree")
        return sample.shape[1]


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
    Case(
        "implicit-lattice",
        "Implicit BCC lattice",
        "bcc_lattice",
        0.12,
        "Eight capsule struts join the cell center to its corners. Benchmark-defined implicit union, sampled through Microgen's public from_field API.",
        "min distance to segments [0, (±L/2, ±L/2, ±L/2)] ≤ 0.12 L",
    ),
    Case(
        "spinodoid",
        "Seeded spinodoid",
        "spinodoid",
        0,
        "Microgen Spinodoid field, seed 20261010; fixed 16-cell spectral construction, k0=4π/L, bandwidth=3/L, retained-mode threshold 0.05. Only meshing resolution varies.",
        "Σ Aᵢ cos(2π mᵢ·(x/L + 1/2) + φᵢ) ≥ 0; shared signed field = tanh(−Σ)",
    ),
]

PROFILES = {
    "meshers-fast": {"optimize_passes": 0},
    "meshers-balanced": {"optimize_passes": 4},
    "meshers-quality": {"optimize_passes": 12},
}
DEFAULT_WORKFLOWS = ["vtk", "mmg", "mmg-microgen", *PROFILES]
WORKFLOWS = [
    "vtk",
    "mmg",
    "mmg-microgen",
    "vtk-tetra",
    "mmg-tetra",
    "vtk-tetra-clean",
    "mmg-tetra-clean",
    *PROFILES,
]
