"""Use Microgen's TPMS callables and a frozen Microgen spectral spinodoid.

The spinodoid spectrum is generated once from Microgen's source, independently
of meshing resolution. Its coefficients are copied into each prepared run.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
from microgen_benchmark.cases import Case

TPMS = ["gyroid", "schwarz_p", "schwarz_d", "neovius", "schoen_iwp", "schoen_frd",
        "fischer_koch_s", "pmy", "lidinoid", "split_p"]


@lru_cache(maxsize=1)
def spectral_data():
    with np.load(Path(__file__).parent / "spinodoid-spectrum.npz") as archive:
        return {name: archive[name].copy() for name in archive.files}


class MicrogenCase(Case):
    @property
    def implicit(self):
        return self.function == "spinodoid_spectral"

    def field(self, x, y, z):
        if self.implicit:
            spectrum = spectral_data()
            arrays = np.broadcast_arrays(x, y, z)
            coordinates = [a.ravel() + self.scale / 2 for a in arrays]
            values = np.empty(coordinates[0].size)
            for start in range(0, len(values), 4096):
                sl = slice(start, start + 4096)
                points = np.stack([a[sl] / self.scale for a in coordinates], axis=1)
                angles = 2 * np.pi * (points @ spectrum["modes"].T) + spectrum["phases"]
                values[sl] = np.cos(angles) @ spectrum["amplitudes"] - spectrum["threshold"]
            # Microgen Spinodoid defines positive raw field as solid. Match its
            # public Shape negative-inside convention; tanh leaves the zero set.
            return np.tanh(-values.reshape(arrays[0].shape))
        from microgen.shape import surface_functions

        raw = getattr(surface_functions, self.function)(2 * np.pi * x / self.scale,
                                                      2 * np.pi * y / self.scale,
                                                      2 * np.pi * z / self.scale)
        return raw / (self.width(x, y, z) / 2)


def make_case(name, definition):
    cls = MicrogenCase if definition.get("microgen_callable") else Case
    return cls(name, name, definition["function"], definition["offset"],
               "Periodic unit cell", "Frozen shared Microgen field")
