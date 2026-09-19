"""A varying-thickness gyroid, without periodic constraints."""

import meshers
import numpy as np


def graded_gyroid(x, y, z):
    x_phase, y_phase, z_phase = (2 * np.pi * v for v in (x, y, z))
    g = (
        np.sin(x_phase) * np.cos(y_phase)
        + np.sin(y_phase) * np.cos(z_phase)
        + np.sin(z_phase) * np.cos(x_phase)
    )
    thickness = 0.35 + 0.15 * x
    return g / thickness


mesh = meshers.generate(
    graded_gyroid, cells=16, band=(-1.0, 1.0), geometry_tolerance=0.1, optimize_passes=1
)
mesh.write_vtkhdf("graded.vtkhdf")
print(mesh.diagnostics)
