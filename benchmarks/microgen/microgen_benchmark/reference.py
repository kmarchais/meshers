"""Fixed sampling diagnostics independent of the evaluated tetrahedral meshes."""

import numpy as np
from scipy.stats import qmc

from .preview import surface_preview


def target(case):
    sample = (qmc.Sobol(3, scramble=True, seed=20261010).random_base2(18) - 0.5) * case.scale
    material = case.material(*sample.T)
    volume_fraction = float(material.mean())
    coarse_fraction = float(material[: 2**17].mean())
    record = {
        "target_volume_fraction": volume_fraction,
        "target_volume_sample_count": len(sample),
        "target_volume_half_sample_delta": abs(volume_fraction - coarse_fraction),
        "target_preview_intervals": 32,
    }
    try:
        shape = case.shape(32)
        surface = case.volume_grid(shape, 32).extract_surface().triangulate()
        faces = surface.faces.reshape(-1, 4)[:, 1:]
        record["target_preview"] = surface_preview(np.asarray(surface.points), faces, case.scale)
    except Exception as exc:
        record["target_preview_error"] = f"{type(exc).__name__}: {exc}"
    return record
