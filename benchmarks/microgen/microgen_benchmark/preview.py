"""Display-only surface simplification; never used by validation or FEA."""

import numpy as np


def surface_preview(points, faces, scale, budget=1600):
    import pyvista as pv

    cells = np.column_stack((np.full(len(faces), 3), faces)).ravel()
    surface = pv.PolyData(points / scale, cells)
    if len(faces) > budget:
        surface = surface.decimate_pro(
            1 - budget / len(faces),
            preserve_topology=True,
            boundary_vertex_deletion=False,
        )
    triangles = np.asarray(surface.points)[surface.faces.reshape(-1, 4)[:, 1:]]
    return triangles.round(6).tolist()
