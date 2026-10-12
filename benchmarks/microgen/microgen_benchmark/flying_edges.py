"""Explicit direct-surface VTK control; failures retain the common geometry gates."""

import time


def generate(case, n):
    """Extract the zero band, cap the physical box, then conform cap edges."""
    import numpy as np
    import pyvista as pv
    from vtkmodules.vtkCommonDataModel import vtkPlane, vtkPlaneCollection
    from vtkmodules.vtkFiltersGeneral import vtkClipClosedSurface

    start = time.perf_counter()
    h = case.scale / n
    grid = pv.ImageData(
        dimensions=(n + 3,) * 3, spacing=(h,) * 3, origin=(-case.scale / 2 - h,) * 3
    )
    x = np.asarray(grid.points)
    field = case.field(*x.T)
    value = (
        field - case.band[1]
        if case.implicit
        else np.maximum(case.band[0] - field, field - case.band[1])
    )
    halo = np.any(np.abs(x / case.scale) > 0.5 + 1e-10, axis=1)
    value[halo] = 1
    grid["signed_band"] = value
    poly = grid.contour([0], scalars="signed_band", method="flying_edges", compute_normals=False)
    poly.points = np.asarray(poly.points, dtype=np.float64)
    poly = poly.clean(
        tolerance=case.scale * 1e-12,
        absolute=True,
        lines_to_points=False,
        polys_to_lines=False,
        strips_to_polys=False,
    )
    poly = poly.compute_normals(
        consistent_normals=True, auto_orient_normals=True, split_vertices=False
    )
    planes = vtkPlaneCollection()
    for axis in range(3):
        for side in [-1, 1]:
            normal = np.zeros(3)
            normal[axis] = -side
            origin = np.zeros(3)
            origin[axis] = side * case.scale / 2
            plane = vtkPlane()
            plane.SetNormal(normal)
            plane.SetOrigin(origin)
            planes.AddItem(plane)
    clip = vtkClipClosedSurface()
    clip.SetInputData(poly)
    clip.SetClippingPlanes(planes)
    clip.SetTolerance(case.scale * 1e-10)
    clip.GenerateFacesOn()
    clip.Update()
    poly = (
        pv.wrap(clip.GetOutput())
        .triangulate()
        .clean(
            tolerance=case.scale * 1e-12,
            absolute=True,
            lines_to_points=False,
            polys_to_lines=False,
            strips_to_polys=False,
        )
    )
    from microgen_benchmark.surface_cleanup import conform_cap_edges

    points = np.asarray(poly.points).copy()
    faces, cleanup = conform_cap_edges(points, poly.faces.reshape(-1, 4)[:, 1:], case.scale)
    return points, faces, time.perf_counter() - start, cleanup
