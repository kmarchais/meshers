"""Explicit public backends. No silent algorithm or quality fallback."""

import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np

from .cases import PROFILES


def tetrahedra(grid):
    import pyvista as pv

    types = np.asarray(grid.celltypes)
    if not np.all(types == pv.CellType.TETRA):
        raise ValueError("Backend returned non-linear-tetrahedral cells")
    return np.asarray(grid.points).copy(), np.asarray(grid.cells_dict[pv.CellType.TETRA]).copy()


def generate(case, workflow, cells, config):
    import microgen  # Import cost is excluded equally for every workflow.
    import pyvista as pv
    from vtkmodules.vtkCommonCore import vtkSMPTools

    if workflow == "mmg":
        import mmgpy
    elif workflow == "mmg-microgen":
        from microgen.remesh import remesh_keeping_boundaries_for_fem

        # mmgpy installs Microgen's conventional command alias beside Python.
        os.environ["PATH"] = str(Path(sys.executable).parent) + os.pathsep + os.environ["PATH"]
        if not shutil.which("mmg3d_O3"):
            raise ImportError("mmgpy's mmg3d_O3 command alias is unavailable")
    elif workflow in PROFILES:
        import meshers
    vtkSMPTools.Initialize(config["threads"])
    if not hasattr(microgen, "Tpms"):
        raise ImportError("Microgen does not provide Tpms")
    start = time.perf_counter()
    shape = case.shape(cells)
    setup_s = time.perf_counter() - start
    start = time.perf_counter()
    if workflow in {"vtk", "mmg", "mmg-microgen"}:
        grid = shape.grid_sheet.triangulate()
        extraction_s = time.perf_counter() - start
        if workflow == "mmg":
            # Standard MMG3D adaptation of Microgen's VTK tetrahedra.
            # The default protects physical box faces while allowing MMG to
            # improve the curved surface. A fully fixed surface is optional.
            native = mmgpy.from_pyvista(grid)
            triangles = np.asarray(native.get_triangles())
            coordinates = np.asarray(grid.points)[triangles] / case.scale
            caps = np.any(
                np.all(np.isclose(coordinates, 0.5, atol=1e-9), axis=1)
                | np.all(np.isclose(coordinates, -0.5, atol=1e-9), axis=1),
                axis=1,
            )
            native.set_required_triangles(np.flatnonzero(caps).astype(np.int32))
            native.remesh(
                hsiz=case.scale / cells,
                hausd=config["geometry_tolerance"] * case.scale,
                hgrad=1.3,
                nosurf=config["mmg_boundary"] == "fixed",
                verbose=-1,
            )
            grid = mmgpy.to_pyvista(native)
        elif workflow == "mmg-microgen":
            grid = remesh_keeping_boundaries_for_fem(
                grid,
                periodic=False,
                hsiz=case.scale / cells,
                hausd=config["geometry_tolerance"] * case.scale,
                hgrad=1.3,
            )
        points, tets = tetrahedra(grid)
    elif workflow in PROFILES:
        native = meshers.generate(
            case.field,
            bounds=(-case.scale / 2, case.scale / 2) * 3,
            cells=cells,
            band=(-1, 1),
            minimum_quality=config["minimum_quality"],
            geometry_tolerance=config["geometry_tolerance"] * case.scale,
            max_tetrahedra=config["max_tetrahedra"],
            threads=config["threads"],
            **PROFILES[workflow],
        )
        points, tets = np.asarray(native.points).copy(), np.asarray(native.tetrahedra).copy()
        extraction_s = None
    else:
        raise ValueError(f"Unknown workflow {workflow}")
    generation_s = time.perf_counter() - start
    # Check the physical field outside the generation timer. Microgen stores
    # the raw periodic field, before its sheet half-offset is applied.
    sample = np.asarray(shape.grid.points)
    expected = case.field(*sample.T) * 0.5 * case.width(*sample.T)
    raw = np.asarray(shape.grid.point_data["surface"])
    if not np.allclose(expected, raw, rtol=1e-12, atol=1e-12):
        raise ValueError("Case formula and Microgen's sampled surface field disagree")
    offset = np.broadcast_to(np.asarray(shape.offset).ravel(order="F"), (len(sample),))
    if not np.allclose(offset, case.width(*sample.T), rtol=1e-12, atol=1e-12):
        raise ValueError("Actual Microgen offset differs from the declared graded field")
    return (
        points,
        tets,
        {
            "setup_s": setup_s,
            "generation_s": generation_s,
            "mesh_s": setup_s + generation_s,
            "vtk_extraction_s": extraction_s,
            "field_samples_checked": len(sample),
            "vtk_version": pv.vtk_version_info,
        },
    )
