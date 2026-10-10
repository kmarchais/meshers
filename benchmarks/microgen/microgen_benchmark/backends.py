"""Explicit public backends. No silent algorithm or quality fallback."""

import os
import shutil
import sys
import time
from contextlib import nullcontext
from pathlib import Path
from tempfile import NamedTemporaryFile
from unittest.mock import patch

import numpy as np

from .cases import PROFILES


def reopenable_temporary_file(*args, **kwargs):
    """Python 3.12 permits Windows reopening while retaining context cleanup."""
    kwargs["delete_on_close"] = False
    return NamedTemporaryFile(*args, **kwargs)


def tetrahedra(grid):
    import pyvista as pv

    types = np.asarray(grid.celltypes)
    if not np.all(types == pv.CellType.TETRA):
        raise ValueError("Backend returned non-linear-tetrahedral cells")
    return np.asarray(grid.points).copy(), np.asarray(grid.cells_dict[pv.CellType.TETRA]).copy()


def generate(case, workflow, cells, config, progress=None):
    callback = progress or (lambda **values: None)
    reporting_s = 0.0

    def progress(**values):
        nonlocal reporting_s
        start = time.perf_counter()
        callback(**values)
        reporting_s += time.perf_counter() - start

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
    progress(stage="Microgen field setup", failure_owner="Microgen")
    start = time.perf_counter()
    shape = case.shape(cells)
    setup_s = time.perf_counter() - start
    progress(setup_s=setup_s)
    start = time.perf_counter()
    reporting_start = reporting_s
    if workflow in {"vtk", "mmg", "mmg-microgen"}:
        progress(stage="Microgen VTK extraction", failure_owner="Microgen / VTK")
        grid = case.volume_grid(shape, cells).triangulate()
        extraction_s = time.perf_counter() - start - (reporting_s - reporting_start)
        if workflow == "mmg":
            progress(stage="MMG binding remeshing", failure_owner="MMG binding pipeline")
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
            progress(stage="Microgen MMG remeshing", failure_owner="Microgen / MMG pipeline")
            # Microgen reopens a delete=True temporary file before its context
            # closes. POSIX permits this; Windows needs delete_on_close=False.
            # Scope the file-lifecycle compatibility to this isolated worker.
            context = (
                patch("microgen.remesh.NamedTemporaryFile", reopenable_temporary_file)
                if os.name == "nt"
                else nullcontext()
            )
            with context:
                grid = remesh_keeping_boundaries_for_fem(
                    grid,
                    periodic=False,
                    hsiz=case.scale / cells,
                    hausd=config["geometry_tolerance"] * case.scale,
                    hgrad=1.3,
                )
        points, tets = tetrahedra(grid)
    elif workflow in PROFILES:
        progress(stage="Meshers native generation", failure_owner="Meshers")
        native = meshers.generate(
            case.field,
            bounds=(-case.scale / 2, case.scale / 2) * 3,
            cells=cells,
            band=case.band,
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
    generation_s = time.perf_counter() - start - (reporting_s - reporting_start)
    progress(mesh_s=setup_s + generation_s, generation_s=generation_s)
    # Check the physical field outside the generation timer. Microgen stores
    # the raw periodic field, before its sheet half-offset is applied.
    progress(stage="shared field preflight", failure_owner="Benchmark field adapter")
    sample_count = case.verify_field(shape, cells)
    return (
        points,
        tets,
        {
            "setup_s": setup_s,
            "generation_s": generation_s,
            "mesh_s": setup_s + generation_s,
            "vtk_extraction_s": extraction_s,
            "field_samples_checked": sample_count,
            "vtk_version": pv.vtk_version_info,
        },
    )
