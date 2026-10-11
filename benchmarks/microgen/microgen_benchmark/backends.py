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


def clean_tetrahedra(points, tets, scale):
    """Declared floating-point cleanup for an additional VTK control only."""
    normalized = np.asarray(points, float) / scale
    _, representatives, inverse = np.unique(
        np.round(normalized / 1e-12), axis=0, return_index=True, return_inverse=True
    )
    points, tets = np.asarray(points)[representatives], inverse[tets]
    v = points[tets] / scale
    det = np.linalg.det(v[:, 1:] - v[:, :1])
    keep = np.abs(det) > 1e-15
    return (
        points,
        tets[keep],
        {
            "cleanup_merge_tolerance_over_L": 1e-12,
            "cleanup_det_tolerance_over_L3": 1e-15,
            "cleanup_merged_nodes": len(normalized) - len(points),
            "cleanup_discarded_cells": int((~keep).sum()),
            "cleanup_discarded_abs_volume_over_L3": float(np.abs(det[~keep]).sum() / 6),
        },
    )


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

    if workflow in {"mmg", "mmg-tetra", "mmg-tetra-clean"}:
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
    tetra_seed = workflow in {"vtk-tetra", "mmg-tetra", "vtk-tetra-clean", "mmg-tetra-clean"}
    explicit_field = workflow in PROFILES or tetra_seed
    shape = None if explicit_field else case.shape(cells)
    if explicit_field and case.function == "spinodoid":
        from .cases import spinodoid_source

        spinodoid_source(case.scale)
    setup_s = time.perf_counter() - start
    progress(setup_s=setup_s)
    start = time.perf_counter()
    reporting_start = reporting_s
    cleanup = {}
    native_diagnostics = {}
    if workflow in {"vtk", "mmg", "mmg-microgen"} or tetra_seed:
        progress(stage="Microgen VTK extraction", failure_owner="Microgen / VTK")
        if tetra_seed:
            # Triangulate the background first, then clip linear tetrahedra.
            # This is a distinct control, not a repair of Microgen's output.
            grid = (
                pv.ImageData(
                    dimensions=(cells + 1,) * 3,
                    spacing=(case.scale / cells,) * 3,
                    origin=(-case.scale / 2,) * 3,
                )
                .cast_to_unstructured_grid()
                .triangulate()
            )
            grid.point_data["benchmark_field"] = case.field(*np.asarray(grid.points).T)
            grid = grid.clip_scalar(scalars="benchmark_field", value=case.band[0], invert=False)
            grid = grid.clip_scalar(
                scalars="benchmark_field", value=case.band[1], invert=True
            ).triangulate()
        else:
            grid = case.volume_grid(shape, cells).triangulate()
        if workflow.endswith("-clean"):
            seed_points, seed_tets = tetrahedra(grid)
            seed_points, seed_tets, cleanup = clean_tetrahedra(seed_points, seed_tets, case.scale)
            grid = pv.UnstructuredGrid(
                np.c_[np.full(len(seed_tets), 4), seed_tets].ravel(),
                np.full(len(seed_tets), pv.CellType.TETRA),
                seed_points,
            )
        extraction_s = time.perf_counter() - start - (reporting_s - reporting_start)
        if workflow in {"mmg", "mmg-tetra", "mmg-tetra-clean"}:
            progress(stage="MMG binding remeshing", failure_owner="MMG binding pipeline")
            # Standard MMG3D adaptation of Microgen's VTK tetrahedra.
            # The default protects physical box faces while allowing MMG to
            # improve the curved surface. A fully fixed surface is optional.
            native = mmgpy.from_pyvista(grid)
            triangles = np.asarray(native.get_triangles())
            coordinates = np.asarray(grid.points)[triangles] / case.scale
            caps = np.any(
                np.all(np.isclose(coordinates, 0.5, atol=1e-9, rtol=0), axis=1)
                | np.all(np.isclose(coordinates, -0.5, atol=1e-9, rtol=0), axis=1),
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
        native_diagnostics = {
            "meshers_evaluator": native.diagnostics.get("evaluator"),
            "meshers_callback_calls": native.diagnostics.get("callback_calls"),
            "meshers_field_preparation_s": native.diagnostics.get("preparation_seconds"),
        }
        extraction_s = None
    else:
        raise ValueError(f"Unknown workflow {workflow}")
    generation_s = time.perf_counter() - start - (reporting_s - reporting_start)
    progress(mesh_s=setup_s + generation_s, generation_s=generation_s)
    # Check the physical field outside the generation timer. Microgen stores
    # the raw periodic field, before its sheet half-offset is applied.
    progress(stage="shared field preflight", failure_owner="Benchmark field adapter")
    sample_count = case.verify_field(shape if shape is not None else case.shape(cells), cells)
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
            **cleanup,
            **native_diagnostics,
        },
    )


def generate_surface(case, workflow, cells, config, progress=None):
    """Time the public surface path, or the full volume-to-surface cost."""
    if workflow == "vtk-flying-edges":
        import pyvista as pv
        from vtkmodules.vtkCommonCore import vtkSMPTools

        from .flying_edges import generate as flying_edges

        vtkSMPTools.Initialize(config["threads"])
        if progress:
            progress(
                stage="Flying Edges surface extraction", failure_owner="VTK Flying Edges control"
            )
        start = time.perf_counter()
        if case.function == "spinodoid":
            from .cases import spinodoid_source

            spinodoid_source(case.scale)
        setup_s = time.perf_counter() - start
        points, faces, generation_s, cleanup = flying_edges(case, cells)
        checked = case.verify_field(case.shape(cells), cells)
        return (
            points,
            faces,
            {
                "setup_s": setup_s,
                "generation_s": generation_s,
                "mesh_s": setup_s + generation_s,
                "field_samples_checked": checked,
                "vtk_version": pv.vtk_version_info,
                "generated_tetrahedra": 0,
                "surface_pipeline": "Flying Edges, positive halo, six planar caps, conforming cap edges",
                "surface_grid_dimensions": [cells + 3] * 3,
                **cleanup,
            },
        )
    if workflow != "vtk":
        from .geometry import boundary

        points, tets, timing = generate(case, workflow, cells, config, progress)
        start = time.perf_counter()
        # Orientation normalization is necessary before oriented face extraction.
        from .geometry import canonical_mesh

        points, tets, _, _ = canonical_mesh(points, tets)
        faces = boundary(points, tets)
        timing["surface_extraction_s"] = time.perf_counter() - start
        timing["mesh_s"] += timing["surface_extraction_s"]
        timing["surface_pipeline"] = "volume generation plus boundary extraction"
        timing["generated_tetrahedra"] = len(tets)
        return points, faces, timing
    import pyvista as pv
    from vtkmodules.vtkCommonCore import vtkSMPTools

    vtkSMPTools.Initialize(config["threads"])
    if progress:
        progress(stage="Microgen surface field setup", failure_owner="Microgen / VTK")
    start = time.perf_counter()
    shape = case.shape(cells)
    setup_s = time.perf_counter() - start
    start = time.perf_counter()
    poly = shape.generate_surface_mesh(resolution=cells + 1).triangulate()
    faces = poly.faces.reshape(-1, 4)[:, 1:]
    points = np.asarray(poly.points).copy()
    generation_s = time.perf_counter() - start
    checked = case.verify_field(shape, cells)
    return (
        points,
        faces.copy(),
        {
            "setup_s": setup_s,
            "generation_s": generation_s,
            "mesh_s": setup_s + generation_s,
            "field_samples_checked": checked,
            "vtk_version": pv.vtk_version_info,
            "generated_tetrahedra": 0,
            "surface_pipeline": "Microgen public surface generator",
        },
    )
