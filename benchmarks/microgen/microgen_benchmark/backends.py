"""Explicit public backends. No silent algorithm or quality fallback."""

import hashlib
import importlib.metadata
import os
import shutil
import time
from contextlib import ExitStack, nullcontext
from pathlib import Path
from tempfile import NamedTemporaryFile
from unittest.mock import patch

import numpy as np

from .cases import MICROGEN_MMG_WORKFLOWS, PROFILES


def reopenable_temporary_file(*args, **kwargs):
    """Python 3.12 permits Windows reopening while retaining context cleanup."""
    kwargs["delete_on_close"] = False
    return NamedTemporaryFile(*args, **kwargs)


def configure_mmg_cli(directory=None):
    """Use a native CLI executable; the wheel only supplies the default binary."""
    packaged = directory is None
    directory = Path(directory or importlib.metadata.distribution("mmgpy").locate_file("mmgpy/bin"))
    executable = directory / ("mmg3d_O3.exe" if os.name == "nt" else "mmg3d_O3")
    if not executable.is_file():
        raise FileNotFoundError(f"Native MMG command is unavailable: {executable}")
    if packaged and os.name != "nt" and not os.access(executable, os.X_OK):
        # Wheel extraction can lose executable bits; the package's launcher
        # restores these too. Do not change permissions on a user installation.
        executable.chmod(executable.stat().st_mode | 0o111)
    os.environ["PATH"] = str(directory.resolve()) + os.pathsep + os.environ["PATH"]
    resolved = shutil.which("mmg3d_O3")
    if resolved is None or not Path(resolved).samefile(executable):
        raise RuntimeError("Microgen's MMG command does not resolve to the selected executable")
    with executable.open("rb") as stream:
        sha = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"mmg_cli_executable": str(executable.resolve()), "mmg_cli_sha256": sha}


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

    cli_metadata = {}
    if workflow in {"mmg", "mmg-tetra", "mmg-tetra-clean", "mmg-tetra-tagged"}:
        import mmgpy
    elif workflow in MICROGEN_MMG_WORKFLOWS:
        from microgen.remesh import remesh_keeping_boundaries_for_fem

        cli_metadata = configure_mmg_cli(config.get("mmg_bin_directory"))
        progress(**cli_metadata)
    elif workflow in PROFILES:
        import meshers
    vtkSMPTools.Initialize(config["threads"])
    mmg_hausd = config.get("mmg_hausd", config["geometry_tolerance"])
    if workflow in MICROGEN_MMG_WORKFLOWS:
        cli_metadata["mmg_first_pass_hausd_over_L"] = (
            mmg_hausd if workflow == "mmg-microgen-first-distance" else None
        )
    if workflow.startswith("mmg"):
        progress(mmg_hausd_over_L=mmg_hausd, mmg_hsiz_over_L=1 / cells, mmg_hgrad=1.3)
    if not hasattr(microgen, "Tpms"):
        raise ImportError("Microgen does not provide Tpms")
    progress(stage="Microgen field setup", failure_owner="Microgen")
    start = time.perf_counter()
    tetra_seed = workflow in {
        "vtk-tetra",
        "mmg-tetra",
        "vtk-tetra-clean",
        "mmg-tetra-clean",
        "mmg-tetra-tagged",
    }
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
    if workflow in {"vtk", "mmg", *MICROGEN_MMG_WORKFLOWS} or tetra_seed:
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
        if workflow.endswith("-clean") or workflow == "mmg-tetra-tagged":
            seed_points, seed_tets = tetrahedra(grid)
            seed_points, seed_tets, cleanup = clean_tetrahedra(seed_points, seed_tets, case.scale)
            grid = pv.UnstructuredGrid(
                np.c_[np.full(len(seed_tets), 4), seed_tets].ravel(),
                np.full(len(seed_tets), pv.CellType.TETRA),
                seed_points,
            )
        extraction_s = time.perf_counter() - start - (reporting_s - reporting_start)
        if workflow in {"mmg", "mmg-tetra", "mmg-tetra-clean", "mmg-tetra-tagged"}:
            progress(stage="MMG binding remeshing", failure_owner="MMG binding pipeline")
            # Standard MMG3D adaptation of Microgen's VTK tetrahedra.
            # The default protects physical box faces while allowing MMG to
            # improve the curved surface. A fully fixed surface is optional.
            # Volume-only conversion has no triangle entities until MMG's
            # analysis. Supply actual boundary triangles before attaching
            # references or required flags; an empty list silently protects none.
            poly = grid.extract_surface(
                pass_pointid=True, algorithm="dataset_surface"
            ).triangulate()
            original = np.asarray(poly.point_data["vtkOriginalPointIds"])
            triangles = original[poly.faces.reshape(-1, 4)[:, 1:]].astype(np.int32)
            mixed = pv.UnstructuredGrid(
                np.concatenate([grid.cells, np.c_[np.full(len(triangles), 3), triangles].ravel()]),
                np.concatenate(
                    [
                        grid.celltypes,
                        np.full(len(triangles), pv.CellType.TRIANGLE, dtype=np.uint8),
                    ]
                ),
                grid.points,
            )
            # Mixed input allocates triangle storage in MMG before setting it.
            native = mmgpy.from_pyvista(mixed)
            if len(native.get_triangles()) != len(triangles) or not len(triangles):
                raise ValueError("MMG boundary triangle entities were not supplied")
            coordinates = np.asarray(grid.points)[triangles] / case.scale
            caps = np.any(
                np.all(np.isclose(coordinates, 0.5, atol=1e-9, rtol=0), axis=1)
                | np.all(np.isclose(coordinates, -0.5, atol=1e-9, rtol=0), axis=1),
                axis=1,
            )
            cleanup["mmg_boundary_triangles_supplied"] = len(triangles)
            cleanup["mmg_cap_triangles_supplied"] = int(caps.sum())
            if workflow == "mmg-tetra-tagged":
                # Distinct boundary patches preserve their interfaces while
                # allowing MMG to retriangulate the planar cap interiors.
                samples = np.asarray(grid.points)[triangles].mean(axis=1)
                value = case.field(*samples.T)
                refs = np.where(
                    np.abs(value - case.band[0]) < np.abs(value - case.band[1]), 1, 2
                ).astype(np.int64)
                for axis in range(3):
                    for side, reference in [(-0.5, 3 + 2 * axis), (0.5, 4 + 2 * axis)]:
                        cap_patch = np.all(
                            np.isclose(coordinates[:, :, axis], side, atol=1e-9, rtol=0),
                            axis=1,
                        )
                        refs[cap_patch] = reference
                native.set_triangles(triangles, refs)
                cleanup.update(
                    mmg_boundary_policy="tagged wall and box patches",
                    mmg_boundary_patch_counts={
                        str(ref): int((refs == ref).sum()) for ref in np.unique(refs)
                    },
                    mmg_required_cap_triangles=0,
                )
            else:
                native.set_required_triangles(np.flatnonzero(caps).astype(np.int32))
                cleanup.update(
                    mmg_boundary_policy="required box cap triangles",
                    mmg_required_cap_triangles=int(caps.sum()),
                )
            native.remesh(
                hsiz=case.scale / cells,
                hausd=mmg_hausd * case.scale,
                hgrad=1.3,
                nosurf=config["mmg_boundary"] == "fixed",
                verbose=-1,
            )
            grid = mmgpy.to_pyvista(native)
        elif workflow in MICROGEN_MMG_WORKFLOWS:
            progress(stage="Microgen MMG remeshing", failure_owner="Microgen / MMG pipeline")
            # Microgen reopens a delete=True temporary file before its context
            # closes. POSIX permits this; Windows needs delete_on_close=False.
            # Scope the file-lifecycle compatibility to this isolated worker.
            context = (
                patch("microgen.remesh.NamedTemporaryFile", reopenable_temporary_file)
                if os.name == "nt"
                else nullcontext()
            )
            with ExitStack() as contexts:
                contexts.enter_context(context)
                if workflow == "mmg-microgen-first-distance":
                    from microgen.external import Mmg

                    original_mmg = Mmg.mmg3d

                    def distance_control(**parameters):
                        if parameters.get("nofem"):
                            parameters["hausd"] = mmg_hausd * case.scale
                        return original_mmg(**parameters)

                    contexts.enter_context(patch.object(Mmg, "mmg3d", distance_control))
                grid = remesh_keeping_boundaries_for_fem(
                    grid,
                    periodic=config.get("periodic", False),
                    hsiz=case.scale / cells,
                    hausd=mmg_hausd * case.scale,
                    hgrad=1.3,
                )
        points, tets = tetrahedra(grid)
    elif workflow in PROFILES:
        native_band = (
            None
            if case.implicit and config.get("meshers_solid_region", "band") == "negative"
            else case.band
        )
        region_metadata = {
            "meshers_region": "negative" if native_band is None else "band",
            "meshers_band": native_band,
        }
        progress(stage="Meshers native generation", failure_owner="Meshers", **region_metadata)
        native = meshers.generate(
            case.field,
            bounds=(-case.scale / 2, case.scale / 2) * 3,
            cells=cells,
            band=native_band,
            minimum_quality=config["minimum_quality"],
            geometry_tolerance=config["geometry_tolerance"] * case.scale,
            max_tetrahedra=config["max_tetrahedra"],
            threads=config["threads"],
            periodic=(config.get("periodic", False),) * 3,
            **PROFILES[workflow],
        )
        points, tets = np.asarray(native.points).copy(), np.asarray(native.tetrahedra).copy()
        native_diagnostics = {
            **region_metadata,
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
            **cli_metadata,
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
