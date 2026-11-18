"""Measure one live Schwarz P workflow at a declared effective-stiffness target."""

import argparse
import datetime as dt
import hashlib
import inspect
import json
import os
import shutil
import time
from contextlib import ExitStack
from importlib.metadata import distribution, version
from pathlib import Path
from tempfile import NamedTemporaryFile
from types import SimpleNamespace
from unittest.mock import patch

import fedoo as fd
import microgen
import microgen.remesh
import numpy as np
import pyvista as pv
from fedoo_homogenization import homogenize
from microgen.shape import surface_functions
from periodic_geometry import BASELINE, audit
from release_wheel_comparison import generate, hardware_snapshot
from scipy.linalg import eigvalsh
from scipy.spatial import cKDTree
from vtkmodules.vtkCommonCore import vtkSMPTools

REFERENCE = BASELINE.with_name("schwarz-reference.json")
PROBE = np.array([-0.364, -0.347, -0.313])
PACKAGES = {
    "meshers": "0.2.1",
    "microgen": "2.0.0b1",
    "mmgpy": "0.17.1",
    "fedoo": "1.0.0b2",
    "numpy": "2.5.3",
    "scipy": "1.18.1",
    "pyvista": "0.49.0",
    "vtk": "9.7.1",
}
SOURCE_HASHES = {
    "microgen/remesh.py": "b93c8b8150c1c40ccf5312d57785331d52f114d2217a8b59296b7484cf6f2612",
    "microgen/shape/tpms.py": "1c2876a8f65441e6d7a0f21397a3b781b6f03a338f244c7db5c462c13348f955",
}


def digest(path):
    """Return the SHA-256 of an input or exported mesh."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reopenable(*args, **kwargs):
    """Return a Windows-reopenable file with context cleanup."""
    return NamedTemporaryFile(*args, **kwargs, delete_on_close=False)


def live_microgen(route, cells, folder, row, budget):
    """Write the volume produced by the pinned Microgen VTK and MMG implementations.

    Returns:
        Path to the exported tetrahedra.

    Raises:
        RuntimeError: A source, MMG executable or mesh size fails its precondition.
    """
    sources = {
        "microgen/shape/tpms.py": digest(inspect.getfile(microgen.Tpms)),
        "microgen/remesh.py": digest(microgen.remesh.__file__),
    }
    row["source_sha256"] = sources
    if sources != SOURCE_HASHES:
        raise RuntimeError("Microgen source differs from the frozen VTK route")
    folder_bin = Path(distribution("mmgpy").locate_file("mmgpy/bin")).resolve()
    os.environ["PATH"] = str(folder_bin) + os.pathsep + os.environ["PATH"]
    executable = shutil.which("mmg3d_O3")
    if executable is None or Path(executable).resolve().parent != folder_bin:
        raise RuntimeError("The mmgpy wheel's MMG executable is unavailable")
    row["mmg_sha256"] = digest(executable)
    start = time.perf_counter()
    settings = {"resolution": cells} if cells else {}
    shape = microgen.Tpms(surface_functions.schwarz_p, offset=0.0, **settings)
    seed = shape.generate_volume_mesh(type_part="lower skeletal").triangulate()
    row.update(
        resolution=int(shape.resolution),
        seed_seconds=time.perf_counter() - start,
        seed_nodes=seed.n_points,
        seed_cells=seed.n_cells,
    )
    if max(seed.n_points, seed.n_cells) > budget:
        raise RuntimeError("VTK seed exceeds the mesh budget")
    options = {"hausd": 0.001, "hmax": 0.06, "hgrad": 1.3} if route == "microgen-tuned" else {}
    original = microgen.Mmg.mmg3d
    row["mmg_calls"] = []

    def measured(**params):
        options_used = {k: v for k, v in params.items() if k not in {"input", "output"}}
        row["mmg_calls"].append(options_used)
        begin = time.perf_counter()
        try:
            return original(**params)
        finally:
            options_used["seconds"] = time.perf_counter() - begin

    with ExitStack() as stack:
        if os.name == "nt":
            stack.enter_context(patch.object(microgen.remesh, "NamedTemporaryFile", reopenable))
        stack.enter_context(patch.object(microgen.Mmg, "mmg3d", measured))
        grid = microgen.remesh.remesh_keeping_boundaries_for_fem(seed, periodic=True, **options)
    path = folder / "volume.vtu"
    grid.save(path)
    return path


def solver_input(grid, quality):
    """Return Tet4 connectivity and low/high node correspondences for fedoo.

    Raises:
        ValueError: Opposite faces lack a bijective node map.
    """
    points = np.asarray(grid.points, dtype=float)
    boundary = grid.extract_surface(pass_pointid=True, algorithm="dataset_surface")
    surface = boundary.point_data["vtkOriginalPointIds"][boundary.regular_faces]
    pairs = []
    for axis in range(3):
        low, high = [np.flatnonzero(abs(points[:, axis] - side) < 1e-12) for side in (-0.5, 0.5)]
        other = [a for a in range(3) if a != axis]
        distance, index = cKDTree(points[high][:, other]).query(points[low][:, other])
        if len(low) != len(high) or np.max(distance) >= 1e-12 or len(np.unique(index)) != len(high):
            raise ValueError("Could not reconstruct a bijective periodic node map")
        pairs.append(np.column_stack([low, high[index]]))
    return SimpleNamespace(
        points=points + 0.5,
        tetrahedra=grid.cells_dict[pv.CellType.TETRA],
        surface=surface,
        periodic_pairs=pairs,
        diagnostics={"minimum_mmg_quality": quality},
    )


def run(route, cells, folder, row):
    """Generate, audit and solve one mesh, preserving the first failed stage.

    Raises:
        ValueError: The fixed interior stress probe is missing.
        RuntimeError: A package version or VTK backend differs from the study.
    """
    if row["versions"] != PACKAGES:
        raise RuntimeError("Install the pinned comparison packages")
    vtkSMPTools.SetBackend("Sequential")
    if vtkSMPTools.GetBackend() != "Sequential":
        raise RuntimeError("A single-threaded VTK backend is required")
    baseline = json.loads(BASELINE.read_text())
    start = time.perf_counter()
    if route == "native":
        cfg = dict(baseline["routes"]["schwarz_native"]["config"], cells=cells)
        row["config"] = cfg
        mesh = generate("schwarz_native", cfg, baseline["limits"])
        path = folder / "volume.vtkhdf"
        mesh.write_vtkhdf(path)
    else:
        path = live_microgen(route, cells, folder, row, baseline["limits"]["max_nodes_or_elements"])
    grid = pv.read(path)
    row["generation_and_io_seconds"] = time.perf_counter() - start
    row["mesh_sha256"] = digest(path)
    start = time.perf_counter()
    geometry = audit(grid, "schwarz", baseline)
    row.update(geometry=geometry, audit_seconds=time.perf_counter() - start)
    if not geometry["sampled_admission"]:
        row["status"] = "geometry-rejected"
        return
    probe = int(grid.find_containing_cell(PROBE))
    if probe < 0:
        raise ValueError("The fixed interior stress probe is missing")
    row["stress_probe_xyz"] = PROBE.tolist()
    mesh = solver_input(grid, geometry["minimum_quality"])
    assembly_seconds = 0.0
    original_assembly = fd.Assembly.assemble_global_mat

    def timed_assembly(assembly, *args, **kwargs):
        nonlocal assembly_seconds
        begin = time.perf_counter()
        try:
            return original_assembly(assembly, *args, **kwargs)
        finally:
            assembly_seconds += time.perf_counter() - begin

    with patch.object(fd.Assembly, "assemble_global_mat", timed_assembly):
        fea = homogenize(mesh, stress_probe_cell=probe)
    row["assembly_seconds"] = assembly_seconds
    row["fea"] = fea
    ref = np.array(json.loads(REFERENCE.read_text())["stiffness"])
    tensor = np.array(fea["stiffness_engineering_voigt"])
    delta = (tensor + tensor.T - ref - ref.T) / 2
    discrepancy = float(abs(eigvalsh(delta, (ref + ref.T) / 2)).max())
    row["relative_energy_discrepancy"] = discrepancy
    row["status"] = "within-reference-target" if discrepancy <= 0.01 else "outside-reference-target"


def main():
    """Save a timed workflow record including rejected meshes and solver failures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("route", choices=("native", "microgen-default", "microgen-tuned"))
    parser.add_argument("output", type=Path)
    parser.add_argument("--cells", type=int, help="VTK resolution or native background cells")
    args = parser.parse_args()
    if args.cells is not None and args.cells < 2:
        parser.error("--cells must be at least 2")
    if args.route != "microgen-default" and args.cells is None:
        parser.error("Specify --cells for this route")
    args.output.mkdir()
    hardware = hardware_snapshot()
    row = {
        "recorded_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "route": args.route,
        "cells": args.cells,
        "platform": hardware["platform"],
        "hardware": hardware,
        "reference_sha256": digest(REFERENCE),
        "versions": {pkg: version(pkg) for pkg in PACKAGES},
    }
    start = time.perf_counter()
    try:
        run(args.route, args.cells, args.output, row)
    except Exception as exc:
        row.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        row["workflow_seconds"] = time.perf_counter() - start
        (args.output / "result.json").write_text(
            json.dumps(row, indent=2, allow_nan=False) + "\n", encoding="utf8"
        )


if __name__ == "__main__":
    main()
