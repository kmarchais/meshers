# Copyright (c) 2026 Kevin Marchais

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
from dataclasses import dataclass
from functools import partial
from importlib.metadata import distribution, version
from pathlib import Path
from tempfile import NamedTemporaryFile
from types import SimpleNamespace
from typing import NotRequired, TypedDict, Unpack
from unittest.mock import patch

import fedoo as fd
import numpy as np
import pyvista as pv
from fedoo_homogenization import AuditedCG, homogenize
from historical_schwarz_repair import repaired_grid
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


class _MmgPass(TypedDict):
    """Keywords supplied by Microgen's pinned two-pass remesher."""

    input: str
    output: str
    nofem: NotRequired[bool]
    ls: NotRequired[bool]
    nr: NotRequired[bool]
    hausd: NotRequired[float | None]
    hgrad: NotRequired[float | None]
    hmax: NotRequired[float | None]
    hmin: NotRequired[float | None]
    hsiz: NotRequired[float | None]


def _mmg_executable():
    folder_bin = Path(distribution("mmgpy").locate_file("mmgpy/bin")).resolve()
    os.environ["PATH"] = str(folder_bin) + os.pathsep + os.environ["PATH"]
    executable = shutil.which("mmg3d_O3")
    if executable is None or Path(executable).resolve().parent != folder_bin:
        msg = "The mmgpy wheel's MMG executable is unavailable"
        raise RuntimeError(msg)
    return executable


@dataclass(frozen=True)
class SchwarzSettings:
    """Identify the mesh route, resolution and physical settings for one trial."""

    route: str
    cells: int | None
    native_sizing: dict[str, float] | None = None
    shifted: bool = False
    repair_hmax: float = 0.1


def live_microgen(trial: SchwarzSettings, folder, row, budget, *, observe=None):
    """Write the volume produced by the pinned Microgen VTK and MMG implementations.

    Returns:
        Path to the exported tetrahedra.

    Raises:
        RuntimeError: A source, MMG executable or mesh size fails its precondition.
    """
    # Only the historical comparison needs the VTK and MMG backend from Microgen.
    import microgen.remesh  # ruff: ignore[import-outside-top-level]
    from microgen.shape import surface_functions  # ruff: ignore[import-outside-top-level]

    sources = {
        "microgen/shape/tpms.py": digest(inspect.getfile(microgen.Tpms)),
        "microgen/remesh.py": digest(microgen.remesh.__file__),
    }
    row["source_sha256"] = sources
    if sources != SOURCE_HASHES:
        msg = "Microgen source differs from the frozen VTK route"
        raise RuntimeError(msg)
    row["mmg_sha256"] = digest(_mmg_executable())
    start = time.perf_counter()
    settings = {"resolution": trial.cells} if trial.cells else {}
    shape = microgen.Tpms(surface_functions.schwarz_p, offset=0.0, **settings)
    seed = shape.generate_volume_mesh(type_part="lower skeletal").triangulate()
    row.update(
        resolution=int(shape.resolution),
        seed_seconds=time.perf_counter() - start,
        seed_nodes=seed.n_points,
        seed_cells=seed.n_cells,
    )
    if max(seed.n_points, seed.n_cells) > budget:
        msg = "VTK seed exceeds the mesh budget"
        raise RuntimeError(msg)
    if observe is not None:
        observe("seed", seed)
    options = (
        {"hausd": 0.001, "hmax": 0.06, "hgrad": 1.3}
        if trial.route in {"microgen-tuned", "single-pass-repair"}
        else {}
    )
    single_pass = trial.route in {"single-pass-repair", "microgen-repaired"}
    original = microgen.Mmg.mmg3d
    row["mmg_calls"] = []

    def measured(**params: Unpack[_MmgPass]):
        if single_pass and params.get("ls"):
            shutil.copyfile(params["input"], params["output"])
            return None
        if single_pass:
            params.update(hausd=0.0005, hmax=trial.repair_hmax, hgrad=1.3)
        options_used = {k: v for k, v in params.items() if k not in {"input", "output"}}
        row["mmg_calls"].append(options_used)
        begin = time.perf_counter()
        try:
            result = original(**params)
        finally:
            options_used["seconds"] = time.perf_counter() - begin
        if observe is not None:
            observe(f"pass-{len(row['mmg_calls'])}", pv.UnstructuredGrid(params["output"]))
        return result

    with ExitStack() as stack:
        if os.name == "nt":
            reopenable = partial(NamedTemporaryFile, delete_on_close=False)
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
            msg = "Could not reconstruct a bijective periodic node map"
            raise ValueError(msg)
        pairs.append(np.column_stack([low, high[index]]))
    return SimpleNamespace(
        points=points + 0.5,
        tetrahedra=grid.cells_dict[pv.CellType.TETRA],
        surface=surface,
        periodic_pairs=pairs,
        diagnostics={"minimum_mmg_quality": quality},
    )


def run(trial: SchwarzSettings, folder, row, *, solver_factory=AuditedCG):
    """Generate, audit and solve one mesh, preserving the first failed stage.

    Raises:
        RuntimeError: A package version or VTK backend differs from the study.
    """
    if trial.route != "native-candidate" and row["versions"] != PACKAGES:
        msg = "Install the pinned comparison packages"
        raise RuntimeError(msg)
    vtkSMPTools.SetBackend("Sequential")
    if vtkSMPTools.GetBackend() != "Sequential":
        msg = "A single-threaded VTK backend is required"
        raise RuntimeError(msg)
    baseline = json.loads(BASELINE.read_text())
    start = time.perf_counter()
    if trial.route in {"native", "native-candidate"}:
        cfg = dict(baseline["routes"]["schwarz_native"]["config"], cells=trial.cells)
        cfg.update(trial.native_sizing or {})
        if trial.shifted:
            cfg["phase_x"] = 0.125
        row["config"] = cfg
        mesh = generate("schwarz_native", cfg, baseline["limits"])
        path = folder / "volume.vtkhdf"
        mesh.write_vtkhdf(path)
    else:
        historical = trial.route == "microgen-repaired"
        path = live_microgen(
            trial,
            folder,
            row,
            baseline["limits"]["max_nodes_or_elements"],
        )
        if historical:
            source = pv.read(path)
            if max(source.n_points, source.n_cells) > baseline["limits"]["max_nodes_or_elements"]:
                msg = "Single-pass MMG output exceeds the mesh budget before repair"
                raise RuntimeError(msg)
            repair_start = time.perf_counter()
            repaired = repaired_grid(source, row)
            row["repair_seconds"] = time.perf_counter() - repair_start
            path = folder / "repaired.vtu"
            repaired.save(path)
    grid = pv.read(path)
    row["generation_and_io_seconds"] = time.perf_counter() - start
    row["mesh_sha256"] = digest(path)
    start = time.perf_counter()
    geometry = audit(grid, "schwarz-shifted" if trial.shifted else "schwarz", baseline)
    row.update(geometry=geometry, audit_seconds=time.perf_counter() - start)
    if not geometry["sampled_admission"]:
        row["status"] = "geometry-rejected"
        return
    solve_admitted(grid, geometry, row, solver_factory=solver_factory)


def solve_admitted(grid, geometry, row, *, solver_factory=AuditedCG):
    """Solve six strain loads and their stress-controlled checks on an admitted grid.

    Raises:
        ValueError: The fixed stress probe lies outside the returned mesh.
    """
    probe = int(grid.find_containing_cell(PROBE))
    if probe < 0:
        msg = "The fixed interior stress probe is missing"
        raise ValueError(msg)
    row["stress_probe_xyz"] = PROBE.tolist()
    mesh = solver_input(grid, geometry["minimum_quality"])
    mesh_hash = hashlib.sha256()
    mesh_hash.update(np.ascontiguousarray(mesh.points, dtype="<f8").tobytes())
    mesh_hash.update(np.ascontiguousarray(mesh.tetrahedra, dtype="<i8").tobytes())
    row["solver_mesh_sha256"] = mesh_hash.hexdigest()
    assembly_seconds = 0.0
    original_assembly = fd.Assembly.assemble_global_mat

    def timed_assembly(assembly: fd.Assembly, compute: str = "all") -> None:
        nonlocal assembly_seconds
        begin = time.perf_counter()
        try:
            return original_assembly(assembly, compute)
        finally:
            assembly_seconds += time.perf_counter() - begin

    with patch.object(fd.Assembly, "assemble_global_mat", timed_assembly):
        fea = homogenize(mesh, stress_probe_cell=probe, solver_factory=solver_factory)
    row["assembly_seconds"] = assembly_seconds
    row["fea"] = fea
    ref = np.array(json.loads(REFERENCE.read_text())["stiffness"])
    tensor = np.array(fea["stiffness_engineering_voigt"])
    delta = (tensor + tensor.T - ref - ref.T) / 2
    discrepancy = float(abs(eigvalsh(delta, (ref + ref.T) / 2)).max())
    row["relative_energy_discrepancy"] = discrepancy
    row["status"] = "within-reference-target" if discrepancy <= 0.01 else "outside-reference-target"


def main():
    """Save a workflow record.

    Returns:
        A failing exit code if the candidate misses geometry or stiffness admission.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "route",
        choices=("native", "native-candidate", "microgen-default", "microgen-tuned"),
    )
    parser.add_argument("output", type=Path)
    parser.add_argument("--cells", type=int, help="VTK resolution or native background cells")
    args = parser.parse_args()
    candidate = args.route == "native-candidate"
    if candidate:
        if args.cells is not None:
            parser.error("native-candidate uses the fixed 24-cell refinement control")
        args.cells = 24
    if args.cells is not None and args.cells < 2:
        parser.error("--cells must be at least 2")
    if args.route != "microgen-default" and args.cells is None:
        parser.error("Specify --cells for this route")
    args.output.mkdir()
    hardware = hardware_snapshot()
    row = {
        "recorded_at": dt.datetime.now(dt.UTC).isoformat(),
        "route": args.route,
        "cells": args.cells,
        "platform": hardware["platform"],
        "hardware": hardware,
        "reference_sha256": digest(REFERENCE),
        "versions": {
            pkg: version(pkg)
            for pkg in PACKAGES
            if not candidate or pkg not in {"mmgpy", "microgen"}
        },
    }
    start = time.perf_counter()
    try:
        settings = {"max_length": 0.032, "adapt_tolerance": 0.0003} if candidate else None
        run(SchwarzSettings(args.route, args.cells, native_sizing=settings), args.output, row)
    except Exception as exc:
        row.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        row["workflow_seconds"] = time.perf_counter() - start
        (args.output / "result.json").write_text(
            json.dumps(row, indent=2, allow_nan=False) + "\n",
            encoding="utf8",
        )

    return int(candidate and row["status"] != "within-reference-target")


if __name__ == "__main__":
    raise SystemExit(main())
