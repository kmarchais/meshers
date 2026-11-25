# Copyright (c) 2026 Kevin Marchais

"""Check the built Meshers wheel against live Microgen F-reps and audit the VTU readback."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import tempfile
from importlib import import_module
from importlib.metadata import distribution, version
from pathlib import Path
from typing import TYPE_CHECKING

import meshers
import numpy as np
import pyvista as pv
from microgen_frep_audit import QUALITY, TOLERANCE, VOLUME_LIMIT, audit
from microgen_frep_cases import CASES, specimen
from microgen_spatial_audit import audit as spatial_audit
from microgen_spatial_audit import missed_checks
from microgen_spatial_cases import BOUNDS, SHIFT
from microgen_spatial_cases import specimen as infill

if TYPE_CHECKING:
    from periodic_bench import Coordinates

LOGGER = logging.getLogger(__name__)
REGRESSIONS = ("repeat", "bored-infill", "graded-cube", "graded-cube-translated")
BUDGET = 250_000
OPTIONS = {
    "geometry_tolerance": TOLERANCE,
    "minimum_quality": QUALITY,
    "max_tetrahedra": BUDGET,
    "threads": 1,
    "optimize_passes": 4,
}


def readback(grid: pv.UnstructuredGrid) -> pv.UnstructuredGrid:
    """Return cells, coordinates and boundary metadata read from the exported VTU."""
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "result.vtu"
        grid.save(path)
        return pv.read(path)


def chart(x: Coordinates) -> Coordinates:
    """Return physical cube coordinates with finer spacing near each cap."""
    return 0.65 * x + 0.175 * np.sin(np.pi * x)


def measure(name: str) -> dict:
    """Return measurements and missed physical limits after a VTU round trip."""
    if name in CASES:
        case = specimen(name)
        options = {
            **OPTIONS,
            "cells": 12,
            "adapt": False,
            "periodic": (case.periodic,) * 3,
            "geometry_tolerance_mode": "absolute",
            **case.options,
        }
        grid = case.shape.generate_volume_mesh(backend="meshers", **options)
        metrics, misses = audit(case, readback(grid))
    else:
        options = {**OPTIONS, "cells": 24, "adapt": True, "geometry_tolerance_mode": "absolute"}
        if name == "bored-infill":
            spatial = infill(name)
            grid = spatial.grid(options)
        else:
            shift = np.asarray(SHIFT)
            if name == "graded-cube-translated":
                shift += np.array([1 / 32, -3 / 64, 1 / 64])
            spatial = infill("distance-infill", tuple(shift))
            # The chart can lower quality on the final physical coordinates.
            options = {**OPTIONS, "cells": 32, "minimum_quality": 0.16, "bounds": BOUNDS}

            def mapped(x: Coordinates, y: Coordinates, z: Coordinates) -> Coordinates:
                return spatial.level(chart(x), chart(y), chart(z))

            mesh = meshers.generate(mapped, periodic=(False,) * 3, **options)
            grid = pv.UnstructuredGrid(
                {pv.CellType.TETRA: mesh.tetrahedra},
                spatial.world(chart(mesh.points)),
            )
            grid.field_data["meshers_boundary_faces"] = mesh.surface.copy()
            grid.field_data["meshers_boundary_tags"] = mesh.boundary_tags.copy()
            options["phase_shift"] = shift.tolist()
        metrics = spatial_audit(spatial, readback(grid))
        misses = missed_checks(metrics)
    if metrics["tetrahedra"] > BUDGET:
        misses.append("tetrahedron_budget")
    return {
        "status": "missed" if misses else "admitted",
        "missed": misses,
        "metrics": metrics,
        "options": options,
    }


def main() -> int:
    """Write one record per attempted specimen, including rejections and audit errors.

    Returns:
        Nonzero for regressions, missing results or a broken audit.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=["regressions", "inventory"], default="regressions")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    names = CASES if args.suite == "inventory" else REGRESSIONS
    run = {
        "suite": args.suite,
        "planned": names,
        "packages": {
            pkg: version(pkg) for pkg in ("microgen", "meshers", "numpy", "pyvista", "scipy", "vtk")
        },
        "microgen_source": json.loads(
            distribution("microgen").read_text("direct_url.json") or "{}",
        ).get("vcs_info"),
        "native_sha256": hashlib.sha256(
            Path(import_module("meshers._meshers").__file__).read_bytes(),
        ).hexdigest(),
        "geometry_tolerance": TOLERANCE,
        "physical_quality_min": QUALITY,
        "relative_volume_limit": VOLUME_LIMIT,
        "tetrahedron_budget": BUDGET,
        "records": [],
    }
    failed = False
    for name in names:
        try:
            record = measure(name)
        except meshers.MeshingError as error:
            record = {"status": "rejected", "reason": str(error)}
        except Exception as error:
            LOGGER.exception("Could not complete %s", name)
            record = {"status": "incomplete", "reason": str(error)}
        run["records"].append({"case": name, **record})
        args.output.write_text(json.dumps(run, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        LOGGER.info(
            "%s: %s %s",
            name,
            record["status"],
            record.get("missed", record.get("reason", "")),
        )
        if record["status"] == "incomplete" or (
            name in REGRESSIONS and record["status"] != "admitted"
        ):
            failed = True
    return int(failed)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())
