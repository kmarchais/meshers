"""One-process TPMS scaling measurement. Run once per size for peak memory."""

import argparse
import ctypes
import json
import os
import time

import meshers
import numpy as np


def split_p(x, y, z):
    x, y, z = (2 * np.pi * v for v in (x, y, z))
    return (
        1.1
        * (
            np.sin(2 * x) * np.cos(y) * np.sin(z)
            + np.sin(2 * y) * np.cos(z) * np.sin(x)
            + np.sin(2 * z) * np.cos(x) * np.sin(y)
        )
        - 0.2
        * (
            np.cos(2 * x) * np.cos(2 * y)
            + np.cos(2 * y) * np.cos(2 * z)
            + np.cos(2 * z) * np.cos(2 * x)
        )
        - 0.4 * (np.cos(2 * x) + np.cos(2 * y) + np.cos(2 * z))
    )


def gyroid(x, y, z):
    x, y, z = (2 * np.pi * v for v in (x, y, z))
    return np.sin(x) * np.cos(y) + np.sin(y) * np.cos(z) + np.sin(z) * np.cos(x)


def peak_megabytes():
    class MemoryCounters(ctypes.Structure):
        _fields_ = [
            ("size", ctypes.c_ulong),
            ("page_fault_count", ctypes.c_ulong),
            ("peak_working_set", ctypes.c_size_t),
            ("working_set", ctypes.c_size_t),
            ("quota_peak_paged_pool", ctypes.c_size_t),
            ("quota_paged_pool", ctypes.c_size_t),
            ("quota_peak_nonpaged_pool", ctypes.c_size_t),
            ("quota_nonpaged_pool", ctypes.c_size_t),
            ("pagefile", ctypes.c_size_t),
            ("peak_pagefile", ctypes.c_size_t),
        ]

    counters = MemoryCounters()
    counters.size = ctypes.sizeof(counters)
    ctypes.windll.kernel32.GetCurrentProcess.restype = ctypes.c_void_p
    ctypes.windll.psapi.GetProcessMemoryInfo.argtypes = (
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_ulong,
    )
    process = ctypes.windll.kernel32.GetCurrentProcess()
    if not ctypes.windll.psapi.GetProcessMemoryInfo(
        process, ctypes.byref(counters), counters.size
    ):
        raise OSError("GetProcessMemoryInfo failed")
    return round(counters.peak_working_set / 2**20, 1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=("surface", "volume"))
    parser.add_argument("geometry", choices=("gyroid", "split_p"))
    parser.add_argument("--unit-cells", type=int, default=1)
    parser.add_argument("--resolution", type=int, default=16)
    parser.add_argument("--passes", type=int, default=10)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--tile", action="store_true")
    args = parser.parse_args()
    extent = (1 if args.tile else args.unit_cells) / 2
    bounds = (-extent, extent) * 3
    cells = args.resolution * (1 if args.tile else args.unit_cells)
    field = meshers.compile_field(globals()[args.geometry])
    start = time.perf_counter()
    if args.kind == "surface":
        result = meshers.generate_surface(
            field,
            bounds=bounds,
            cells=cells,
            band=(-0.25, 0.25),
            periodic=(True, True, True),
            polish_passes=args.passes,
        )
        output = {
            "points": len(result.points),
            "elements": len(result.triangles),
            "minimum_angle": result.diagnostics["minimum_angle_degrees"],
            "actual_cells": result.diagnostics["background_cells"],
            "native_seconds": result.diagnostics["seconds"],
        }
    else:
        result = meshers.generate(
            field,
            bounds=bounds,
            cells=cells,
            band=(-0.25, 0.25),
            periodic=(True, True, True),
            geometry_tolerance=0.01,
            minimum_quality=0,
            max_tetrahedra=10_000_000,
            optimize_passes=args.passes,
            threads=args.threads,
        )
        output = {
            "points": len(result.points),
            "elements": len(result.tetrahedra),
            "diagnostics": result.diagnostics,
        }
    generation_seconds = time.perf_counter() - start
    if args.tile:
        result = meshers.tile_periodic(
            result, bounds=bounds, repeats=(args.unit_cells,) * 3
        )
        output["points"] = len(result.points)
        output["elements"] = (
            len(result.triangles) if args.kind == "surface" else len(result.tetrahedra)
        )
        output["tiling_seconds"] = round(
            time.perf_counter() - start - generation_seconds, 3
        )
        if args.kind == "volume":
            output["diagnostics"] = result.diagnostics
    output.update(
        kind=args.kind,
        geometry=args.geometry,
        unit_cells=args.unit_cells,
        resolution=args.resolution,
        requested_cells=args.resolution * args.unit_cells,
        tiled=args.tile,
        generation_seconds=round(generation_seconds, 3),
        passes=args.passes,
        threads=args.threads,
        total_seconds=round(time.perf_counter() - start, 3),
        peak_megabytes=peak_megabytes(),
        pid=os.getpid(),
    )
    print(json.dumps(output, default=str), flush=True)


if __name__ == "__main__":
    main()
