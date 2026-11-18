"""Replay native Tet4 cases against frozen Microgen-study records on a release wheel."""

import argparse
import datetime as dt
import hashlib
import json
import sys
import time
from dataclasses import replace
from importlib import import_module
from importlib.metadata import version
from pathlib import Path

import meshers
import numpy as np
from fedoo_homogenization import element_geometry, homogenize
from scipy.linalg import eigvalsh

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from benchmark_hardware import hardware_snapshot

BASELINE = Path(__file__).with_name("release-wheel") / "baseline.json"


def generate(case, cfg, limits):
    """Generate the frozen periodic specimen in its original coordinates.

    Returns:
        Native linear tetrahedral mesh.
    """

    def gyroid(x, y, z):
        a, b, c = 2 * np.pi * x, 2 * np.pi * y, 2 * np.pi * z
        return np.sin(a) * np.cos(b) + np.sin(b) * np.cos(c) + np.sin(c) * np.cos(a)

    def schwarz(x, y, z):
        return (
            np.cos(2 * np.pi * (x + cfg.get("phase_x", 0)))
            + np.cos(2 * np.pi * y)
            + np.cos(2 * np.pi * z)
        )

    return meshers.generate(
        gyroid if case.startswith("gyroid") else schwarz,
        band=(-0.55, 0.55) if case.startswith("gyroid") else None,
        bounds=(-0.5, 0.5) * 3,
        cells=cfg["cells"],
        adapt=cfg["adapt"],
        order=cfg.get("order", 1),
        periodic=(True,) * 3,
        threads=1,
        optimize_passes=cfg.get("optimization_passes", 1),
        geometry_tolerance=cfg["tolerance"],
        minimum_quality=limits["minimum_quality"],
        min_length=cfg["min_length"],
        max_length=cfg.get("max_length", 0.1),
        gradation=cfg.get("gradation", 0.3),
        adapt_iterations=4,
        adapt_tolerance=cfg.get("adapt_tolerance", 0.002),
        max_tetrahedra=limits["max_nodes_or_elements"],
    )


def solve_case(case, baseline, row):
    """Generate and solve one case, updating the record after each stage.

    Raises:
        RuntimeError: A mesh misses the count or volume limit.
    """
    route = baseline["routes"][case]
    limits = baseline["limits"]
    start = time.perf_counter()
    mesh = generate(case, route["config"], limits)
    row["generation_seconds"] = time.perf_counter() - start
    volume = float(element_geometry(mesh)[0].sum())
    key = case.split("_")[0]
    relative = volume / baseline["reference_volumes"][key] - 1
    row.update(
        nodes=len(mesh.points),
        elements=len(mesh.tetrahedra),
        volume=volume,
        relative_volume_error=relative,
    )
    if (
        max(row["nodes"], row["elements"]) > limits["max_nodes_or_elements"]
        or abs(relative) > limits["relative_volume_error"]
    ):
        raise RuntimeError("Count or volume outside the frozen study limits")
    # The existing solver verifies all three translated caps in [0,1]^3.
    shifted = replace(mesh, points=mesh.points + 0.5)
    result = homogenize(shifted)
    row["fea"] = result
    if "archived_stiffness" in route:
        ref = np.asarray(route["archived_stiffness"])
        tensor = np.asarray(result["stiffness_engineering_voigt"])
        diff = tensor - ref
        row["worst_relative_energy_change_from_archive"] = float(
            np.abs(eigvalsh((diff + diff.T) / 2, (ref + ref.T) / 2)).max()
        )
    row["status"] = "solved"


def main():
    """Record one fresh generation and six-load study, preserving rejected cases.

    Raises:
        RuntimeError: The installed version differs from the pinned release.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=["gyroid_uniform", "gyroid_adaptive", "schwarz_native"])
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    if version("meshers") != baseline["release"]:
        raise RuntimeError("Install the specified release wheel in an isolated environment")
    route = baseline["routes"][args.case]
    extension = Path(import_module("meshers._meshers").__file__)
    hardware = hardware_snapshot()
    row = {
        "recorded_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "case": args.case,
        "config": route["config"],
        "versions": {name: version(name) for name in ("meshers", "fedoo", "numpy", "scipy")},
        "baseline_sha256": hashlib.sha256(BASELINE.read_bytes()).hexdigest(),
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "platform": hardware["platform"],
        "hardware": hardware,
        "release_wheel_comparative_admission": "pending independent two-sided distance audit",
    }
    with args.output.open("x", encoding="utf-8") as output:
        try:
            solve_case(args.case, baseline, row)
        except Exception as exc:
            row.update(status="failed", error=str(exc))
            raise
        finally:
            json.dump(row, output, indent=2, allow_nan=False)
            output.write("\n")


if __name__ == "__main__":
    main()
