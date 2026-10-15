"""Interleaved wall-clock timing at matched, independently measured geometry error.

Meant for the pueue `gpu` group (one job at a time): timings taken next to other
heavy work are not valid. Each configuration is warmed up once in this process, then
measured `--repeats` times in rotation with the others. Mesh counts and the
independent Hausdorff distance are recorded once per configuration. The FEA block
times mesh generation plus the six-load homogenization for a 1% target.

Run: uv run --no-sync python experiments/adaptive/timing.py [--threads N] [--fea]
`--threads` also caps the surface pipeline's rayon pool. The pip VTK wheel runs its
SMP filters sequentially; the backend is recorded with the results.
"""

import argparse
import json
import os
import platform
import statistics
import sys
import time
from pathlib import Path

if "--threads" in sys.argv:  # must precede the meshers import that starts rayon
    os.environ["RAYON_NUM_THREADS"] = sys.argv[sys.argv.index("--threads") + 1]

import numpy as np

sys.path.insert(0, str(Path.home() / "bin"))
from pueue_progress import report, result  # no-ops outside pueue

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# fea puts the frozen extension named by MESHERS_BUILD first on sys.path; timing_job.sh
# builds and freezes HEAD so the timed code is the committed code.
import fea  # noqa: E402  (frozen build path, periodic-fea modules)
import meshers  # noqa: E402
import scorecard  # noqa: E402

BOX = (-0.5, 0.5) * 3


# Flying Edges grids bracketing H = 1e-3 (independent scorecard, see RESULTS.md).
FLYING_EDGES = {"sphere": (28, 32), "ellipsoid": (64, 80), "gyroid": (40, 48), "schwarz-p": (36, 40)}
MODES = {  # label: (balanced, fast, iterations)
    "adaptive": (False, False, 5),
    "balanced": (True, False, 5),
    "fast": (False, True, 2),
    "fast-balanced": (True, True, 2),
}


def surface_configs(threads):
    out = {}
    for name, grids in FLYING_EDGES.items():
        case = scorecard.CASES[name]
        for n in grids:
            out[f"{name}/flying-edges-{n}"] = (case, lambda c=case, n=n: scorecard.flying_edges(c, n))
        for label, mode in MODES.items():
            out[f"{name}/{label}-1e-3"] = (case, lambda c=case, m=mode: _surface(c, 24, *m))
    return out


def _surface(case, cells, balanced, fast, iterations):
    mesh = meshers.generate_surface(
        case.function,
        bounds=BOX,
        cells=cells,
        band=case.band,
        geometry_tolerance=1e-3,
        adaptive=True,
        iterations=iterations,
        balanced=balanced,
        fast=fast,
    )
    return mesh.points, mesh.triangles


def volume_configs(threads):
    case = scorecard.CASES["gyroid"]
    common = dict(bounds=BOX, band=case.band, periodic=(True,) * 3, threads=threads)

    def lattice():
        return meshers.generate(
            case.function, cells=32, geometry_tolerance=1.0, optimize_passes=1, **common
        )

    def adapted():
        return meshers.generate(
            case.function,
            cells=24,
            geometry_tolerance=1e-3,
            optimize_passes=1,
            adapt=True,
            max_length=0.1,
            **common,
        )

    def mmg():
        seed = meshers.generate(
            case.function,
            cells=24,
            geometry_tolerance=1.0,
            optimize_passes=1,
            **{**common, "periodic": (False,) * 3},
        )
        return scorecard.mmg_adapt(case, seed, 1e-3, 0.1)

    return {
        "gyroid-volume/lattice-32": (case, lattice),
        "gyroid-volume/adapt-1e-3": (case, adapted),
        "gyroid-volume/mmg3d-1e-3-nonperiodic": (case, mmg),
    }


class Progress:
    """Fraction over the whole job: `--part i/n` places this process's share."""

    def __init__(self, part, total):
        index, parts = part
        self.base, self.span, self.total = index / parts, 1 / parts, total
        self.done, self.start = 0, time.perf_counter()

    def step(self, message):
        self.done += 1
        rate = (time.perf_counter() - self.start) / self.done
        fraction = self.base + self.span * self.done / self.total
        report(fraction, message, eta_s=rate * (self.total - self.done))


def measure(configs, repeats, progress, tag):
    rows = {}
    for label, (case, make) in configs.items():
        mesh = make()  # warm-up and geometry record
        if isinstance(mesh, tuple):
            row = scorecard.score(case, *mesh)
        else:
            row = scorecard.score_volume(case, mesh)
        rows[label] = {"metrics": row, "seconds": []}
        progress.step(f"{tag}: scored {label}")
    for repeat in range(repeats):
        for label, (_, make) in configs.items():
            start = time.perf_counter()
            make()
            rows[label]["seconds"].append(time.perf_counter() - start)
        progress.step(f"{tag}: repeat {repeat + 1}/{repeats}")
    for row in rows.values():
        row["median_s"] = statistics.median(row["seconds"])
    return rows


def fea_block(threads):
    """Mesh + homogenization to a 1% stiffness target: linear vs curved quadratic."""
    import p2
    from homogenization import homogenize

    out = {}
    start = time.perf_counter()
    linear = fea.generate("lattice:cells=48")
    generation = time.perf_counter() - start
    shifted = fea.Shifted(
        linear.points + 0.5,
        linear.tetrahedra,
        linear.surface,
        linear.periodic_pairs,
        linear.diagnostics,
    )
    result = homogenize(shifted)
    error = fea.tensor_difference(np.asarray(result["stiffness_engineering_voigt"]), fea.reference())
    out["linear-48"] = {
        "nodes": len(linear.points),
        "elements": len(linear.tetrahedra),
        "generation_s": generation,
        "fea_s": result["seconds"],
        "error": error["worst_relative_strain_energy_difference"],
    }
    start = time.perf_counter()
    quadratic = meshers.generate(
        scorecard.gyroid,
        bounds=BOX,
        band=(-1.0, 1.0),
        periodic=(True,) * 3,
        cells=9,
        geometry_tolerance=1.0,
        optimize_passes=1,
        threads=threads,
        order=2,
    )
    generation = time.perf_counter() - start
    start = time.perf_counter()
    stiffness = p2.homogenize_p2(quadratic.points, quadratic.tetrahedra)
    solve = time.perf_counter() - start
    error = fea.tensor_difference((stiffness + stiffness.T) / 2, fea.reference())
    out["quadratic-9"] = {
        "nodes": len(quadratic.points),
        "elements": len(quadratic.tetrahedra),
        "generation_s": generation,
        "fea_s": solve,
        "error": error["worst_relative_strain_energy_difference"],
    }
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--fea", action="store_true")
    parser.add_argument("--part", default="0/1", help="i/n: this run's share of the progress bar")
    args = parser.parse_args()
    part = tuple(int(v) for v in args.part.split("/"))
    surface, volume = surface_configs(args.threads), volume_configs(args.threads)
    progress = Progress(part, len(surface) + len(volume) + 2 * args.repeats)
    report(part[0] / part[1], f"threads {args.threads}: starting")
    from vtkmodules.vtkCommonCore import vtkSMPTools

    record = {
        "vtk_smp_backend": vtkSMPTools.GetBackend(),
        "vtk_smp_threads": vtkSMPTools.GetEstimatedNumberOfThreads(),
        "meshers": meshers.__file__,
        "threads": args.threads,
        "repeats": args.repeats,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "surface": measure(surface, args.repeats, progress, f"threads {args.threads} surface"),
        "volume": measure(volume, args.repeats, progress, f"threads {args.threads} volume"),
    }
    if args.fea:
        record["fea"] = fea_block(args.threads)
    out = HERE / f"timing-threads{args.threads}.json"
    out.write_text(json.dumps(record, indent=1), encoding="utf-8")
    for block in ("surface", "volume"):
        for label, row in record[block].items():
            m = row["metrics"]
            count = m.get("tetrahedra", m["triangles"])
            result(
                f"t{args.threads} {label}",
                elements=count,
                hausdorff=m["hausdorff"],
                median={"value": row["median_s"], "unit": "s"},
                jsonl=HERE / "timing-results.jsonl",
            )
    for label, row in record.get("fea", {}).items():
        print(label, row)


if __name__ == "__main__":
    main()
