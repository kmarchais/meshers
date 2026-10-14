"""Homogenized-stiffness error per DOF for the periodic Gyroid sheet (raw band +-0.55).

Reuses the audited fedoo homogenization of the periodic-fea experiment (P1 tet4,
six strain loads, single-thread CG). The reference is the quadratic extrapolation of
the archived native64 and native80 tensors; measured convergence order there is ~2.1.
Run: uv run --no-sync python experiments/adaptive/fea.py <config> [...]
"""

import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PERIODIC = HERE.parents[2] / "periodic-microgen/experiments/periodic-fea"
sys.path.insert(0, str(PERIODIC))
sys.path.insert(0, str(HERE))
# Experiments import a frozen copy of the extension so rebuilds cannot collide.
BUILD = HERE / "builds" / __import__("os").environ.get("MESHERS_BUILD", "dfd0e339")
if BUILD.exists():
    sys.path.insert(0, str(BUILD))

import meshers  # noqa: E402
import scorecard  # noqa: E402
from homogenization import homogenize  # noqa: E402
from periodic import tensor_difference  # noqa: E402

RUNS = PERIODIC / "runs"
ARCHIVED = {
    64: RUNS / "reference-and-seed-v1/02-periodic-gyroid-meshers-pass1-64/result.json",
    80: RUNS / "reference-80-v1/00-periodic-gyroid-meshers-pass1-80/result.json",
}


def archived(cells):
    data = json.loads(ARCHIVED[cells].read_text())
    return np.asarray(data.get("homogenization", data)["stiffness_engineering_voigt"])


def reference():
    c64, c80 = archived(64), archived(80)
    return (80**2 * c80 - 64**2 * c64) / (80**2 - 64**2)


@dataclass(frozen=True)
class Shifted:
    """Mesh in the [0,1]^3 frame expected by the homogenization audit."""

    points: np.ndarray
    tetrahedra: np.ndarray
    surface: np.ndarray
    periodic_pairs: tuple
    diagnostics: dict


def generate(config):
    kind, _, rest = config.partition(":")
    args = dict(item.split("=") for item in rest.split(",") if item)
    common = dict(
        bounds=(-0.5, 0.5) * 3,
        band=(-1.0, 1.0),
        periodic=(True, True, True),
        optimize_passes=int(args.get("passes", 1)),
        threads=1,
    )
    if kind == "lattice":
        return meshers.generate(
            scorecard.gyroid, cells=int(args["cells"]), geometry_tolerance=1.0, **common
        )
    if kind == "adapt":
        return meshers.generate(
            scorecard.gyroid,
            cells=int(args.get("cells", 24)),
            geometry_tolerance=float(args["tol"]),
            adapt=True,
            max_length=float(args.get("hmax", 0.1)),
            layers=float(args.get("layers", 0)),
            gradation=float(args.get("hgrad", 0.3)),
            **common,
        )
    raise ValueError(kind)


def main():
    ref = reference()
    out = HERE / "fea.jsonl"
    for config in sys.argv[1:]:
        start = time.perf_counter()
        mesh = generate(config)
        generation = time.perf_counter() - start
        shifted = Shifted(
            mesh.points + 0.5,
            mesh.tetrahedra,
            mesh.surface,
            mesh.periodic_pairs,
            mesh.diagnostics,
        )
        result = homogenize(shifted)
        error = tensor_difference(np.asarray(result["stiffness_engineering_voigt"]), ref)
        q, dihedral, _ = scorecard.tet_metrics(mesh.points, mesh.tetrahedra)
        row = {
            "config": config,
            "nodes": len(mesh.points),
            "tetrahedra": len(mesh.tetrahedra),
            "q_min": float(q.min()),
            "dihedral_min_deg": float(dihedral.min()),
            "worst_energy_error": error["worst_relative_strain_energy_difference"],
            "eigenvalues": error["relative_energy_difference_eigenvalues"],
            "solid_volume_fraction": result["solid_volume_fraction"],
            "generation_s_unreliable": generation,
            "fea_s_unreliable": result["seconds"],
        }
        with out.open("a") as log:
            log.write(json.dumps(row) + "\n")
        print(
            f"{config:40s} nodes={row['nodes']:7d} tets={row['tetrahedra']:8d} "
            f"err={row['worst_energy_error']:.4%} qmin={row['q_min']:.3f} "
            f"fea={row['fea_s_unreliable']:.0f}s",
            flush=True,
        )


if __name__ == "__main__":
    main()
