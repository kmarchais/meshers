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


def balance(mesh):
    """Move wall vertices outward by half the mean centroid sag of their wall faces.

    Box coordinates stay fixed and periodic copies share one offset, so the cell stays
    periodic. Linear boundary chords then deviate to both sides of the wall.
    """
    points = np.array(mesh.points)
    surface = np.asarray(mesh.surface)
    corners = points[surface]
    cap = np.zeros(len(surface), dtype=bool)
    for a in range(3):
        for side in (-0.5, 0.5):
            cap |= np.all(np.isclose(corners[:, :, a], side, atol=1e-12), axis=1)
    wall = surface[~cap]
    centroid = points[wall].mean(axis=1)
    value = scorecard.gyroid(*centroid.T)
    level = np.where(value > 0, 1.0, -1.0)
    grad = scorecard.gradient(scorecard.gyroid, centroid)
    # Signed outward distance: outward is +grad on the upper wall, -grad on the lower.
    signed = (value - level) / np.linalg.norm(grad, axis=1) * level
    total = np.zeros(len(points))
    count = np.zeros(len(points))
    np.add.at(total, wall.ravel(), np.repeat(signed, 3))
    np.add.at(count, wall.ravel(), 1)
    # Periodic copies form one class; pool their faces.
    parent = np.arange(len(points))
    for pairs in mesh.periodic_pairs:
        for a, b in pairs:
            ra, rb = parent[a], parent[b]
            while parent[ra] != ra:
                ra = parent[ra]
            while parent[rb] != rb:
                rb = parent[rb]
            parent[max(ra, rb)] = min(ra, rb)
    root = np.array([i for i in range(len(points))])
    for i in range(len(points)):
        r = i
        while parent[r] != r:
            r = parent[r]
        root[i] = r
    class_total = np.zeros(len(points))
    class_count = np.zeros(len(points))
    np.add.at(class_total, root, total)
    np.add.at(class_count, root, count)
    moving = class_count[root] > 0
    offset = np.zeros(len(points))
    offset[moving] = -0.5 * class_total[root][moving] / class_count[root][moving]
    v = points[moving]
    vv = scorecard.gyroid(*v.T)
    lv = np.where(vv > 0, 1.0, -1.0)
    normal = scorecard.gradient(scorecard.gyroid, v) * lv[:, None]
    locked = np.isclose(v, -0.5, atol=1e-12) | np.isclose(v, 0.5, atol=1e-12)
    normal[locked] = 0
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    points[moving] = v + offset[moving][:, None] * normal
    tets = np.asarray(mesh.tetrahedra)
    p = points[tets]
    det = np.einsum("ij,ij->i", p[:, 1] - p[:, 0], np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))
    assert (det > 0).all(), "balancing inverted an element"
    return Shifted(points, tets, mesh.surface, mesh.periodic_pairs, mesh.diagnostics)


def main():
    ref = reference()
    out = HERE / "fea.jsonl"
    for config in sys.argv[1:]:
        start = time.perf_counter()
        balanced = config.startswith("balanced-")
        mesh = generate(config.removeprefix("balanced-"))
        if balanced:
            mesh = balance(mesh)
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
