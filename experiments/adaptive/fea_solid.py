"""Homogenized-stiffness error per DOF for a thick periodic solid: Schwarz P network.

Solid: cos(2πx) + cos(2πy) + cos(2πz) < 0 in [-0.5, 0.5]^3, about 50% solid with thick
ligaments, E=2.5, ν=0.25. Linear meshes use the audited fedoo homogenization;
quadratic meshes come from the native `order=2` output (fedoo tet10). The reference
is the p=4 extrapolation of the quadratic 16- and 20-cell lattice tensors.
Run: uv run --no-sync python experiments/adaptive/fea_solid.py <config> [...]
Configs: lattice:cells=N[,order=2] or adapt:tol=T[,hmax=H,cells=N,order=2]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import fea  # noqa: E402  (frozen build path, periodic-fea modules)
import meshers  # noqa: E402
import p2  # noqa: E402
from homogenization import homogenize  # noqa: E402
from periodic import tensor_difference  # noqa: E402

TAU = 2 * np.pi
OUT = HERE / "fea_solid.jsonl"


def schwarz_p(x, y, z):
    return np.cos(TAU * x) + np.cos(TAU * y) + np.cos(TAU * z)


def generate(config):
    kind, _, rest = config.partition(":")
    args = dict(item.split("=") for item in rest.split(",") if item)
    order = int(args.get("order", 1))
    common = dict(
        bounds=(-0.5, 0.5) * 3,
        periodic=(True, True, True),
        optimize_passes=1,
        threads=1,
        order=order,
    )
    if kind == "lattice":
        return order, meshers.generate(
            schwarz_p, cells=int(args["cells"]), geometry_tolerance=1.0, **common
        )
    return order, meshers.generate(
        schwarz_p,
        cells=int(args.get("cells", 16)),
        geometry_tolerance=float(args["tol"]),
        adapt=True,
        max_length=float(args.get("hmax", 0.2)),
        **common,
    )


def stiffness(order, mesh):
    if order == 2:
        c = p2.homogenize_p2(mesh.points, mesh.tetrahedra)
        return (c + c.T) / 2
    shifted = fea.Shifted(
        mesh.points + 0.5, mesh.tetrahedra, mesh.surface, mesh.periodic_pairs, mesh.diagnostics
    )
    return np.asarray(homogenize(shifted)["stiffness_engineering_voigt"])


def reference():
    rows = {r["config"]: np.array(r["stiffness"]) for r in map(json.loads, OUT.open())}
    c16, c20 = rows["lattice:cells=16,order=2"], rows["lattice:cells=20,order=2"]
    return (20**4 * c20 - 16**4 * c16) / (20**4 - 16**4)


def main():
    for config in sys.argv[1:]:
        if config == "report":
            ref = reference()
            for r in map(json.loads, OUT.open()):
                e = tensor_difference(np.array(r["stiffness"]), ref)
                print(
                    f"{r['config']:34s} nodes={r['nodes']:7d} elements={r['elements']:7d} "
                    f"err={e['worst_relative_strain_energy_difference']:.4%}"
                )
            continue
        start = time.perf_counter()
        order, mesh = generate(config)
        generation = time.perf_counter() - start
        start = time.perf_counter()
        c = stiffness(order, mesh)
        row = {
            "config": config,
            "nodes": len(mesh.points),
            "elements": len(mesh.tetrahedra),
            "stiffness": c.tolist(),
            "generation_s_unreliable": generation,
            "fea_s_unreliable": time.perf_counter() - start,
        }
        with OUT.open("a") as log:
            log.write(json.dumps(row) + "\n")
        print(config, row["nodes"], row["elements"], f"{row['fea_s_unreliable']:.0f}s", flush=True)


if __name__ == "__main__":
    main()
