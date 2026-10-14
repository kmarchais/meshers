"""Quadratic (tet10) periodic homogenization of the Gyroid sheet with curved walls.

Mid-edge nodes on implicit walls are Newton-projected onto the exact level set, so
the boundary is quadratic; interior and box-plane mid-edge nodes stay straight.
Curved elements are checked for a positive Jacobian at every quadrature point and
offending edges fall back to straight midpoints. Reference as in fea.py.
Run: uv run --no-sync python experiments/adaptive/p2.py <config> [...]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import fea  # noqa: E402  (sets the periodic-fea path)
import fedoo as fd  # noqa: E402
import scorecard  # noqa: E402
from homogenization import AuditedCG  # noqa: E402
from periodic import tensor_difference  # noqa: E402

LO, HI = -0.5, 0.5
EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]  # VTK / fedoo tet10 order


def project(points, locked, levels, iterations=30):
    """Newton projection of points onto wall levels with locked box coordinates."""
    q = points.copy()
    for _ in range(iterations):
        r = scorecard.gyroid(*q.T) - levels
        g = scorecard.gradient(scorecard.gyroid, q)
        g[locked] = 0
        step = (r / np.maximum((g * g).sum(1), 1e-30))[:, None] * g
        q -= step
        if np.abs(step).max() < 1e-14:
            break
    return q


def to_tet10(mesh, curve=True):
    """Add mid-edge nodes; wall edges are curved onto the exact surface."""
    points, tets = np.asarray(mesh.points), np.asarray(mesh.tetrahedra)
    pairs = np.concatenate([tets[:, list(e)] for e in EDGES])
    keys = np.sort(pairs, axis=1)
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    mid = points[unique].mean(axis=1)
    # Wall edges: both ends on a wall face of the boundary.
    surface = np.asarray(mesh.surface)
    corners = points[surface]
    cap = np.zeros(len(surface), dtype=bool)
    for a in range(3):
        for side in (LO, HI):
            cap |= np.all(np.isclose(corners[:, :, a], side, atol=1e-12), axis=1)
    wall = surface[~cap]
    wall_edges = np.unique(
        np.sort(np.concatenate([wall[:, [0, 1]], wall[:, [1, 2]], wall[:, [2, 0]]]), axis=1), axis=0
    )
    is_wall = np.zeros(len(unique), dtype=bool)
    lookup = {tuple(e): i for i, e in enumerate(unique)}
    for e in wall_edges:
        is_wall[lookup[tuple(e)]] = True
    # Box coordinates shared by both ends stay locked; canonicalize hi to lo so that
    # opposite periodic faces get bitwise identical in-plane coordinates.
    ends = points[unique[is_wall]]
    locked = np.isclose(ends[:, 0], ends[:, 1], atol=0) & (
        np.isclose(ends[:, 0], LO, atol=1e-12) | np.isclose(ends[:, 0], HI, atol=1e-12)
    )
    start = mid[is_wall].copy()
    high = locked & np.isclose(start, HI, atol=1e-12)
    start[high] = LO
    levels = np.where(scorecard.gyroid(*start.T) > 0, 1.0, -1.0)
    curved = project(start, locked, levels)
    curved[high] = HI
    if curve:
        mid[is_wall] = curved
    straight = points[unique].mean(axis=1)
    curved_index = np.flatnonzero(is_wall)
    canonical = np.round(np.where(np.isclose(straight, HI, atol=1e-12), LO, straight), 11)
    tet10 = np.concatenate([tets, len(points) + inverse.reshape(6, -1).T], axis=1)
    # Straighten curved edges of invalid elements, together with their periodic copies.
    for _ in range(10):
        nodes = np.concatenate([points, mid])
        bad = ~jacobian_ok(nodes, tet10)
        if not bad.any():
            break
        edges = np.unique(tet10[bad, 4:].ravel() - len(points))
        keys = {tuple(k) for k in canonical[edges]}
        flat = np.array([tuple(k) in keys for k in canonical[curved_index]])
        mid[curved_index[flat]] = straight[curved_index[flat]]
    nodes = np.concatenate([points, mid])
    return nodes, tet10, int(is_wall.sum())


def jacobian_ok(nodes, tet10):
    """Positive Jacobian at fedoo's 15 Gauss points and the vertices."""
    element = fd.lib_elements.element_list.get_element("tet10")(15)
    xi = np.concatenate([element.xi_pg, element.xi_nd[:4]])
    ok = np.ones(len(tet10), dtype=bool)
    x = nodes[tet10]
    for derivative in element.shape_function_derivative(xi):
        jac = np.einsum("an,tnb->tab", derivative, x)
        ok &= np.linalg.det(jac) > 0
    return ok


def homogenize_p2(nodes, tet10):
    fd.ModelingSpace("3D")
    mesh = fd.Mesh(nodes + 0.5, tet10.copy(), "tet10")
    assembly = fd.Assembly.create(
        fd.weakform.StressEquilibrium(fd.constitutivelaw.ElasticIsotrop(2.5, 0.25)), mesh
    )
    solver = AuditedCG()
    stiffness = fd.homogen.get_homogenized_stiffness(
        assembly, meshperio=True, solver=solver, rigid_body_constraint="pin"
    )
    return np.asarray(stiffness)


def main():
    ref = fea.reference()
    out = HERE / "p2.jsonl"
    for config in sys.argv[1:]:
        straight = config.startswith("straight-")
        mesh = fea.generate(config.removeprefix("straight-"))
        nodes, tet10, curved = to_tet10(mesh, curve=not straight)
        ok = jacobian_ok(nodes, tet10)
        start = time.perf_counter()
        stiffness = homogenize_p2(nodes, tet10)
        seconds = time.perf_counter() - start
        error = tensor_difference((stiffness + stiffness.T) / 2, ref)
        row = {
            "config": config,
            "order": 2,
            "vertices": len(mesh.points),
            "nodes": len(nodes),
            "tetrahedra": len(tet10),
            "curved_edges": curved,
            "invalid_jacobian_elements": int((~ok).sum()),
            "worst_energy_error": error["worst_relative_strain_energy_difference"],
            "eigenvalues": error["relative_energy_difference_eigenvalues"],
            "stiffness": stiffness.tolist(),
            "fea_s_unreliable": seconds,
        }
        with out.open("a") as log:
            log.write(json.dumps(row) + "\n")
        print(
            f"P2 {config:36s} nodes={row['nodes']:7d} tets={row['tetrahedra']:7d} "
            f"err={row['worst_energy_error']:.4%} bad_jac={row['invalid_jacobian_elements']} "
            f"fea={seconds:.0f}s",
            flush=True,
        )


if __name__ == "__main__":
    main()
