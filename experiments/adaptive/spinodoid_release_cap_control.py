"""Bounded two-vertex feasibility control on the remaining frozen cap defect.

Diagnostic only. Run inside macbook_stage.py; no repaired mesh is promoted.
"""
import argparse
import json
from pathlib import Path

import microgen_workflow as w


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args()
    a.run.mkdir(parents=True, exist_ok=False)
    np, _, meshers, *_ = w.imports()
    import comparison_fixtures as fixtures
    import scorecard
    source = w.HERE/"runs/macbook-spinodoid-release-isotropic-capture-v1/capture/spinodoid/02-optimized.npz"
    arrays = np.load(source)
    points, cells, surface = [arrays[k] for k in ("points", "tetrahedra", "surface")]
    compiled = meshers.compile_field(fixtures.case("spinodoid", w.HERE/"runs/macbook-comparison-fixtures-v1").function)
    coordinates = points.copy()
    coordinates[coordinates == .5] = -.5
    _, classes = np.unique(np.round(coordinates*1e11).astype(np.int64), axis=0, return_inverse=True)
    mask = np.zeros(len(points), dtype=np.uint8)
    for face in surface:
        xyz = points[face]
        label = 1
        for axis in range(3):
            if (xyz[:, axis] == -.5).all() or (xyz[:, axis] == .5).all():
                label = 1 << (2+2*axis)
        mask[face] |= label
    vertex = int(np.argmin(np.linalg.norm(points-[-.41158049093388716, -.5, .21607614855859372], axis=1)))
    plane = int(mask[vertex]) & ~1
    assert plane == 16
    incident = np.flatnonzero(np.isin(classes[cells], classes[vertex]).any(axis=1))
    neighbors = np.unique(cells[incident])
    neighbors = neighbors[np.isin(mask[neighbors], [plane, plane | 1]) & (classes[neighbors] != classes[vertex])]
    neighbors = sorted(set(map(int, neighbors)), key=lambda i:float(np.linalg.norm(coordinates[i]-coordinates[vertex])))
    neighbors = list({int(classes[i]):i for i in reversed(neighbors)}.values())

    def project(point):
        q = point.copy()
        for _ in range(30):
            residual = float(compiled(*q))
            gradient = np.array(compiled.gradient(*q));gradient[1] = 0
            norm2 = float(gradient@gradient)
            if abs(residual)/max(norm2**.5, 1e-300) <= 1e-13:
                return q
            q -= np.clip(residual/norm2, -.05/max(norm2**.5, 1e-300), .05/max(norm2**.5, 1e-300))*gradient
        raise AssertionError("Projection did not converge")

    target = project(coordinates[vertex])
    trials = []
    for neighbor in neighbors:
        support = np.flatnonzero(np.isin(classes[cells], [classes[vertex], classes[neighbor]]).any(axis=1))
        xyz = points[cells[support]].copy()
        vmask = classes[cells[support]] == classes[vertex]
        wmask = classes[cells[support]] == classes[neighbor]
        origin = coordinates[neighbor]
        incident_w = cells[np.isin(classes[cells], classes[neighbor]).any(axis=1)]
        lengths = [np.linalg.norm(points[incident_w[:, i]]-points[incident_w[:, j]], axis=1).max()
            for i in range(4) for j in range(i+1, 4)]
        radius = .25*max(lengths)
        best = dict(neighbor=neighbor, wall=bool(mask[neighbor]&1), minimum_quality=0., feasible=False)
        for halving in range(7):
            for dx, dz in [(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)]:
                candidate = origin+radius*2**(-halving)*np.array([dx,0,dz])
                if mask[neighbor]&1:
                    candidate = project(candidate)
                elif float(compiled(*candidate)) > 0:
                    continue
                if abs(candidate).max() > .5:
                    continue
                moved = xyz+vmask[:,:,None]*(target-coordinates[vertex])+wmask[:,:,None]*(candidate-origin)
                flat = moved.reshape(-1,3)
                quality, _, volume = scorecard.tet_metrics(flat, np.arange(len(flat)).reshape(-1,4))
                minimum = float(quality.min()) if volume.min() > 0 else 0.
                if minimum > best["minimum_quality"]:
                    best.update(minimum_quality=minimum, feasible=minimum>=.1, target=target.tolist(),
                        neighbor_target=candidate.tolist(), halving=halving, direction=[dx,dz])
        trials.append(best)
    w.save(a.run/"result.json", dict(complete=True, source=str(source), source_sha256=w.digest(source),
        vertex=vertex, target=target.tolist(), trials=trials, native_admitted=False,
        scope="Feasible two-vertex witnesses only; exact field and held cap plane; quality checked over all incident periodic copies; no timing or final geometry certification"))
    print("Cap-only best", max(t["minimum_quality"] for t in trials if not t["wall"]),
        "curve best", max(t["minimum_quality"] for t in trials if t["wall"]), flush=True)


if __name__ == "__main__":
    main()
