"""Local Jacobian feasibility of moving saved PMY midpoint witnesses to the wall.

Each trial starts from the same frozen mesh and moves one complete periodic
midpoint class. It is a diagnosis, not a new mesh or an admission decision.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys

import microgen_workflow as w


def worker(a):
    import numpy as np
    from fedoo.lib_elements.tetrahedron import Tet10
    from quadrature_audit import duffy_rule
    from periodic_trace import periodic_pairs
    from pmy_neck_control import raw_gradient
    import tpms_family_fixtures as fixtures

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    source = json.loads((a.source/"result.json").read_text())
    observed = json.loads((a.observer/"result.json").read_text())
    assert source["mesh_sha256"] == observed["mesh_sha256"] == w.digest(a.source/"mesh.npz")
    case = fixtures.case("tpms-pmy-sheet", a.fixtures)
    with np.load(a.source/"mesh.npz") as archive:
        points, cells, surface = (archive[key] for key in ("points", "tetrahedra", "surface"))
    pairs, _, _ = periodic_pairs(points, surface, declared_bounds=([-.5]*3, [.5]*3))
    parent = np.arange(len(points))

    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i

    for axis in pairs:
        for i, j in axis:
            i, j = root(i), root(j)
            parent[max(i, j)] = min(i, j)
    roots = np.array([root(i) for i in range(len(points))])
    endpoints = {}
    for tet in cells:
        for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))):
            pair = tuple(sorted((int(tet[i]), int(tet[j]))))
            old = endpoints.setdefault(int(tet[k+4]), pair)
            assert old == pair
    corners = np.unique(cells[:, :4])
    element = Tet10(15)
    native_bary = np.vstack((np.eye(4),
        [(np.eye(4)[i]+np.eye(4)[j])/2 for i in range(4) for j in range(i+1, 4)],
        [(1-np.eye(4)[i])/3 for i in range(4)], [[.25]*4],
        [[b if j != i else c for j in range(4)] for c, b in ((.7, .1), (.1, .3)) for i in range(4)]))
    # Symmetric native set is invariant under barycentric convention changes.
    probes = np.vstack((native_bary[:, :3], duffy_rule(5)[0]))
    derivatives = element.shape_function_derivative(probes)

    def measures(x):
        det = np.linalg.det(np.einsum("san,tnb->tsab", derivatives, x))
        low, high = det.min(1), det.max(1)
        ratio = low/high
        return dict(minimum_determinant=float(low.min()), minimum_ratio=float(ratio.min()),
                    nonpositive_elements=int((low <= 0).sum()),
                    ratio_failures=int((ratio <= .001).sum()),
                    sampled_valid=bool((low > 0).all() and (ratio > .001).all()))

    linear = np.vstack((element.xi_nd[:4], element.xi_nd[4:]))
    control = measures(linear[None])
    assert abs(control["minimum_determinant"]-1) < 1e-13 and control["sampled_valid"]
    selected = [row for row in observed["witnesses"] if row["kind"] == "edge-midpoint"]
    assert len(selected) <= 8
    result = dict(complete=False, control=control, rows=[], no_meshing=True,
        no_admission_change=True, source_mesh_sha256=source["mesh_sha256"],
        scope="Independent one-class trials from frozen mesh. Saved witnesses guide diagnosis only. No global candidate, boundary guarantee, closest-point certificate or stiffness claim.")
    with Job(parts={str(i): 1 for i in range(len(selected))}) as job:
        for index, witness in enumerate(selected):
            with job.part(str(index)) as part:
                node = witness["mesh_node"]
                members = np.flatnonzero(roots == root(node))
                assert not np.isin(members, corners).any()
                incident = np.flatnonzero(np.isin(cells[:, 4:], members).any(1))
                xyz = points[cells[incident]].copy()
                start = points[node]
                ends = points[list(endpoints[node])]
                midpoint = ends.mean(0)
                corner_values = case.function(*ends.T)
                corner_levels = np.asarray(case.levels)[np.argmin(abs(corner_values[:, None]-case.levels), axis=1)]
                target = float(corner_levels[0])
                row = dict(node=node, periodic_members=members.tolist(), incident_elements=incident.tolist(),
                    endpoint_nodes=list(endpoints[node]), endpoint_field_values=corner_values.tolist(),
                    endpoint_levels=corner_levels.tolist(), same_endpoint_level=bool((corner_levels == target).all()),
                    straight_midpoint=midpoint.tolist(), current=start.tolist(),
                    current_field_value=float(case.function(*start)),
                    current_displacement_from_chord=float(np.linalg.norm(start-midpoint)),
                    edge_length=float(np.linalg.norm(ends[1]-ends[0])), baseline=measures(xyz), trials=[])
                targets = {"saved_standard_foot": witness["standard_foot"]}
                feasible = witness["closest_found_feasible"]
                if feasible is not None and feasible["level"] == target:
                    targets["saved_feasible_foot"] = feasible["foot"]
                # Newton part of the native projection, without claiming an exact
                # replay of Rust arithmetic or its directional fallback.
                projected = midpoint.copy()
                locked = np.isclose(abs(midpoint), .5, atol=1e-12, rtol=0)
                converged = False
                for _ in range(30):
                    residual = float(case.function(*projected)-target)
                    gradient = raw_gradient(projected[None])[0]
                    gradient[locked] = 0
                    length = float(np.linalg.norm(gradient))
                    if length <= 1e-150:
                        break
                    if abs(residual)/length <= 1e-13:
                        converged = True
                        break
                    projected -= np.clip(residual/length**2, -.05/length, .05/length)*gradient
                row["native_style_newton"] = dict(converged=converged, point=projected.tolist(),
                    inside_box=bool((abs(projected) <= .5).all()),
                    chord_to_target=float(np.linalg.norm(projected-midpoint)))
                if converged and (abs(projected) <= .5).all():
                    targets["native_style_newton"] = projected
                for name, candidate in targets.items():
                    delta = np.asarray(candidate)-start
                    assert (abs(delta[locked]) < 1e-12).all()
                    for fraction in (1., .5, .25, .125):
                        trial = xyz.copy()
                        for member in members:
                            trial[cells[incident] == member] += fraction*delta
                        row["trials"].append(dict(target=name, fraction=fraction,
                            displacement=float(np.linalg.norm(fraction*delta)), **measures(trial)))
                np.savez_compressed(a.run/f"witness-{index}.npz", points=xyz, cells=cells[incident],
                    members=members, probes=probes, endpoint_nodes=np.asarray(endpoints[node]))
                result["rows"].append(row)
                part.update(1, "Checked independent midpoint motions", incident_tets=len(incident),
                    valid_trials=sum(trial["sampled_valid"] for trial in row["trials"]))
                w.save(a.run/"result.json", result)
    result.update(complete=True, full_wall_projection_valid=sum(any(trial["fraction"] == 1 and
        trial["sampled_valid"] for trial in row["trials"]) for row in result["rows"]),
        all_endpoints_same_level=all(row["same_endpoint_level"] for row in result["rows"]))
    w.save(a.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-rollback-native-v1/edge")
    parser.add_argument("--observer", type=Path, default=w.HERE/"runs/macbook-pmy-saved776-forward-v1")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = parser.parse_args()
    for name in ("run", "source", "observer", "fixtures"):
        setattr(a, name, getattr(a, name).resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    paths = [Path(__file__), Path(w.__file__), a.source/"result.json", a.source/"mesh.npz",
        a.observer/"result.json", *list((a.fixtures/"tpms-pmy-sheet").iterdir())]
    paths += [w.HERE/name for name in ("pmy_neck_control.py", "quadrature_audit.py", "periodic_trace.py",
                                      "tpms_family_fixtures.py", "job_parts.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in paths}
    for index, path in enumerate(paths):
        if path.suffix == ".py":
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=180, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, no_meshing=True))
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    from job_parts import worker_part
    with Job(parts={"curvature": 1}) as job:
        with worker_part(job, "curvature") as part:
            record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
                a.run/"worker.log", part, timeout=180)
        w.save(a.run/"worker.json", record)
        assert record["returncode"] == 0, "Read worker.log"
        assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
        row = json.loads((a.run/"result.json").read_text())
        job.result("curvature-feasibility", witnesses=len(row["rows"]),
            full_wall_projection_valid=row["full_wall_projection_valid"],
            all_endpoints_same_level=row["all_endpoints_same_level"])


if __name__ == "__main__":
    main()
