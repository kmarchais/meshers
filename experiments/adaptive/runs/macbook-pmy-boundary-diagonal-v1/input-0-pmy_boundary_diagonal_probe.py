"""Test boundary-diagonal alternatives in eight frozen PMY two-tet edge stars.

Independent local feasibility only. Preserve every outer quadratic interface,
but allow the two physical wall faces to change. No global candidate is made.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys

import microgen_workflow as w


def worker(a):
    import numpy as np
    import tpms_family_fixtures as fixtures
    from pmy_neck_control import raw_gradient
    from quadratic_edge_removal_probe import boundary, measure, oriented
    from nested_refinement import EDGES, positive_quality
    from fedoo.lib_elements.tetrahedron import Tet10
    from quadrature_audit import duffy_rule

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    saved = json.loads((a.source/"result.json").read_text())
    selected = json.loads((a.selection/"result.json").read_text())
    assert selected["complete"] and selected["source_mesh_sha256"] == saved["mesh_sha256"] == w.digest(a.source/"mesh.npz")
    with np.load(a.source/"mesh.npz") as f:
        points, cells, surface = (f[key] for key in ("points", "tetrahedra", "surface"))
    case = fixtures.case("tpms-pmy-sheet", a.fixtures)
    midpoint = {}
    for cell in cells:
        for k, (i, j) in enumerate(EDGES):
            key = tuple(sorted((int(cell[i]), int(cell[j]))))
            assert midpoint.setdefault(key, int(cell[k+4])) == cell[k+4]
    bary = np.array([[i/13, j/13, 1-(i+j)/13] for i in range(14) for j in range(14-i)])
    weights = w.tri_values(bary)
    element = Tet10(15)
    probes = np.vstack((element.xi_nd, element.xi_pg, [[.25]*3], duffy_rule(5)[0]))
    derivative = element.shape_function_derivative(probes)

    def metrics(x):
        det = np.linalg.det(np.einsum("san,tnb->tsab", derivative, x))
        low, high = det.min(1), det.max(1)
        return dict(minimum_determinant=float(low.min()), minimum_ratio=float((low/high).min()),
            positive_samples=bool((low > 0).all()), ratio_guard=bool((low/high > .001).all()))

    def field_measure(coordinates, faces, level):
        samples = np.einsum("sn,fnd->fsd", weights, coordinates[faces]).reshape(-1, 3)
        residual = abs(case.function(*samples.T)-level)
        gradient = np.linalg.norm(raw_gradient(samples), axis=1)
        return dict(maximum_field_residual=float(residual.max()),
            maximum_normalized_residual=float((residual/np.maximum(gradient, 1e-30)).max()),
            minimum_gradient_norm=float(gradient.min()), samples=len(samples))

    def face_nodes(nodes, mids):
        return tuple(nodes)+tuple(mids[tuple(sorted((nodes[i], nodes[j])))] for i, j in ((0, 1), (1, 2), (2, 0)))

    controls = metrics(element.xi_nd[None])
    assert abs(controls["minimum_determinant"]-1) < 1e-13
    result = dict(complete=False, rows=[], control=controls, no_global_meshing=True, local_connectivity_trials=True, no_global_candidate=True,
        no_admission_change=True, selection_source=str(a.selection),
        scope="Two boundary faces change; corners, counts and outer quadratic interface stay fixed. Saved failure witnesses select this diagnostic only. Field residual is not a certified distance.")
    with Job(parts={str(i): 1 for i in range(len(selected["rows"]))}) as job:
        for index, witness in enumerate(selected["rows"]):
            with job.part(str(index)) as part:
                node = witness["node"]
                u, v = witness["endpoint_nodes"]
                incident = np.flatnonzero((cells[:, 4:] == node).any(1))
                assert len(incident) == 2
                old = cells[incident]
                shared = set(old[0, :4]) & set(old[1, :4])
                assert len(shared) == 3 and u in shared and v in shared
                center = int((shared-{u, v}).pop())
                left = int((set(old[0, :4])-shared).pop())
                right = int((set(old[1, :4])-shared).pop())
                wall_ids = np.flatnonzero(((surface[:, :3] == u).any(1)) & ((surface[:, :3] == v).any(1)))
                assert len(wall_ids) == 2
                assert {tuple(sorted(face[:3])) for face in surface[wall_ids]} == {
                    tuple(sorted((u, v, left))), tuple(sorted((u, v, right)))}
                vertices = [u, v, center, left, right]
                assert np.max(abs(points[vertices])) < .5-1e-10, "This diagnostic excludes cap classes"
                values = case.function(*points[[u, v, left, right]].T)
                levels = np.asarray(case.levels)[np.argmin(abs(values[:, None]-case.levels), axis=1)]
                assert (levels == levels[0]).all()
                level = float(levels[0])
                row = dict(node=node, old_edge=[u, v], proposed_edge=[left, right], center=center,
                    elements=incident.tolist(), old_corner_q_min=float(positive_quality(points[old[:, :4]]).min()),
                    old_curved=metrics(points[old]), old_boundary=field_measure(points, surface[wall_ids], level), trials=[])
                key = tuple(sorted((left, right)))
                if key in midpoint:
                    row["rejection"] = "Alternative edge already exists elsewhere"
                else:
                    raw = np.array([[left, right, center, u], [left, right, center, v]])
                    sign = np.linalg.det(points[raw[:, 1:]]-points[raw[:, :1]])
                    row["common_face_signed_determinants"] = sign.tolist()
                    if sign.prod() >= 0:
                        row["rejection"] = "Alternative two tetrahedra overlap or degenerate"
                    else:
                        corner = np.array([oriented(points, tet) for tet in raw])
                        q = positive_quality(points[corner])
                        row.update(new_corner_q_min=float(q.min()), corner_quality_passed=bool((q >= .1).all()))
                        mids = dict(midpoint); mids[key] = node
                        new = np.array([list(tet)+[mids[tuple(sorted((int(tet[i]), int(tet[j]))))] for i, j in EDGES] for tet in corner])
                        walls = np.array([face_nodes(tuple(sorted(face)), mids) for face in ((u, left, right), (v, left, right))])
                        old_boundary = boundary(old)
                        old_walls = {face_nodes(tuple(sorted(face[:3])), midpoint) for face in surface[wall_ids]}
                        new_walls = set(map(tuple, walls))
                        assert old_boundary-old_walls == boundary(new)-new_walls
                        row["outer_quadratic_interface_identical"] = True
                        straight = points[[left, right]].mean(0)
                        projected = straight.copy(); converged = False
                        for _ in range(30):
                            residual = float(case.function(*projected)-level)
                            gradient = raw_gradient(projected[None])[0]
                            length = float(np.linalg.norm(gradient))
                            if length <= 1e-150:
                                break
                            if abs(residual)/length <= 1e-13:
                                converged = True
                                break
                            projected -= np.clip(residual/length**2, -.05/length, .05/length)*gradient
                        row.update(projection_converged=converged, projected_midpoint=projected.tolist(),
                            straight_midpoint=straight.tolist())
                        fractions = (1., .5, .25, 0.) if converged and (abs(projected) < .5).all() else (0.,)
                        for fraction in fractions:
                            trial = points.copy()
                            trial[node] = straight if fraction == 0 else straight+fraction*(projected-straight)
                            curved = metrics(trial[new])
                            field = field_measure(trial, walls, level)
                            _, volume = measure(trial[new])
                            _, before_volume = measure(points[old])
                            row["trials"].append(dict(fraction=fraction, curved=curved, boundary=field,
                                local_volume_change=float(volume.sum()-before_volume.sum()),
                                sampled_quality_and_field_improved=bool((q >= .1).all() and curved["positive_samples"] and curved["ratio_guard"] and
                                    field["maximum_normalized_residual"] < row["old_boundary"]["maximum_normalized_residual"])))
                        np.savez_compressed(a.run/f"star-{index}.npz", original_coordinates=points[old],
                            candidate_corner_nodes=corner, candidate_tets=new, candidate_wall=walls,
                            moved_node=np.asarray(node), projected_midpoint=projected, original_wall=surface[wall_ids])
                result["rows"].append(row)
                part.update(1, "Boundary diagonal checked", corner_quality=row.get("new_corner_q_min"),
                    feasible=sum(t["sampled_quality_and_field_improved"] for t in row["trials"]))
                w.save(a.run/"result.json", result)
    result.update(complete=True, locally_feasible_count=sum(any(t["sampled_quality_and_field_improved"]
        for t in row["trials"]) for row in result["rows"]))
    w.save(a.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-pmy-rollback-native-v1/edge")
    parser.add_argument("--selection", type=Path, default=w.HERE/"runs/macbook-pmy-curvature-feasibility-v1")
    parser.add_argument("--fixtures", type=Path, default=w.HERE/"runs/macbook-tpms-family-fixtures-v1")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = parser.parse_args()
    for key, value in vars(a).items():
        if isinstance(value, Path):
            setattr(a, key, value.resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    paths = [Path(__file__), Path(w.__file__), a.source/"result.json", a.source/"mesh.npz", a.selection/"result.json",
             *list((a.fixtures/"tpms-pmy-sheet").iterdir())]
    paths += [w.HERE/name for name in ("pmy_neck_control.py", "quadratic_edge_removal_probe.py", "nested_refinement.py",
        "quadrature_audit.py", "tpms_family_fixtures.py", "job_parts.py")]
    hashes = {str(path.resolve()): w.digest(path) for path in paths}
    for index, path in enumerate(paths):
        if path.suffix == ".py":
            shutil.copyfile(path, a.run/f"input-{index}-{path.name}")
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=180, memory_bytes=w.MEMORY, numerical_thread_policy=w.NUMERICAL_THREAD_POLICY))
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    from job_parts import worker_part
    with Job(parts={"diagonals": 1}) as job:
        with worker_part(job, "diagonals") as part:
            receipt = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
                a.run/"worker.log", part, timeout=180)
        w.save(a.run/"worker.json", receipt)
        assert receipt["returncode"] == 0, "Read worker.log"
        assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
        row = json.loads((a.run/"result.json").read_text())
        job.result("boundary-diagonal-probe", inspected=len(row["rows"]), locally_feasible=row["locally_feasible_count"])


if __name__ == "__main__":
    main()
