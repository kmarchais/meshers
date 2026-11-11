"""Bounded saved-mesh PMY diagnosis. No generation, mesh edits, or gate changes.

The driver limits the worker to 180 seconds and 8 GiB. The worker reports real
stages through Job. Run through qrun with --progress report and --check-start.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
RUNS = w.HERE / "runs"
GATE = 0.002
WITNESSES = 8


def basis(bary):
    import numpy as np

    a, b, c = np.asarray(bary).T
    return np.array([a * (2 * a - 1), b * (2 * b - 1), c * (2 * c - 1),
                     4 * a * b, 4 * b * c, 4 * c * a]).T


def bezier(nodes):
    """Quadratic Lagrange nodes to the six triangular Bernstein controls."""
    import numpy as np

    result = np.array(nodes, copy=True)
    for k, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
        result[..., 3 + k, :] = 2 * nodes[..., 3 + k, :] - (nodes[..., a, :] + nodes[..., b, :]) / 2
    return result


def child_maps():
    import numpy as np

    v = np.eye(3)
    a, b, c = (v[[i, j]].mean(0) for i, j in ((0, 1), (1, 2), (2, 0)))
    maps = []
    for tri in (np.array([v[0], a, c]), np.array([a, v[1], b]),
                np.array([c, b, v[2]]), np.array([a, b, c])):
        nodes = np.vstack((tri, tri[[0, 1]].mean(0), tri[[1, 2]].mean(0), tri[[2, 0]].mean(0)))
        x, y, z = nodes.T
        maps.append(bezier(np.array([x*x, y*y, z*z, 2*x*y, 2*y*z, 2*z*x]).T))
    return np.array(maps)


def exclude_ball(controls, query, radius=GATE, depth_limit=10, box_limit=50000):
    """Exclude a ball from every complete Tri6 patch by convex-hull subdivision.

    A positive result covers all faces, including caps, without a tessellation
    assumption. Bounds have 1e-12 padding; this is floating-point evidence,
    not an interval-arithmetic certificate. A negative result is inconclusive.
    """
    import numpy as np

    transforms = child_maps()
    active = controls
    discarded_min = float("inf")
    visited = 0
    for depth in range(depth_limit + 1):
        lo, hi = active.min(axis=1) - 1e-12, active.max(axis=1) + 1e-12
        delta = np.maximum(np.maximum(lo - query, query - hi), 0)
        lower = np.linalg.norm(delta, axis=1)
        visited += len(active)
        keep = lower <= radius
        if (~keep).any():
            discarded_min = min(discarded_min, float(lower[~keep].min()))
        if not keep.any():
            return dict(excluded=True, lower_bound=discarded_min, radius=radius,
                        depth=depth, visited_boxes=visited, unresolved_boxes=0)
        active = active[keep]
        if depth == depth_limit or visited + 4 * len(active) > box_limit:
            return dict(excluded=False, lower_bound=min(discarded_min, float(lower[keep].min())),
                        radius=radius, depth=depth, visited_boxes=visited, unresolved_boxes=len(active))
        active = np.einsum("cij,fjk->fcik", transforms, active).reshape(-1, 6, 3)
    raise AssertionError("Unreachable subdivision exit")


def newton(case, probe):
    """Replay the historic observer's feasible, but not closest, Newton foot."""
    import numpy as np
    from scorecard import gradient

    levels = np.array([v for v in case.levels if np.isfinite(v)])
    level = levels[np.argmin(abs(levels - case.function(*probe)))]
    locked = np.isclose(abs(probe), 0.5, atol=1e-12, rtol=0)
    foot = probe.copy()
    for _ in range(40):
        g = gradient(case.function, foot[None])[0]
        g[locked] = 0
        step = (case.function(*foot) - level) * g / max(float(g @ g), 1e-30)
        foot -= step
        if np.linalg.norm(step) < 1e-12:
            break
    residual = abs(float(case.function(*foot) - level)) / max(float(np.linalg.norm(gradient(case.function, foot[None]))), 1e-30)
    return foot, dict(distance=float(np.linalg.norm(foot - probe)), level=float(level),
                     residual_distance=residual, feasible=bool(residual <= 1e-9 and np.max(abs(foot)) <= 0.5 + 1e-9))


def closest_feasible(case, probe, starts, levels=None):
    """Bounded multistart constrained distance; each feasible result is an upper bound."""
    import numpy as np
    from scipy.optimize import minimize
    from scorecard import gradient

    levels = [v for v in case.levels if np.isfinite(v)] if levels is None else levels
    locked = np.isclose(abs(probe), 0.5, atol=1e-12, rtol=0)
    lo, hi = np.full(3, -0.5), np.full(3, 0.5)
    lo[locked] = hi[locked] = probe[locked]
    candidates = []
    scale = 0.02
    for level in levels:
        for start in starts:
            start = np.clip(start, lo, hi)
            result = minimize(lambda x: float(x @ x), (start - probe) / scale,
                              jac=lambda x: 2*x, method="SLSQP",
                              bounds=list(zip((lo - probe) / scale, (hi - probe) / scale)),
                              constraints={"type": "eq",
                                           "fun": lambda x: float(case.function(*(probe + scale*x)) - level) / scale,
                                           "jac": lambda x: gradient(case.function, (probe + scale*x)[None])[0]},
                              options={"maxiter": 100, "ftol": 1e-12})
            foot = probe + scale * result.x
            residual = abs(float(case.function(*foot) - level)) / max(float(np.linalg.norm(gradient(case.function, foot[None]))), 1e-30)
            if np.isfinite(foot).all() and np.isfinite(residual) and residual < 1e-9 and (foot >= lo - 1e-10).all() and (foot <= hi + 1e-10).all():
                candidates.append(dict(foot=foot.tolist(), distance=float(np.linalg.norm(foot - probe)),
                                       level=float(level), residual_distance=residual,
                                       optimizer_success=bool(result.success)))
    return min(candidates, key=lambda x: x["distance"]) if candidates else None


def nearest_patch(nodes, query):
    import numpy as np
    from scipy.optimize import minimize

    def point(uv):
        return basis([[uv[0], uv[1], 1 - uv.sum()]])[0] @ nodes

    candidates = []
    for start in ([1/3, 1/3], [0., 0.], [1., 0.], [0., 1.]):
        result = minimize(lambda uv: float(np.sum((point(uv) - query)**2)), start,
                          method="SLSQP", bounds=[(0, 1), (0, 1)],
                          constraints={"type": "ineq", "fun": lambda uv: 1 - uv.sum()},
                          options={"maxiter": 100, "ftol": 1e-16})
        if min(result.x) >= -1e-10 and sum(result.x) <= 1 + 1e-10:
            foot = point(result.x)
            candidates.append(dict(foot=foot.tolist(), distance=float(np.linalg.norm(foot - query)),
                                   barycentric=[float(result.x[0]), float(result.x[1]), float(1-result.x.sum())]))
    return min(candidates, key=lambda x: x["distance"])


def locator(nodes, divisions=4):
    """A coarse locator selects witnesses only. Exclusion uses exact patches."""
    import numpy as np
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

    ij = [(i, j) for i in range(divisions + 1) for j in range(divisions + 1 - i)]
    lookup = {v: k for k, v in enumerate(ij)}
    bary = np.array([[i/divisions, j/divisions, 1-(i+j)/divisions] for i, j in ij])
    faces = []
    for i, j in ij:
        if i+j < divisions:
            faces.append([lookup[i, j], lookup[i+1, j], lookup[i, j+1]])
        if i+j < divisions-1:
            faces.append([lookup[i+1, j], lookup[i+1, j+1], lookup[i, j+1]])
    vertices = np.einsum("sn,fnd->fsd", basis(bary), nodes).reshape(-1, 3)
    faces = (np.array(faces)[None] + np.arange(len(nodes))[:, None, None]*len(bary)).reshape(-1, 3)
    poly = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(vertices, deep=True))
    poly.SetPoints(points)
    cells = vtk.vtkCellArray()
    cells.SetCells(len(faces), numpy_to_vtkIdTypeArray(np.column_stack((np.full(len(faces), 3), faces)).astype(np.int64).ravel(), deep=True))
    poly.SetPolys(cells)
    tree = vtk.vtkStaticCellLocator()
    tree.SetDataSet(poly)
    tree.BuildLocator()
    return tree


def controls():
    import numpy as np
    from scorecard import Case

    tri = np.array([[-.4, -.4, 0.], [.4, -.4, 0.], [0., .4, 0.]])
    nodes = np.vstack((tri, tri[[0, 1]].mean(0), tri[[1, 2]].mean(0), tri[[2, 0]].mean(0)))
    assert exclude_ball(bezier(nodes)[None], np.array([0., 0., .01]))["excluded"]
    assert not exclude_ball(bezier(nodes)[None], np.array([0., 0., .001]), depth_limit=5)["excluded"]
    # Midpoint elevation makes the patch pass through the query. Corner-only
    # distance would falsely exclude it, so this exercises the curved bounds.
    curved = nodes.copy()
    curved[3:, 2] = .03
    query = basis([[1/3]*3])[0] @ curved
    assert not exclude_ball(bezier(curved)[None], query, depth_limit=5)["excluded"]
    plane = Case("plane", lambda x, y, z: z, (-float("inf"), 0.))
    p = np.array([.1, .1, .01])
    np.testing.assert_allclose(closest_feasible(plane, p, [p])["distance"], .01, atol=1e-9)
    roots = Case("three-planes", lambda x, y, z: x*(x-.3)*(x+.3), (-float("inf"), 0.))
    p = np.array([.165, 0., 0.])
    _, old = newton(roots, p)
    best = closest_feasible(roots, p, [p, np.array([.3, 0., 0.]), np.zeros(3)])
    np.testing.assert_allclose(best["distance"], .135, atol=1e-8)
    assert old["feasible"] and old["distance"] > 3*best["distance"]
    return dict(plane_far_excluded=True, plane_near_retained=True, curved_patch_retained=True,
                analytic_plane_distance=True, analytic_newton_branch_jump_detected=True)


def patch_trial(case, points, cells, face):
    """One full midpoint projection on a copied local star; never a candidate mesh."""
    import numpy as np
    from scorecard import gradient
    from quadratic_edge_removal_probe import measure

    xyz = points[face]
    level = case.levels[int(np.argmin(abs(np.array(case.levels) - case.function(*xyz[0]))))]
    target = xyz[3:].copy()
    for i, p in enumerate(target.copy()):
        best = closest_feasible(case, p, [p, *xyz[:3]], [level])
        if best is None:
            return dict(completed=False, reason="No feasible midpoint projection")
        target[i] = best["foot"]
    affected = np.flatnonzero(np.isin(cells[:, 4:], face[3:]).any(axis=1))
    original = points[cells[affected]]
    trial = original.copy()
    for node, p in zip(face[3:], target, strict=True):
        trial[cells[affected] == node] = p
    before, _ = measure(original)
    after, _ = measure(trial)
    normals = gradient(case.function, xyz[:3])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    return dict(completed=True, affected_tets=len(affected),
                minimum_original_jacobian=float(before.min()), minimum_trial_jacobian=float(after.min()),
                invalid_trial_tets=int((after <= 0).sum()),
                midpoint_displacements=np.linalg.norm(target-xyz[3:], axis=1).tolist(),
                midpoint_target=target.tolist(),
                corner_normal_min_cosine=float((normals @ normals.T).min()),
                node_field_values=np.asarray(case.function(*xyz.T)).tolist(),
                scope="Local sensitivity only. Periodic copies and global gates are not evaluated; no mesh is published.")


def worker(a):
    import numpy as np
    import vtk
    from microgen.shape import surface_functions
    from scorecard import gradient
    import tpms_family_fixtures as fixtures

    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    vtk.vtkSMPTools.SetBackend("Sequential")
    assert vtk.vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    source_result = json.loads((a.source / "result.json").read_text())
    historical = json.loads(a.history.read_text())
    detail = next(r for r in historical["rows"] if r["case"] == NAME and r["arm"] == "meshers")
    adaptive = json.loads(a.adaptive.read_text())
    frozen = json.loads((a.fixtures / "result.json").read_text())
    specimen_files = list((a.fixtures / NAME).iterdir())
    frozen_paths = {Path(p).resolve(): h for p, h in frozen["inputs"].items()}
    assert frozen["passed"]
    for path in specimen_files:
        assert w.digest(path) == frozen_paths[path.resolve()], f"Changed fixture: {path}"
    mesh_path = a.source / "mesh.npz"
    assert w.digest(mesh_path) == source_result["mesh_sha256"]
    assert w.digest(mesh_path) == historical["inputs"][str(a.source.relative_to(w.ROOT) / "mesh.npz")]
    files = [Path(__file__), Path(w.__file__), Path(fixtures.__file__), Path(surface_functions.__file__),
             mesh_path, a.source / "result.json", a.history, a.adaptive, a.fixtures / "result.json", *specimen_files]
    files += [w.HERE / n for n in ("scorecard.py", "quadratic_edge_removal_probe.py", "quadrature_audit.py", "nested_refinement.py")]
    files += [w.ROOT / "crates/meshers-core/src" / n for n in ("projection.rs", "tet_adapt.rs", "quadratic.rs", "quadratic/fit.rs")]
    hashes = {str(p.resolve()): w.digest(p) for p in files}
    w.save(a.run / "manifest.json", dict(inputs=hashes, numerical_threads=w.NUMERICAL_THREAD_POLICY,
           gate_distance=GATE, gate_nodes=250000, gate_tets=250000, timeout_s=180, memory_bytes=w.MEMORY,
           protocol="Saved coarse PMY mesh only. Adaptive attempt has no saved mesh because the cap check preceded saving.",
           witness_budget=WITNESSES, no_generation=True, no_admission_claim=True))
    case = fixtures.case(NAME, a.fixtures)
    with np.load(mesh_path) as mesh:
        points, cells, surface = (mesh[k] for k in ("points", "tetrahedra", "surface"))
    assert cells.shape[1] == 10 and surface.shape[1] == 6
    with np.load(a.fixtures / NAME / (NAME + "-samples.npz")) as samples:
        queries = np.concatenate((samples["walls"], samples["caps"]))
        wall_count = len(samples["walls"])
    nodes = points[surface]
    result = dict(source_geometry=source_result["geometry"], adaptive_attempt=dict(
        nodes=adaptive["nodes"], elements=adaptive["elements"], error=adaptive["error"],
        diagnostics=adaptive["diagnostics"], mesh_available=(a.adaptive.parent / "mesh.npz").exists()),
        gates_unchanged=True, candidate_generated=False, inputs=hashes)
    selected_faces = []
    with Job(parts={"controls": 1, "forward": 3, "reverse": 5, "local-stars": 3}) as job:
        with job.part("controls") as part:
            result["controls"] = controls()
            part.update(1, "Analytic controls passed", controls=len(result["controls"]))
        job.result("controls", passed=True)
        with job.part("forward") as part:
            rows = []
            witnesses = detail["details"]["worst_boundary_samples"]
            for i, old in enumerate(witnesses):
                face = int(old["face"])
                p = np.array(old["probe"])
                np.testing.assert_allclose(basis([__import__("scorecard").BARY[old["probe_index"]]])[0] @ nodes[face], p, atol=1e-13, rtol=0)
                foot, observation = newton(case, p)
                np.testing.assert_allclose(observation["distance"], old["distance"], atol=1e-7, rtol=0)
                best = closest_feasible(case, p, [p, foot, *nodes[face, :3]])
                assert best is not None
                first_order = abs(float(case.function(*p) - observation["level"])) / float(np.linalg.norm(gradient(case.function, p[None])))
                rows.append(dict(face=face, probe=p.tolist(), historic=observation, closest_feasible=best,
                                 forward_overstatement_at_least=observation["distance"]-best["distance"],
                                 first_order_distance=first_order))
                selected_faces.append(face)
                part.update((i+1)/len(witnesses), "Comparing fixed forward witnesses", newton=observation["distance"], feasible=best["distance"])
            result["forward_witnesses"] = rows
        job.result("forward", largest_overstatement=max(r["forward_overstatement_at_least"] for r in rows))
        w.save(a.run / "partial.json", result)
        with job.part("reverse") as part:
            tree = locator(nodes)
            distances = np.empty(len(queries))
            face_ids = np.empty(len(queries), dtype=int)
            closest = [0., 0., 0.]
            cell, sub, d2 = vtk.reference(0), vtk.reference(0), vtk.reference(0.)
            for i, query in enumerate(queries):
                tree.FindClosestPoint(query, closest, cell, sub, d2)
                distances[i] = np.sqrt(float(d2))
                face_ids[i] = int(cell)//16
                if i % 512 == 0:
                    part.update(.45*(i+1)/len(queries), "Selecting reverse witnesses", sampled_max=float(distances[:i+1].max()))
            controls_all = bezier(nodes)
            selected = np.argsort(distances)[-WITNESSES:][::-1]
            rows = []
            for i, index in enumerate(selected):
                p = queries[index]
                face = int(face_ids[index])
                exclusion = exclude_ball(controls_all, p)
                nearest = nearest_patch(nodes[face], p)
                level = min(case.levels, key=lambda v: abs(float(case.function(*p))-v))
                residual = abs(float(case.function(*p))-level) / max(float(np.linalg.norm(gradient(case.function, p[None]))), 1e-30)
                assert index >= wall_count or residual < 1e-9
                rows.append(dict(query_index=int(index), kind="wall" if index < wall_count else "cap",
                                 probe=p.tolist(), target_residual_distance=residual,
                                 tessellation4_distance=float(distances[index]), nearest_tessellation_face=face,
                                 exact_patch_feasible_upper_bound=nearest, complete_surface_ball_exclusion=exclusion))
                selected_faces.append(face)
                part.update(.45+.55*(i+1)/len(selected), "Bounding complete curved surface distance", excluded=int(exclusion["excluded"]), lower=exclusion["lower_bound"])
            result["reverse_witnesses"] = rows
            result["reverse_selection"] = dict(queries=len(queries), subdivisions=4,
                                               max_distance=float(distances.max()), above_gate=int((distances > GATE).sum()))
        job.result("reverse", proven_gaps=sum(r["complete_surface_ball_exclusion"]["excluded"] for r in rows))
        w.save(a.run / "partial.json", result)
        with job.part("local-stars") as part:
            rows = []
            faces = list(dict.fromkeys(selected_faces))[:WITNESSES]
            for i, face in enumerate(faces):
                trial = patch_trial(case, points, cells, surface[face])
                rows.append(dict(face=face, **trial))
                part.update((i+1)/len(faces), "Testing full wall projection on local stars", invalid_tets=trial.get("invalid_trial_tets"))
            result["local_star_trials"] = rows
        job.result("local-stars", blocked_faces=sum(r.get("invalid_trial_tets", 0) > 0 for r in rows))
    assert all(w.digest(Path(p)) == h for p, h in hashes.items()), "An input changed during diagnosis"
    result.update(complete=True, input_hashes_verified=True,
                  interpretation_limits="Feasible forward feet prove overstatement only. Complete Bezier exclusion establishes geometric gaps at frozen samples with floating-point padding. Local star trials are not admitted meshes. Global topological equivalence is not certified.")
    w.save(a.run / "result.json", result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--source", type=Path, default=RUNS / "macbook-tpms-family-remaining-v1" / (NAME + "-native-fit"))
    p.add_argument("--fixtures", type=Path, default=RUNS / "macbook-tpms-family-fixtures-v1")
    p.add_argument("--history", type=Path, default=RUNS / "macbook-tpms-hard-failures-v1/result.json")
    p.add_argument("--adaptive", type=Path, default=RUNS / "macbook-tpms-remaining-adaptation-v1" / (NAME + "-native-adaptive/result.json"))
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    for key in ("run", "source", "fixtures", "history", "adaptive"):
        setattr(a, key, getattr(a, key).resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run / "input-pmy_failure_probe.py")
    # The worker's Job owns the progress report. This adapter only lets the
    # existing bounded runner enforce its deadline and process-tree RSS limit.
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
                       a.run / "worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=180)
    w.save(a.run / "worker.json", record)
    assert record["returncode"] == 0, "Probe failed; read worker.log and worker.json"
    assert json.loads((a.run / "result.json").read_text())["complete"]


if __name__ == "__main__":
    main()
