"""Localize crossed-cylinder errors without changing the saved meshes.

This is a diagnostic of the fixed CAD specimen, not a meshing timing. The
selection of crossing faces uses the two original analytic cylinder branches;
the frozen reverse samples are observers only.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import cad_union_comparison as base
import microgen_workflow as w


RUNS = w.HERE / "runs"
NATIVE = RUNS / "macbook-cad-union-pair-v1/native/native-tet4"
ACTUAL = RUNS / "macbook-cad-union-mmg-tolerance-v1/actual/actual-cad-gmsh-public-mmg"


def wall_mask(points, surface):
    import numpy as np

    triangles = points[surface[:, :3]]
    cap = np.any(np.all(abs(triangles + .5) < 1e-10, axis=1) |
                 np.all(abs(triangles - .5) < 1e-10, axis=1), axis=1)
    return ~cap


def branch_value(points):
    """Positive means the exposed wall belongs to the X-axis cylinder."""
    return points[..., 0]**2 - points[..., 1]**2


def face_groups(points, surface):
    import numpy as np

    wall = wall_mask(points, surface)
    branch = branch_value(points[surface[:, :3]])
    crossing = wall & (branch.min(1) < -1e-12) & (branch.max(1) > 1e-12)
    ridge_touch = wall & ~crossing & (abs(branch).min(1) <= 1e-12)
    return dict(crossing=crossing, ridge_touch=ridge_touch,
                smooth=wall & ~crossing & ~ridge_touch, cap=~wall)


def ridge_target(points):
    """An exact point on the intersection, using a deterministic polar map.

    This is a feasible target, not a claim of the nearest point on the ridge.
    """
    import numpy as np

    points = np.asarray(points)
    a = .5 * (abs(points[..., 0]) + abs(points[..., 1]))
    scale = base.fixture.RADIUS / np.maximum(np.hypot(a, points[..., 2]), 1e-30)
    radius = a * scale
    sign = np.where(points[..., :2] < 0, -1., 1.)
    return np.concatenate((sign * radius[..., None],
                           (points[..., 2] * scale)[..., None]), axis=-1)


def describe(source, part):
    import numpy as np
    import scorecard
    from periodic_trace import periodic_pairs

    saved = json.loads((source / "result.json").read_text())
    mesh_file = source / "mesh.npz"
    assert w.digest(mesh_file) == saved["mesh_sha256"]
    with np.load(mesh_file) as mesh:
        points, cells, surface = (mesh[key] for key in ("points", "tetrahedra", "surface"))
    assert cells.shape[1] == 4 and surface.shape[1] == 3
    assert max(len(points), len(cells)) <= w.CAP
    assert saved["qualified"], "Probe only admitted source candidates"
    _, _, periodic = periodic_pairs(points, surface, declared_bounds=([-.5]*3, [.5]*3))
    original = w.specimen
    w.specimen = lambda name: (base.target(), None, None) if name == base.NAME else original(name)
    try:
        measured = w.geometry(source, base.NAME, points, cells, surface,
                              detail_count=32, box_projection_fallback=True)
    finally:
        w.specimen = original
    assert measured["geometry_admitted"]
    for key in ("sampled_two_sided_max", "relative_volume_error", "q_min_corner"):
        np.testing.assert_allclose(measured[key], saved["geometry"][key], atol=1e-12, rtol=1e-10)
    part.update(.4, "Reproduced saved admission", tets=len(cells), nodes=len(points))
    groups = face_groups(points, surface)
    wall = wall_mask(points, surface)
    wall_ids = np.flatnonzero(wall)
    probes = np.einsum("sn,fnd->fsd", scorecard.BARY, points[surface[wall]])
    foot, converged = scorecard.project(base.target(), probes.reshape(-1, 3))
    distances = np.linalg.norm(foot - probes.reshape(-1, 3), axis=1).reshape(probes.shape[:2])
    converged = converged.reshape(probes.shape[:2])
    triangles = points[surface]
    areas = .5*np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                     triangles[:, 2]-triangles[:, 0]), axis=1)
    description = {}
    for name, mask in groups.items():
        selected = mask[wall_ids]
        values = distances[selected]
        description[name] = dict(faces=int(mask.sum()), area=float(areas[mask].sum()),
            diagnostic_forward_max=float(values.max(initial=0)),
            converged_samples=int(converged[selected].sum()), sample_count=int(values.size))
    nodes = np.unique(surface[wall])
    residual = abs(base.fixture.periodic_field(*points[nodes].T))
    description["wall_corners"] = dict(count=len(nodes), maximum_field_residual=float(residual.max()),
        off_wall_above_1e9=int((residual > 1e-9).sum()),
        worst_node=int(nodes[residual.argmax()]), worst_point=points[nodes[residual.argmax()]].tolist())
    for row in measured["worst_boundary_samples"]:
        row["face_group"] = next(name for name, mask in groups.items() if mask[row["face"]])
    part.update(.65, "Classified branch-crossing and smooth faces", crossing_faces=int(groups["crossing"].sum()))
    with np.load(source / (base.NAME + "-samples.npz")) as samples:
        reverse = {}
        for name in ("walls", "caps", "ridges", "rims"):
            values = scorecard.distance_to_mesh(points, surface, samples[name])
            worst = int(values.argmax())
            reverse[name] = dict(count=len(values), maximum=float(values.max()),
                p99=float(np.quantile(values, .99)), worst_point=samples[name][worst].tolist())
    part.update(1, "Localized frozen ridge, rim and wall witnesses",
                ridge_error=reverse["ridges"]["maximum"], rim_error=reverse["rims"]["maximum"])
    return dict(source=str(source), input_mesh_sha256=w.digest(mesh_file), nodes=len(points),
        elements=len(cells), measured=measured, periodic=periodic,
        face_groups=description, reverse_groups=reverse,
        projection_scope="Per-face grouping uses the existing unconstrained Newton diagnostic; shared admission above uses the standard box-aware observer")


def worker(args):
    w.imports()
    import numpy as np
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import Job

    sources = {"native": args.native, "actual_modified_mmg": args.actual}
    paths = {Path(__file__), Path(base.__file__), Path(base.fixture.__file__), Path(w.__file__)}
    paths.update(w.HERE / name for name in ("scorecard.py", "periodic_trace.py", "box_projection.py"))
    for source in sources.values():
        paths.update(source / name for name in ("mesh.npz", "result.json",
            base.NAME + "-samples.npz", base.NAME + "-specimen.json"))
    hashes = {str(p.resolve()): w.digest(p) for p in sorted(paths)}
    for source in sources.values():
        for suffix in ("-samples.npz", "-specimen.json"):
            assert w.digest(source / (base.NAME+suffix)) == w.digest(args.native / (base.NAME+suffix))
    for p in sorted(paths):
        if p.suffix == ".py":
            shutil.copyfile(p, args.run / ("input-"+p.name))
    result = dict(inputs=hashes, scope="Saved Tet4 boundary localization, no geometry change or timing claim",
        hypotheses=["Boundary faces bridge the two exposed cylinder branches and omit the intersection ridge",
                    "Ordinary single-cylinder chord error dominates",
                    "Wall corners remain off the analytic walls"], rows={})
    w.save(args.run / "manifest.json", dict(inputs=hashes, worker_timeout_s=180,
        memory_bytes=w.MEMORY, numerical_thread_policy=w.NUMERICAL_THREAD_POLICY,
        node_cap=w.CAP, tetrahedron_cap=w.CAP, hypotheses=result["hypotheses"]))
    with Job(parts={name: 1 for name in sources}) as job:
        for name, source in sources.items():
            with job.part(name) as part:
                part.update(0, "Loading immutable source mesh", nodes=0)
                result["rows"][name] = describe(source, part)
            row = result["rows"][name]
            job.result(name, ridge_error=row["reverse_groups"]["ridges"]["maximum"],
                crossing_faces=row["face_groups"]["crossing"]["faces"],
                worst_boundary=row["measured"]["sampled_two_sided_max"])
            w.save(args.run / "result.json", result)
    # These controls check the branch classifier and feasible ridge map itself.
    probes = np.array([[.12, .10, .05], [-.09, .10, -.08], [.01, -.02, .14]])
    target = ridge_target(probes)
    np.testing.assert_allclose(target[:, 0]**2+target[:, 2]**2, base.fixture.RADIUS**2, atol=1e-16)
    np.testing.assert_allclose(target[:, 1]**2+target[:, 2]**2, base.fixture.RADIUS**2, atol=1e-16)
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    result.update(complete=True, input_hashes_verified=True, original_meshes_unchanged=True,
                  ridge_projection_feasible_control=True)
    w.save(args.run / "result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--native", type=Path, default=NATIVE)
    parser.add_argument("--actual", type=Path, default=ACTUAL)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "native", "actual"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run / "worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=180)
    w.save(args.run / "worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    sys.path.insert(0, str(Path.home() / "bin"))
    from pueue_progress import result
    for name, row in json.loads((args.run / "result.json").read_text())["rows"].items():
        result(name, ridge_error=row["reverse_groups"]["ridges"]["maximum"],
               crossing_faces=row["face_groups"]["crossing"]["faces"],
               worst_boundary=row["measured"]["sampled_two_sided_max"])


if __name__ == "__main__":
    main()
