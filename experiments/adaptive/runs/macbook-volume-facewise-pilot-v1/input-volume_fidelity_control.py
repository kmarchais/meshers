"""Bounded signed-wall-bias control on saved native Tet10 meshes.

The controller reads the implicit field, never a reference volume. A first-order
boundary integral estimates signed volume bias. One bounded normal correction
tries to remove that bias without increasing its own surface-residual maximum.
The fixed specimen volume and independent boundary grids are held-out observers.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
import time

import microgen_workflow as w
import numpy as np

import comparison_fixtures as fixtures
from comparison_mesh_audit import observe
from periodic_trace import periodic_pairs
from quadratic_edge_removal_probe import measure
import scorecard


def barycentric(degree):
    return np.array([(i, j, degree-i-j) for i in range(degree+1)
                     for j in range(degree+1-i)], dtype=float)/degree


def tri_basis(bary):
    a, b, c = np.asarray(bary).T
    values = np.array([a*(2*a-1), b*(2*b-1), c*(2*c-1), 4*a*b, 4*b*c, 4*c*a]).T
    du = np.array([1-4*a, 4*b-1, np.zeros_like(a), 4*(a-b), 4*c, -4*c]).T
    dv = np.array([1-4*a, np.zeros_like(a), 4*c-1, -4*b, 4*b, 4*(a-c)]).T
    return values, du, dv


def quadrature(degree):
    x, weights = np.polynomial.legendre.leggauss(degree)
    x, weights = (x+1)/2, weights/2
    u, v = np.meshgrid(x, x, indexing="ij")
    wu, wv = np.meshgrid(weights, weights, indexing="ij")
    return np.column_stack((1-u.ravel()-(1-u.ravel())*v.ravel(),
                            u.ravel(), (1-u.ravel())*v.ravel())), (wu*wv*(1-u)).ravel()


def wall_data(target, points, surface):
    corners = points[surface[:, :3]]
    caps = np.any(np.all(abs(corners+.5) < 1e-10, axis=1)
                  | np.all(abs(corners-.5) < 1e-10, axis=1), axis=1)
    wall = surface[~caps]
    center = points[wall[:, :3]].mean(1)
    levels = np.array([v for v in target.levels if np.isfinite(v)])
    selected = levels[np.argmin(abs(target.function(*center.T)[:, None]-levels), axis=1)]
    sign = np.where(selected == target.levels[1], 1., -1.)
    return wall, selected, sign


def samples(target, points, wall, levels, sign, bary, weights=None):
    basis, du, dv = tri_basis(bary)
    coordinates = points[wall]
    positions = np.einsum("sn,fnd->fsd", basis, coordinates)
    flat = positions.reshape(-1, 3)
    gradients = scorecard.gradient(target.function, flat).reshape(positions.shape)
    norm = np.linalg.norm(gradients, axis=2)
    assert np.all(norm > 1e-10), "Stationary wall sample requires a separate control"
    normals = gradients/norm[..., None]
    residual = sign[:, None]*(target.function(*flat.T).reshape(norm.shape)-levels[:, None])/norm
    result = dict(residual=residual, normals=normals*sign[:, None, None], basis=basis,
                  positions=positions)
    if weights is not None:
        cross = np.cross(np.einsum("sn,fnd->fsd", du, coordinates),
                         np.einsum("sn,fnd->fsd", dv, coordinates))
        area = np.linalg.norm(cross, axis=2)
        result.update(weights=area*weights[None], cross=cross,
                      signed_bias=float(np.sum(residual*area*weights[None])),
                      area=float(np.sum(area*weights[None])))
    return result


def integral(target, points, wall, levels, sign, degree):
    bary, weights = quadrature(degree)
    values = samples(target, points, wall, levels, sign, bary, weights)
    return dict(gauss_degree=degree, signed_bias=values["signed_bias"], area=values["area"],
                residual_max=float(abs(values["residual"]).max()))


def correct(target, points, cells, surface, classes, original_points, *, facewise=False):
    """One field-only linear correction, with local bounds and eight halvings."""
    wall, levels, sign = wall_data(target, points, surface)
    mids = np.unique(wall[:, 3:])
    roots, inverse = np.unique(classes[mids], return_inverse=True)
    columns = np.full(len(points), -1, dtype=int)
    columns[mids] = inverse
    locked = abs(abs(points)-.5) < 1e-10
    raw = np.zeros_like(points)
    raw[mids] = scorecard.gradient(target.function, points[mids])
    raw[locked] = 0
    grouped = np.zeros_like(points)
    np.add.at(grouped, classes, raw)
    direction = grouped[classes]
    direction /= np.maximum(np.linalg.norm(direction, axis=1), 1e-30)[:, None]
    edges = wall[:, np.array([[0, 1], [1, 2], [2, 0]])]
    lengths = np.linalg.norm(points[edges[:, :, 1]]-points[edges[:, :, 0]], axis=2)
    limits = np.full(len(roots), np.inf)
    np.minimum.at(limits, columns[wall[:, 3:]].ravel(), .02*lengths.ravel())
    original_limits = limits.copy()
    already_used = np.zeros(len(roots))
    np.maximum.at(already_used, inverse, np.linalg.norm(points[mids]-original_points[mids], axis=1))
    # Triangle inequality preserves the original 2% budget even though this
    # field-gradient direction can differ from the first fit's direction.
    limits = np.maximum(limits-already_used, 0.)
    bary, weights = quadrature(6)
    data = samples(target, points, wall, levels, sign, bary, weights)
    before6 = integral(target, points, wall, levels, sign, 6)
    before8 = integral(target, points, wall, levels, sign, 8)
    before_grid = samples(target, points, wall, levels, sign, barycentric(6))
    maximum = float(abs(before_grid["residual"]).max())
    face_maximum = np.max(abs(before_grid["residual"]), axis=1)
    coeff = np.einsum("fsd,fed->fse", data["normals"], direction[wall[:, 3:]])
    coeff *= data["basis"][None, :, 3:]
    columns3 = np.broadcast_to(columns[wall[:, None, 3:]], coeff.shape)
    gradient = np.bincount(columns3.ravel(),
                          weights=(coeff*data["weights"][..., None]).ravel(), minlength=len(roots))
    diagonal = np.bincount(columns3.ravel(),
                          weights=(coeff**2*data["weights"][..., None]).ravel(), minlength=len(roots))
    diagonal = np.maximum(diagonal, diagonal.max()*1e-12)
    denominator = float(np.sum(gradient**2/diagonal))
    step = -data["signed_bias"]*(gradient/diagonal)/max(denominator, 1e-30)

    # Give each of three midpoint contributions a third of the remaining
    # linearized distance budget. A dominant-error face cannot be worsened
    # by compensating changes elsewhere, even when its corners are fixed.
    grid_coeff = np.einsum("fsd,fed->fse", before_grid["normals"], direction[wall[:, 3:]])
    grid_coeff *= before_grid["basis"][None, :, 3:]
    lower, upper = -limits.copy(), limits.copy()
    gc = np.broadcast_to(columns[wall[:, None, 3:]], grid_coeff.shape)
    nonzero = abs(grid_coeff) > 1e-14
    denominator_grid = np.where(nonzero, grid_coeff, 1)
    allowance = face_maximum[:, None] if facewise else maximum
    lo = (-allowance-before_grid["residual"])[..., None]/(3*denominator_grid)
    hi = (allowance-before_grid["residual"])[..., None]/(3*denominator_grid)
    np.maximum.at(lower, gc[nonzero], np.minimum(lo, hi)[nonzero])
    np.minimum.at(upper, gc[nonzero], np.maximum(lo, hi)[nonzero])
    assert np.all(lower <= 1e-12) and np.all(upper >= -1e-12)
    step = np.clip(step, lower, upper)
    touched = np.zeros(len(points), dtype=bool)
    touched[mids[abs(step[inverse]) > 0]] = True
    affected = np.any(touched[cells], axis=1)
    records = []
    moved = points.copy()
    quadrature_stable = abs(before6["signed_bias"]-before8["signed_bias"]) <= .1*abs(before8["signed_bias"])
    if quadrature_stable and affected.any():
        for fraction in (1., .5, .25, .125, .0625, .03125, .015625, .0078125):
            trial = points.copy()
            trial[mids] += direction[mids]*(fraction*step[inverse])[:, None]
            after = integral(target, trial, wall, levels, sign, 6)
            grid = samples(target, trial, wall, levels, sign, barycentric(6))
            grid_max = float(abs(grid["residual"]).max())
            face_increase = float(np.max(np.max(abs(grid["residual"]), axis=1)-face_maximum))
            distance_ok = face_increase <= 1e-12 if facewise else grid_max <= maximum+1e-12
            minimum = float(measure(trial[cells[affected]])[0].min())
            accepted = bool(abs(after["signed_bias"]) < abs(before6["signed_bias"])
                            and distance_ok and minimum > 0)
            records.append(dict(fraction=fraction, signed_bias=after["signed_bias"],
                                grid_max=grid_max, minimum_sampled_jacobian=minimum,
                                maximum_facewise_grid_increase=face_increase,
                                accepted=accepted))
            if accepted:
                moved = trial
                break
    np.testing.assert_array_equal(points[np.unique(cells[:, :4])], moved[np.unique(cells[:, :4])])
    np.testing.assert_array_equal(points[locked], moved[locked])
    return moved, dict(before6=before6, before8=before8,
                       after6=integral(target, moved, wall, levels, sign, 6),
                       after8=integral(target, moved, wall, levels, sign, 8),
                       quadrature_stable=bool(quadrature_stable), grid_max_before=maximum,
                       distance_guard="per-face degree6 maximum" if facewise else "global degree6 maximum",
                       active_groups=int(np.count_nonzero(step)), attempts=records,
                       maximum_move=float(np.linalg.norm(moved-points, axis=1).max()),
                       affected_elements=int(affected.sum()),
                       total_original_budget_ratio=float(np.max(np.linalg.norm(moved[mids]-original_points[mids], axis=1)/original_limits[inverse])),
                       reference_volume_used=False, benchmark_samples_used=False,
                       displacement_limit_shortest_edge_fraction=.02)


def controls():
    from scipy.spatial import ConvexHull

    points = np.vstack((np.eye(3), -np.eye(3)))
    faces = ConvexHull(points).simplices
    for _ in range(2):
        nodes = points.tolist()
        edge_ids = {}
        refined = []
        for a, b, c in faces:
            mids = []
            for u, v in ((a, b), (b, c), (c, a)):
                edge = tuple(sorted((int(u), int(v))))
                if edge not in edge_ids:
                    edge_ids[edge] = len(nodes)
                    vector = points[u]+points[v]
                    nodes.append((vector/np.linalg.norm(vector)).tolist())
                mids.append(edge_ids[edge])
            x, y, z = mids
            refined.extend(((a, x, z), (x, b, y), (z, y, c), (x, y, z)))
        points, faces = np.array(nodes), np.array(refined)
    nodes, edge_ids, surface = points.tolist(), {}, []
    for face in faces:
        if np.dot(np.cross(points[face[1]]-points[face[0]], points[face[2]]-points[face[0]]), points[face].mean(0)) < 0:
            face = face[[0, 2, 1]]
        mids = []
        for u, v in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge = tuple(sorted((int(u), int(v))))
            if edge not in edge_ids:
                edge_ids[edge] = len(nodes)
                vector = points[u]+points[v]
                nodes.append((vector/np.linalg.norm(vector)).tolist())
            mids.append(edge_ids[edge])
        surface.append(list(face)+mids)
    points, surface = np.array(nodes)*.3, np.array(surface)
    sphere = scorecard.Case("analytic-sphere", lambda x, y, z: x*x+y*y+z*z-.09, (-np.inf, 0.))
    bary, weights = quadrature(8)
    data = samples(sphere, points, surface, np.zeros(len(surface)), np.ones(len(surface)), bary, weights)
    volume = float(np.sum(np.sum(data["positions"]*data["cross"], axis=2)*weights)/3)
    exact = 4*np.pi*.3**3/3
    error = volume-exact
    assert data["signed_bias"] < 0 and error < 0
    assert abs(data["signed_bias"]-error) < .02*abs(error)
    perm = surface[:, [1, 2, 0, 4, 5, 3]]
    permuted = samples(sphere, points, perm, np.zeros(len(surface)), np.ones(len(surface)), bary, weights)
    np.testing.assert_allclose(permuted["signed_bias"], data["signed_bias"], atol=1e-12, rtol=0)
    flat = points[surface[0]].copy()
    flat[:, 2] = .2
    plane = scorecard.Case("analytic-plane", lambda x, y, z: z-.2, (-np.inf, 0.))
    planar = samples(plane, flat, np.arange(6)[None], np.array([0.]), np.array([1.]), bary, weights)
    assert abs(planar["signed_bias"]) < 1e-15
    return dict(sphere_exact_volume=exact, sphere_mesh_volume=volume,
                sphere_actual_signed_error=error, sphere_predicted_signed_bias=data["signed_bias"],
                sphere_relative_prediction_error=abs(data["signed_bias"]-error)/abs(error),
                plane_bias_zero=True, cyclic_face_order_invariant=True)


def worker(a):
    w.imports()
    selection_file = w.ROOT/"handoff/macbook/face-fit-selection-20261110.json"
    selection = json.loads(selection_file.read_text())
    entry = next(e for e in selection["cases"] if e["case"] == a.case)
    source = a.source/f"{a.case}-repeat-0-fit"
    baseline_source = a.source/f"{a.case}-repeat-0-baseline"
    saved = json.loads((source/"result.json").read_text())
    paths = [Path(__file__), selection_file, source/"result.json", source/"mesh.npz",
             baseline_source/"result.json", baseline_source/"mesh.npz",
             source/(a.case+"-samples.npz"), source/(a.case+"-specimen.json"),
             w.HERE/"comparison_fixtures.py", w.HERE/"comparison_mesh_audit.py",
             w.HERE/"quadratic_edge_removal_probe.py", w.HERE/"scorecard.py",
             w.HERE/"microgen_workflow.py", w.HERE/"periodic_trace.py",
             w.HERE/"box_projection.py", w.HERE/"facet_incidence_audit.py"]
    if a.case.startswith("spinodoid"):
        paths.append(fixtures.coefficients_path(a.case, Path(entry["fixtures"])))
    if a.case == "graded-gyroid-density-linear":
        paths.append(Path(entry["fixtures"])/"linear-density-calibration.json")
    hashes = {str(p.resolve()): w.digest(p) for p in paths}
    assert w.digest(source/"mesh.npz") == saved["mesh_sha256"]
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(source/(a.case+suffix), a.run/(a.case+suffix))
    with np.load(source/"mesh.npz") as m:
        points, cells, surface = (m[k].copy() for k in ("points", "tetrahedra", "surface"))
    baseline_result = json.loads((baseline_source/"result.json").read_text())
    assert w.digest(baseline_source/"mesh.npz") == baseline_result["mesh_sha256"]
    with np.load(baseline_source/"mesh.npz") as m:
        original_points = m["points"].copy()
        np.testing.assert_array_equal(cells, m["tetrahedra"])
        np.testing.assert_array_equal(surface, m["surface"])
    assert cells.shape[1] == 10 and max(len(points), len(cells)) <= w.CAP
    target = fixtures.case(a.case, Path(entry["fixtures"]))
    required = fixtures.periodic(a.case)
    classes = periodic_pairs(points, surface)[1] if required else np.arange(len(points))
    control = controls()
    started = time.perf_counter()
    moved, correction = correct(target, points, cells, surface, classes, original_points,
                                 facewise=a.facewise)
    assert correction["total_original_budget_ratio"] <= 1+1e-10
    seconds = time.perf_counter()-started
    np.savez_compressed(a.run/"mesh.npz", points=moved, tetrahedra=cells, surface=surface)
    observations = [observe(a.run, a.case, target, p, cells, surface, required,
                            periodic_auditor=periodic_pairs) for p in (points, moved)]
    before, after = observations
    for key in ("volume", "sampled_two_sided_max", "relative_volume_error"):
        np.testing.assert_allclose(before[0][key], saved["geometry"][key], atol=1e-12, rtol=1e-12)
    np.testing.assert_array_equal(list(before[2].values()), list(after[2].values()))
    original_bary = scorecard.BARY
    scorecard.BARY = np.unique(np.vstack((original_bary, barycentric(13))), axis=0)
    try:
        dense = [observe(a.run, a.case, target, p, cells, surface, required,
                         periodic_auditor=periodic_pairs)[0] for p in (points, moved)]
    finally:
        scorecard.BARY = original_bary
    reference = json.loads((a.run/(a.case+"-specimen.json")).read_text())
    exact = target.volume is not None
    report = dict(case=a.case, source=str(source), inputs=hashes, controls=control,
                  nodes=len(points), elements=len(cells), correction_seconds=seconds,
                  correction=correction,
                  before=dict(geometry=before[0], periodic=before[1], quality=before[2]),
                  after=dict(geometry=after[0], periodic=after[1], quality=after[2]),
                  dense_before=dense[0], dense_after=dense[1],
                  reference=dict(volume=reference["volume"], analytic_exact=exact,
                                 scramble_spread=reference.get("volume_scramble_spread"),
                                 sobol_estimates=reference.get("volume_sobol_estimates"),
                                 uncertainty_kind="Analytic exact value" if exact else "Scramble spread, not a rigorous error bound",
                                 original_native_fit_volume_change=before[0]["volume"]-baseline_result["geometry"]["volume"],
                                 predicted_volume_from_before6=before[0]["volume"]-correction["before6"]["signed_bias"],
                                 predicted_volume_from_before8=before[0]["volume"]-correction["before8"]["signed_bias"]),
                  counts_corners_and_corner_quality_unchanged=True, cap_coordinates_unchanged=True,
                  admitted=after[0]["geometry_admitted"], dense_admitted=dense[1]["geometry_admitted"],
                  independent_standard_boundary_nonworse=after[0]["sampled_two_sided_max"] <= before[0]["sampled_two_sided_max"]+1e-12,
                  independent_dense_boundary_nonworse=dense[1]["sampled_two_sided_max"] <= dense[0]["sampled_two_sided_max"]+1e-12,
                  mesh_sha256=w.digest(a.run/"mesh.npz"),
                  scope="One saved-mesh field-only correction. First-order bias estimate, not an exact volume certificate. No generation-time comparison or FE solve.")
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    report["inputs_unchanged"] = True
    w.save(a.run/"result.json", report)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--source", type=Path, default=w.HERE/"runs/macbook-native-face-fit-nine-v1")
    p.add_argument("--cases", default="schwarz-d-solid,gyroid")
    p.add_argument("--facewise", action="store_true", help="Use each face's original degree6 distance envelope")
    p.add_argument("--worker", action="store_true")
    p.add_argument("--case")
    a = p.parse_args()
    a.run, a.source = a.run.resolve(), a.source.resolve()
    if a.worker:
        return worker(a)
    from job_parts import worker_part
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run/"input-volume_fidelity_control.py")
    names = a.cases.split(",")
    w.save(a.run/"manifest.json", dict(cases=names, source=str(a.source),
           numerical_threads=w.NUMERICAL_THREAD_POLICY, memory_bytes=w.MEMORY,
           worker_timeout_s=300, nodes_and_tets_cap=w.CAP,
           facewise_distance_guard=a.facewise,
           method="One bounded signed wall residual correction. Reference volume is held out.",
           no_native_changes=True))
    rows = []
    with Job(parts=names) as job:
        for name in names:
            directory = a.run/name
            directory.mkdir()
            with worker_part(job, name) as part:
                record = w.bounded([sys.executable, str(Path(__file__).resolve()), "--worker",
                                    "--case", name, "--run", str(directory), "--source", str(a.source),
                                    *(["--facewise"] if a.facewise else [])],
                                   directory/"worker.log", part, timeout=300)
                w.save(directory/"worker.json", record)
                result = json.loads((directory/"result.json").read_text()) if (directory/"result.json").exists() else {}
                rows.append(dict(case=name, worker=record, directory=str(directory),
                                 complete=bool(result), admitted=result.get("admitted")))
                job.result(name, admitted=result.get("admitted"),
                           correction_s=result.get("correction_seconds"),
                           returncode=record["returncode"])
                w.save(a.run/"summary.json", dict(rows=rows, complete=len(rows)==len(names)))
                assert record["returncode"] == 0, "Read worker.log and worker.json"


if __name__ == "__main__":
    main()
