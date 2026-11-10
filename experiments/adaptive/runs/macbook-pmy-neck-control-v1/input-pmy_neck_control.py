"""One diagnostic control: insert a field-proved void path into the fixed BCC grid.

The cells16 grid, its phase, all original vertices, the physical field, and the
frozen specimen stay fixed. A stationary point is found from raw-field starts
in one declared lattice cell. Local stellar subdivisions make an outside path
an edge chain before clipping. This is a topology control, not an admitted mesh.
The driver enforces 180 seconds and 8 GiB; submit it through qrun.
"""

import argparse
from itertools import product
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
CELL = (1, 12, 1)  # Frozen diagnostic cell, not selected from admission samples.
CELLS = 16
RUNS = w.HERE / "runs"
EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
FACES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))


def raw_gradient(points):
    import numpy as np

    x, y, z = (2*np.pi*np.asarray(points)).T
    return 2*np.pi*np.array([
        -2*np.sin(x)*np.cos(y)*np.cos(z) + 2*np.cos(2*x)*np.sin(y) + np.cos(x)*np.sin(2*z),
        -2*np.cos(x)*np.sin(y)*np.cos(z) + np.sin(2*x)*np.cos(y) + 2*np.cos(2*y)*np.sin(z),
        -2*np.cos(x)*np.cos(y)*np.sin(z) + 2*np.sin(x)*np.cos(2*z) + np.sin(2*y)*np.cos(z),
    ]).T


def critical_search(case):
    import numpy as np
    from scipy.optimize import root

    lo = -.5 + np.array(CELL)/CELLS
    hi = lo + 1/CELLS
    found = []
    for fraction in product((.25, .5, .75), repeat=3):
        start = lo + np.array(fraction)/CELLS
        solve = root(lambda p: raw_gradient(p[None])[0], start)
        p = solve.x
        if (not np.isfinite(p).all() or np.linalg.norm(raw_gradient(p[None])) > 1e-8
                or (p < lo).any() or (p > hi).any()
                or any(np.linalg.norm(p-r["point"]) < 1e-7 for r in found)):
            continue
        h = 1e-6
        hessian = np.column_stack([(raw_gradient((p+np.eye(3)[a]*h)[None])[0]
                                   - raw_gradient((p-np.eye(3)[a]*h)[None])[0])/(2*h) for a in range(3)])
        values, vectors = np.linalg.eigh((hessian+hessian.T)/2)
        found.append(dict(point=p.tolist(), value=float(case.function(*p)),
                          gradient_norm=float(np.linalg.norm(raw_gradient(p[None]))),
                          hessian_eigenvalues=values.tolist(), positive_axis=vectors[:, -1].tolist()))
    selected = [r for r in found if r["value"] > case.levels[1]
                and r["hessian_eigenvalues"][1] < 0 < r["hessian_eigenvalues"][2]]
    assert len(selected) == 1, "The declared cell must contain one upper-void neck"
    return selected[0], found


def barycentric(points, cells, query):
    import numpy as np

    xyz = points[cells]
    b = np.linalg.solve((xyz[:, 1:]-xyz[:, :1]).transpose(0, 2, 1),
                        (query-xyz[:, 0])[..., None])[..., 0]
    return np.column_stack((1-b.sum(1), b))


def containing(points, cells, query, tolerance=1e-10):
    import numpy as np

    xyz = points[cells]
    selected = np.flatnonzero(((query >= xyz.min(1)-tolerance) & (query <= xyz.max(1)+tolerance)).all(1))
    if not len(selected):
        return selected, np.empty((0, 4))
    b = barycentric(points, cells[selected], query)
    inside = b.min(1) >= -tolerance
    return selected[inside], b[inside]


def oriented(points, cells):
    import numpy as np

    cells = np.asarray(cells, dtype=np.int64).copy()
    d = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
    assert np.min(abs(d)) > 1e-20, "Degenerate local tetrahedron"
    cells[d < 0, 1], cells[d < 0, 2] = cells[d < 0, 2].copy(), cells[d < 0, 1].copy()
    return cells


def insert(points, cells, query):
    """Stellar subdivision of the point's complete carrier star."""
    import numpy as np

    distances = np.linalg.norm(points-query, axis=1)
    if distances.min() < 1e-11:
        return points, cells, int(distances.argmin())
    hits, b = containing(points, cells, query, 1e-8)
    assert len(hits), "Insertion escaped the fixed background patch"
    carrier = cells[hits[0]][b[0] > 1e-8]
    assert 2 <= len(carrier) <= 4
    incident = np.flatnonzero(np.isin(cells, carrier).sum(1) == len(carrier))
    new_id = len(points)
    result = []
    for tet in cells[incident]:
        for node in carrier:
            result.append(np.where(tet == node, new_id, tet))
    points = np.vstack((points, query))
    cells = np.vstack((np.delete(cells, incident, axis=0), oriented(points, result)))
    assert len(points) <= 256 and len(cells) <= 2000, "Local insertion budget exceeded"
    return points, cells, new_id


def void_margin(case, start, end):
    """A conservative minimum from linear interpolation and an analytic f'' bound.

    PMY has |H_ii| <= 7*(2*pi)^2 and |H_ij| <= 4*(2*pi)^2.
    Linear interpolation differs from f by at most max|f''|*dt^2/8.
    """
    import numpy as np

    delta = end-start
    n = 256
    probes = start + np.linspace(0, 1, n+1)[:, None]*delta
    second = (2*np.pi)**2*(7*float(delta @ delta)
                         + 8*sum(abs(delta[a]*delta[b]) for a, b in ((0, 1), (1, 2), (2, 0))))
    observed = float(np.min(case.function(*probes.T)-case.levels[1]))
    return dict(sampled_minimum=observed, second_derivative_bound=second,
                conservative_minimum=observed-second/(8*n*n)-1e-12)


def endpoints(case, points, critical, axis):
    import numpy as np

    delta = points-critical
    distance = np.linalg.norm(delta, axis=1)
    projection = delta @ axis
    outside = case.function(*points.T) > case.levels[1] + 1e-6
    selected = []
    for sign in (-1, 1):
        candidates = np.flatnonzero(outside & (sign*projection > .01) & (distance < .12))
        for index in candidates[np.argsort(distance[candidates], kind="stable")]:
            margin = void_margin(case, critical, points[index])
            if margin["conservative_minimum"] > 0:
                selected.append(dict(node=int(index), point=points[index].tolist(), **margin))
                break
        else:
            raise AssertionError("No field-proved original void endpoint")
    return selected


def insert_segment(points, cells, start_id, target):
    """Walk a segment, inserting each crossed face and connecting consecutive points."""
    import numpy as np

    start = points[start_id].copy()
    delta = target-start
    along = 0.
    chain = [start_id]
    for _ in range(48):
        after = min(1., along + 1e-7)
        hits, _ = containing(points, cells, start+after*delta, 1e-10)
        assert len(hits), "Void path left the fixed local patch"
        tet = cells[hits[0]:hits[0]+1]
        b0 = barycentric(points, tet, start)[0]
        b1 = barycentric(points, tet, target)[0]
        slope = b1-b0
        times = [-b0[j]/slope[j] for j in range(4) if slope[j] < -1e-12
                 and -b0[j]/slope[j] > along+1e-8]
        stop = min([1., *times])
        assert stop > along+1e-9
        query = target if stop >= 1-1e-9 else start+stop*delta
        points, cells, new_id = insert(points, cells, query)
        edge = (np.isin(cells, [chain[-1], new_id]).sum(1) == 2).any()
        assert edge, "Stellar walk did not leave an explicit edge path"
        chain.append(new_id)
        if stop >= 1-1e-9:
            return points, cells, chain
        along = stop
    raise AssertionError("Void path exceeded 48 insertions")


def clip(case, points, cells):
    """Small replay of native edge-root clipping and pulling tetrahedralization.

    Roots use the unchanged physical function. Only sign-changing edges are
    cut, matching the production Builder's sampling limitation. No smoothing,
    snapping, quality repair, or curvature is included in this stage control.
    """
    import numpy as np
    from scipy.optimize import brentq

    p = [v.copy() for v in points]
    values = list(map(float, case.function(*points.T)))
    cache = {}

    def rank(node):
        return tuple(np.rint(np.array(p[node])*1e11).astype(np.int64))

    def cut(a, b, level):
        if abs(values[a]-level) < 1e-12:
            return a
        if abs(values[b]-level) < 1e-12:
            return b
        a, b = sorted((a, b))
        key = (a, b, level)
        if key not in cache:
            root = brentq(lambda t: float(case.function(*(p[a]+t*(p[b]-p[a])))-level), 0, 1, xtol=1e-14)
            cache[key] = len(p)
            p.append(p[a]+root*(p[b]-p[a]))
            values.append(level)
        return cache[key]

    def halfspace(faces, level, sign):
        out, edges = [], []
        for face in faces:
            poly, cuts = [], []
            for a, b in zip(face, face[1:]+face[:1], strict=True):
                ia, ib = sign*(values[a]-level) <= 1e-12, sign*(values[b]-level) <= 1e-12
                if ia:
                    poly.append(a)
                if ia != ib:
                    node = cut(a, b, level)
                    poly.append(node)
                    cuts.append(node)
            poly = [node for i, node in enumerate(poly) if not i or poly[i-1] != node]
            if len(poly) > 1 and poly[0] == poly[-1]:
                poly.pop()
            if len(poly) >= 3:
                out.append(poly)
            cuts = sorted(set(cuts))
            if len(cuts) == 2:
                edges.append(cuts)
        if len(edges) >= 3:
            ring, previous = [edges[0][0]], -1
            for _ in edges:
                following = next((b if a == ring[-1] else a for a, b in edges
                                  if (a == ring[-1] and b != previous) or (b == ring[-1] and a != previous)), None)
                if following is None or following == ring[0]:
                    break
                previous = ring[-1]
                ring.append(following)
            if len(ring) >= 3:
                out.append(ring)
        return out

    result = []
    for tet in cells:
        faces = [list(map(int, tet[list(face)])) for face in FACES]
        for level, sign in ((case.levels[1], 1), (case.levels[0], -1)):
            if np.isfinite(level):
                faces = halfspace(faces, level, sign)
        if not faces:
            continue
        anchor = min({i for face in faces for i in face}, key=rank)
        for face in faces:
            if anchor in face:
                continue
            pivot = min(range(len(face)), key=lambda j: rank(face[j]))
            face = face[pivot:]+face[:pivot]
            for j in range(1, len(face)-1):
                candidate = [anchor, face[0], face[j], face[j+1]]
                xyz = np.asarray([p[i] for i in candidate])
                if abs(np.linalg.det(xyz[1:]-xyz[:1])) > 1e-18:
                    result.append(candidate)
    p = np.asarray(p)
    cells = oriented(p, result)
    used, inverse = np.unique(cells, return_inverse=True)
    return p[used], inverse.reshape(-1, 4)


def measure(points, cells):
    import numpy as np
    from bistellar_probe import oriented_boundary
    from mmg_quality_localization import quality

    q = quality(points, cells)
    d = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
    boundary = oriented_boundary(cells)
    counts = np.unique(np.sort(cells[:, np.asarray(FACES)].reshape(-1, 3), axis=1), axis=0, return_counts=True)[1]
    return dict(nodes=len(points), tets=len(cells), minimum_quality=float(q.min()),
                below_quality_gate=int((q < .1).sum()), minimum_determinant=float(d.min()), volume=float(d.sum()/6),
                oriented_faces_consistent=boundary is not None, nonmanifold_faces=int((counts > 2).sum()),
                within_caps=bool(max(len(points), len(cells)) <= w.CAP))


def saved_occupancy(source, query):
    import numpy as np
    from scipy.optimize import root

    with np.load(source/"mesh.npz") as mesh:
        points, cells = mesh["points"], mesh["tetrahedra"]
    assert cells.shape[1] == 10
    selected = np.flatnonzero(np.linalg.norm(points[cells[:, :4]].mean(1)-query, axis=1) < .08)
    for index in selected:
        coordinates = points[cells[index]]
        initial = barycentric(points, cells[index:index+1, :4], query)[0, 1:]

        def mapped(uv):
            b = np.r_[1-uv.sum(), uv]
            weights = np.r_[b*(2*b-1), 4*np.prod(b[np.asarray(EDGES)], axis=1)]
            return weights @ coordinates

        solve = root(lambda uv: mapped(uv)-query, initial)
        b = np.r_[1-solve.x.sum(), solve.x]
        residual = float(np.linalg.norm(mapped(solve.x)-query))
        if b.min() > 1e-8 and residual < 1e-10:
            return dict(element=int(index), barycentric=b.tolist(), map_residual=residual,
                        point_in_material=True, corner_coordinates=coordinates[:4].tolist())
    return dict(point_in_material=False, searched_tets=len(selected),
                scope="Bounded inverse-map search; absence is not a global exclusion proof")


def controls():
    import numpy as np
    from scorecard import Case
    from bistellar_probe import oriented_boundary

    p = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    t = np.arange(4)[None]
    q = p.mean(0)
    refined, subtets, _ = insert(p, t, q)
    assert oriented_boundary(t) == oriented_boundary(subtets)
    np.testing.assert_allclose(measure(p, t)["volume"], measure(refined, subtets)["volume"], atol=1e-14)
    plane = Case("plane", lambda x, y, z: x+y+z-.5, (-np.inf, 0.))
    cp, ct = clip(plane, p, t)
    np.testing.assert_allclose(measure(cp, ct)["volume"], 1/48, atol=1e-14)
    hidden = Case("hidden-void", lambda x, y, z: .02-(x-.25)**2-(y-.25)**2-(z-.25)**2, (-np.inf, 0.))
    cp, ct = clip(hidden, p, t)
    assert len(containing(cp, ct, q)[0])
    cp, ct = clip(hidden, refined, subtets)
    assert not len(containing(cp, ct, q)[0])
    return dict(stellar_boundary_unchanged=True, stellar_volume_unchanged=True,
                analytic_plane_clip=True, hidden_void_recovered=True)


def worker(a):
    import numpy as np
    import meshers
    from microgen.shape import surface_functions
    from scorecard import gradient
    from bistellar_probe import oriented_boundary
    import tpms_family_fixtures as fixtures

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    fixture = a.fixtures/NAME
    protocol = json.loads((a.fixtures/"protocol.json").read_text())
    frozen = json.loads((a.fixtures/"result.json").read_text())
    known = {Path(p).resolve(): h for p, h in frozen["inputs"].items()}
    files = [Path(__file__), Path(w.__file__), Path(fixtures.__file__), Path(surface_functions.__file__),
             Path(meshers._meshers.__file__), a.source/"mesh.npz", a.source/"result.json",
             a.fixtures/"protocol.json", a.fixtures/"result.json", *fixture.iterdir()]
    files += [w.HERE/n for n in ("scorecard.py", "bistellar_probe.py", "mmg_quality_localization.py")]
    assert all(w.digest(p) == known[p.resolve()] for p in fixture.iterdir())
    expected = next(h for p, h in protocol["source_hashes"].items() if Path(p).name == "surface_functions.py")
    assert w.digest(Path(surface_functions.__file__)) == expected
    inputs = {str(p.resolve()): w.digest(p) for p in files}
    case = fixtures.case(NAME, a.fixtures)
    results = dict(inputs=inputs, cell=CELL, cells=CELLS, unchanged_physical_specimen=True,
                   unchanged_base_grid=True, unchanged_admission_gates=True, admitted=False,
                   scope="Local preclip topology experiment. No global candidate or FEA result.")
    w.save(a.run/"manifest.json", dict(inputs=inputs, timeout_s=180, memory_bytes=w.MEMORY,
           numerical_threads=w.NUMERICAL_THREAD_POLICY, cell=CELL, witness_selection="27 raw-field root starts in declared cell; no admission samples"))
    with Job(parts={"controls": 1, "critical": 1, "background": 2, "path-insertion": 2, "clipping": 3}) as job:
        with job.part("controls") as part:
            results["controls"] = controls()
            part.update(1, "Analytic clipping and subdivision controls passed", controls=4)
        job.result("controls", passed=True)
        with job.part("critical") as part:
            probe = -.5+np.random.default_rng(47).random((256, 3))
            np.testing.assert_allclose(raw_gradient(probe), gradient(case.function, probe), atol=2e-7, rtol=1e-7)
            critical, all_critical = critical_search(case)
            point, axis = np.array(critical["point"]), np.array(critical["positive_axis"])
            results.update(critical=critical, stationary_points=all_critical)
            results["saved_curved_occupancy"] = saved_occupancy(a.source, point)
            assert results["saved_curved_occupancy"]["point_in_material"], "Frozen saved-neck failure did not replay"
            part.update(1, "Found upper-void saddle from field", level_gap=critical["value"]-case.levels[1])
        job.result("critical", level_gap=critical["value"]-case.levels[1])
        with job.part("background") as part:
            # A full solid extracts the production BCC lattice without cutting.
            background = meshers.generate(lambda x, y, z: 0*x-1, bounds=w.BOUNDS, cells=CELLS,
                snap=0., optimize_passes=0, adapt=False, order=1, minimum_quality=.1,
                geometry_tolerance=.002, periodic=(True,)*3, threads=1, max_tetrahedra=w.CAP)
            points, cells = np.asarray(background.points), np.asarray(background.tetrahedra)
            assert max(len(points), len(cells)) <= w.CAP
            lo = -.5+(np.array(CELL)-1)/CELLS
            hi = -.5+(np.array(CELL)+2)/CELLS
            centers = points[cells].mean(1)
            cells = cells[((centers >= lo) & (centers <= hi)).all(1)]
            used, inverse = np.unique(cells, return_inverse=True)
            points, cells = points[used], inverse.reshape(-1, 4)
            assert len(points) < 256 and len(cells) < 2000
            original_points, original_cells = points.copy(), cells.copy()
            hits, bary = containing(points, cells, point)
            assert len(hits)
            values = case.function(*points.T)
            results["containing_raw_tetrahedra"] = [dict(cell=int(t), corners=points[cells[t]].tolist(),
                field_values=values[cells[t]].tolist(), barycentric=b.tolist(),
                interpolated_value=float(b @ values[cells[t]]),
                all_corners_material=bool(((values[cells[t]] >= case.levels[0]) & (values[cells[t]] <= case.levels[1])).all()))
                for t, b in zip(hits, bary, strict=True)]
            ends = endpoints(case, points, point, axis)
            results["void_path_endpoints"] = ends
            results["background_before"] = measure(points, cells)
            np.savez_compressed(a.run/"background-before.npz", points=points, tetrahedra=cells)
            part.update(1, "Located raw lattice and proved physical void path", local_tets=len(cells), local_nodes=len(points))
        job.result("background", tets=len(cells), nodes=len(points))
        with job.part("path-insertion") as part:
            points, cells, first = insert_segment(points, cells, ends[0]["node"], point)
            part.update(.5, "Inserted first half of void path", inserted_nodes=len(points)-len(original_points))
            points, cells, second = insert_segment(points, cells, first[-1], np.array(ends[1]["point"]))
            chain = first+second[1:]
            np.testing.assert_array_equal(points[:len(original_points)], original_points)
            assert oriented_boundary(cells) == oriented_boundary(original_cells)
            before, after = measure(original_points, original_cells), measure(points, cells)
            np.testing.assert_allclose(after["volume"], before["volume"], rtol=1e-12, atol=1e-14)
            assert np.min(case.function(*points[chain].T)) > case.levels[1]
            for x, y in zip(chain, chain[1:]):
                assert (np.isin(cells, [x, y]).sum(1) == 2).any()
            results.update(background_after=after, original_nodes_bitwise_unchanged=True,
                           outer_facets_identical=True, inserted_nodes=len(points)-len(original_points),
                           void_edge_path=chain)
            np.savez_compressed(a.run/"background-after.npz", points=points, tetrahedra=cells, void_edge_path=chain)
            part.update(1, "Void path is a conforming edge chain", inserted_nodes=results["inserted_nodes"], minimum_q=after["minimum_quality"])
        job.result("path-insertion", inserted_nodes=results["inserted_nodes"], minimum_q=after["minimum_quality"])
        w.save(a.run/"partial.json", results)
        with job.part("clipping") as part:
            samples = np.vstack([np.array(end["point"])+np.linspace(0, 1, 129)[:, None]*(point-np.array(end["point"])) for end in ends])
            path_checks = {}
            for i, (label, p, t) in enumerate((("before", original_points, original_cells), ("after", points, cells))):
                cp, ct = clip(case, p, t)
                row = measure(cp, ct)
                occupied = [j for j, query in enumerate(samples) if len(containing(cp, ct, query)[0])]
                path_checks[label] = dict(material_samples=len(occupied), sample_indices=occupied,
                                          critical_point_in_material=bool(len(containing(cp, ct, point)[0])))
                results["clipped_"+label] = row
                np.savez_compressed(a.run/("clipped-"+label+".npz"), points=cp, tetrahedra=ct)
                part.update((i+1)/2, "Clipped unchanged field on local background", tets=len(ct), wrong_void_samples=len(occupied), minimum_q=row["minimum_quality"])
            results["void_path_checks"] = path_checks
            results["raw_clipping_already_misses_void"] = path_checks["before"]["material_samples"] > 0
            results["local_topology_recovered"] = bool(path_checks["before"]["material_samples"] > 0
                and path_checks["after"]["material_samples"] == 0
                and results["clipped_after"]["oriented_faces_consistent"]
                and results["clipped_after"]["nonmanifold_faces"] == 0)
        job.result("clipping", topology_recovered=results["local_topology_recovered"], minimum_q=results["clipped_after"]["minimum_quality"])
    assert all(w.digest(Path(p)) == h for p, h in inputs.items())
    results.update(complete=True, input_hashes_verified=True,
                   limits="Local edge-root clipping replay, not an integrated native change. Root-sampled path exclusion is finite; the physical void path has a separate analytic interpolation bound. Quality, distance, volume, periodicity and count gates still apply to any future global candidate.")
    w.save(a.run/"result.json", results)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--fixtures", type=Path, default=RUNS/"macbook-tpms-family-fixtures-v1")
    p.add_argument("--source", type=Path, default=RUNS/"macbook-tpms-family-remaining-v1"/(NAME+"-native-fit"))
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    for name in ("run", "fixtures", "source"):
        setattr(a, name, getattr(a, name).resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run/"input-pmy_neck_control.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        a.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=180)
    w.save(a.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"


if __name__ == "__main__":
    main()
