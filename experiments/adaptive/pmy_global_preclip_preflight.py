"""Global one-round BCC classification preflight, without clipping or curving.

All original nodes remain fixed. Physical-cap face decisions and new nodes are
paired exactly. Each insertion visits only active descendants of its original
incident tetrahedra. Strict live caps remain 250,000 nodes and tetrahedra.
"""

import argparse
from itertools import permutations, product
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w
import pmy_barycenter_control as local
import pmy_neck_control as neck


NAME = "tpms-pmy-sheet"
RUNS = w.HERE/"runs"
CELLS = 16


class BudgetExceeded(RuntimeError):
    pass


def inventory(cells):
    import numpy as np

    raw = cells[:, np.asarray(neck.FACES)].reshape(-1, 3)
    faces, first, inverse, count = np.unique(np.sort(raw, axis=1), axis=0,
                                            return_index=True, return_inverse=True, return_counts=True)
    assert count.max() <= 2, "Nonmanifold background face"
    parents = np.full((len(faces), 2), -1, dtype=np.int64)
    parent = np.repeat(np.arange(len(cells)), 4)
    slot = (parent != first[inverse]//4).astype(np.int8)
    parents[inverse, slot] = parent
    inversions = ((raw[:, 0] > raw[:, 1]).astype(np.int8)
                  + (raw[:, 0] > raw[:, 2]).astype(np.int8)
                  + (raw[:, 1] > raw[:, 2]).astype(np.int8))
    sums = np.bincount(inverse, weights=1-2*(inversions % 2), minlength=len(faces))
    assert (sums[count == 2] == 0).all(), "Incoherent shared background faces"
    assert (abs(sums[count == 1]) == 1).all()
    return faces, count, parents, raw[first[count == 1]]


def paired_selection(points, faces, count, requested, centers):
    """Synchronize opposite physical triangles; transverse coordinates stay exact."""
    import numpy as np

    selected = requested.copy()
    centers = centers.copy()
    groups = {}
    for index in np.flatnonzero(count == 1):
        xyz = points[faces[index]]
        planes = [(axis, side) for axis in range(3) for side in (-.5, .5)
                  if (xyz[:, axis] == side).all()]
        assert len(planes) == 1, "Boundary face is not on one physical cell plane"
        axis, side = planes[0]
        transverse = [i for i in range(3) if i != axis]
        key = (axis, tuple(sorted(tuple(p[transverse]) for p in xyz)))
        groups.setdefault(key, []).append((int(index), side))
    pairs, forced = [], []
    for (axis, _), members in sorted(groups.items()):
        assert len(members) == 2 and sorted(side for _, side in members) == [-.5, .5]
        lo, hi = [index for index, _ in sorted(members, key=lambda item: item[1])]
        pairs.append((lo, hi))
        if selected[lo] or selected[hi]:
            forced.extend(i for i in (lo, hi) if not selected[i])
            selected[[lo, hi]] = True
            centers[lo, axis] = -.5
            centers[hi] = centers[lo]
            centers[hi, axis] = .5
        np.testing.assert_array_equal(centers[hi]-centers[lo], np.eye(3)[axis])
    operations = [[int(i)] for i in np.flatnonzero(selected & (count == 2))]
    operations += [[lo, hi] for lo, hi in pairs if selected[lo]]
    operations.sort(key=min)
    return selected, centers, operations, pairs, forced


def selection(case, points, cells, faces, count):
    import numpy as np

    values = np.asarray(case.function(*points.T))
    face_centers = points[faces].mean(1)
    tet_centers = points[cells].mean(1)
    face_actual = np.asarray(case.function(*face_centers.T))
    tet_actual = np.asarray(case.function(*tet_centers.T))
    face_affine = values[faces].mean(1)
    tet_affine = values[cells].mean(1)
    raw_faces = local.classify(face_actual, case.levels) != local.classify(face_affine, case.levels)
    selected_tets = local.classify(tet_actual, case.levels) != local.classify(tet_affine, case.levels)
    selected_faces, face_centers, groups, pairs, forced = paired_selection(points, faces, count, raw_faces, face_centers)
    return dict(values=values, face_centers=face_centers, tet_centers=tet_centers,
                face_actual=face_actual, face_affine=face_affine, tet_actual=tet_actual, tet_affine=tet_affine,
                requested_faces=raw_faces, selected_faces=selected_faces, selected_tets=selected_tets,
                groups=groups, pairs=pairs, forced=forced)


def subdivide(points, cells, faces, parents, plan, progress=None):
    import numpy as np

    original_nodes = len(points)
    nodes = points.tolist()
    descendants = [[tuple(map(int, tet))] for tet in cells]
    live = len(cells)
    selected_tets = np.flatnonzero(plan["selected_tets"])
    total = len(plan["groups"])+len(selected_tets)
    completed = 0
    visits = 0
    maximum_descendants = 1

    def check(added_nodes, added_tets):
        if len(nodes)+added_nodes > w.CAP or live+added_tets > w.CAP:
            raise BudgetExceeded(f"Strict background cap: proposed {len(nodes)+added_nodes} nodes / {live+added_tets} tetrahedra")

    for group in plan["groups"]:
        updates = {}
        additions = []
        for face_index in group:
            carrier = set(map(int, faces[face_index]))
            new_id = len(nodes)+len(additions)
            additions.append(plan["face_centers"][face_index].tolist())
            for parent in parents[face_index]:
                if parent < 0:
                    continue
                old = updates.get(int(parent), descendants[parent])
                visits += len(old)
                matched = [i for i, tet in enumerate(old) if carrier.issubset(tet)]
                assert len(matched) == 1, "Original face was lost in its descendant star"
                index = matched[0]
                tet = old[index]
                replacements = [tuple(new_id if node == vertex else node for node in tet)
                                for vertex in sorted(carrier)]
                updates[int(parent)] = old[:index]+old[index+1:]+replacements
        growth = sum(len(new)-len(descendants[parent]) for parent, new in updates.items())
        check(len(additions), growth)
        nodes.extend(additions)
        for parent, new in updates.items():
            descendants[parent] = new
            maximum_descendants = max(maximum_descendants, len(new))
        live += growth
        completed += 1
        if progress is not None:
            progress(completed/max(total, 1), len(nodes), live)

    for parent in selected_tets:
        old = descendants[parent]
        visits += len(old)
        coordinates = np.asarray([[nodes[node] for node in tet] for tet in old])
        query = plan["tet_centers"][parent]
        b = np.linalg.solve((coordinates[:, 1:]-coordinates[:, :1]).transpose(0, 2, 1),
                            (query-coordinates[:, 0])[..., None])[..., 0]
        bary = np.column_stack((1-b.sum(1), b))
        hits = np.flatnonzero(bary.min(1) >= -1e-9)
        assert len(hits), "Tet center escaped its own original descendants"
        carrier = {old[hits[0]][j] for j in range(4) if bary[hits[0], j] > 1e-9}
        assert 2 <= len(carrier) <= 4
        new_id = len(nodes)
        replacement = []
        for tet in old:
            if carrier.issubset(tet):
                replacement.extend(tuple(new_id if node == vertex else node for node in tet)
                                   for vertex in sorted(carrier))
            else:
                replacement.append(tet)
        growth = len(replacement)-len(old)
        check(1, growth)
        nodes.append(query.tolist())
        descendants[parent] = replacement
        maximum_descendants = max(maximum_descendants, len(replacement))
        live += growth
        completed += 1
        if progress is not None:
            progress(completed/max(total, 1), len(nodes), live)
    result = np.array(nodes)
    np.testing.assert_array_equal(result[:original_nodes], points)
    final_cells = np.array([tet for children in descendants for tet in children], dtype=np.int64)
    owner = np.repeat(np.arange(len(cells)), [len(children) for children in descendants])
    assert len(final_cells) == live
    assert len(result) == original_nodes+int(plan["selected_faces"].sum())+len(selected_tets)
    return result, final_cells, owner, descendants, dict(
        original_nodes_bitwise_unchanged=True, local_descendant_visits=visits,
        maximum_descendants_per_original_tet=maximum_descendants,
        global_searches_during_insertion=0, peak_live_nodes=len(result), peak_live_tets=live)


def metrics(points, cells):
    import numpy as np
    from mmg_quality_localization import quality
    from periodic_trace import periodic_pairs

    faces, count, _, boundary = inventory(cells)
    q = quality(points, cells)
    det = np.linalg.det(points[cells[:, 1:]]-points[cells[:, :1]])
    assert det.min() > 0 and len(np.unique(cells)) == len(points)
    assert (points >= -.5).all() and (points <= .5).all()
    _, _, periodic = periodic_pairs(points, boundary, declared_bounds=([-.5]*3, [.5]*3))
    assert periodic["periodic_coordinate_errors"] == [0., 0., 0.]
    xyz = points[boundary]
    areas = np.linalg.norm(np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0]), axis=1)/2
    caps = []
    for axis in range(3):
        for side in (-.5, .5):
            mask = (xyz[:, :, axis] == side).all(1)
            area = float(areas[mask].sum())
            np.testing.assert_allclose(area, 1., atol=1e-12, rtol=0)
            caps.append(dict(axis=axis, side=side, faces=int(mask.sum()), area=area))
    np.testing.assert_allclose(det.sum()/6, 1., atol=1e-12, rtol=0)
    return dict(nodes=len(points), tets=len(cells), boundary_faces=len(boundary), unique_faces=len(faces),
                minimum_quality=float(q.min()), below_quality_gate=int((q < .1).sum()),
                quality_floor_met=bool(q.min() >= .1), minimum_determinant=float(det.min()),
                volume=float(det.sum()/6), coherent_oriented_faces=True, nonmanifold_faces=0,
                exact_periodic_nodes_and_facets=True, periodic=periodic, caps=caps,
                within_caps=bool(max(len(points), len(cells)) <= w.CAP)), det, boundary


def interpolate(points, cells, values, query):
    import numpy as np

    b = neck.barycentric(points, cells, query)
    selected = b.min(1) >= -1e-9
    assert selected.any(), "Observer point is outside selected descendants"
    estimates = np.einsum("ij,ij->i", b[selected], values[cells[selected]])
    assert np.ptp(estimates) < 1e-9, "Field interpolation is discontinuous across shared faces"
    return float(estimates.mean())


def observe(source, case, before, after, descendants):
    import numpy as np
    from scipy.spatial import cKDTree

    frozen = json.loads((source/"result.json").read_text())
    critical = np.array(frozen["critical"]["point"])
    ends = frozen["void_path_endpoints"]
    queries = np.vstack([np.array(end["point"])+np.linspace(0, 1, 129)[:, None]*(critical-np.array(end["point"])) for end in ends])
    margins = [neck.void_margin(case, critical, np.array(end["point"])) for end in ends]
    assert all(row["conservative_minimum"] > 0 for row in margins)
    original_points, original_cells = before
    new_points = after
    centers = original_points[original_cells].mean(1)
    radius = float(np.linalg.norm(original_points[original_cells]-centers[:, None], axis=2).max())+1e-12
    tree = cKDTree(centers)
    old_values, new_values = case.function(*original_points.T), case.function(*new_points.T)
    observed = {"baseline": [], "one_round": []}
    for query in queries:
        candidates = np.array(tree.query_ball_point(query, radius), dtype=np.int64)
        bary = neck.barycentric(original_points, original_cells[candidates], query)
        parents = candidates[bary.min(1) >= -1e-9]
        assert len(parents)
        children = np.array([tet for parent in parents for tet in descendants[parent]], dtype=np.int64)
        observed["baseline"].append(interpolate(original_points, original_cells[parents], old_values, query))
        observed["one_round"].append(interpolate(new_points, children, new_values, query))
    result = {}
    for name, values in observed.items():
        classes = local.classify(values, case.levels)
        result[name] = dict(samples=len(queries), falsely_material=int((classes == 0).sum()),
                            wrongly_below_lower_wall=int((classes == -1).sum()),
                            true_upper_void_class=int((classes == 1).sum()),
                            critical_affine_value=values[128], critical_affine_class=int(classes[128]))
    result.update(physical_void_margins=margins,
                  scope="Affine background classification only. These counts are not clipped-mesh occupancy or admission evidence.")
    return result, dict(queries=queries, actual=np.asarray(case.function(*queries.T)), **{k: np.asarray(v) for k, v in observed.items()})


def controls():
    import numpy as np
    from scorecard import Case

    bits = list(product((0, 1), repeat=3))
    lookup = {b: i for i, b in enumerate(bits)}
    points = np.array(bits, dtype=float)-.5
    cells = []
    for order in permutations(range(3)):
        p = [0, 0, 0]
        tet = [lookup[tuple(p)]]
        for axis in order:
            p[axis] = 1
            tet.append(lookup[tuple(p)])
        cells.append(tet)
    cells = neck.oriented(points, cells)
    faces, count, parents, _ = inventory(cells)
    case = Case("periodic-cosines", lambda x, y, z: np.cos(2*np.pi*x)+np.cos(2*np.pi*y)+np.cos(2*np.pi*z), (-.4, .4))
    plan = selection(case, points, cells, faces, count)
    mask = np.zeros(len(faces), dtype=bool)
    lo, hi = plan["pairs"][0]
    mask[lo] = True
    paired, _, _, _, forced = paired_selection(points, faces, count, mask, points[faces].mean(1))
    assert paired[lo] and paired[hi] and forced == [hi]
    refined, tets, owners, _, stats = subdivide(points, cells, faces, parents, plan)
    before, determinants, _ = metrics(points, cells)
    after, refined_determinants, _ = metrics(refined, tets)
    np.testing.assert_allclose(np.bincount(owners, weights=refined_determinants), determinants, atol=1e-14, rtol=0)
    assert stats["global_searches_during_insertion"] == 0 and after["volume"] == before["volume"]
    return dict(asymmetric_periodic_selection_paired=True, exact_periodic_new_nodes_and_faces=True,
                coherent_descendant_subdivision=True, per_parent_volume_preserved=True,
                original_nodes_unchanged=True, cap_area_preserved=True)


def worker(a):
    import numpy as np
    import meshers
    from microgen.shape import surface_functions
    import tpms_family_fixtures as fixtures

    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    frozen = json.loads((a.fixtures/"result.json").read_text())
    protocol = json.loads((a.fixtures/"protocol.json").read_text())
    known = {Path(p).resolve(): h for p, h in frozen["inputs"].items()}
    fixture = a.fixtures/NAME
    assert all(w.digest(p) == known[p.resolve()] for p in fixture.iterdir())
    expected = next(h for p, h in protocol["source_hashes"].items() if Path(p).name == "surface_functions.py")
    assert w.digest(Path(surface_functions.__file__)) == expected
    files = [Path(__file__), Path(w.__file__), Path(neck.__file__), Path(local.__file__), Path(fixtures.__file__),
             Path(surface_functions.__file__), Path(meshers._meshers.__file__), Path(meshers.__file__),
             a.fixtures/"result.json", a.fixtures/"protocol.json", a.observers/"result.json", *fixture.iterdir()]
    files += [w.HERE/n for n in ("scorecard.py", "periodic_trace.py", "mmg_quality_localization.py")]
    hashes = {str(p.resolve()): w.digest(p) for p in files}
    case = fixtures.case(NAME, a.fixtures)
    result = dict(inputs=hashes, cells=CELLS, rounds=1, periodic=[True]*3,
                  physical_specimen_unchanged=True, original_grid_and_phase_unchanged=True,
                  admission_gates_unchanged=True, admitted=False, clipping_performed=False, solve_performed=False)
    w.save(a.run/"manifest.json", dict(inputs=hashes, timeout_s=180, memory_bytes=w.MEMORY,
           numerical_threads=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
           selection="All original face/tet center three-way classification mismatches; exact paired physical boundary selection"))
    with Job(parts={"controls": 1, "background-selection": 2, "subdivision": 4, "geometry-audit": 3, "observers": 1}) as job:
        with job.part("controls") as part:
            result["controls"] = controls()
            part.update(1, "Periodic subdivision controls passed", controls=len(result["controls"]))
        job.result("controls", passed=True)
        with job.part("background-selection") as part:
            mesh = meshers.generate(lambda x, y, z: 0*x-1, bounds=w.BOUNDS, cells=CELLS,
                snap=0., optimize_passes=0, adapt=False, order=1, minimum_quality=.1,
                geometry_tolerance=.002, periodic=(True,)*3, threads=1, max_tetrahedra=w.CAP)
            original_points, original_cells = np.asarray(mesh.points), np.asarray(mesh.tetrahedra)
            assert max(len(original_points), len(original_cells)) <= w.CAP
            faces, count, parents, _ = inventory(original_cells)
            plan = selection(case, original_points, original_cells, faces, count)
            result["selection"] = dict(face_samples=len(faces), tet_samples=len(original_cells),
                requested_faces=int(plan["requested_faces"].sum()), selected_faces=int(plan["selected_faces"].sum()),
                selected_internal_faces=int((plan["selected_faces"] & (count == 2)).sum()),
                selected_physical_boundary_faces=int((plan["selected_faces"] & (count == 1)).sum()),
                selected_tets=int(plan["selected_tets"].sum()), paired_boundary_classes=len(plan["pairs"]),
                forced_pair_members=plan["forced"], frozen_before_insertion=True)
            np.savez_compressed(a.run/"selection.npz", faces=faces, face_counts=count, face_parents=parents,
                **{k: v for k, v in plan.items() if isinstance(v, np.ndarray)})
            np.savez_compressed(a.run/"background-before.npz", points=original_points, tetrahedra=original_cells)
            result["selection_sha256"] = w.digest(a.run/"selection.npz")
            part.update(1, "Frozen global classification plan", faces=result["selection"]["selected_faces"], tets=result["selection"]["selected_tets"])
        job.result("background-selection", selected_faces=result["selection"]["selected_faces"], selected_tets=result["selection"]["selected_tets"])
        with job.part("subdivision") as part:
            try:
                points, cells, owner, descendants, stats = subdivide(original_points, original_cells, faces, parents, plan,
                    progress=lambda fraction, nodes, tets: part.update(fraction, "Updating original-tet descendant stars", nodes=nodes, tets=tets))
            except BudgetExceeded as error:
                result.update(complete=True, status="resource-rejected", error=str(error), preflight_passed=False)
                w.save(a.run/"result.json", result)
                return
            result["subdivision"] = stats
            np.savez_compressed(a.run/"background-one-round.npz", points=points, tetrahedra=cells, original_tet=owner)
            part.update(1, "Completed strict-cap one-round subdivision", nodes=len(points), tets=len(cells))
        job.result("subdivision", nodes=len(points), tets=len(cells), maximum_descendants=stats["maximum_descendants_per_original_tet"])
        with job.part("geometry-audit") as part:
            before, old_det, _ = metrics(original_points, original_cells)
            part.update(.4, "Audited original BCC geometry", minimum_q=before["minimum_quality"])
            after, new_det, surface = metrics(points, cells)
            per_parent = np.bincount(owner, weights=new_det, minlength=len(original_cells))
            relative = float(np.max(abs(per_parent-old_det)/old_det))
            assert relative < 1e-10
            result.update(background=dict(baseline=before, one_round=after), maximum_parent_volume_error=relative,
                          preflight_passed=bool(after["quality_floor_met"] and after["within_caps"]))
            np.savez_compressed(a.run/"boundary-one-round.npz", surface=surface)
            values = case.function(*points.T)
            result["field_node_class_counts"] = {str(c): int((local.classify(values, case.levels) == c).sum()) for c in (-1, 0, 1)}
            part.update(1, "Audited exact periodic boundary and every parent volume", minimum_q=after["minimum_quality"], below_q=after["below_quality_gate"])
        job.result("geometry-audit", minimum_q=after["minimum_quality"], periodic_error=max(after["periodic"]["periodic_coordinate_errors"]), parent_volume_error=relative)
        w.save(a.run/"partial.json", result)
        with job.part("observers") as part:
            observation, arrays = observe(a.observers, case, (original_points, original_cells), points, descendants)
            result["path_affine_observer"] = observation
            np.savez_compressed(a.run/"path-observer.npz", **arrays)
            part.update(1, "Observed preclip classifications after all insertion decisions", before=observation["baseline"]["falsely_material"], after=observation["one_round"]["falsely_material"])
        job.result("observers", falsely_material=observation["one_round"]["falsely_material"])
    assert w.digest(a.run/"selection.npz") == result["selection_sha256"]
    assert all(w.digest(Path(p)) == h for p, h in hashes.items())
    result.update(complete=True, status="preflight-complete", input_hashes_verified=True,
                  limits="Background feasibility only, before original-node snapping and root clipping. No clipped quality, physical boundary accuracy, admission, or solve is established. Finite barycenter tests can miss other features.")
    w.save(a.run/"result.json", result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--fixtures", type=Path, default=RUNS/"macbook-tpms-family-fixtures-v1")
    p.add_argument("--observers", type=Path, default=RUNS/"macbook-pmy-neck-control-v1")
    p.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    a = p.parse_args()
    for name in ("run", "fixtures", "observers"):
        setattr(a, name, getattr(a, name).resolve())
    if a.worker:
        return worker(a)
    a.run.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, a.run/"input-pmy_global_preclip_preflight.py")
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        a.run/"worker.log", SimpleNamespace(update=lambda *args, **kwargs: None), timeout=180)
    w.save(a.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"


if __name__ == "__main__":
    main()
