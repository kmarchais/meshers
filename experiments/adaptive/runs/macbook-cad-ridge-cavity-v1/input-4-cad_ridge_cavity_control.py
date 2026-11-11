"""Two conforming edge-star split orders in one saved CAD ridge cavity.

The frozen worst witness selects a diagnostic cavity only. Splitting first
reproduces each parent Tet10 map exactly. New ridge corners and wall midpoints
then use the analytic branch constraints, without observer-driven optimization.
No production change, parameter sweep, FEA or generation timing is claimed.
"""

import argparse
from collections import Counter
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w
import cad_intersection_probe as probe
import cad_intersection_curved_control as curved
import cad_union_comparison as base


SOURCE = w.HERE / "runs/macbook-cad-intersection-curved-v1"


def tet_shape(bary):
    import numpy as np
    from nested_refinement import EDGES

    bary = np.asarray(bary)
    return np.concatenate((bary*(2*bary-1), 4*bary[..., EDGES[:, 0]]*bary[..., EDGES[:, 1]]), axis=-1)


def midpoint_map(cells):
    from nested_refinement import EDGES

    mapping = {}
    for tet in cells:
        for index, (a, b) in enumerate(EDGES):
            edge = tuple(sorted((int(tet[a]), int(tet[b]))))
            mid = int(tet[4+index])
            assert mapping.setdefault(edge, mid) == mid, "Quadratic midpoint crack"
    return mapping


def boundary(cells):
    import numpy as np
    from nested_refinement import FACES

    mids = midpoint_map(cells)
    oriented = cells[:, FACES].reshape(-1, 3)
    _, first, counts = np.unique(np.sort(oriented, axis=1), axis=0, return_index=True, return_counts=True)
    assert counts.max() <= 2, "Nonmanifold tetrahedron faces"
    faces = oriented[first[counts == 1]]
    return np.column_stack((faces, np.asarray([[mids[tuple(sorted((int(a), int(b))))]
        for a, b in zip(face, face[[1, 2, 0]], strict=True)] for face in faces], dtype=int)))


def boundary_keys(faces):
    return {tuple(sorted(map(int, face[:3]))): tuple(sorted(map(int, face))) for face in faces}


def split_parameter(points, edge, mid, plane_sign):
    import numpy as np

    a, b = map(int, edge)
    values = points[[a, mid, b], 0] - plane_sign*points[[a, mid, b], 1]
    va, vm, vb = values
    coefficients = [2*(va+vb-2*vm), 4*vm-3*va-vb, va]
    roots = np.roots(coefficients if abs(coefficients[0]) > 1e-16 else coefficients[1:])
    eligible = [float(root.real) for root in roots if abs(root.imag) < 1e-10 and 1e-8 < root.real < 1-1e-8]
    assert len(eligible) == 1, f"Expected one exposed ridge crossing, found {eligible}"
    t = eligible[0]
    weights = np.array([(1-t)*(1-2*t), 4*t*(1-t), t*(2*t-1)])
    point = weights @ points[[a, mid, b]]
    assert abs(point[0]-plane_sign*point[1]) < 1e-12
    return t, point


def split_star(points, cells, owners, edge, parameter):
    """Conforming subdivision, with every new node evaluated in its parent map."""
    import numpy as np
    from nested_refinement import EDGES
    from quadratic_edge_removal_probe import measure

    a, b = map(int, edge)
    selected = np.flatnonzero((cells[:, :4] == a).any(1) & (cells[:, :4] == b).any(1))
    assert len(selected), "Selected edge has no incident tetrahedra"
    mids = midpoint_map(cells)
    old_mid = mids[tuple(sorted((a, b)))]
    t = parameter
    position = np.array([(1-t)*(1-2*t), 4*t*(1-t), t*(2*t-1)]) @ points[[a, old_mid, b]]
    node = len(points)
    nodes = points.tolist()+[position.tolist()]
    children, child_owners = [], []
    reproduction = 0.
    evaluation = np.vstack((np.eye(4), np.full((1, 4), .25), np.array([[.1, .2, .3, .4], [.4, .1, .3, .2]])))
    for parent in selected:
        tet = cells[parent]
        corner = tet[:4]
        ia, ib = int(np.flatnonzero(corner == a)[0]), int(np.flatnonzero(corner == b)[0])
        inserted_bary = (1-t)*np.eye(4)[ia]+t*np.eye(4)[ib]
        for replaced in (ia, ib):
            corners = corner.copy()
            corners[replaced] = node
            mapping = np.eye(4)
            mapping[replaced] = inserted_bary
            child = list(map(int, corners))
            for u, v in EDGES:
                key = tuple(sorted((int(corners[u]), int(corners[v]))))
                expected = tet_shape((mapping[u]+mapping[v])/2) @ points[tet]
                if key not in mids:
                    mids[key] = len(nodes)
                    nodes.append(expected.tolist())
                else:
                    np.testing.assert_allclose(nodes[mids[key]], expected, atol=2e-14, rtol=0)
                child.append(mids[key])
            xyz = np.asarray([nodes[index] for index in child])
            error = abs(tet_shape(evaluation) @ xyz-tet_shape(evaluation @ mapping) @ points[tet]).max()
            reproduction = max(reproduction, float(error))
            children.append(child)
            child_owners.append(int(owners[parent]))
    nodes = np.asarray(nodes)
    children = np.asarray(children)
    _, before_volume = measure(points[cells[selected]], positive_rule=False)
    _, after_volume = measure(nodes[children], positive_rule=False)
    per_parent = after_volume.reshape(-1, 2).sum(1)
    relative = float(np.max(abs(per_parent-before_volume)/abs(before_volume)))
    assert reproduction < 1e-12 and relative < 1e-10
    keep = np.ones(len(cells), dtype=bool)
    keep[selected] = False
    combined = np.vstack((cells[keep], children))
    new_owners = np.r_[owners[keep], child_owners]
    assert max(len(nodes), len(combined)) <= w.CAP
    return nodes, combined, new_owners, node, dict(edge=[a, b], parameter=float(t),
        original_edge_midpoint=old_mid, new_corner=node, incident_tets=len(selected),
        original_parents=sorted(set(map(int, owners[selected]))),
        map_reproduction_max=reproduction, maximum_parent_volume_relative_error=relative)


def nearest_straddling_face(points, surface, query):
    import numpy as np
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

    chosen = np.flatnonzero(probe.face_groups(points, surface)["crossing"])
    vertices, triangles = w.tessellate(points, surface[chosen], 16)
    poly = vtk.vtkPolyData()
    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(vertices), deep=True))
    poly.SetPoints(pts)
    faces = vtk.vtkCellArray()
    flat = np.c_[np.full(len(triangles), 3), triangles].astype(np.int64).ravel()
    faces.SetCells(len(triangles), numpy_to_vtkIdTypeArray(flat, deep=True))
    poly.SetPolys(faces)
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(poly)
    locator.BuildLocator()
    point = [0., 0., 0.]
    cell, sub, d2 = vtk.reference(0), vtk.reference(0), vtk.reference(0.)
    locator.FindClosestPoint(query, point, cell, sub, d2)
    face = int(chosen[int(cell)//256])
    return face, dict(face=face, query=list(map(float, query)),
        closest_tessellated_point=point, distance=float(np.sqrt(d2)), subdivisions=16,
        search_scope="All original branch-straddling Tri6 boundary faces")


def branch_projection(points, cells, surface, new_corners, original_point_count):
    """Move new corners affinely into their stars, then fit new wall edges."""
    import numpy as np

    moved = points.copy()
    mids = midpoint_map(cells)
    targets = probe.ridge_target(points[new_corners])
    displacements = {node: target-points[node] for node, target in zip(new_corners, targets, strict=True)}
    for node, delta in displacements.items():
        moved[node] += delta
    for (a, b), mid in mids.items():
        moved[mid] += .5*(displacements.get(a, 0)+displacements.get(b, 0))
    # Preserve all old nodes, including the complete outer cavity interface.
    np.testing.assert_array_equal(moved[:original_point_count], points[:original_point_count])
    boundary_edges = set()
    crossing = []
    for index in np.flatnonzero(probe.wall_mask(moved, surface)):
        face = surface[index]
        signs = probe.branch_value(moved[face[:3]])
        if signs.min() < -1e-12 and signs.max() > 1e-12:
            crossing.append(int(index))
        for k in range(3):
            key = tuple(sorted((int(face[k]), int(face[(k+1)%3]))))
            boundary_edges.add(key)
    rows = []
    for edge in sorted(boundary_edges):
        mid = mids[edge]
        if mid < original_point_count:
            continue
        signs = probe.branch_value(moved[list(edge)])
        assert not (signs.min() < -1e-12 and signs.max() > 1e-12), "A new boundary edge still straddles the ridge"
        branch = 0 if signs.max() > 1e-12 else 1 if signs.min() < -1e-12 else 2
        row = dict(edge=list(edge), midpoint=mid, branch=branch)
        candidate = moved[mid].copy()
        if branch == 2:
            candidate = probe.ridge_target(candidate[None])[0]
            row["constraint"] = "both-cylinder-walls"
        else:
            axes = [1, 2] if branch == 0 else [0, 2]
            candidate[axes] *= base.fixture.RADIUS/np.linalg.norm(candidate[axes])
            row["constraint"] = "x-axis-cylinder" if branch == 0 else "y-axis-cylinder"
        fields = np.array([np.hypot(candidate[1], candidate[2]), np.hypot(candidate[0], candidate[2])])-base.fixture.RADIUS
        row.update(displacement=float(np.linalg.norm(candidate-moved[mid])), fields=fields.tolist(),
                   exposed=bool(fields.min() >= -1e-12 and abs(fields.min()) < 1e-10))
        # Keep a failed candidate intact for diagnosis; do not silently switch branches.
        moved[mid] = candidate
        rows.append(row)
    return moved, dict(corner_displacements={str(k): v.tolist() for k, v in displacements.items()},
        new_wall_midpoints=rows, all_new_wall_midpoints_exposed=all(row["exposed"] for row in rows),
        original_nodes_bitwise_unchanged=bool(np.array_equal(moved[:original_point_count], points[:original_point_count])),
        remaining_crossing_faces=crossing)


def compact(points, cells, surface):
    import numpy as np

    used = np.unique(cells)
    assert np.isin(surface, used).all()
    mapping = np.full(len(points), -1, dtype=int)
    mapping[used] = np.arange(len(used))
    return points[used], mapping[cells], mapping[surface], mapping


def candidate(original_points, original_cells, original_surface, edge_order, plane_sign, query, output, part):
    import numpy as np
    from nested_refinement import positive_quality
    from quadratic_edge_removal_probe import measure
    from facet_incidence_audit import audit
    from comparison_mesh_audit import quadratic_conformity
    from periodic_trace import periodic_pairs

    points, cells = original_points.copy(), original_cells.copy()
    owners = np.arange(len(cells))
    operations, new_corners = [], []
    row = dict(edge_order=[list(edge) for edge in edge_order], accepted=False, rejection_reasons=[])
    try:
        for index, edge in enumerate(edge_order):
            mid = midpoint_map(cells)[tuple(edge)]
            parameter, _ = split_parameter(points, edge, mid, plane_sign)
            points, cells, owners, node, operation = split_star(points, cells, owners, edge, parameter)
            operations.append(operation)
            new_corners.append(node)
            part.update(.15+.15*index, "Subdivided complete curved edge star", incident_tets=operation["incident_tets"])
        parents = np.unique([parent for op in operations for parent in op["original_parents"]])
        local = np.isin(owners, parents)
        original_cavity_nodes = np.unique(original_cells[parents])
        cap_clearance = float(.5-abs(original_points[original_cavity_nodes]).max())
        assert cap_clearance > 1e-10, "Cavity reaches a periodic cap; isolated operation rejected"
        surface = boundary(cells)
        before_interface = boundary_keys(boundary(original_cells[parents]))
        before_physical = boundary_keys(original_surface)
        before_interface = {key: value for key, value in before_interface.items() if key not in before_physical}
        after_interface = boundary_keys(boundary(cells[local]))
        after_physical = boundary_keys(surface)
        after_interface = {key: value for key, value in after_interface.items() if key not in after_physical}
        assert before_interface == after_interface, "Complete quadratic cavity interface differs"
        interface_nodes = np.unique([node for face in before_interface.values() for node in face])
        np.testing.assert_array_equal(points[interface_nodes], original_points[interface_nodes])
        before_deformation_quality = positive_quality(points[cells[local, :4]])
        points, projection = branch_projection(points, cells, surface, new_corners, len(original_points))
        np.testing.assert_array_equal(points[interface_nodes], original_points[interface_nodes])
        quality = positive_quality(points[cells[local, :4]])
        minimum, volumes = measure(points[cells[local]])
        old_minimum, old_volumes = measure(original_points[original_cells[parents]])
        original_q = positive_quality(original_points[original_cells[parents, :4]])
        original_global_q = float(positive_quality(original_points[original_cells[:, :4]]).min())
        new_ridge = tuple(sorted(new_corners))
        assert new_ridge in midpoint_map(cells), "The new intersection is not a mesh edge"
        incident_wall = [i for i, face in enumerate(surface) if set(new_ridge) <= set(face[:3])]
        assert len(incident_wall) == 2, "Ridge edge must have two exposed boundary faces"
        local_boundary_keys = set(boundary_keys(boundary(cells[local])))
        local_wall_ids = np.array([i for i, face in enumerate(surface)
            if tuple(sorted(map(int, face[:3]))) in local_boundary_keys], dtype=int)
        local_crossing = [i for i in projection["remaining_crossing_faces"] if i in set(local_wall_ids)]
        from scorecard import distance_to_mesh
        local_vertices, local_faces = w.tessellate(points, surface[local_wall_ids], 16)
        local_witness = float(distance_to_mesh(local_vertices, local_faces, query[None])[0])
        if quality.min() < .1:
            row["rejection_reasons"].append("corner-quality-below-0.1")
        if minimum.min() <= 0:
            row["rejection_reasons"].append("nonpositive-sampled-curved-Jacobian")
        if not projection["all_new_wall_midpoints_exposed"]:
            row["rejection_reasons"].append("new-wall-midpoint-is-hidden-inside-another-union-branch")
        compact_points, compact_cells, compact_surface, mapping = compact(points, cells, surface)
        row["facet_audit"] = audit(compact_points, compact_cells, compact_surface)
        row["quadratic_conformity"] = quadratic_conformity(compact_points, compact_cells, compact_surface)
        assert row["quadratic_conformity"]["passed"]
        _, _, row["periodic"] = periodic_pairs(compact_points, compact_surface, declared_bounds=([-.5]*3, [.5]*3))
        if not row["facet_audit"]["oriented_facet_consistency_passed"]:
            row["rejection_reasons"].append("facet-incidence")
        if max(len(compact_points), len(compact_cells)) > w.CAP:
            row["rejection_reasons"].append("node-or-tetrahedron-cap")
        np.savez_compressed(output/"mesh.npz", points=compact_points, tetrahedra=compact_cells, surface=compact_surface)
        np.savez_compressed(output/"cavity.npz", original_parent_ids=parents, original_cells=original_cells[parents],
            after_cells=cells[local], before_points=original_points, after_points=points,
            interface_nodes=interface_nodes, new_ridge_nodes=new_corners)
        qkeys = ("min", "p01", "p10", "median", "p90", "p99", "max")
        quantiles = [0, .01, .1, .5, .9, .99, 1]
        row.update(operations=operations, projection=projection, original_parent_count=len(parents),
            final_local_tetrahedra=int(local.sum()), nodes=len(compact_points), elements=len(compact_cells),
            original_global_minimum_q=original_global_q, original_cavity_minimum_q=float(original_q.min()),
            subdivided_before_motion_minimum_q=float(before_deformation_quality.min()),
            final_cavity_minimum_q=float(quality.min()), original_cavity_minimum_jacobian=float(old_minimum.min()),
            final_cavity_minimum_jacobian=float(minimum.min()),
            cavity_volume_change=float(volumes.sum()-old_volumes.sum()),
            original_global_quality=dict(zip(qkeys, map(float, np.quantile(positive_quality(original_points[original_cells[:, :4]]), quantiles)))),
            candidate_global_quality=dict(zip(qkeys, map(float, np.quantile(positive_quality(compact_points[compact_cells[:, :4]]), quantiles)))),
            periodic_cap_clearance=cap_clearance, interface_faces=len(before_interface),
            complete_quadratic_interface_unchanged=True, interface_nodes_bitwise_unchanged=True,
            ridge_edge_original_ids=list(new_ridge), ridge_edge_compacted_ids=mapping[list(new_ridge)].tolist(),
            ridge_boundary_faces=incident_wall, local_remaining_crossing_faces=local_crossing,
            frozen_witness_distance_to_local_surface=local_witness,
            local_distance_scope="Upper bound on distance to the complete candidate surface; frozen witness does not select the candidate",
            mesh_sha256=w.digest(output/"mesh.npz"), accepted=not row["rejection_reasons"])
        part.update(1, "Checked local cavity and full conformity", admitted_locally=row["accepted"],
            minimum_q=float(quality.min()), minimum_jacobian=float(minimum.min()))
    except Exception as error:
        import traceback
        row.update(error_type=type(error).__name__, error=str(error), traceback=traceback.format_exc(), operations=operations)
        row["rejection_reasons"].append("cavity-control-exception")
        np.savez_compressed(output/"partial.npz", points=points, tetrahedra=cells, owners=owners)
        part.update(1, "Retained rejected cavity construction", failed=True)
    w.save(output/"result.json", row)
    return row


def controls():
    import numpy as np
    from microgen_tet10_control import straight_elevation
    from nested_refinement import FACES

    points = np.array([[0., 0., 0.], [1., 0., 0.], [.5, .8, 0.], [.5, .3, .9]])
    tets = np.array([[0, 1, 2, 3]])
    lifted = straight_elevation(points, tets, tets[:, FACES].reshape(-1, 3))
    pp, cc, _, _, stats = split_star(lifted["points"], lifted["tetrahedra"], np.array([0]), (0, 1), .3)
    assert len(cc) == 2 and stats["map_reproduction_max"] < 1e-12
    from facet_incidence_audit import audit
    pp, cc, ss, _ = compact(pp, cc, boundary(cc))
    assert audit(pp, cc, ss)["oriented_facet_consistency_passed"]
    return dict(arbitrary_fraction_exact_parent_map=True, affine_volume_preserved=True,
                conforming_oriented_edge_star=True)


def worker(args):
    np, _, meshers, *_ = w.imports()
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetBackend() == "Sequential" and vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    previous = json.loads((args.source/"result.json").read_text())
    source = args.source/"native-curved-tet10"
    mesh = source/"mesh.npz"
    assert previous["native"]["qualified"] and w.digest(mesh) == previous["native"]["mesh_sha256"]
    with np.load(mesh) as data:
        points, cells, surface = (data[key] for key in ("points", "tetrahedra", "surface"))
    paths = {Path(__file__), Path(w.__file__), Path(base.__file__), Path(probe.__file__), Path(curved.__file__),
        Path(base.fixture.__file__), args.source/"result.json", mesh, source/"result.json"}
    paths.update(w.HERE/name for name in ("nested_refinement.py", "scorecard.py", "periodic_trace.py", "comparison_mesh_audit.py",
        "quadratic_edge_removal_probe.py", "quadrature_audit.py", "facet_incidence_audit.py", "microgen_tet10_control.py", "box_projection.py"))
    paths.update(source/(base.NAME+suffix) for suffix in ("-samples.npz", "-specimen.json"))
    hashes = {str(p.resolve()): w.digest(p) for p in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(source/(base.NAME+suffix), args.run/(base.NAME+suffix))
    result = dict(inputs=hashes, source=str(source), complete=False, candidates=[],
        production_change=False, solve_performed=False, generation_timing_claim=False,
        vtk_backend=vtkSMPTools.GetBackend(), numerical_threads=1,
        scope="One worst-witness diagnostic cavity; prospective production selection must use branch identities")
    w.save(args.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=300, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        hypotheses=["Conforming edge-star splits expose an explicit shared ridge edge",
            "Unequal root split positions can violate the unchanged quality floor",
            "Consistent corner/midside motion is needed to retain positive curved Jacobians"],
        candidate_count=2, selection="The two split orders on the nearest branch-straddling face to the frozen worst ridge witness"))
    with Job(parts={"controls": 1, "cavity": 1, "split-order-0": 2, "split-order-1": 2, "global-audit": 4}) as job:
        with job.part("controls") as part:
            part.update(0, "Checking exact-map edge subdivision", controls=0)
            result["controls"] = controls()
            part.update(1, "Conforming subdivision controls passed", controls=len(result["controls"]))
        job.result("controls", passed=True)
        with job.part("cavity") as part:
            query = np.asarray(previous["localization"]["reverse_groups"]["ridges"]["worst_point"])
            part.update(0, "Locating the saved worst ridge cavity",
                        crossing_faces=int(probe.face_groups(points, surface)["crossing"].sum()))
            face, selection = nearest_straddling_face(points, surface, query)
            corners = surface[face, :3]
            branch = probe.branch_value(points)
            edges = [tuple(sorted(map(int, edge))) for edge in zip(corners, corners[[1, 2, 0]], strict=True)
                     if branch[edge[0]]*branch[edge[1]] < -1e-24]
            assert len(edges) == 2, "Selected face does not have exactly two branch-crossing edges"
            edges = sorted(edges)
            plane_sign = 1 if query[0]*query[1] > 0 else -1
            result["selection"] = dict(**selection, face_nodes=surface[face].tolist(),
                corner_coordinates=points[corners].tolist(), branch_values=branch[corners].tolist(),
                edges=[list(edge) for edge in edges], plane_sign=plane_sign,
                observer_used_for_localization_only=True)
            part.update(1, "Selected one straddling face and its edge stars", edges=len(edges))
        job.result("cavity", face=face, ridge_distance=selection["distance"])
        for index, order in enumerate((edges, edges[::-1])):
            directory = args.run/f"split-order-{index}"
            directory.mkdir()
            with job.part(f"split-order-{index}") as part:
                part.update(0, "Preparing conforming local candidate", order=index)
                row = candidate(points, cells, surface, order, plane_sign, query, directory, part)
            row["directory"] = str(directory)
            result["candidates"].append(row)
            job.result(f"split-order-{index}", locally_accepted=row["accepted"],
                minimum_q=row.get("final_cavity_minimum_q"), minimum_jacobian=row.get("final_cavity_minimum_jacobian"),
                rejection="; ".join(row["rejection_reasons"]))
            w.save(args.run/"result.json", result)
        with job.part("global-audit") as part:
            feasible = [row for row in result["candidates"] if row["accepted"]]
            if feasible:
                chosen = max(feasible, key=lambda row: row["final_cavity_minimum_q"])
                with np.load(Path(chosen["directory"])/"mesh.npz") as data:
                    arrays = tuple(data[key] for key in ("points", "tetrahedra", "surface"))
                result["selected_by_local_quality"] = chosen["directory"]
                result["full_admission"] = curved.archive(args.run, "selected-candidate", arrays, part)
                # Keep the witness observation separate from candidate selection.
                from scorecard import distance_to_mesh
                vertices, faces = w.tessellate(arrays[0], arrays[2], 16)
                result["frozen_worst_witness_after"] = float(distance_to_mesh(vertices, faces, query[None])[0])
                part.update(1, "Completed independent global admission", admitted=result["full_admission"]["qualified"])
            else:
                part.update(1, "No locally feasible candidate; failures retained", feasible=0)
                result["full_admission"] = dict(qualified=False, status="no-local-candidate")
        job.result("global-audit", admitted=result["full_admission"]["qualified"],
                   witness_distance=result.get("frozen_worst_witness_after"))
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    result.update(complete=True, input_hashes_verified=True, original_mesh_and_samples_unchanged=True)
    w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.run, args.source = args.run.resolve(), args.source.resolve()
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import report, result
    report(0, "Starting bounded CAD ridge cavity control", eta_s=300, memory_ceiling_gib=w.MEMORY/1024**3)
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=300)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    row = json.loads((args.run/"result.json").read_text())
    for index, candidate_row in enumerate(row["candidates"]):
        result(f"split-order-{index}", accepted=candidate_row["accepted"],
               q_min=candidate_row.get("final_cavity_minimum_q"),
               minimum_jacobian=candidate_row.get("final_cavity_minimum_jacobian"),
               rejection="; ".join(candidate_row["rejection_reasons"]))
    result("cavity-control", full_admission=row["full_admission"]["qualified"],
           witness_after=row.get("frozen_worst_witness_after"))


if __name__ == "__main__":
    main()
