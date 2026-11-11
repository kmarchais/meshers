"""One simultaneous curved replay of an existing three-corner ridge chain.

The saved Tet4 control supplies provenance for a possible chain, not proof
that its isolated motion is feasible. Targets use the analytic branch
intersection. Complete corner stars, interface, Jacobians and admission are
rechecked on the untouched curved Tet10 baseline.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import microgen_workflow as w
import cad_intersection_probe as probe
import cad_intersection_curved_control as curved
import cad_ridge_cavity_control as cavity
import cad_union_comparison as base


SOURCE = w.HERE/"runs/macbook-cad-intersection-curved-v1"
PROVENANCE = w.HERE/"runs/macbook-cad-intersection-ridge-v1"
CHAIN = (3184, 4797, 9548)


def move_chain(points, cells, surface, selected):
    import numpy as np

    moved = points.copy()
    mids = cavity.midpoint_map(cells)
    target = probe.ridge_target(points[selected])
    displacement = {int(node): point-points[node] for node, point in zip(selected, target, strict=True)}
    for node, delta in displacement.items():
        moved[node] += delta
    affected_mids = []
    for (a, b), mid in mids.items():
        if a in displacement or b in displacement:
            moved[mid] += .5*(displacement.get(a, 0)+displacement.get(b, 0))
            affected_mids.append(mid)
    walls = surface[probe.wall_mask(points, surface)]
    wall_edges = {tuple(sorted(map(int, pair))) for face in walls[:, :3]
                  for pair in zip(face, face[[1, 2, 0]], strict=True)}
    projected = []
    for edge in sorted(wall_edges):
        mid = mids[edge]
        if mid not in affected_mids:
            continue
        signs = probe.branch_value(moved[list(edge)])
        assert not (signs.min() < -1e-12 and signs.max() > 1e-12), "Moved wall edge crosses unresolved branches"
        mask = 1 if signs.max() > 1e-12 else 2 if signs.min() < -1e-12 else 3
        point = moved[mid].copy()
        if mask == 3:
            point = probe.ridge_target(point[None])[0]
        else:
            axes = [1, 2] if mask == 1 else [0, 2]
            point[axes] *= base.fixture.RADIUS/np.linalg.norm(point[axes])
        fields = np.array([np.hypot(point[1], point[2]), np.hypot(point[0], point[2])])-base.fixture.RADIUS
        projected.append(dict(edge=list(edge), midpoint=mid, constraint_mask=mask,
            additional_projection_displacement=float(np.linalg.norm(point-moved[mid])),
            branch_values=fields.tolist(), exposed=bool(fields.min() >= -1e-12 and abs(fields.min()) < 1e-10)))
        moved[mid] = point
    return moved, dict(corner_displacements={str(k): v.tolist() for k, v in displacement.items()},
        affected_midpoints=sorted(affected_mids), projected_wall_midpoints=projected,
        all_projected_midpoints_exposed=all(row["exposed"] for row in projected),
        maximum_corner_displacement=max(float(np.linalg.norm(v)) for v in displacement.values()),
        maximum_midpoint_displacement=float(np.linalg.norm(moved[affected_mids]-points[affected_mids], axis=1).max()))


def worker(args):
    np, *_ = w.imports()
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    from nested_refinement import positive_quality
    from quadratic_edge_removal_probe import measure
    from comparison_mesh_audit import quadratic_conformity
    from facet_incidence_audit import audit
    from periodic_trace import periodic_pairs
    from scorecard import distance_to_mesh
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    source = args.source/"native-curved-tet10"
    previous = json.loads((args.source/"result.json").read_text())
    provenance = json.loads((args.provenance/"result.json").read_text())
    assert previous["native"]["qualified"] and provenance["candidate"]["qualified"]
    assert w.digest(source/"mesh.npz") == previous["native"]["mesh_sha256"]
    assert w.digest(args.provenance/"ridge-control/mesh.npz") == provenance["candidate"]["mesh_sha256"]
    with np.load(source/"mesh.npz") as data:
        points, cells, surface = (data[key] for key in ("points", "tetrahedra", "surface"))
    paths = {Path(__file__), Path(w.__file__), Path(base.__file__), Path(probe.__file__), Path(curved.__file__),
        Path(cavity.__file__), Path(base.fixture.__file__), args.source/"result.json", source/"result.json",
        source/"mesh.npz", args.provenance/"result.json", args.provenance/"ridge-control/mesh.npz"}
    paths.update(w.HERE/name for name in ("nested_refinement.py", "scorecard.py", "periodic_trace.py", "comparison_mesh_audit.py",
        "quadratic_edge_removal_probe.py", "quadrature_audit.py", "facet_incidence_audit.py", "box_projection.py"))
    paths.update(source/(base.NAME+suffix) for suffix in ("-samples.npz", "-specimen.json"))
    hashes = {str(path.resolve()): w.digest(path) for path in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    for suffix in ("-samples.npz", "-specimen.json"):
        shutil.copyfile(source/(base.NAME+suffix), args.run/(base.NAME+suffix))
    result = dict(inputs=hashes, source=str(source), provenance=str(args.provenance), chain=list(CHAIN),
        complete=False, numerical_threads=1, vtk_backend=vtkSMPTools.GetBackend(),
        admitted=False, production_change=False, solve_performed=False, generation_timing_claim=False,
        scope="One simultaneous isolated curved chain motion, not a replay of all177 sequential Tet4 moves")
    w.save(args.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=240, memory_bytes=w.MEMORY,
        numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
        hypothesis="Reuse the existing ridge chain to avoid near-endpoint sliver insertion, then apply branch-aware quadratic geometry",
        fixed_chain=list(CHAIN), selection="Existing two-edge chain identified in saved747; no observer-distance optimization"))
    selected = np.asarray(CHAIN)
    with Job(parts={"provenance": 1, "motion": 1, "local-audit": 2, "global-admission": 4}) as job:
        with job.part("provenance") as part:
            part.update(0, "Checking immutable chain provenance", corners=len(selected))
            attempts = [next(row for row in provenance["repair"]["attempts"] if row["node"] == node) for node in CHAIN]
            assert all(row["status"] == "accepted" for row in attempts)
            np.testing.assert_array_equal(probe.ridge_target(points[selected]), np.asarray([row["point"] for row in attempts]))
            result["prior_sequential_corner_records"] = attempts
            result["prior_records_are_not_isolated_feasibility_proof"] = True
            local = np.isin(cells[:, :4], selected).any(1)
            local_nodes = np.unique(cells[local])
            clearance = float(.5-abs(points[local_nodes]).max())
            assert clearance > 1e-10, "Complete stars reach periodic caps"
            original_interface = cavity.boundary_keys(cavity.boundary(cells[local]))
            physical = cavity.boundary_keys(surface)
            interface = {key: value for key, value in original_interface.items() if key not in physical}
            interface_nodes = np.unique([node for face in interface.values() for node in face])
            assert not np.isin(interface_nodes, selected).any(), "Selected corner lies on cavity interface"
            result.update(cavity_tetrahedra=int(local.sum()), cavity_nodes=len(local_nodes),
                cavity_interface_faces=len(interface), periodic_cap_clearance=clearance)
            part.update(1, "Complete stars are interior to the periodic cell", tets=int(local.sum()), cap_clearance=clearance)
        job.result("provenance", cavity_tets=int(local.sum()), chain_corners=len(selected))
        with job.part("motion") as part:
            part.update(0, "Applying simultaneous chain and midside motion", corners=len(selected))
            moved, result["motion"] = move_chain(points, cells, surface, selected)
            np.savez_compressed(args.run/"candidate.npz", points=moved, tetrahedra=cells, surface=surface)
            result["candidate_sha256"] = w.digest(args.run/"candidate.npz")
            part.update(1, "Candidate archived before admission", moved_mids=len(result["motion"]["affected_midpoints"]))
        job.result("motion", moved_corners=len(selected), mids=len(result["motion"]["affected_midpoints"]))
        w.save(args.run/"result.json", result)
        with job.part("local-audit") as part:
            part.update(0, "Checking interface, quality and positive125", tets=int(local.sum()))
            np.testing.assert_array_equal(moved[interface_nodes], points[interface_nodes])
            outside = np.setdiff1d(np.arange(len(points)), local_nodes)
            np.testing.assert_array_equal(moved[outside], points[outside])
            before_q = positive_quality(points[cells[:, :4]])
            after_q = positive_quality(moved[cells[:, :4]])
            old_jac, old_vol = measure(points[cells[local]])
            new_jac, new_vol = measure(moved[cells[local]])
            facets = audit(moved, cells, surface)
            topology = quadratic_conformity(moved, cells, surface)
            assert topology["passed"]
            _, _, periodic = periodic_pairs(moved, surface, declared_bounds=([-.5]*3, [.5]*3))
            ridge_faces = []
            for a, b in zip(CHAIN[:-1], CHAIN[1:], strict=True):
                faces = np.flatnonzero((surface[:, :3] == a).any(1) & (surface[:, :3] == b).any(1))
                assert len(faces) == 2
                signs = probe.branch_value(moved[surface[faces, :3]])
                branch_ids = [1 if values.max() > 1e-12 else 2 if values.min() < -1e-12 else 3 for values in signs]
                ridge_faces.append(dict(edge=[a, b], boundary_faces=faces.tolist(), branches=branch_ids,
                    opposite_exposed_branches=sorted(branch_ids) == [1, 2]))
            quantiles = [0, .01, .1, .5, .9, .99, 1]
            keys = ("min", "p01", "p10", "median", "p90", "p99", "max")
            reasons = []
            if after_q.min() < .1:
                reasons.append("corner-quality-below-0.1")
            if new_jac.min() <= 0:
                reasons.append("nonpositive-sampled-curved-Jacobian")
            if not result["motion"]["all_projected_midpoints_exposed"]:
                reasons.append("hidden-union-wall-midpoint")
            if not all(edge["opposite_exposed_branches"] for edge in ridge_faces):
                reasons.append("ridge-does-not-separate-two-branches")
            if not facets["oriented_facet_consistency_passed"]:
                reasons.append("facet-incidence")
            if max(len(moved), len(cells)) > w.CAP:
                reasons.append("node-or-tetrahedron-cap")
            local_physical = np.array([i for i, face in enumerate(surface)
                if tuple(sorted(map(int, face[:3]))) in original_interface], dtype=int)
            query = np.asarray(previous["localization"]["reverse_groups"]["ridges"]["worst_point"])
            vertices, faces = w.tessellate(moved, surface[local_physical], 16)
            witness_distance = float(distance_to_mesh(vertices, faces, query[None])[0])
            result["local"] = dict(feasible=not reasons, rejection_reasons=reasons,
                original_global_quality=dict(zip(keys, map(float, np.quantile(before_q, quantiles)))),
                candidate_global_quality=dict(zip(keys, map(float, np.quantile(after_q, quantiles)))),
                original_local_minimum_q=float(before_q[local].min()), candidate_local_minimum_q=float(after_q[local].min()),
                original_local_minimum_jacobian=float(old_jac.min()), candidate_local_minimum_jacobian=float(new_jac.min()),
                local_volume_change=float(new_vol.sum()-old_vol.sum()), facet_audit=facets, quadratic_conformity=topology,
                periodic=periodic, ridge_edges=ridge_faces, interface_nodes_bitwise_unchanged=True,
                outside_cavity_nodes_bitwise_unchanged=True, counts_and_connectivity_unchanged=True,
                nodes=len(points), elements=len(cells), frozen_witness_distance_to_local_surface=witness_distance,
                witness_scope="Frozen observer only; local distance upper-bounds distance to the whole candidate surface")
            np.savez_compressed(args.run/"cavity.npz", parent_indices=np.flatnonzero(local), interface_nodes=interface_nodes,
                old_points=points[local_nodes], moved_points=moved[local_nodes], original_node_ids=local_nodes,
                cells=cells[local], before_quality=before_q[local], after_quality=after_q[local],
                before_jacobian=old_jac, after_jacobian=new_jac)
            part.update(1, "Finished complete affected-star checks", feasible=not reasons,
                q_min=float(after_q.min()), minimum_jacobian=float(new_jac.min()))
        job.result("local-audit", feasible=not reasons, q_min=float(after_q.min()),
            minimum_jacobian=float(new_jac.min()), witness_distance=witness_distance)
        w.save(args.run/"result.json", result)
        with job.part("global-admission") as part:
            if not reasons:
                result["full_admission"] = curved.archive(args.run, "chain-candidate", (moved, cells, surface), part)
                result["admitted"] = result["full_admission"]["qualified"]
            else:
                result["full_admission"] = dict(qualified=False, status="local-rejection", reasons=reasons)
                part.update(1, "Rejected candidate retained; no global success claim", feasible=False)
        job.result("global-admission", admitted=result["admitted"],
            boundary=result["full_admission"].get("geometry", {}).get("sampled_two_sided_max"))
    assert all(w.digest(Path(path)) == digest for path, digest in hashes.items())
    result.update(complete=True, input_hashes_verified=True, original_meshes_and_samples_unchanged=True)
    w.save(args.run/"result.json", result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--provenance", type=Path, default=PROVENANCE)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for name in ("run", "source", "provenance"):
        setattr(args, name, getattr(args, name).resolve())
    if args.worker:
        return worker(args)
    args.run.mkdir(parents=True, exist_ok=False)
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import report, result
    report(0, "Starting bounded curved ridge-chain control", eta_s=240, memory_ceiling_gib=w.MEMORY/1024**3)
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=240)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    row = json.loads((args.run/"result.json").read_text())
    result("curved-ridge-chain", locally_feasible=row["local"]["feasible"], admitted=row["admitted"],
        q_min=row["local"]["candidate_global_quality"]["min"],
        minimum_jacobian=row["local"]["candidate_local_minimum_jacobian"],
        local_witness=row["local"]["frozen_witness_distance_to_local_surface"],
        rejection="; ".join(row["local"]["rejection_reasons"]))


if __name__ == "__main__":
    main()
