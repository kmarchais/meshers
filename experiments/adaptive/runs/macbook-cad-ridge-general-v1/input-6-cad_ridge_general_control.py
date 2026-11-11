"""Frozen field-derived disjoint ridge chains under two quality policies.

At most sixteen three-corner chains are selected from the unchanged CAD mesh.
Selection uses branch identity, analytic joint distance and original curved
field residuals. The two policies see exactly the same proposed motions. No
frozen surface query or target volume is used to select or accept a motion.
"""

import argparse
import json
from pathlib import Path
import shutil
import sys
import traceback
from types import SimpleNamespace

import microgen_workflow as w
import cad_intersection_probe as probe
import cad_intersection_curved_control as curved
import cad_ridge_chain_control as chain
import cad_ridge_cavity_control as cavity
import cad_union_comparison as base


SOURCE = w.HERE/"runs/macbook-cad-intersection-curved-v1/native-curved-tet10"
MAX_CHAINS = 16
ROUND_OFF = 1e-13


def proposal_plan(points, cells, surface):
    import numpy as np
    import scorecard

    wall_ids = np.flatnonzero(probe.wall_mask(points, surface))
    wall = surface[wall_ids]
    wall_nodes = np.unique(wall[:, :3])
    branches = probe.branch_value(points)
    edge_faces = {}
    corner_faces = {}
    for face_id in wall_ids:
        corners = surface[face_id, :3]
        for node in corners:
            corner_faces.setdefault(int(node), []).append(int(face_id))
        for a, b in zip(corners, corners[[1, 2, 0]], strict=True):
            edge_faces.setdefault(tuple(sorted((int(a), int(b)))), []).append(int(face_id))
    edges = np.array(sorted(edge_faces))
    lengths = np.linalg.norm(points[edges[:, 0]]-points[edges[:, 1]], axis=1)
    shortest = np.full(len(points), np.inf)
    np.minimum.at(shortest, edges[:, 0], lengths)
    np.minimum.at(shortest, edges[:, 1], lengths)
    targets = probe.ridge_target(points[wall_nodes])
    displacement = np.linalg.norm(targets-points[wall_nodes], axis=1)
    ratios = displacement/shortest[wall_nodes]
    eligible = set(map(int, wall_nodes[ratios <= .25]))
    ratio_by_node = dict(zip(map(int, wall_nodes), map(float, ratios), strict=True))
    graph = {}
    edge_rejections = {"ineligible-endpoint": 0, "not-two-wall-faces": 0, "not-opposite-branches": 0}
    for edge, faces in sorted(edge_faces.items()):
        if not set(edge) <= eligible:
            edge_rejections["ineligible-endpoint"] += 1
            continue
        if len(faces) != 2:
            edge_rejections["not-two-wall-faces"] += 1
            continue
        third = [next(int(node) for node in surface[face, :3] if node not in edge) for face in faces]
        if branches[third[0]]*branches[third[1]] >= -1e-24:
            edge_rejections["not-opposite-branches"] += 1
            continue
        a, b = edge
        graph.setdefault(a, []).append(b)
        graph.setdefault(b, []).append(a)
    # This proxy is derived only from the unchanged field and mesh faces.
    probes = np.einsum("sn,fnd->fsd", w.tri_values(scorecard.BARY, 2), points[wall])
    residual = abs(base.fixture.periodic_field(*probes.reshape(-1, 3).T)).reshape(len(wall), -1).max(1)
    error_by_face = dict(zip(map(int, wall_ids), map(float, residual), strict=True))
    parents_by_corner = {}
    for parent, tet in enumerate(cells[:, :4]):
        for node in tet:
            parents_by_corner.setdefault(int(node), []).append(parent)
    proposals = []
    structural_skips = []
    for center, neighbors in sorted(graph.items()):
        if len(neighbors) != 2:
            structural_skips.append(dict(center=center, reason="ridge-graph-degree-not-two", degree=len(neighbors)))
            continue
        selected = [min(neighbors), center, max(neighbors)]
        parents = sorted({p for node in selected for p in parents_by_corner[node]})
        nodes = np.unique(cells[parents])
        clearance = float(.5-abs(points[nodes]).max())
        if clearance <= 1e-10:
            structural_skips.append(dict(center=center, chain=selected,
                reason="periodic-cap-pairing-not-implemented-for-this-control", cap_clearance=clearance,
                classification="unimplemented coverage, not a repaired failure"))
            continue
        incident_faces = sorted({face for node in selected for face in corner_faces[node]})
        score = max(error_by_face[face] for face in incident_faces)
        proposals.append(dict(chain=selected, center=center, rank_field_residual=score,
            maximum_corner_displacement_ratio=max(ratio_by_node[node] for node in selected),
            parents=parents, complete_star_nodes=nodes.tolist(), wall_faces=incident_faces,
            cap_clearance=clearance))
    proposals.sort(key=lambda row: (-row["rank_field_residual"], row["center"]))
    frozen, occupied = [], set()
    for proposal in proposals:
        if len(frozen) >= MAX_CHAINS:
            proposal["selection"] = "bounded-proposal-limit"
        elif occupied.intersection(proposal["complete_star_nodes"]):
            proposal["selection"] = "overlapping-complete-star-node-including-midpoints"
        else:
            proposal["selection"] = "selected"
            proposal["proposal_index"] = len(frozen)
            occupied.update(proposal["complete_star_nodes"])
            frozen.append(proposal)
    return dict(eligible_wall_corners=len(eligible), wall_corner_count=len(wall_nodes),
        graph_nodes=len(graph), graph_edges=sum(map(len, graph.values()))//2,
        edge_rejections=edge_rejections, structural_skips=structural_skips,
        ranked_proposals=proposals, frozen=frozen, maximum_chains=MAX_CHAINS,
        maximum_corner_displacement_shortest_wall_edge=.25,
        selected_complete_stars_include_all_quadratic_nodes=True,
        rank="Descending original Tri6 max absolute implicit field residual over incident wall faces; center ID tie-break",
        no_frozen_queries_or_volume_targets_used=True)


def quality_summary(values):
    import numpy as np

    return dict(zip(("min", "p01", "p10", "median", "p90", "p99", "max"),
                    map(float, np.quantile(values, [0, .01, .1, .5, .9, .99, 1]))))


def evaluate(proposal, points, cells, surface, original_quality, physical_faces, output):
    import numpy as np
    from nested_refinement import positive_quality
    from quadratic_edge_removal_probe import measure

    row = dict(proposal_index=proposal["proposal_index"], chain=proposal["chain"],
        accepted_minimum=False, accepted_sorted=False, common_rejections=[],
        minimum_rejections=[], sorted_rejections=[])
    updates = None
    try:
        selected = np.asarray(proposal["chain"])
        parents = np.asarray(proposal["parents"])
        nodes = np.asarray(proposal["complete_star_nodes"])
        moved, motion = chain.move_chain(points, cells, surface, selected, ridge_fallback=True)
        # Save every proposed local candidate before assessing it.
        np.savez_compressed(output/"candidate-cavity.npz", nodes=nodes, points_before=points[nodes],
            points_after=moved[nodes], parents=parents, cells=cells[parents])
        row["candidate_sha256"] = w.digest(output/"candidate-cavity.npz")
        changed = np.flatnonzero(np.any(moved != points, axis=1))
        assert np.isin(changed, nodes).all(), "Motion escapes its complete affected stars"
        interface = cavity.boundary_keys(cavity.boundary(cells[parents]))
        interface = {key: value for key, value in interface.items() if key not in physical_faces}
        interface_nodes = np.unique([node for face in interface.values() for node in face])
        np.testing.assert_array_equal(moved[interface_nodes], points[interface_nodes])
        old_q = original_quality[parents]
        new_q = positive_quality(moved[cells[parents, :4]])
        old_jac, old_vol = measure(points[cells[parents]])
        new_jac, new_vol = measure(moved[cells[parents]])
        sorted_change = np.sort(new_q)-np.sort(old_q)
        np.savez_compressed(output/"local-metrics.npz", original_quality=old_q, candidate_quality=new_q,
            sorted_quality_change=sorted_change, original_minimum_jacobian=old_jac,
            candidate_minimum_jacobian=new_jac)
        ridge_edges = []
        for a, b in zip(proposal["chain"][:-1], proposal["chain"][1:], strict=True):
            face_ids = np.flatnonzero((surface[:, :3] == a).any(1) & (surface[:, :3] == b).any(1))
            if len(face_ids) != 2:
                row["common_rejections"].append("ridge-edge-does-not-have-two-boundary-faces")
                continue
            values = probe.branch_value(moved[surface[face_ids, :3]])
            labels = [1 if v.max() > 1e-12 else 2 if v.min() < -1e-12 else 3 for v in values]
            opposite = sorted(labels) == [1, 2]
            ridge_edges.append(dict(edge=[a, b], faces=face_ids.tolist(), labels=labels, opposite_branches=opposite))
            if not opposite:
                row["common_rejections"].append("ridge-does-not-separate-two-exposed-branches")
        if not motion["all_projected_midpoints_exposed"]:
            row["common_rejections"].append("hidden-projected-wall-midpoint")
        if new_jac.min() <= 0:
            row["common_rejections"].append("nonpositive-positive125-Jacobian")
        if new_q.min() < .1:
            row["common_rejections"].append("corner-quality-below-0.1")
        if new_q.min() < original_quality.min()-ROUND_OFF:
            row["minimum_rejections"].append("below-original-global-minimum-q")
        if sorted_change.min() < -ROUND_OFF:
            row["sorted_rejections"].append("sorted-local-quality-distribution-regression")
        row["accepted_minimum"] = not row["common_rejections"] and not row["minimum_rejections"]
        row["accepted_sorted"] = row["accepted_minimum"] and not row["sorted_rejections"]
        row.update(motion=motion, ridge_edges=ridge_edges, changed_nodes=changed.tolist(),
            cavity_tetrahedra=len(parents), interface_nodes_bitwise_unchanged=True,
            original_local_quality=quality_summary(old_q), candidate_local_quality=quality_summary(new_q),
            minimum_sorted_quality_change=float(sorted_change.min()),
            regressed_sorted_quality_entries=int((sorted_change < -ROUND_OFF).sum()),
            original_local_minimum_jacobian=float(old_jac.min()), candidate_local_minimum_jacobian=float(new_jac.min()),
            local_volume_change=float(new_vol.sum()-old_vol.sum()))
        updates = (changed, moved[changed])
    except Exception as error:
        row.update(error_type=type(error).__name__, error=str(error), traceback=traceback.format_exc())
        row["common_rejections"].append("candidate-construction-or-interface-error")
    w.save(output/"result.json", row)
    return row, updates


def worker(args):
    np, *_ = w.imports()
    from vtkmodules.vtkCommonCore import vtkSMPTools
    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    from nested_refinement import positive_quality
    from periodic_trace import periodic_pairs
    from comparison_mesh_audit import quadratic_conformity
    from facet_incidence_audit import audit
    sys.path.insert(0, str(Path.home()/"bin"))
    from pueue_progress import Job

    source = args.source
    saved = json.loads((source/"result.json").read_text())
    assert saved["qualified"] and w.digest(source/"mesh.npz") == saved["mesh_sha256"]
    with np.load(source/"mesh.npz") as data:
        points, cells, surface = (data[key] for key in ("points", "tetrahedra", "surface"))
    assert cells.shape[1] == 10 and surface.shape[1] == 6
    original_quality = positive_quality(points[cells[:, :4]])
    np.testing.assert_allclose(original_quality.min(), saved["geometry"]["q_min_corner"], atol=1e-14, rtol=0)
    paths = {Path(__file__), Path(w.__file__), Path(base.__file__), Path(base.fixture.__file__), Path(probe.__file__),
        Path(curved.__file__), Path(cavity.__file__), Path(chain.__file__), source/"mesh.npz", source/"result.json"}
    paths.update(w.HERE/name for name in ("scorecard.py", "nested_refinement.py", "quadratic_edge_removal_probe.py",
        "quadrature_audit.py", "periodic_trace.py", "comparison_mesh_audit.py", "facet_incidence_audit.py", "box_projection.py"))
    hashes = {str(path.resolve()): w.digest(path) for path in sorted(paths)}
    for index, path in enumerate(sorted(paths)):
        if path.suffix == ".py":
            shutil.copyfile(path, args.run/f"input-{index}-{path.name}")
    result = dict(inputs=hashes, source=str(source), complete=False, candidates=[], policies={},
        production_change=False, solve_performed=False, generation_timing_claim=False,
        vtk_backend=vtkSMPTools.GetBackend(), numerical_threads=1,
        scope="At most16 field-selected disjoint chains, same frozen proposals for global-minimum and sorted-local quality guards")
    physical_faces = cavity.boundary_keys(surface)
    outputs = {name: points.copy() for name in ("minimum", "sorted")}
    accepted = {name: [] for name in outputs}
    with Job(parts={"freeze-proposals": 1, "local-trials": 4, "minimum-audit": 3, "sorted-audit": 3}) as job:
        with job.part("freeze-proposals") as part:
            part.update(0, "Deriving ridge graph from original fields and mesh", nodes=len(points), tets=len(cells))
            plan = proposal_plan(points, cells, surface)
            w.save(args.run/"frozen-proposals.json", plan)
            result["proposal_sha256"] = w.digest(args.run/"frozen-proposals.json")
            result["selection"] = {key: value for key, value in plan.items() if key not in ("ranked_proposals", "frozen")}
            result["frozen_proposal_count"] = len(plan["frozen"])
            # Only after ranking is frozen are observer inputs copied and hashed.
            for suffix in ("-samples.npz", "-specimen.json"):
                path = source/(base.NAME+suffix)
                hashes[str(path.resolve())] = w.digest(path)
                shutil.copyfile(path, args.run/path.name)
            result["observer_inputs_accessed_only_after_proposal_freeze"] = True
            w.save(args.run/"manifest.json", dict(inputs=hashes, worker_timeout_s=300, memory_bytes=w.MEMORY,
                numerical_thread_policy=w.NUMERICAL_THREAD_POLICY, node_cap=w.CAP, tetrahedron_cap=w.CAP,
                proposal_sha256=result["proposal_sha256"], maximum_proposals=MAX_CHAINS,
                sorted_comparison_roundoff=ROUND_OFF, policies={
                    "minimum": "Original global corner q minimum, exposure, positive125 and complete-star interface guards",
                    "sorted": "All minimum-policy guards plus every sorted local corner q no worse within roundoff"},
                cap_policy="Skip cap-touching stars as unimplemented paired-operation coverage; no periodicity relaxation",
                observer_inputs_accessed_only_after_proposal_freeze=True))
            part.update(1, "Frozen disjoint chains before observer input access", chains=len(plan["frozen"]),
                        eligible_corners=plan["eligible_wall_corners"])
        job.result("freeze-proposals", chains=len(plan["frozen"]), graph_edges=plan["graph_edges"])
        with job.part("local-trials") as part:
            occupied = set()
            for index, proposal in enumerate(plan["frozen"]):
                assert not occupied.intersection(proposal["complete_star_nodes"])
                occupied.update(proposal["complete_star_nodes"])
                directory = args.run/f"proposal-{index:02d}"
                directory.mkdir()
                row, updates = evaluate(proposal, points, cells, surface, original_quality, physical_faces, directory)
                result["candidates"].append(row)
                for name in outputs:
                    if row["accepted_"+name]:
                        nodes, moved = updates
                        outputs[name][nodes] = moved
                        accepted[name].append(index)
                part.update((index+1)/max(1, len(plan["frozen"])), "Testing identical motions under both quality guards",
                    completed=index+1, minimum_accepted=len(accepted["minimum"]), sorted_accepted=len(accepted["sorted"]))
                w.save(args.run/"result.json", result)
            if not plan["frozen"]:
                part.update(1, "No eligible disjoint ridge chains", accepted=0)
        job.result("local-trials", minimum_accepted=len(accepted["minimum"]), sorted_accepted=len(accepted["sorted"]))
        for name, moved in outputs.items():
            with job.part(name+"-audit") as part:
                part.update(0, "Checking the composed full mesh", accepted=len(accepted[name]), nodes=len(points))
                q = positive_quality(moved[cells[:, :4]])
                delta = np.sort(q)-np.sort(original_quality)
                np.savez_compressed(args.run/(name+"-quality.npz"), original_sorted=np.sort(original_quality),
                    candidate_sorted=np.sort(q), difference=delta)
                strict = bool((delta >= -ROUND_OFF).all())
                if name == "sorted":
                    assert strict, "Disjoint local sorted dominance did not compose globally"
                assert q.min() >= original_quality.min()-ROUND_OFF
                assert len(moved) == len(points) and max(len(moved), len(cells)) <= w.CAP
                assert quadratic_conformity(moved, cells, surface)["passed"]
                assert audit(moved, cells, surface)["oriented_facet_consistency_passed"]
                _, _, periodic = periodic_pairs(moved, surface, declared_bounds=([-.5]*3, [.5]*3))
                record = curved.archive(args.run, name+"-candidate", (moved, cells, surface), part)
                record.update(accepted_proposals=accepted[name], unchanged_baseline=bool(np.array_equal(moved, points)),
                    full_sorted_quality_nondecreasing=strict, minimum_sorted_quality_change=float(delta.min()),
                    regressed_sorted_quality_entries=int((delta < -ROUND_OFF).sum()),
                    original_quality=quality_summary(original_quality), candidate_quality=quality_summary(q),
                    fixed_counts=True, complete_star_nodes_disjoint=True, periodic_composition_check=periodic)
                geometry = record.get("geometry", {})
                record["sampled_boundary_improved"] = bool(geometry.get("sampled_two_sided_max", float("inf")) < saved["geometry"]["sampled_two_sided_max"])
                result["policies"][name] = record
            job.result(name+"-audit", admitted=record["qualified"], accepted=len(accepted[name]),
                boundary=geometry.get("sampled_two_sided_max"), sorted_quality_nondecreasing=strict)
            w.save(args.run/"result.json", result)
    assert w.digest(args.run/"frozen-proposals.json") == result["proposal_sha256"]
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
    report(0, "Starting bounded general ridge-chain control", eta_s=300, memory_ceiling_gib=w.MEMORY/1024**3)
    record = w.bounded([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
        args.run/"worker.log", SimpleNamespace(update=lambda *a, **k: None), timeout=300)
    w.save(args.run/"worker.json", record)
    assert record["returncode"] == 0, "Read worker.log and worker.json"
    rows = json.loads((args.run/"result.json").read_text())
    for name, row in rows["policies"].items():
        geometry = row.get("geometry", {})
        result(name, accepted=len(row["accepted_proposals"]), admitted=row["qualified"],
            boundary=geometry.get("sampled_two_sided_max"), volume_error=geometry.get("relative_volume_error"),
            full_sorted_quality_nondecreasing=row["full_sorted_quality_nondecreasing"],
            regressed_sorted_quality_entries=row["regressed_sorted_quality_entries"],
            unchanged_baseline=row["unchanged_baseline"])


if __name__ == "__main__":
    main()
