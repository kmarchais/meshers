"""Saved-mesh forward observer localization; never mesh generation or admission."""

import json
from pathlib import Path

import microgen_workflow as w


NAME = "tpms-pmy-sheet"
GATE = .002
GLOBAL_WORST = 8


def taylor_ball(case, query):
    """PMY Fourier modes give ||H|| <= (6+5+5+5)*(2*pi)^2 everywhere."""
    import numpy as np
    from pmy_neck_control import raw_gradient

    value = float(case.function(*query))
    gradient = float(np.linalg.norm(raw_gradient(query[None])[0]))
    hessian_bound = 21*(2*np.pi)**2
    excursion = gradient*GATE+.5*hessian_bound*GATE**2+1e-12
    excluded = [abs(value-level) > excursion for level in case.levels]
    outside = value < case.levels[0] or value > case.levels[1]
    cap_distance = float(np.min(.5-abs(query)))
    return dict(radius=GATE, field_value=value, gradient_norm=gradient, hessian_operator_bound=hessian_bound,
        maximum_field_excursion=excursion, wall_levels_excluded=bool(all(excluded)),
        material_boundary_ball_excluded=bool(all(excluded) and (outside or cap_distance > GATE+1e-12)),
        distance_to_nearest_box_plane=cap_distance,
        scope="Analytic global PMY Hessian bound with floating-point padding; not interval arithmetic. Cap material is conservatively included.")


def run(a):
    import numpy as np
    from scipy.spatial import cKDTree
    from vtkmodules.vtkCommonCore import vtkSMPTools
    import scorecard
    import pmy_failure_probe as probe
    import tpms_family_fixtures as fixtures
    from box_projection import constrained_foot

    assert vtkSMPTools.SetBackend("Sequential")
    vtkSMPTools.Initialize(1)
    assert vtkSMPTools.GetEstimatedNumberOfThreads() == 1
    directory = a.run/"candidate"
    row = json.loads((directory/"result.json").read_text())
    assert w.digest(directory/"mesh.npz") == row["mesh_sha256"]
    case = fixtures.case(NAME, a.fixtures)
    with np.load(directory/"mesh.npz") as mesh:
        points, cells, surface = mesh["points"], mesh["tetrahedra"], mesh["surface"]
    assert surface.shape[1] == 6 and cells.shape[1] == 10
    controls = probe.controls()
    assert taylor_ball(case, np.zeros(3))["material_boundary_ball_excluded"]
    corners = points[surface[:, :3]]
    caps = np.any(np.all(np.isclose(corners, -.5, atol=1e-12, rtol=0), axis=1) |
                  np.all(np.isclose(corners, .5, atol=1e-12, rtol=0), axis=1), axis=1)
    wall_ids = np.flatnonzero(~caps)
    wall = surface[wall_ids]
    weights = w.tri_values(scorecard.BARY)
    samples = np.einsum("sn,fnd->fsd", weights, points[wall]).reshape(-1, 3)
    locked = np.isclose(abs(samples), .5, atol=1e-12, rtol=0)
    foot = samples.copy()
    values = case.function(*samples.T)
    levels_array = np.asarray(case.levels)
    levels = levels_array[np.argmin(abs(values[:, None]-levels_array), axis=1)]
    for _ in range(40):
        gradient = scorecard.gradient(case.function, foot)
        gradient[locked] = 0
        remaining = case.function(*foot.T)-levels
        step = remaining[:, None]*gradient/np.maximum((gradient**2).sum(1), 1e-30)[:, None]
        foot -= step
        if np.max(np.linalg.norm(step, axis=1)) < 1e-12:
            break
    residual = abs(case.function(*foot.T)-levels)/np.maximum(np.linalg.norm(scorecard.gradient(case.function, foot), axis=1), 1e-30)
    failed = (residual > 1e-9) | np.any(abs(foot) > .5+1e-9, axis=1)
    raw_foot = foot.copy()
    fallback = []
    for index in np.flatnonzero(failed):
        alternative, info = constrained_foot(case.function, samples[index], levels[index], locked[index], foot[index])
        fallback.append(dict(sample_index=int(index), **info))
        if info["feasible_target"]:
            foot[index] = alternative
    residual = abs(case.function(*foot.T)-levels)/np.maximum(np.linalg.norm(scorecard.gradient(case.function, foot), axis=1), 1e-30)
    failed = (~np.isfinite(foot).all(1)) | (~np.isfinite(residual)) | (residual > 1e-9) | np.any(abs(foot) > .5+1e-9, axis=1)
    distance = np.linalg.norm(foot-samples, axis=1)
    recorded = row.get("geometry", {}).get("sampled_forward_max")
    if recorded is not None:
        np.testing.assert_allclose(distance.max(), recorded, atol=1e-12, rtol=0)
    kinds = []
    for bary in scorecard.BARY:
        kinds.append("corner" if (bary == 1).any() else "edge-midpoint" if (bary == .5).sum() == 2
                     else "edge-interior" if (bary == 0).any() else "face-interior")
    sample_kinds = np.tile(np.asarray(kinds), len(wall))
    category_summary = {kind: dict(samples=int((sample_kinds == kind).sum()),
        maximum=float(distance[sample_kinds == kind].max()), above_gate=int((distance[sample_kinds == kind] > GATE).sum())) for kind in sorted(set(kinds))}
    # Global worst distinct sample coordinates, plus the worst of each category.
    selected, seen = [], set()
    for index in np.argsort(-distance, kind="stable"):
        key = tuple(samples[index])
        if key not in seen:
            selected.append(int(index))
            seen.add(key)
        if len(selected) == GLOBAL_WORST:
            break
    for kind in sorted(set(kinds)):
        ids = np.flatnonzero(sample_kinds == kind)
        index = int(ids[np.argmax(distance[ids])])
        if index not in selected:
            selected.append(index)
    np.savez_compressed(directory/"forward-samples.npz", points=samples, raw_newton_feet=raw_foot,
        standard_feet=foot, standard_distances=distance, failed=failed, wall_face_ids=wall_ids,
        barycentric=scorecard.BARY, sample_kind=np.asarray(kinds), selected_indices=np.asarray(selected),
        source_mesh_sha256=np.asarray(row["mesh_sha256"]))
    trace = json.loads((directory/"native-snap-trace.json").read_text())
    restored = sorted({node for state in trace for node in state["restored_vertices"]})
    parsed = json.loads((directory/"native-input-bits.json").read_text())
    background_bits = np.asarray(parsed["point_bits"], dtype=np.uint64)
    background = background_bits.view(np.float64)
    restored_tree = cKDTree(background[restored]) if restored else None
    coordinate_lookup = {}
    for index, bits in enumerate(background_bits):
        coordinate_lookup.setdefault(tuple(map(int, bits)), []).append(index)
    rows = []
    for index in selected:
        wall_index, probe_index = divmod(index, len(weights))
        face_index = int(wall_ids[wall_index])
        nodes = surface[face_index]
        xyz, query = points[nodes], samples[index]
        starts = [query, foot[index], *xyz[:3], xyz[3:].mean(0)]
        feasible = probe.closest_feasible(case, query, starts)
        # The already validated standard foot is itself a feasible upper bound.
        if not failed[index] and (feasible is None or distance[index] < feasible["distance"]):
            feasible = dict(foot=foot[index].tolist(), distance=float(distance[index]), level=float(levels[index]),
                residual_distance=float(residual[index]), optimizer_success=None, source="standard-feasible-foot")
        node_associations = []
        for node in nodes:
            matches = coordinate_lookup.get(tuple(map(int, points[node].view(np.uint64))), [])
            node_associations.append(dict(mesh_node=int(node), exact_coordinate_coincidences=matches,
                coincides_with_restored_background_vertex=any(i in restored for i in matches)))
        nearest = None
        if restored_tree is not None:
            distances, locations = restored_tree.query(np.vstack((query, xyz)))
            nearest = dict(probe_distance=float(distances[0]), nearest_restored_background_vertex=int(restored[int(locations[0])]),
                face_node_distances=distances[1:].tolist())
        node_weights = weights[probe_index]
        mesh_node = int(nodes[np.argmax(node_weights)]) if np.count_nonzero(node_weights) == 1 and node_weights.max() == 1 else None
        rows.append(dict(sample_index=index, face=face_index, probe_index=probe_index, kind=kinds[probe_index],
            barycentric=scorecard.BARY[probe_index].tolist(), probe=query.tolist(), mesh_node=mesh_node,
            node_ids=nodes.tolist(), coordinates=xyz.tolist(), standard_foot=foot[index].tolist(),
            standard_distance=float(distance[index]), standard_failed=bool(failed[index]), locked=locked[index].tolist(),
            closest_found_feasible=feasible,
            forward_overstatement_at_least=float(distance[index]-feasible["distance"]) if feasible is not None else None,
            taylor_ball=taylor_ball(case, query), midpoint_displacements=np.linalg.norm(
                xyz[3:]-xyz[np.asarray(((0, 1), (1, 2), (2, 0)))].mean(1), axis=1).tolist(),
            background_coordinate_associations=node_associations, nearest_restored_spatial_association=nearest,
            association_scope="Coordinate coincidence and spatial proximity only; original source IDs are not retained through clipping/optimization, so lineage and causal association remain unknown"))
    w.save(directory/"forward-localization.json", dict(status="observed", complete=True, controls=controls,
        mesh_sha256=row["mesh_sha256"], raw_standard_forward_max=float(distance.max()), recorded_forward_max=recorded,
        standard_observer_reproduced=recorded is not None, category_summary=category_summary,
        sample_archive_sha256=w.digest(directory/"forward-samples.npz"), samples=len(samples), fallback=fallback,
        witness_budget="8 distinct global worst coordinates plus worst sample per category, maximum12", witnesses=rows,
        summary=dict(witnesses=len(rows), material_boundary_balls_excluded=sum(r["taylor_ball"]["material_boundary_ball_excluded"] for r in rows),
            feasible_within_gate=sum(r["closest_found_feasible"] is not None and r["closest_found_feasible"]["distance"] <= GATE for r in rows),
            maximum_proven_overstatement=max((r["forward_overstatement_at_least"] or 0) for r in rows)),
        no_admission_change=True, no_meshing=True,
        limits="Multistart feet are feasible upper bounds, not globally closest certificates. Raw geometry audits stay unchanged. Taylor exclusions use a global analytic Hessian bound and floating-point padding. Spatial associations do not prove rollback caused an error."))
