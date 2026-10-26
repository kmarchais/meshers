"""Geometry-only split selection with quality-aware longer-edge fallback."""

import numpy as np
from nested_refinement import EDGES, child_maps, positive_quality
from crease_directional_refinement import edge_metric_squared
from longest_edge_closure import closure_scores

QUALITY_POLICY = dict(child_quality_minimum=.1,continuation_depth=1,
    continuation='A locally longest child edge admits both grandchildren at the same quality floor',
    continuation_required_only_if_sampled_child_metric_unmet=True,
    fallback='Strictly longer incident edges in complete periodic orbits',
    fallback_leaves_checked_by_same_child_rule=True)


def select_with_quality_closure(points, cells, surface, scores, tensor, *, periodic=True):
    from fedoo.lib_elements.tetrahedron import Tet10

    pairs = np.sort(cells[:,EDGES].reshape(-1,2),axis=1)
    edges, inverse = np.unique(pairs,axis=0,return_inverse=True)
    local_edges = inverse.reshape(-1,6)
    classes = np.arange(len(points))
    if periodic:
        from periodic import periodic_pairs
        _, classes, _ = periodic_pairs(points,surface)
    keys, orbit = np.unique(np.sort(classes[edges],axis=1),axis=0,return_inverse=True)
    lookup = {tuple(edge):i for i,edge in enumerate(edges)}
    weight = np.zeros(len(keys))
    for edge,value in scores.items():
        assert np.isfinite(value) and value >= 0
        group = orbit[lookup[edge]]
        weight[group] = max(weight[group],value)
    requested = weight > 0
    readiness = np.full(len(keys),-1,dtype=np.int8)
    element = Tet10()
    barycentric = np.vstack((np.eye(4),np.eye(4)[EDGES].mean(axis=1)))
    weights = [element.shape_function((barycentric@mapping)@element.xi_nd[:4])
        for offset in range(6) for mapping in child_maps(offset)]
    weights = np.array(weights).reshape(6,2,10,10)
    examined = 0

    def ready_orbits(marked):
        nonlocal examined
        unknown = marked & (readiness < 0)
        readiness[unknown] = 1
        parents, offsets = np.nonzero(unknown[orbit[local_edges]])
        examined += len(parents)
        for offset in range(6):
            selected = parents[offsets == offset]
            for start in range(0,len(selected),256):
                take = selected[start:start+256]
                children = np.einsum('can,bnd->bcad',weights[offset],points[cells[take]])
                child_points = children.reshape(-1,10,3)
                immediate = positive_quality(child_points[:,:4]) >= .1+1e-12
                flat = child_points.reshape(-1,3)
                child_cells = np.arange(len(flat)).reshape(-1,10)
                metric = edge_metric_squared(flat,child_cells[:,EDGES].reshape(-1,2),
                    child_cells[:,4:].ravel(),tensor).reshape(-1,6)
                needs_more = (metric > (4/3)**2).any(axis=1)
                allowed = np.ones((len(child_points),6),dtype=bool)
                for second, endpoints in enumerate(EDGES):
                    for endpoint in endpoints:
                        xyz = child_points[:,:4].copy()
                        xyz[:,endpoint] = child_points[:,4+second]
                        allowed[:,second] &= positive_quality(xyz) >= .1+1e-12
                length2 = ((child_points[:,EDGES[:,1]]-child_points[:,EDGES[:,0]])**2).sum(axis=2)
                longest = length2 >= length2.max(axis=1)[:,None]*(1-1e-12)
                can_continue = (allowed & longest).any(axis=1)
                ready = (immediate & (~needs_more | can_continue)).reshape(-1,2).all(axis=1)
                np.minimum.at(readiness,orbit[local_edges[take,offset]],ready.astype(np.int8))
        return readiness == 1

    direct = requested & ready_orbits(requested)
    deferred = requested & ~direct
    result = {tuple(map(int,edge)):float(weight[group])
        for edge,group in zip(edges,orbit,strict=True) if direct[group]}
    fallback, closure = {}, None
    if deferred.any():
        blocked = {tuple(map(int,edge)):float(weight[group])
            for edge,group in zip(edges,orbit,strict=True) if deferred[group]}
        fallback, closure = closure_scores(points,cells,surface,blocked,periodic=periodic)
        proposed = np.zeros(len(keys),dtype=bool)
        for edge in fallback:
            proposed[orbit[lookup[edge]]] = True
        ready = ready_orbits(proposed)
        for edge,value in fallback.items():
            if ready[orbit[lookup[edge]]]:
                result[edge] = max(result.get(edge,0.),value)
        rejected = int((proposed & ~ready).sum())
    else:
        rejected = 0
    return result, dict(requested_orbits=int(requested.sum()),
        direct_requested_orbits=int(direct.sum()),deferred_requested_orbits=int(deferred.sum()),
        fallback_orbits_rejected_by_continuation=rejected,
        examined_parent_edge_pairs=examined,longest_edge_fallback=closure,
        selection='Accept requested splits only if children meet the metric or retain a legal longest-edge split; otherwise propagate to longer incident edges',
        fe_fields_used=False,quality_floor_relaxed=False)
