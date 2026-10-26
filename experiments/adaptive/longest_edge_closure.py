"""Propagate requested splits to longer incident edges before bisection."""

import numpy as np
from nested_refinement import EDGES


def closure_scores(points, cells, surface, scores, *, periodic=True):
    pairs = np.sort(cells[:, EDGES].reshape(-1, 2), axis=1)
    edges, inverse = np.unique(pairs, axis=0, return_inverse=True)
    classes = np.arange(len(points))
    if periodic:
        from periodic import periodic_pairs
        _, classes, _ = periodic_pairs(points, surface)
    keys, orbit = np.unique(np.sort(classes[edges], axis=1), axis=0, return_inverse=True)
    groups = [[] for _ in keys]
    for edge, group in enumerate(orbit):
        groups[group].append(edge)
    length2 = ((points[edges[:, 1]]-points[edges[:, 0]])**2).sum(axis=1)
    group_length = np.zeros(len(groups))
    np.maximum.at(group_length, orbit, length2)
    local_edges = inverse.reshape(-1, 6)
    order = np.argsort(inverse, kind='stable')
    bounds = np.r_[0, np.cumsum(np.bincount(inverse, minlength=len(edges)))]
    ids = {tuple(edge): i for i, edge in enumerate(edges)}
    weight = np.zeros(len(groups))
    for edge, value in scores.items():
        assert np.isfinite(value) and value >= 0
        group = orbit[ids[edge]]
        weight[group] = max(weight[group], value)
    requested = set(np.flatnonzero(weight > 0).tolist())
    pending, dependencies = list(requested), {}
    while pending:
        group = pending.pop()
        if group in dependencies:
            continue
        star = np.unique(np.concatenate([order[bounds[e]:bounds[e+1]]//6 for e in groups[group]]))
        local = local_edges[star]
        # Each dependency is strictly longer, so this graph cannot cycle.
        largest = local[np.arange(len(star)), length2[local].argmax(axis=1)]
        needs = largest[length2[largest] > group_length[group]*(1+1e-12)]
        children = set(orbit[needs].tolist())
        assert all(group_length[child] > group_length[group]*(1+1e-12) for child in children)
        dependencies[group] = children
        pending.extend(children-dependencies.keys())
    # Max avoids counting the same request repeatedly along converging paths.
    for group in sorted(dependencies, key=lambda g: group_length[g]):
        for child in dependencies[group]:
            weight[child] = max(weight[child], weight[group])
    leaves = [group for group, children in dependencies.items() if not children]
    result = {tuple(map(int, edges[e])): float(weight[group])
        for group in leaves for e in groups[group]}
    return result, dict(requested_orbits=len(requested), closure_orbits=len(dependencies),
        locally_longest_leaf_orbits=len(leaves),
        requested_already_locally_longest=sum(not dependencies[g] for g in requested),
        selection='Propagate requests to strictly longer incident edges; bisect locally longest leaves first',
        fe_fields_used=False, quality_floor_relaxed=False)
