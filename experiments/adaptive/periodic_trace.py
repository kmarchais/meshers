"""Independent material-cap and complete node-class periodicity audit.

Keep the historical all-plane observer unchanged. A wall can touch a box
corner whose opposite image has no volume support. Every cap node still
requires a bijective match; all present images, including contacts, must
belong to one connected translational constraint class.
"""
import numpy as np
from scipy.spatial import cKDTree


def periodic_pairs(points, surface, tolerance=1e-10, *, require_facets=True):
    points, surface = np.asarray(points), np.asarray(surface)
    if surface.ndim != 2 or surface.shape[1] not in (3, 6):
        raise ValueError('Expected Tri3 or Tri6 surface')
    lower, upper = points.min(axis=0), points.max(axis=0)
    if not np.allclose(upper-lower, 1, atol=tolerance, rtol=0):
        raise ValueError('Periodic specimen does not span the declared unit cube')
    parent = np.arange(len(points))

    def root(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    pairs, errors, contacts, facet_checks = [], [], [], []
    for axis in range(3):
        masks = [np.all(abs(points[surface[:, :3], axis]-value) <= tolerance, axis=1)
            for value in (lower[axis], upper[axis])]
        caps = [surface[mask] for mask in masks]
        for cap, value in zip(caps, (lower[axis], upper[axis]), strict=True):
            if len(cap) and not np.all(abs(points[cap, axis]-value) <= tolerance):
                raise ValueError(f'Periodic cap node leaves axis {axis} plane')
        low, high = [np.flatnonzero(abs(points[:, axis]-value) <= tolerance)
            for value in (lower[axis], upper[axis])]
        if not len(low) or not len(high):
            raise ValueError(f'Missing opposite periodic face on axis {axis}')
        shift = np.eye(3)[axis]
        candidates = cKDTree(points[low]).query_ball_point(points[high]-shift, tolerance)
        if any(len(match) > 1 for match in candidates):
            raise ValueError(f'Ambiguous periodic nodes on axis {axis}')
        selected = np.array([i for i, match in enumerate(candidates) if match], dtype=int)
        mapped = low[np.array([candidates[i][0] for i in selected], dtype=int)]
        matched_high = high[selected]
        if len(np.unique(mapped)) != len(mapped):
            raise ValueError(f'Periodic pairing is not bijective on axis {axis}')
        missing = [np.setdiff1d(low, mapped), np.setdiff1d(high, matched_high)]
        if any(np.isin(nodes, cap).any() for nodes, cap in zip(missing, caps, strict=True)):
            raise ValueError(f'Unmatched periodic cap node on axis {axis}')
        mapping = np.arange(len(points)); mapping[matched_high] = mapped
        faces = [{tuple(sorted(row)) for row in mapping[cap]} for cap in caps]
        matched_facets=bool(faces[0]) and faces[0]==faces[1] and len(caps[0])==len(caps[1])
        facet_checks.append(matched_facets)
        if require_facets and not matched_facets:
            raise ValueError(f'Opposite-face triangulations differ on axis {axis}')
        axis_pairs = np.column_stack((mapped, matched_high))
        errors.append(float(np.max(abs(points[matched_high]-points[mapped]-shift), initial=0)))
        contacts.append([nodes.tolist() for nodes in missing])
        for lo, hi in axis_pairs:
            parent[root(hi)] = root(lo)
        pairs.append(axis_pairs)
    classes = np.array([root(i) for i in range(len(points))])
    boundary = np.any((abs(points-lower) <= tolerance) | (abs(points-upper) <= tolerance), axis=1)
    ids = np.flatnonzero(boundary)
    canonical = points[ids].copy()
    for axis in range(3):
        canonical[abs(canonical[:, axis]-upper[axis]) <= tolerance, axis] -= 1
    aliases = cKDTree(canonical).query_pairs(tolerance, output_type='ndarray')
    if len(aliases) and np.any(classes[ids[aliases[:, 0]]] != classes[ids[aliases[:, 1]]]):
        raise ValueError('Disconnected periodic edge/corner node class')
    return tuple(pairs), classes, dict(periodic_pair_counts=[len(p) for p in pairs],
        periodic_coordinate_errors=errors, periodic_pair_tolerance=tolerance,
        periodic_equivalence_classes=len(np.unique(classes)), periodic_triangles_match=all(facet_checks),
        full_facet_check_required=require_facets,
        complete_existing_node_classes=True, unmatched_noncap_contact_nodes=contacts,
        periodic_gate='material-cap-and-complete-node-classes-v1')
