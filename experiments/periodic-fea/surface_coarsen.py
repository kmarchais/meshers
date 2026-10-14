"""Conservative surface-edge contraction prototype; no meshing speed claim.

Boundary-coned volume and surface links follow the manifold conditions in
Dey et al., Topology Preserving Edge Contraction. Independent global validation
remains required. Only disjoint vertex cavities are changed in this pass.
"""
from itertools import combinations
import time

import numpy as np


def link(facets, simplex):
    simplex = set(simplex)
    result = set()
    for facet in facets:
        vertices = set(map(int, facet))
        if simplex <= vertices:
            remaining = sorted(vertices - simplex)
            for size in range(len(remaining) + 1):
                result.update(combinations(remaining, size))
    return result


def contractible(volume, surface, a, b):
    # Cone the boundary to a dummy vertex for the volume condition. The
    # boundary itself also has to satisfy its own link condition.
    augmented = list(volume) + [tuple(map(int, f)) + (-1,) for f in surface]
    return all(link(facets, [a]) & link(facets, [b]) == link(facets, [a, b])
               for facets in [augmented, surface])


def incidence(cells, count):
    vertices = cells.ravel()
    order = np.argsort(vertices, kind="stable")
    starts = np.r_[0, np.cumsum(np.bincount(vertices, minlength=count))]
    owners = np.repeat(np.arange(len(cells)), cells.shape[1])
    return lambda v: owners[order[starts[v]:starts[v + 1]]]


def coarsen(points, tets, surface, node_budget, max_edge_length, minimum_quality=.1,
            *, protected=None, levels=None, geometry_check=None, checkpoint=lambda **_: None):
    """Contract short wall edges to existing endpoints, preserving protected vertices.

    Returns retained original vertex IDs as well as compacted arrays. A zero
    level marks ineligible vertices. geometry_check tests newly changed wall
    triangles. Caps are protected automatically in the unit-cell prototype.
    """
    count = len(points)
    caps = np.any(np.isclose(np.abs(points), .5, atol=1e-9, rtol=0), axis=1)
    protected = caps if protected is None else caps | np.asarray(protected, bool)
    levels = np.ones(count, int) if levels is None else np.asarray(levels)
    edges = np.unique(np.sort(np.concatenate([surface[:, [0, 1]], surface[:, [1, 2]],
                                             surface[:, [2, 0]]]), axis=1), axis=0)
    length2 = np.sum((points[edges[:, 0]] - points[edges[:, 1]])**2, axis=1)
    eligible = (~protected[edges].any(axis=1) & (levels[edges[:, 0]] != 0)
                & (levels[edges[:, 0]] == levels[edges[:, 1]])
                & (length2 > 1e-20) & (length2 <= max_edge_length**2))
    candidates = np.flatnonzero(eligible)
    tie = np.random.default_rng(20261014).permutation(len(candidates))
    candidates = candidates[np.lexsort((tie, length2[candidates]))]
    tet_incident, face_incident = incidence(tets, count), incidence(surface, count)
    changed = tets.copy()
    alive = np.ones(len(tets), bool)
    blocked = np.zeros(count, bool)
    removed = []
    rejected = {"link": 0, "quality": 0, "geometry": 0}
    tested = 0
    considered = 0
    last_report = 0.
    for index, ei in enumerate(candidates):
        if len(removed) >= node_budget:
            break
        considered = index + 1
        now = time.monotonic()
        if now - last_report >= 1:
            checkpoint(fraction=index / max(len(candidates), 1), considered=index,
                       candidates=len(candidates), removed=len(removed), tested=tested)
            last_report = now
        a, b = map(int, edges[ei])
        if blocked[[a, b]].any():
            continue
        star = np.union1d(tet_incident(a), tet_incident(b))
        vertices = np.unique(tets[star])
        if blocked[vertices].any():
            continue
        faces = np.union1d(face_incident(a), face_incident(b))
        tested += 1
        if not contractible(tets[star], surface[faces], a, b):
            rejected["link"] += 1
            continue
        accepted = False
        for remove, keep in [(b, a), (a, b)]:
            incident = tet_incident(remove)
            old = tets[incident]
            vanished = np.any(old == keep, axis=1)
            replacement = old[~vanished].copy()
            replacement[replacement == remove] = keep
            if not len(replacement):
                continue
            xyz = points[replacement]
            determinant = np.linalg.det(xyz[:, 1:] - xyz[:, :1])
            edge_sum = sum(np.sum((xyz[:, i] - xyz[:, j])**2, axis=1)
                           for i, j in combinations(range(4), 2))
            quality = np.sqrt(432 * determinant**2 / edge_sum**3)
            if not np.isfinite(quality).all() or determinant.min() <= 0 or quality.min() < minimum_quality:
                rejected["quality"] += 1
                continue
            old_faces = surface[face_incident(remove)]
            new_faces = old_faces[~np.any(old_faces == keep, axis=1)].copy()
            new_faces[new_faces == remove] = keep
            if geometry_check is not None and not geometry_check(points[new_faces]):
                rejected["geometry"] += 1
                continue
            changed[incident[~vanished]] = replacement
            alive[incident[vanished]] = False
            blocked[vertices] = True
            removed.append(remove)
            accepted = True
            break
        if accepted:
            assert not protected[removed[-1]]
    survivors = changed[alive]
    retained, inverse = np.unique(survivors, return_inverse=True)
    assert count - len(retained) == len(removed), "Contraction removed an unintended vertex"
    assert not np.intersect1d(retained, removed).size
    checkpoint(fraction=1., considered=considered, candidates=len(candidates),
               removed=len(removed), tested=tested)
    return points[retained].copy(), inverse.reshape(survivors.shape), retained, {
        "removed_nodes": len(removed), "eligible_edges": len(candidates),
        "tested_cavities": tested, "rejections": rejected,
        "scope": "One disjoint surface-edge contraction pass on an archived mesh",
    }
