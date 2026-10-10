"""Independent linear-tetrahedron validation and sampled geometry diagnostics."""

import itertools

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

EDGES = list(itertools.combinations(range(4), 2))
FACES = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])


def canonical_mesh(points, tets):
    p, t = np.asarray(points, float).copy(), np.asarray(tets, np.int64).copy()
    if p.ndim != 2 or p.shape[1] != 3 or t.ndim != 2 or t.shape[1] != 4:
        raise ValueError("Expected points (N,3) and linear tetrahedra (M,4)")
    if not len(t) or not np.isfinite(p).all() or t.min() < 0 or t.max() >= len(p):
        raise ValueError("Empty or malformed tetrahedral mesh")
    used, inverse = np.unique(t, return_inverse=True)
    p, t = p[used], inverse.reshape(t.shape)
    det = np.linalg.det(p[t[:, 1:]] - p[t[:, :1]])
    flipped = det < 0
    t[flipped, 1], t[flipped, 2] = t[flipped, 2].copy(), t[flipped, 1].copy()
    return p, t, int(flipped.sum()), len(points) - len(p)


def boundary(points, tets):
    faces = tets[:, FACES].reshape(-1, 3)
    _, inverse, counts = np.unique(
        np.sort(faces, axis=1), axis=0, return_inverse=True, return_counts=True
    )
    if np.any(counts > 2):
        raise ValueError("Nonmanifold tetrahedral face incidence")
    if len(np.unique(np.sort(tets, axis=1), axis=0)) != len(tets):
        raise ValueError("Duplicate tetrahedra")
    order = np.argsort(inverse, kind="stable")
    starts = np.r_[0, np.cumsum(counts)[:-1]]
    double = starts[counts == 2]
    links_i, links_j = order[double] // 4, order[double + 1] // 4
    graph = sparse.coo_matrix(
        (np.ones(len(double)), (links_i, links_j)), shape=(len(tets), len(tets))
    )
    _, components = connected_components(graph, directed=False)
    owners = np.stack([tets.ravel(), np.repeat(components, 4)], axis=1)
    unique_owners = np.unique(owners, axis=0)
    _, ownership_count = np.unique(unique_owners[:, 0], return_counts=True)
    if np.any(ownership_count != 1):
        raise ValueError("Material components touch only at an edge or vertex")
    if len(double):
        i, j = order[double], order[double + 1]
        shared = points[faces[i]]
        n = np.cross(shared[:, 1] - shared[:, 0], shared[:, 2] - shared[:, 0])
        centers = points[tets].mean(axis=1)
        left = np.einsum("ij,ij->i", n, centers[i // 4] - shared[:, 0])
        right = np.einsum("ij,ij->i", n, centers[j // 4] - shared[:, 0])
        if np.any(left * right >= 0):
            raise ValueError("Adjacent tetrahedra lie on the same side of a shared face")
    surface = faces[counts[inverse] == 1]
    edges = np.sort(
        np.concatenate([surface[:, [0, 1]], surface[:, [1, 2]], surface[:, [2, 0]]]), axis=1
    )
    _, ec = np.unique(edges, axis=0, return_counts=True)
    if np.any(ec != 2):
        raise ValueError("Boundary edges are not a closed two-manifold")
    # A pinched boundary vertex can have two closed edge fans. Edge counts
    # alone miss it, so require one connected circular link at every vertex.
    links = {}
    for a, b, c in surface.tolist():
        for vertex, left, right in [(a, b, c), (b, c, a), (c, a, b)]:
            link = links.setdefault(vertex, {})
            link.setdefault(left, set()).add(right)
            link.setdefault(right, set()).add(left)
    for link in links.values():
        if any(len(neighbors) != 2 for neighbors in link.values()):
            raise ValueError("Boundary vertex link is not a circle")
        seen, pending = set(), [next(iter(link))]
        while pending:
            node = pending.pop()
            if node not in seen:
                seen.add(node)
                pending.extend(link[node] - seen)
        if len(seen) != len(link):
            raise ValueError("Boundary vertex has disconnected fans")
    return surface


def overlapping_pair(vertices, max_pairs=2_000_000, tolerance=1e-10):
    """Strict convex-polytope SAT, with a spatial broad phase and a checked cap.

    Touching cells are allowed. This floating-point test cannot certify contacts
    closer than the declared dimensionless tolerance. A cap never means pass.
    """
    lo, hi = vertices.min(axis=1), vertices.max(axis=1)
    centers = (lo + hi) / 2
    radii = np.linalg.norm((hi - lo) / 2, axis=1)
    tree = cKDTree(centers)
    maximum = radii.max()
    checked = 0
    for i, center in enumerate(centers):
        candidates = np.asarray(tree.query_ball_point(center, radii[i] + maximum), dtype=int)
        candidates = candidates[candidates > i]
        if not len(candidates):
            continue
        overlap = np.minimum(hi[i], hi[candidates]) - np.maximum(lo[i], lo[candidates])
        candidates = candidates[np.all(overlap > tolerance, axis=1)]
        checked += len(candidates)
        if checked > max_pairs:
            raise RuntimeError(f"Overlap validation incomplete: pair cap {max_pairs} exceeded")
        for start in range(0, len(candidates), 256):
            ids = candidates[start : start + 256]
            a, b = vertices[i], vertices[ids]
            ea = np.array([a[y] - a[x] for x, y in EDGES])
            eb = np.stack([b[:, y] - b[:, x] for x, y in EDGES], axis=1)
            na = np.cross(a[FACES[:, 1]] - a[FACES[:, 0]], a[FACES[:, 2]] - a[FACES[:, 0]])
            nb = np.cross(
                b[:, FACES[:, 1]] - b[:, FACES[:, 0]], b[:, FACES[:, 2]] - b[:, FACES[:, 0]]
            )
            cross = np.cross(ea[None, :, None, :], eb[:, None, :, :]).reshape(len(ids), 36, 3)
            axes = np.concatenate([np.broadcast_to(na, (len(ids), 4, 3)), nb, cross], axis=1)
            size = np.linalg.norm(axes, axis=2)
            valid = size > 1e-15
            axes /= np.maximum(size[..., None], 1e-30)
            pa = np.einsum("kaj,vj->kav", axes, a)
            pb = np.einsum("kaj,kvj->kav", axes, b)
            interval = np.minimum(pa.max(axis=2), pb.max(axis=2)) - np.maximum(
                pa.min(axis=2), pb.min(axis=2)
            )
            hit = np.all((interval > tolerance) | ~valid, axis=1)
            if hit.any():
                return [int(i), int(ids[np.flatnonzero(hit)[0]])], checked
    return None, checked


def validate(points, tets, case, max_pairs=2_000_000):
    x = points / case.scale
    if np.any(np.abs(x) > 0.5 + 1e-8):
        raise ValueError("Mesh leaves the declared specimen box")
    vertices = x[tets]
    det = np.linalg.det(vertices[:, 1:] - vertices[:, :1])
    if np.any(det <= 1e-15):
        raise ValueError(
            "Degenerate tetrahedron at the declared 1e-15 volume determinant tolerance"
        )
    surface = boundary(x, tets)
    hit, checked = overlapping_pair(vertices, max_pairs)
    if hit:
        raise ValueError(f"Tetrahedral interiors overlap: {hit}")
    edge_sum = sum(np.sum((vertices[:, a] - vertices[:, b]) ** 2, axis=1) for a, b in EDGES)
    q = np.sqrt(432 * det**2 / edge_sum**3)
    metrics = {
        "tetrahedra": len(tets),
        "nodes": len(points),
        "volume": float(det.sum() / 6 * case.scale**3),
        "solid_fraction": float(det.sum() / 6),
        "q_min": float(q.min()),
        "q_p05": float(np.percentile(q, 5)),
        "q_median": float(np.median(q)),
        "q_below_01": int((q < 0.1).sum()),
        "overlap_pairs_checked": checked,
        "geometry_status": "checked-at-declared-tolerances",
    }
    if case.function != "solid":
        triangles = points[surface]
        normalized = triangles / case.scale
        cap = np.any(
            np.all(np.isclose(normalized, 0.5, atol=1e-9), axis=1)
            | np.all(np.isclose(normalized, -0.5, atol=1e-9), axis=1),
            axis=1,
        )
        sample = np.concatenate([triangles[~cap].reshape(-1, 3), triangles[~cap].mean(axis=1)])
        if len(sample):
            value = case.field(*sample.T)
            gradient = []
            for axis in range(3):
                shift = np.zeros(3)
                shift[axis] = case.scale * 1e-6
                gradient.append(
                    (case.field(*(sample + shift).T) - case.field(*(sample - shift).T))
                    / (2e-6 * case.scale)
                )
            size = np.linalg.norm(np.array(gradient).T, axis=1)
            if np.any(size < 1e-12 / case.scale):
                metrics["sampled_geometry_distance"] = None
            else:
                distance = np.abs(np.abs(value) - 1) / size / case.scale
                metrics["sampled_geometry_distance"] = float(distance.max())
    else:
        metrics["sampled_geometry_distance"] = 0.0
    return metrics, surface
