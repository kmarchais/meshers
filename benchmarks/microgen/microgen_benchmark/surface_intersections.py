"""Independent floating-point intersection check for non-adjacent triangles.
Shared-vertex pairs are excluded. This is a tolerance-based check, not exact
arithmetic certification. Run with native triangle JSON paths.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def segment_hits(a, b, tri):
    direction = b - a
    e1 = tri[:, 1] - tri[:, 0]
    e2 = tri[:, 2] - tri[:, 0]
    h = np.cross(direction, e2)
    det = np.sum(e1 * h, axis=1)
    valid = abs(det) > 1e-12 * np.linalg.norm(direction, axis=1) * np.linalg.norm(
        e1, axis=1
    ) * np.linalg.norm(e2, axis=1)
    inv = np.divide(1, det, out=np.zeros_like(det), where=valid)
    s = a - tri[:, 0]
    u = inv * np.sum(s * h, axis=1)
    q = np.cross(s, e1)
    v = inv * np.sum(direction * q, axis=1)
    t = inv * np.sum(e2 * q, axis=1)
    return (
        valid
        & (u >= -1e-10)
        & (v >= -1e-10)
        & (u + v <= 1 + 1e-10)
        & (t >= -1e-10)
        & (t <= 1 + 1e-10)
    )


def cross2(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def coplanar_hits(a, b, normal):
    drop = np.abs(normal).argmax(1)
    axes = np.array([[1, 2], [0, 2], [0, 1]])[drop]
    a = np.take_along_axis(a, axes[:, None, :], axis=2)
    b = np.take_along_axis(b, axes[:, None, :], axis=2)
    hit = np.zeros(len(a), bool)
    for first, second in [(a, b), (b, a)]:
        for vertex in range(3):
            signs = cross2(np.roll(first, -1, axis=1) - first, second[:, vertex, None, :] - first)
            hit |= np.all(signs >= -1e-12, axis=1) | np.all(signs <= 1e-12, axis=1)
    for i in range(3):
        for j in range(3):
            p = a[:, i]
            q = a[:, (i + 1) % 3]
            r = b[:, j]
            s = b[:, (j + 1) % 3]
            hit |= (cross2(q - p, r - p) * cross2(q - p, s - p) < -1e-24) & (
                cross2(s - r, p - r) * cross2(s - r, q - r) < -1e-24
            )
    return hit


def check(p, f, max_pairs=2_000_000):
    v = p[f]
    center = v.mean(1)
    radius = np.linalg.norm(v - center[:, None, :], axis=2).max(1)
    low = v.min(1)
    high = v.max(1)
    tree = cKDTree(center)
    tested = 0
    hits = []
    for start in range(0, len(f), 512):
        stop = min(start + 512, len(f))
        neighbors = tree.query_ball_point(
            center[start:stop], radius[start:stop] + radius.max() + 1e-10
        )
        pairs = np.array(
            [(i + start, j) for i, ids in enumerate(neighbors) for j in ids if j > i + start],
            dtype=int,
        ).reshape(-1, 2)
        if not len(pairs):
            continue
        left, right = pairs.T
        keep = np.all(high[left] >= low[right] - 1e-10, axis=1) & np.all(
            high[right] >= low[left] - 1e-10, axis=1
        )
        keep &= ~np.any(f[left, :, None] == f[right, None, :], axis=(1, 2))
        left = left[keep]
        right = right[keep]
        tested += len(left)
        if tested > max_pairs:
            raise RuntimeError(
                f"Surface intersection validation incomplete: pair cap {max_pairs} exceeded"
            )
        if not len(left):
            continue
        a = v[left]
        b = v[right]
        hit = np.zeros(len(left), bool)
        for j in range(3):
            hit |= segment_hits(a[:, j], a[:, (j + 1) % 3], b) | segment_hits(
                b[:, j], b[:, (j + 1) % 3], a
            )
        normal = np.cross(a[:, 1] - a[:, 0], a[:, 2] - a[:, 0])
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        coplanar = np.all(
            abs(np.sum((b - a[:, 0, None, :]) * normal[:, None, :], axis=2)) < 1e-10, axis=1
        )
        hit[coplanar] |= coplanar_hits(a[coplanar], b[coplanar], normal[coplanar])
        hits.extend(np.column_stack([left[hit], right[hit]]).tolist())
    return dict(
        nonadjacent_candidate_pairs=tested,
        intersecting_pairs=len(hits),
        first_pairs=hits[:20],
        coordinate_tolerance=1e-10,
    )


if __name__ == "__main__":
    for path in sys.argv[1:]:
        d = json.loads(Path(path).read_text())
        result = check(np.array(d["points"]), np.array(d["faces"]))
        print(json.dumps(dict(mesh=path, **result)), flush=True)
        assert not result["intersecting_pairs"]
