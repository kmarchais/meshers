"""Sew collinear cap subdivisions without filling holes or moving vertices."""

import numpy as np
from scipy.spatial import cKDTree


def conform_cap_edges(points, faces, scale, rounds=8):
    x = points / scale
    faces = np.asarray(faces, dtype=np.int64)
    area = np.linalg.norm(
        np.cross(x[faces[:, 1]] - x[faces[:, 0]], x[faces[:, 2]] - x[faces[:, 0]]), axis=1
    )
    discarded = int((area <= 1e-15).sum())
    faces = faces[area > 1e-15]
    splits = 0
    for _ in range(rounds):
        edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
        keys, inverse, counts = np.unique(
            np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True
        )
        open_edges = keys[counts == 1]
        if not len(open_edges):
            break
        vertices = np.unique(open_edges)
        tree = cKDTree(x[vertices])
        replacements = {}
        for index in np.flatnonzero(counts[inverse] == 1):
            face_index = int(index % len(faces))
            if face_index in replacements:
                continue
            a, b = edges[index]
            vector = x[b] - x[a]
            length2 = vector @ vector
            if length2 <= 1e-24:
                continue
            candidates = vertices[
                tree.query_ball_point((x[a] + x[b]) / 2, np.sqrt(length2) / 2 + 1e-10)
            ]
            t = (x[candidates] - x[a]) @ vector / length2
            distance = np.linalg.norm(x[candidates] - x[a] - t[:, None] * vector, axis=1)
            inner = (t > 1e-10) & (t < 1 - 1e-10) & (distance <= 1e-10)
            candidates = candidates[inner][np.argsort(t[inner])]
            if not len(candidates):
                continue
            c = next(i for i in faces[face_index] if i != a and i != b)
            chain = np.r_[a, candidates, b]
            replacements[face_index] = np.c_[chain[:-1], chain[1:], np.full(len(chain) - 1, c)]
            splits += len(candidates)
        if not replacements:
            break
        faces = np.concatenate([replacements.get(i, face[None]) for i, face in enumerate(faces)])
    return faces, {
        "cap_edge_splits": splits,
        "cap_discarded_degenerate_triangles": discarded,
        "cap_stitch_tolerance_over_L": 1e-10,
    }
