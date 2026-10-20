"""Bounded fixed-vertex cavity enumeration, not a general volume mesher."""

from itertools import combinations

import numpy as np

from bistellar_probe import oriented_boundary


def replacement_boundary(points, original, old_faces, new_faces):
    before = oriented_boundary(original)
    old_keys = {tuple(sorted(map(int, face))) for face in old_faces}
    assert before is not None and old_keys <= set(before)
    normal = sum(before[face] * np.cross(points[face[1]] - points[face[0]],
        points[face[2]] - points[face[0]]) for face in old_keys)
    target = {key: value for key, value in before.items() if key not in old_keys}
    for face in new_faces:
        key = tuple(sorted(map(int, face)))
        direction = float(np.dot(normal, np.cross(points[key[1]] - points[key[0]], points[key[2]] - points[key[0]])))
        assert direction != 0 and key not in target
        target[key] = 1 if direction > 0 else -1
    return target


def enumerate_cavity(points, original, target):
    vertices = np.unique(original)
    assert len(vertices) <= 6, "Bound enumeration to six existing vertices"
    candidates = np.array(list(combinations(map(int, vertices), 4)))
    determinants = np.linalg.det(points[candidates[:, 1:]] - points[candidates[:, :1]])
    candidates = candidates[determinants != 0]
    determinants = determinants[determinants != 0]
    negative = determinants < 0
    candidates[negative, 1], candidates[negative, 2] = candidates[negative, 2].copy(), candidates[negative, 1].copy()
    assert len(candidates) <= 15
    boundaries = [oriented_boundary(cell[None]) for cell in candidates]
    keys = sorted(set(target).union(*(set(boundary) for boundary in boundaries)))
    matrix = np.array([[boundary.get(key, 0) for boundary in boundaries] for key in keys], dtype=np.int16)
    desired = np.array([target.get(key, 0) for key in keys])
    masks = ((np.arange(1, 1 << len(candidates), dtype=np.uint16)[:, None] >> np.arange(len(candidates))) & 1).astype(np.int16)
    matching = np.flatnonzero(np.all(masks @ matrix.T == desired, axis=1))
    order = matching[np.argsort(masks[matching].sum(axis=1), kind="stable")]
    old_volume = float(np.linalg.det(points[original[:, 1:]] - points[original[:, :1]]).sum())
    assert old_volume > 0
    report = {"vertices": vertices.tolist(), "candidate_tetrahedra": len(candidates),
        "subsets_checked": len(masks), "oriented_boundary_matches": len(matching), "passed": False}
    for index in order:
        selected = candidates[masks[index].astype(bool)]
        volume = float(np.linalg.det(points[selected[:, 1:]] - points[selected[:, :1]]).sum())
        if abs(volume - old_volume) / old_volume >= 1e-10:
            continue
        faces = np.sort(selected[:, np.array(((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)))].reshape(-1, 3), axis=1)
        _, counts = np.unique(faces, axis=0, return_counts=True)
        if counts.max() > 2:
            continue
        assert oriented_boundary(selected) == target
        report.update(passed=True, alternative_tetrahedra=len(selected),
            relative_volume_difference=abs(volume - old_volume) / old_volume)
        return selected, report
    return None, report
