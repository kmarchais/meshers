"""Independent translational node, face and edge/corner equivalence checks."""

import numpy as np
from scipy.spatial import cKDTree


def periodic_pairs(points, surface, tolerance=1e-10):
    """Match all opposite unit-box faces, including their triangulations.

    Coordinates may be centered or translated. Return pairs and equivalence
    classes; reject unmatched, duplicate and ambiguous matches.
    """
    points = np.asarray(points)
    lower, upper = points.min(axis=0), points.max(axis=0)
    if not np.allclose(upper - lower, 1, atol=tolerance, rtol=0):
        raise ValueError("Periodic specimen does not span the declared unit cube")
    parent = np.arange(len(points))

    def root(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    pairs, errors = [], []
    for axis in range(3):
        low = np.flatnonzero(np.abs(points[:, axis] - lower[axis]) <= tolerance)
        high = np.flatnonzero(np.abs(points[:, axis] - upper[axis]) <= tolerance)
        if not len(low) or len(low) != len(high):
            raise ValueError(f"Opposite-face node counts differ on axis {axis}")
        shift = np.eye(3)[axis]
        matches = cKDTree(points[low]).query_ball_point(points[high] - shift, tolerance)
        if any(len(match) != 1 for match in matches):
            raise ValueError(f"Unmatched or ambiguous periodic nodes on axis {axis}")
        mapped = low[np.array([match[0] for match in matches])]
        if len(np.unique(mapped)) != len(low):
            raise ValueError(f"Periodic pairing is not bijective on axis {axis}")
        axis_pairs = np.column_stack([mapped, high])
        errors.append(float(np.max(np.abs(points[high] - points[mapped] - shift))))
        mapping = np.arange(len(points))
        mapping[high] = mapped
        faces = []
        for value in [lower[axis], upper[axis]]:
            selected = np.all(np.abs(points[surface, axis] - value) <= tolerance, axis=1)
            faces.append({tuple(sorted(row)) for row in mapping[surface[selected]]})
        if not faces[0] or faces[0] != faces[1]:
            raise ValueError(f"Opposite-face triangulations differ on axis {axis}")
        for low_node, high_node in axis_pairs:
            parent[root(high_node)] = root(low_node)
        pairs.append(axis_pairs)
    classes = np.array([root(i) for i in range(len(parent))])
    return (
        tuple(pairs),
        classes,
        {
            "periodic_pair_counts": [len(pair) for pair in pairs],
            "periodic_coordinate_errors": errors,
            "periodic_pair_tolerance": tolerance,
            "periodic_equivalence_classes": len(np.unique(classes)),
            "periodic_triangles_match": True,
        },
    )


def tensor_difference(stiffness, reference):
    """Worst relative energy difference over all macroscopic strain combinations.

    The reference must be positive definite. This is reference disagreement,
    not a continuum-error certificate or a fine/coarse error bound.
    """
    from scipy.linalg import eigvalsh

    stiffness, reference = np.asarray(stiffness), np.asarray(reference)
    if stiffness.shape != (6, 6) or reference.shape != (6, 6):
        raise ValueError("Expected full six-load-case stiffness tensors")
    if not np.isfinite(stiffness).all() or not np.isfinite(reference).all():
        raise ValueError("Nonfinite stiffness tensor")
    if not np.allclose(stiffness, stiffness.T, atol=1e-10, rtol=1e-8):
        raise ValueError("Candidate stiffness is not symmetric")
    if not np.allclose(reference, reference.T, atol=1e-10, rtol=1e-8):
        raise ValueError("Reference stiffness is not symmetric")
    reference = (reference + reference.T) / 2
    if eigvalsh(reference).min() <= 0:
        raise ValueError("Reference stiffness must be positive definite")
    spectrum = eigvalsh((stiffness + stiffness.T) / 2 - reference, reference)
    return {
        "worst_relative_strain_energy_difference": float(np.abs(spectrum).max()),
        "relative_tensor_frobenius_difference": float(
            np.linalg.norm(stiffness - reference) / np.linalg.norm(reference)
        ),
        "relative_energy_difference_eigenvalues": spectrum.tolist(),
    }
