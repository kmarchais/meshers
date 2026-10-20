"""Necessary pairwise lower bound for a common full-stiffness target."""

import numpy as np
from scipy.linalg import eigh, eigvalsh


def symmetric_positive(matrix):
    matrix = np.asarray(matrix, dtype=float)
    assert matrix.shape == (6, 6) and np.isfinite(matrix).all()
    assert np.linalg.norm(matrix - matrix.T) <= 1e-9 * np.linalg.norm(matrix)
    symmetric = (matrix + matrix.T) / 2
    np.linalg.cholesky(symmetric)
    return symmetric


def assess(models, target=.01):
    assert len(models) >= 2 and 0 < target < 1
    matrices = [symmetric_positive(matrix) for matrix in models]
    worst_ratio, witness = 1., None
    for first, matrix in enumerate(matrices):
        for second, reference in enumerate(matrices):
            if first == second:
                continue
            eigenvalues, vectors = eigh(matrix, reference)
            ratio = float(eigenvalues[-1])
            if ratio > worst_ratio:
                worst_ratio = ratio
                witness = first, second, vectors[:, -1]
    if witness is None:
        return {"pairwise_lower_bound": 0., "target": target, "contradiction_passed": False,
            "harmonic_trial_passes_all_models": True, "harmonic_trial_stiffness": matrices[0].tolist(),
            "harmonic_trial_maximum_disagreement": 0., "model_count": len(matrices)}
    first, second, direction = witness
    high = float(direction @ matrices[first] @ direction)
    low = float(direction @ matrices[second] @ direction)
    assert abs(high / low - worst_ratio) < 1e-10 * worst_ratio
    lower_bound = (worst_ratio - 1) / (worst_ratio + 1)
    required_lower, required_upper = (1 - target) * high, (1 + target) * low
    gap = (required_lower - required_upper) / low
    harmonic = 2 * np.linalg.solve(np.linalg.solve(matrices[first], np.eye(6))
        + np.linalg.solve(matrices[second], np.eye(6)), np.eye(6))
    harmonic = symmetric_positive(harmonic)
    disagreements = [float(abs(eigvalsh(harmonic - matrix, matrix)).max()) for matrix in matrices]
    contradiction = gap > 1e-10
    feasible = max(disagreements) <= target - 1e-10
    assert not (contradiction and feasible)
    return {"model_count": len(matrices), "target": target, "pairwise_lower_bound": lower_bound,
        "worst_ordered_pair": [first, second], "maximum_generalized_ratio": worst_ratio,
        "necessary_ratio_limit": (1 + target) / (1 - target), "mixed_strain_direction": direction.tolist(),
        "directional_high_energy": high, "directional_low_energy": low,
        "candidate_energy_required_lower": required_lower, "candidate_energy_required_upper": required_upper,
        "normalized_required_energy_gap": gap, "contradiction_margin": 1e-10,
        "contradiction_passed": contradiction, "harmonic_trial_stiffness": harmonic.tolist(),
        "harmonic_trial_model_disagreements": disagreements,
        "harmonic_trial_maximum_disagreement": max(disagreements), "harmonic_trial_passes_all_models": feasible,
        "scope": "Pairwise necessary lower bound for any symmetric candidate; passing harmonic trial proves only algebraic target feasibility, not a realizable mesh or a continuum reference"}
