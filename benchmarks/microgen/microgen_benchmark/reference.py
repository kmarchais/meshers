"""Fixed sampling diagnostics independent of the evaluated tetrahedral meshes."""

import numpy as np
from scipy.stats import qmc

from .preview import surface_preview


def boundary_samples(case, intervals=32):
    """Bisect grid-edge crossings on the analytic field, plus material box samples.

    This fixed sampling may miss features smaller than its grid spacing. It is
    independent of the generated meshes and is not a Hausdorff certificate.
    """
    axis = np.linspace(-0.5, 0.5, intervals + 1) * case.scale
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
    values = case.field(*grid.reshape(-1, 3).T).reshape(grid.shape[:-1])
    samples = []
    levels = [] if case.function == "solid" else ([0] if case.implicit else case.band)
    for level in levels:
        for dimension in range(3):
            left, right = [slice(None)] * 3, [slice(None)] * 3
            left[dimension], right[dimension] = slice(None, -1), slice(1, None)
            a, b = grid[tuple(left)], grid[tuple(right)]
            va, vb = values[tuple(left)] - level, values[tuple(right)] - level
            crossing = (va * vb < 0) | (va == 0) | (vb == 0)
            a, b, va = a[crossing].copy(), b[crossing].copy(), va[crossing].copy()
            for _ in range(24):
                middle = (a + b) / 2
                vm = case.field(*middle.T) - level
                same = va * vm > 0
                a[same], va[same] = middle[same], vm[same]
                b[~same] = middle[~same]
            samples.append((a + b) / 2)
    for dimension in range(3):
        for end in [0, -1]:
            face = np.take(grid, end, axis=dimension).reshape(-1, 3)
            samples.append(face[case.material(*face.T)])
    return np.unique(np.concatenate(samples), axis=0)


def target(case, output=None):
    sample = (qmc.Sobol(3, scramble=True, seed=20261010).random_base2(18) - 0.5) * case.scale
    material = case.material(*sample.T)
    volume_fraction = float(material.mean())
    coarse_fraction = float(material[: 2**17].mean())
    fractions = [volume_fraction]
    for seed in [20261011, 20261012, 20261013]:
        other = (qmc.Sobol(3, scramble=True, seed=seed).random_base2(18) - 0.5) * case.scale
        fractions.append(float(case.material(*other.T).mean()))
    average = float(np.mean(fractions))
    uncertainty = max(
        abs(np.asarray(fractions) - average).max(), abs(volume_fraction - coarse_fraction)
    )
    record = {
        "target_volume_fraction": average,
        "target_volume_sample_count": 4 * len(sample),
        "target_volume_replicates": fractions,
        "target_volume_uncertainty": float(uncertainty),
        "target_volume_uncertainty_scope": "Empirical scramble spread and nested-sample change, not a statistical confidence bound",
        "target_volume_half_sample_delta": abs(volume_fraction - coarse_fraction),
        "target_preview_intervals": 32,
    }
    if output is not None:
        import hashlib

        path = output / f"reference-{case.id}.npz"
        points = boundary_samples(case)
        np.savez_compressed(path, boundary_points=points)
        record.update(
            target_boundary_artifact=path.name,
            target_boundary_samples=len(points),
            target_boundary_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        )
    try:
        shape = case.shape(32)
        surface = case.volume_grid(shape, 32).extract_surface().triangulate()
        faces = surface.faces.reshape(-1, 4)[:, 1:]
        record["target_preview"] = surface_preview(np.asarray(surface.points), faces, case.scale)
    except Exception as exc:
        record["target_preview_error"] = f"{type(exc).__name__}: {exc}"
    return record
