"""Audit Tet4 meshes of the frozen Gyroid and Schwarz P unit cells in physical space."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import qmc

BASELINE = Path(__file__).with_name("release-wheel") / "baseline.json"
FACES = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])
AXES = "xyz"


def field(points, specimen):
    """Return the reference field values in the centered unit cell."""
    c = np.cos(2 * np.pi * points)
    if specimen == "schwarz-shifted":
        return (c[:, 0] - np.sin(2 * np.pi * points[:, 0])) / np.sqrt(2) + c[:, 1] + c[:, 2]
    if specimen == "schwarz":
        return c.sum(axis=1)
    return (np.sin(2 * np.pi * points) * np.roll(c, -1, axis=1)).sum(axis=1)


def inside(value, specimen):
    """Return the signed implicit residual, negative inside solid material."""
    return abs(value) - 0.55 if specimen == "gyroid" else value


def roots(points, specimen, level, axis):
    """Return the two analytic roots per line and a mask of lines with real roots."""
    probe = points.copy()
    values = []
    for position in (0, 0.25, 0.5):
        probe[:, axis] = position
        values.append(field(probe, specimen))
    offset = (values[0] + values[2]) / 2
    a, b = values[0] - offset, values[1] - offset
    amplitude = np.hypot(a, b)
    ratio = (level - offset) / np.maximum(amplitude, 1e-30)
    valid = (amplitude > 1e-12) & (abs(ratio) <= 1 + 1e-13)
    angle = np.arccos(np.clip(ratio, -1, 1))
    phase = np.arctan2(b, a)
    positions = (phase[:, None] + angle[:, None] * [-1, 1]) / (2 * np.pi)
    positions = (positions + 0.5) % 1 - 0.5
    return positions, valid


def forward(samples, specimen, level, axes):
    """Return a physical upper bound on distance to the bounded level set."""
    best = np.full(len(samples), np.inf)
    for axis in axes:
        positions, valid = roots(samples, specimen, level, axis)
        delta = abs(samples[:, axis, None] - positions)
        # +0.5 and -0.5 are distinct physical points with the same field.
        delta = np.where(
            abs(positions + 0.5) < 1e-14,
            np.minimum(delta, abs(samples[:, axis, None] - 0.5)),
            delta,
        )
        best[valid] = np.minimum(best[valid], delta[valid].min(axis=1))
    return best


def worst(distances, samples):
    """Return sample count and the distance and location of the worst sample."""
    if not len(samples):
        return {"samples": 0, "maximum": None, "location": None}
    idx = np.argmax(distances)
    return {
        "samples": len(samples),
        "maximum": float(distances[idx]) if np.isfinite(distances[idx]) else None,
        "location": samples[idx].tolist(),
    }


def boundary_metrics(points, faces, specimen, net):  # ruff: ignore[complex-structure] - separate wall and cap samples.
    """Return two-sided measurements per wall and cap, plus the cap face masks."""
    levels = (-0.55, 0.55) if specimen == "gyroid" else (0.0,)
    triangles = points[faces]
    caps = {}
    wall = np.ones(len(faces), dtype=bool)
    for axis in range(3):
        for side in (-0.5, 0.5):
            mask = np.all(abs(triangles[:, :, axis] - side) < 1e-10, axis=1)
            caps[axis, side] = mask
            wall &= ~mask
    wall_values = field(triangles.mean(axis=1), specimen)
    metrics = {}

    def measure(name, selection, reference, level, fixed_axis=None):
        faces_here = faces[selection]
        if not len(faces_here):
            metrics[name] = {"missing_faces": True}
            return
        patch = pv.PolyData.from_regular_faces(points, faces_here)
        _, closest = patch.find_closest_cell(reference, return_closest_point=True)
        reverse = worst(np.linalg.norm(reference - closest, axis=1), reference)
        tri = points[faces_here]
        samples = np.concatenate([
            tri.reshape(-1, 3),
            ((tri + np.roll(tri, 1, axis=1)) / 2).reshape(-1, 3),
            tri.mean(axis=1),
        ])
        if fixed_axis is None:
            distance = forward(samples, specimen, level, range(3))
        else:
            values = field(samples, specimen)
            distance = np.zeros(len(samples))
            axes = [a for a in range(3) if a != fixed_axis]
            for limit in levels:
                outside = values > limit if limit == levels[-1] else values < limit
                distance[outside] = forward(samples[outside], specimen, limit, axes)
        metrics[name] = {
            "mesh_to_reference_upper": worst(distance, samples),
            "reference_to_mesh": reverse,
        }

    for level in levels:
        refs = []
        for axis in range(3):
            lines = np.zeros((len(net), 3))
            lines[:, [i for i in range(3) if i != axis]] = net
            coords, valid = roots(lines, specimen, level, axis)
            for column in coords.T:
                lines[:, axis] = column
                refs.append(lines[valid].copy())
        reference = np.concatenate(refs)
        assert np.max(abs(field(reference, specimen) - level)) < 1e-12
        mask = wall if len(levels) == 1 else wall & ((wall_values < 0) == (level < 0))
        measure(f"wall {level:g}", mask, reference, level)
    for (axis, side), mask in caps.items():
        probes = np.full((len(net), 3), side)
        probes[:, [i for i in range(3) if i != axis]] = net
        reference = probes[inside(field(probes, specimen), specimen) <= 0]
        measure(f"{AXES[axis]}{side:+}", mask, reference, None, axis)
    return metrics, caps


def audit(grid, specimen, baseline, power=11):  # ruff: ignore[too-many-locals] - retain measured gate inputs.
    """Return physical-space measurements and failed checks for an exported volume.

    Raises:
        ValueError: The volume is empty or contains cells other than Tet4.
    """
    if set(grid.celltypes) != {pv.CellType.TETRA}:
        raise ValueError("Expected a nonempty linear tetrahedral volume")
    points = np.asarray(grid.points, dtype=float)
    cells = grid.cells_dict[pv.CellType.TETRA]
    verts = points[cells]
    det = np.linalg.det(verts[:, 1:] - verts[:, :1])
    edges = sum(((verts[:, i] - verts[:, j]) ** 2).sum(axis=1) for i in range(4) for j in range(i))
    qualities = np.sqrt(432) * det / np.maximum(edges, 1e-30) ** 1.5
    oriented = cells[:, FACES].reshape(-1, 3)
    keys = np.sort(oriented, axis=1)
    _, ids, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(ids)
    starts = np.cumsum(counts) - counts
    faces = oriented[order[starts[counts == 1]]]
    shared = starts[counts == 2]
    first, second = order[shared], order[shared + 1]
    inversions = sum(
        (oriented[:, i] > oriented[:, j]).astype(int) for i in range(3) for j in range(i + 1, 3)
    )
    graph = coo_matrix(
        (np.ones(len(shared)), (first // 4, second // 4)), shape=(len(cells), len(cells))
    )
    components, labels = connected_components(graph, directed=False)
    volumes = np.bincount(labels, weights=det / 6)
    net = qmc.Sobol(2, scramble=True, seed=119).random_base2(power) - 0.5
    metrics, caps = boundary_metrics(points, faces, specimen, net)
    periodic = {}
    for axis in range(3):

        def trace(side, axis=axis):
            coords = points[faces[caps[axis, side]]][:, :, [i for i in range(3) if i != axis]]
            return sorted(tuple(sorted(map(tuple, t))) for t in np.round(coords, 9))

        periodic[AXES[axis]] = trace(-0.5) == trace(0.5)
    probes = qmc.Sobol(3, scramble=True, seed=971).random_base2(power) - 0.5
    signed = inside(field(probes, specimen), specimen)
    # Exclude the tolerance layer using a global gradient bound of 4*pi*sqrt(3).
    limits = baseline["limits"]
    tol = limits["two_sided_sampled_distance"]
    mask = abs(signed) > 4 * np.pi * np.sqrt(3) * tol
    misses = (grid.find_containing_cell(probes[mask]) >= 0) != (signed[mask] < 0)
    volume = float(volumes.sum())
    failed = [
        name
        for name, row in metrics.items()
        if row.get("missing_faces")
        or any(result["maximum"] is None or result["maximum"] > tol for result in row.values())
    ]
    checks = {
        "periodic_faces": all(periodic.values()),
        "material_witnesses": not misses.any(),
        "signed_orientation": det.min() > 0,
        "quality": qualities.min() >= limits["minimum_quality"],
        "conformity": counts.max() <= 2
        and np.all((inversions[first] - inversions[second]) % 2 == 1),
        "volume": abs(
            volume
            / baseline["reference_volumes"][
                "schwarz" if specimen == "schwarz-shifted" else specimen
            ]
            - 1
        )
        <= limits["relative_volume_error"],
        "budget": max(len(cells), len(points)) <= limits["max_nodes_or_elements"],
        "unit_cell": np.max(abs(points)) <= 0.5 + 1e-10,
    }
    failed.extend(name for name, passed in checks.items() if not passed)
    return {
        "specimen": specimen,
        "nodes": len(points),
        "tetrahedra": len(cells),
        "volume": volume,
        "face_connected_components": int(components),
        "component_signed_volumes": volumes.tolist(),
        "minimum_quality": float(qualities.min()),
        "boundary": metrics,
        "periodic_face_match": periodic,
        "sobol_power": power,
        "material_witnesses": int(mask.sum()),
        "misclassified_witnesses": int(misses.sum()),
        "first_misclassified_witness": probes[mask][misses][0].tolist() if misses.any() else None,
        "failed_checks": failed,
        "sampled_admission": not failed,
    }


def main():
    """Return the admission exit code after recording the exported volume audit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("specimen", choices=("gyroid", "schwarz", "schwarz-shifted"))
    parser.add_argument("mesh", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.output.open("x", encoding="utf8") as stream:
        record = {"mesh_sha256": hashlib.sha256(args.mesh.read_bytes()).hexdigest()}
        try:
            record.update(
                audit(pv.read(args.mesh), args.specimen, json.loads(BASELINE.read_text()))
            )
        except Exception as exc:
            record.update(sampled_admission=False, error=str(exc))
            raise
        finally:
            json.dump(record, stream, indent=2, allow_nan=False)
            stream.write("\n")
    return 0 if record["sampled_admission"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
