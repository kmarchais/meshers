"""Prepare a heuristic scalar size map from a verified saved-field residual map."""

import argparse
import json
import os
from pathlib import Path

for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[name] = "1"
import numpy as np
from run import digest


def size_map(analysis, index, output, cells):
    document = json.loads((analysis / "residual-analysis.json").read_text())
    assert document["state"] == "complete"
    row = document["rows"][index]
    source = Path(row["artifact"])
    assert digest(source / "mesh.npz") == row["mesh_sha256"]
    path = analysis / f"indicators-{index}.npz"
    assert digest(path) == row["indicator_sha256"]
    mesh, values = np.load(source / "mesh.npz"), np.load(path)
    points, tets, score = mesh["points"], mesh["tetrahedra"], values["score"]
    assert len(score) == len(tets) and np.isfinite(score).all() and (score >= 0).all()
    volume = np.linalg.det(points[tets[:, 1:]] - points[tets[:, :1]]) / 6
    assert volume.min() > 0
    # Smooth P1 energy residual scales as k*h^2*volume. For 3D isotropic
    # sizing at a fixed complexity, h is proportional to k^(-1/(3+2)).
    equivalent_edge = (6 * np.sqrt(2) * volume) ** (1 / 3)
    coefficient = score / (volume * equivalent_edge**2)
    classes = np.load(source / "periodic.npz")["classes"]
    nodal_volume = np.bincount(classes[tets].ravel(), weights=np.repeat(volume, 4), minlength=len(points))
    assert np.all(nodal_volume[classes] > 0)
    nodal_coefficient = np.bincount(classes[tets].ravel(), weights=np.repeat(coefficient * volume, 4), minlength=len(points))
    nodal_coefficient = nodal_coefficient[classes] / nodal_volume[classes]
    mean = np.dot(coefficient, volume) / volume.sum()
    base, lower, upper = 1 / cells, 1 / 64, 1 / 24
    if score.sum() < 1e-20:
        sizes = np.full(len(points), base)
    else:
        sizes = np.clip(base * (mean / np.maximum(nodal_coefficient, mean * 1e-12))**.2, lower, upper)
    assert np.isfinite(sizes).all() and (sizes > 0).all()
    assert np.array_equal(sizes, sizes[classes]), "Periodic metric differs between equivalent nodes"
    np.savez_compressed(output, points=points, tetrahedra=tets, sizes=sizes, base_h=base)
    metadata = {"source": str(source), "source_mesh_sha256": row["mesh_sha256"],
                "indicator_sha256": row["indicator_sha256"], "size_map_sha256": digest(output),
                "base_h": base, "clamp": [lower, upper],
                "size_quantiles": np.quantile(sizes, [0, .05, .5, .95, 1]).tolist(),
                "scope": "Offline six-load FEA-guided isotropic heuristic; prior solves are not a one-shot generation cost"}
    output.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata))


def write_mmg_metric(path, input_mesh, points, tets, base_h, expected_hash):
    """Map cached nodal sizes to Microgen's possibly reordered MEDIT vertices."""
    from scipy.spatial import cKDTree

    assert digest(path) == expected_hash
    values = np.load(path)
    assert np.array_equal(values["points"], points) and np.array_equal(values["tetrahedra"], tets)
    assert np.isclose(float(values["base_h"]), base_h, rtol=0, atol=1e-14)
    tokens = input_mesh.read_text().split()
    start = tokens.index("Vertices")
    count = int(tokens[start + 1])
    coordinates = np.asarray(tokens[start + 2:start + 2 + 4 * count], dtype=float).reshape(count, 4)[:, :3]
    distance, order = cKDTree(points).query(coordinates)
    assert count == len(points) and len(np.unique(order)) == count and distance.max() < 1e-12
    sizes = values["sizes"][order]
    assert np.isfinite(sizes).all() and (sizes > 0).all()
    output = input_mesh.with_suffix(".sol")
    with output.open("w") as stream:
        stream.write(f"MeshVersionFormatted 2\nDimension 3\nSolAtVertices\n{count}\n1 1\n")
        np.savetxt(stream, sizes, fmt="%.17g")
        stream.write("End\n")
    return output, {"mmg_metric_h_range": [float(sizes.min()), float(sizes.max())],
                    "mmg_metric_sha256": digest(output), "mmg_size_map_sha256": expected_hash,
                    "mmg_metric_scope": "Offline six-load residual sizing; previous FEA cost excluded from generation timing"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cells", type=int, default=36)
    args = parser.parse_args()
    if not 24 <= args.cells <= 64:
        parser.error("--cells must be within the declared 24..64 sizing clamp")
    size_map(args.analysis.resolve(), args.index, args.output.resolve(), args.cells)
