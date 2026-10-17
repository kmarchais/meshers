"""Audit the tetrahedron-quality ceiling imposed by retained cap triangles."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
sys.path.insert(0, str(Path.home() / "bin"))

import numpy as np
from scipy.spatial import cKDTree
from pueue_progress import Job


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bounds(triangles):
    area = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                   triangles[:, 2] - triangles[:, 0]), axis=1) / 2
    squared_edges = sum(np.sum((triangles[:, first] - triangles[:, second])**2, axis=1)
                        for first, second in ((0, 1), (1, 2), (2, 0)))
    return 4 * np.sqrt(3) * area / squared_edges, squared_edges


def self_check():
    triangles = np.array([[[0, 0, 0], [1, 0, 0], [0.5, np.sqrt(3) / 2, 0]],
                          [[0, 0, 0], [1, 0, 0], [0.03, 0.002, 0]],
                          [[0, 0, 0], [2, 0, 0], [0.7, 0.4, 0]]])
    ceiling, squared_edges = bounds(triangles)
    assert abs(ceiling[0] - 1) < 1e-14
    for factor in (0.1, 0.5, 1, 2, 10):
        apex = triangles.mean(axis=1)
        apex[:, 2] = factor * np.sqrt(2 * squared_edges / 9)
        tetrahedra = np.concatenate((triangles, apex[:, None]), axis=1)
        volume = np.abs(np.linalg.det(tetrahedra[:, 1:] - tetrahedra[:, :1])) / 6
        total_squared = sum(np.sum((tetrahedra[:, first] - tetrahedra[:, second])**2, axis=1)
                            for first, second in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)))
        quality = (12 * (3 * volume)**(2 / 3) / total_squared)**1.5
        assert np.all(quality <= ceiling + 1e-14)
        if factor == 1:
            assert np.allclose(quality, ceiling, atol=1e-14, rtol=1e-13)


def caps(data):
    triangles = data["points"][data["surface"]]
    selected = np.zeros(len(triangles), dtype=bool)
    for axis in range(3):
        for value in (-0.5, 0.5):
            selected |= np.all(abs(triangles[:, :, axis] - value) < 1e-10, axis=1)
    return data["surface"][selected], triangles[selected]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    cases = ("gyroid", "schwarz-solid", "spinodoid", "smooth-lattice")
    rows, hashes = [], {}
    with Job(parts=["analytic-check"] + list(cases)) as job:
        with job.part("analytic-check") as part:
            part.update(0, "Checking the analytic bound and attaining apex", triangles=3)
            self_check()
            job.result("analytic-check", passed=1)
        for name in cases:
            with job.part(name) as part:
                part.update(0, "Reading frozen VTK and MMG caps", meshes=2)
                meshes = {}
                for route in ("vtk-default", "mmg-default"):
                    path = args.source / f"{name}-{route}-0.npz"
                    hashes[str(path)] = digest(path)
                    meshes[route] = np.load(path)
                seed_faces, seed_triangles = caps(meshes["vtk-default"])
                final_faces, final_triangles = caps(meshes["mmg-default"])
                seed_bound, _ = bounds(seed_triangles)
                final_bound, _ = bounds(final_triangles)
                distances, mapped = cKDTree(meshes["mmg-default"]["points"]).query(
                    meshes["vtk-default"]["points"])
                matched = np.all(distances[seed_faces] < 1e-10, axis=1)
                final_set = set(map(tuple, np.sort(final_faces, axis=1)))
                retained = np.array([valid and tuple(sorted(face)) in final_set
                                     for face, valid in zip(mapped[seed_faces], matched, strict=True)])
                row = {"case": name, "seed_caps": len(seed_faces), "mmg_caps": len(final_faces),
                       "seed_minimum_adjacent_tet_quality_ceiling": float(seed_bound.min()),
                       "mmg_minimum_adjacent_tet_quality_ceiling": float(final_bound.min()),
                       "seed_caps_below_quality_gate": int((seed_bound < 0.1).sum()),
                       "mmg_caps_below_quality_gate": int((final_bound < 0.1).sum()),
                       "seed_caps_retained": int(retained.sum()),
                       "retained_caps_below_quality_gate": int((retained & (seed_bound < 0.1)).sum())}
                rows.append(row)
                job.result(name, retained=row["seed_caps_retained"],
                           blocking_caps=row["retained_caps_below_quality_gate"],
                           quality_ceiling=row["mmg_minimum_adjacent_tet_quality_ceiling"])
    assert all(digest(Path(path)) == value for path, value in hashes.items())
    args.output.write_text(json.dumps({"source_sha256": digest(Path(__file__)), "inputs": hashes,
        "inputs_unchanged": True, "quality_gate": 0.1,
        "formula": "q_tet <= 4*sqrt(3)*triangle_area/sum_squared_triangle_edges",
        "scope": "Fixed cap triangles only; not a bound for changed extraction grids or repaired caps.",
        "rows": rows}, indent=2) + "\n")


if __name__ == "__main__":
    main()
