"""Inspect saved meshes for a shape-quality ceiling imposed by fixed cap faces."""

import argparse
import hashlib
import json
from pathlib import Path

import meshio
import numpy as np
from periodic import tensor_difference


def triangle_quality(points, triangles):
    vertices = points[triangles]
    twice_area = np.linalg.norm(
        np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]), axis=1
    )
    edge_sum = sum(
        np.sum((vertices[:, a] - vertices[:, b]) ** 2, axis=1) for a, b in [(0, 1), (1, 2), (2, 0)]
    )
    return 2 * np.sqrt(3) * twice_area / edge_sum


def cap_faces(points, surface):
    vertices = points[surface]
    selected = np.any(
        np.all(np.isclose(vertices, -0.5, atol=1e-10, rtol=0), axis=1)
        | np.all(np.isclose(vertices, 0.5, atol=1e-10, rtol=0), axis=1),
        axis=1,
    )
    return surface[selected]


def face_keys(points, triangles):
    return {
        tuple(sorted(tuple(vertex) for vertex in triangle))
        for triangle in np.round(points[triangles], 12)
    }


def read_seed(path):
    # Meshio cannot read MMG's RequiredTriangles extension. Remove only that
    # metadata section in a separate diagnostic copy, retaining all coordinates
    # and connectivity. The original MMG input remains unchanged.
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    start = next(index for index, line in enumerate(lines) if line.strip() == "RequiredTriangles")
    count = int(lines[start + 1])
    readable = path.with_name("mmg-0-reader-copy.mesh")
    readable.write_text("".join(lines[:start] + lines[start + count + 2 :]), encoding="utf-8")
    return meshio.read(readable), count


def analyze(output):
    document = json.loads((output / "results.json").read_text(encoding="utf-8"))
    records = []
    for row in document["rows"]:
        directory = output / row["artifact"]
        path = directory / "mesh.npz"
        if not path.is_file():
            continue
        mesh = np.load(path)
        points, surface = mesh["points"], mesh["surface"]
        caps = cap_faces(points, surface)
        qualities = triangle_quality(points, caps)
        record = {
            "case": row["case"],
            "workflow": row["workflow"],
            "cells": row["cells"],
            "cap_triangles": len(caps),
            "minimum_cap_triangle_quality": float(qualities.min()),
            "cap_triangles_forcing_tet_quality_below_01": int((qualities < 0.1).sum()),
            "tetrahedra_below_01": row.get("q_below_01"),
            "ceiling_scope": "Every incident tetrahedron has MMG quality at most its fixed face's normalized triangle quality",
        }
        mmg_input = directory / "mmg-0-input.mesh"
        if mmg_input.is_file():
            seed, required_count = read_seed(mmg_input)
            seed_caps = cap_faces(seed.points, seed.cells_dict["triangle"])
            seed_quality = triangle_quality(seed.points, seed_caps)
            output_keys, seed_keys = face_keys(points, caps), face_keys(seed.points, seed_caps)
            record.update(
                seed_cap_triangles=len(seed_caps),
                seed_required_triangles=required_count,
                seed_minimum_cap_triangle_quality=float(seed_quality.min()),
                seed_cap_triangles_forcing_tet_quality_below_01=int((seed_quality < 0.1).sum()),
                cap_faces_retained_at_12_decimal_coordinates=len(output_keys & seed_keys),
                cap_face_coordinates_and_connectivity_retained=output_keys == seed_keys,
            )
        records.append(record)
    tensors = [
        row
        for row in document["rows"]
        if row.get("homogenization_audited") and row["case"] == "periodic-gyroid"
    ]
    differences = []
    for index, candidate in enumerate(tensors):
        for reference in tensors[index + 1 :]:
            differences.append(
                {
                    "candidate": [candidate["workflow"], candidate["cells"]],
                    "reference": [reference["workflow"], reference["cells"]],
                    "scope": "Diagnostic disagreement; neither tensor is a converged continuum reference",
                    **tensor_difference(
                        candidate["homogenization"]["stiffness_engineering_voigt"],
                        reference["homogenization"]["stiffness_engineering_voigt"],
                    ),
                }
            )
    result = {
        "analysis_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "periodic_source_sha256": hashlib.sha256(
            Path(__file__).with_name("periodic.py").read_bytes()
        ).hexdigest(),
        "cap_quality_analysis": records,
        "tensor_disagreements": differences,
        "bound_derivation": (
            "For a base triangle of area A and sum of squared edges S, a tetrahedron's "
            "edge sum is at least 4S/3+3h^2 and determinant is 2Ah. Maximizing "
            "q^2=1728 A^2 h^2/(4S/3+3h^2)^3 gives h^2=2S/9 and "
            "q_max=4 sqrt(3) A/S, the normalized triangle quality."
        ),
    }
    (output / "cap-analysis.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    analyze(parser.parse_args().output.resolve())
