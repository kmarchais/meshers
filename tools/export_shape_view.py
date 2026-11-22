# Copyright (c) 2026 Kevin Marchais

"""Freeze a fixed projection of the quality guide's Meshers 0.2.1 gyroid.

Replay with uv run --no-project --with meshers==0.2.1 tools/export_shape_view.py.
"""

import json
from pathlib import Path

import meshers as mrs
import numpy as np
from numpy.typing import NDArray


def _visible_faces(
    points: NDArray[np.float64],
    cells: NDArray[np.intp],
    eye: NDArray[np.float64],
) -> list[list[int]]:
    # Each projected triangle separates a cell behind it from a cell in front.
    # Draw it when the back cell is selected and the front one is absent or hidden.
    neighbours: dict[tuple[int, ...], list[int]] = {}
    for owner, tet in enumerate(cells):
        for opposite in range(4):
            face = tuple(sorted(int(tet[i]) for i in range(4) if i != opposite))
            a, b, c = points[list(face)]
            normal = np.cross(b - a, c - a)
            in_front = (normal @ eye) * (normal @ (points[tet[opposite]] - a)) > 0
            neighbours.setdefault(face, [-1, -1])[int(in_front)] = owner
    depth = points @ eye
    faces = [[*face, back, front] for face, (back, front) in neighbours.items() if back >= 0]
    faces.sort(key=lambda face: sum(depth[face[:3]]))
    return faces


def main() -> None:
    """Save projected vertices, face neighbours and full-precision cell scores."""
    mesh = mrs.generate(
        "gyroid",
        cells=8,
        band=(-0.5, 0.5),
        periodic=(True,) * 3,
        optimize_passes=0,
        geometry_tolerance=0.05,
    )
    points = mesh.points
    cells = mesh.tetrahedra
    corners = points[cells]
    volume = np.linalg.det(corners[:, 1:] - corners[:, :1]) / 6
    rms = np.sqrt(
        sum(
            np.square(corners[:, i] - corners[:, j]).sum(axis=1) for i in range(4) for j in range(i)
        )
        / 6,
    )
    scores = 6 * np.sqrt(2) * volume / rms**3
    np.testing.assert_allclose(scores.min(), mesh.diagnostics["minimum_mmg_quality"])

    eye = np.array([0.62, -0.76, 0.60])
    eye /= np.linalg.norm(eye)
    right = np.cross([0, 0, 1], eye)
    right /= np.linalg.norm(right)
    screen = np.column_stack((points @ right, -points @ np.cross(eye, right)))
    screen = (screen - screen.min(axis=0)) * (600 / np.ptp(screen, axis=0).max()) + 12

    faces = _visible_faces(points, cells, eye)
    target = Path(__file__).resolve().parents[1] / "docs/site/assets/gyroid-shape.json"
    target.write_text(
        json.dumps(
            {
                "points": screen.round(2).tolist(),
                "scores": scores.tolist(),
                "faces": faces,
            },
            separators=(",", ":"),
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
