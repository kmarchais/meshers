"""Matching box faces and explicit quality, resource and threading controls."""

import meshers
import numpy as np

mesh = meshers.generate(
    "gyroid",
    band=(-0.5, 0.5),
    cells=16,
    periodic=(True, True, True),
    geometry_tolerance=0.1,
    minimum_quality=0.05,
    max_tetrahedra=500_000,
    optimize_passes=4,
    snap=0.2,
    threads=1,
)
for axis, pairs in enumerate(mesh.periodic_pairs):
    shift = np.zeros(3)
    shift[axis] = 1
    assert np.allclose(
        mesh.points[pairs[:, 1]] - mesh.points[pairs[:, 0]], shift, rtol=0, atol=1e-10
    )
    print(axis, len(pairs), "node pairs")
mesh.write_vtkhdf("periodic.vtkhdf")
