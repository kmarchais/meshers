"""A TPMS zero set may pass through background vertices on periodic faces."""

import meshers
import numpy as np
import pytest


def schwarz_d(x, y, z):
    x, y, z = (2 * np.pi * a for a in (x, y, z))
    return (np.sin(x) * np.sin(y) * np.sin(z)
        + np.sin(x) * np.cos(y) * np.cos(z)
        + np.cos(x) * np.sin(y) * np.cos(z)
        + np.cos(x) * np.cos(y) * np.sin(z))


@pytest.mark.parametrize("order,cells", [(1, 24), (2, 16)])
def test_schwarz_d_zero_solid_has_periodic_nodes_and_full_cap_facets(order, cells):
    mesh = meshers.generate(schwarz_d, bounds=(-.5, .5) * 3, cells=cells,
        periodic=(True,) * 3, threads=1, order=order, geometry_tolerance=1.,
        optimize_passes=1, minimum_quality=.1, max_tetrahedra=250000)
    points, surface = np.asarray(mesh.points), np.asarray(mesh.surface)
    assert len(points) <= 250000 and len(mesh.tetrahedra) <= 250000
    for axis, pairs in enumerate(mesh.periodic_pairs):
        assert len(pairs) > 0
        np.testing.assert_allclose(points[pairs[:, 1]] - points[pairs[:, 0]],
            np.broadcast_to(np.eye(3)[axis], (len(pairs), 3)), atol=1e-10)
        caps = []
        for side in (-.5, .5):
            faces = surface[np.all(points[surface[:, :3], axis] == side, axis=1)]
            coordinates = points[faces].copy()
            coordinates[:, :, axis] = 0
            caps.append({tuple(sorted(tuple(p) for p in np.round(face, 10))) for face in coordinates})
        assert caps[0] and caps[0] == caps[1]
