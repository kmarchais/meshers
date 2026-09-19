"""Named boundaries and feature edges, with a rotation-periodic coordinate map."""

import meshers
import numpy as np

constraints = {"right": lambda x, _y, _z: x - 0.501, "back": lambda _x, y, _z: y - 0.751}
straight = meshers.generate_intersection(
    constraints, cells=8, geometry_tolerance=1e-8, periodic=(False, False, True)
)
print(straight.surface_for("right").shape, straight.feature_edges.shape)

radius = 3.0
angle = 1 / radius


def bend(x, y, z):
    return (radius + x) * np.cos(z / radius), (radius + x) * np.sin(z / radius), -y


rotation = np.eye(4)
rotation[:2, :2] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
mesh = meshers.generate_intersection(
    constraints,
    cells=8,
    coordinate_map=bend,
    periodic=(False, False, True),
    periodic_transforms={2: rotation},
    geometry_tolerance=0.02,
    batch_size=128,
)
print(mesh.parameters.shape, mesh.constraint_names, mesh.diagnostics)
mesh.write_vtkhdf("intersection.vtkhdf")
mesh.write_vtkhdf("intersection-surface.vtkhdf", surface=True)
