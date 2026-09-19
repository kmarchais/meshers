"""Built-in TPMS, negative-side solids, arrays and volume/surface export."""

import meshers

for field in ("gyroid", "schwarz_p", "schwarz_d"):
    mesh = meshers.generate(
        field, band=(-0.5, 0.5), cells=16, geometry_tolerance=0.1, optimize_passes=1
    )
    print(field, mesh.points.shape, mesh.tetrahedra.shape, mesh.diagnostics)


def sphere(x, y, z):
    return x * x + y * y + z * z - 0.7**2


mesh = meshers.generate(
    sphere, bounds=(-1, 1, -1, 1, -1, 1), cells=(16, 18, 20), geometry_tolerance=0.05
)
mesh.write_vtkhdf("sphere.vtkhdf")
mesh.write_vtkhdf("sphere-surface.vtkhdf", surface=True)
print(mesh.points.shape, mesh.tetrahedra.shape, mesh.surface.shape)
