# Get started

## Install a wheel

Use Python 3.10+. Wheels target Linux x86-64, Windows x86-64, macOS ARM64
and macOS x86-64.

```sh
python -m pip install "meshers[io]"
# Or: uv pip install "meshers[io]"
```

For meshing without export, install `meshers` without the extra.

`h5py` is needed only for VTKHDF export. NumPy is the only required Python runtime
dependency. Installing from source requires Rust; installing a wheel does not.

## Mesh a sphere

```python
import meshers


def sphere(x, y, z):
    return x * x + y * y + z * z - 0.7**2


mesh = meshers.generate(sphere, bounds=(-1, 1, -1, 1, -1, 1), cells=24, geometry_tolerance=0.03)
mesh.write_vtkhdf("sphere.vtkhdf")
print(mesh.diagnostics)
```

Bounds are `(xmin, xmax, ymin, ymax, zmin, zmax)`. Coordinates and tolerance use
the same physical units. `cells` is the number of lattice intervals per axis,
not a target tetrahedron count. Increase it to resolve smaller features.

The result has `points`, `tetrahedra`, boundary `surface` triangles and
`periodic_pairs`. Surface output is the closed solid boundary, including cube
faces: `mesh.write_vtkhdf("surface.vtkhdf", surface=True)`.

Continue with the [runnable examples](examples.md).
