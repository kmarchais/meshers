# Install and mesh a sphere

<figure class="guide-figure"><img src="assets/sphere.svg" alt="Triangulated sphere surface"><figcaption>The exported surface contains the outer triangles of the tetrahedral volume.</figcaption></figure>

```sh
python -m pip install "meshers[io]"
```

For meshing without file export, install `meshers`. It requires NumPy.

| Source | Python | VTKHDF writer installed by `[io]` |
| --- | --- | --- |
| [PyPI 0.2.1](https://pypi.org/project/meshers/0.2.1/) | 3.10+ | h5py |
| [Development checkout](contributing.md) | 3.11+ | VTK |

These docs follow the development checkout.

```python
import meshers


def sphere(x, y, z):
    return x * x + y * y + z * z - 0.7**2


mesh = meshers.generate(sphere, bounds=(-1, 1, -1, 1, -1, 1), cells=24, geometry_tolerance=0.03)
mesh.write_vtkhdf("sphere.vtkhdf")
mesh.write_vtkhdf("boundary.vtkhdf", surface=True)
```

- **Domain** — three low/high coordinate pairs. The mesher keeps the negative side inside this box.
- **Resolution** — `cells=24` gives 24 background intervals per axis. The final tetrahedron count depends on the cuts.
- **Output** — `points` holds coordinates, `tetrahedra` indexes volume cells and `surface` indexes boundary triangles. Tolerance uses the same length units as the coordinates.

[Try the examples](examples.md) · [Inspect quality and boundary tags](quality.md)
