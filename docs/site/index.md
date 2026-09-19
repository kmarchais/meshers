![Meshers](assets/banner.svg)

# Tetrahedral meshes from implicit geometry

Meshers 0.1 turns bounded implicit solids into CPU float64 tetrahedral meshes.
Use built-in TPMS fields, custom NumPy expressions, or named constraints with a
coordinate map. The Python interface returns NumPy arrays; the Rust core can
also be used directly.

```python
import meshers

mesh = meshers.generate("gyroid", band=(-0.5, 0.5), cells=24, periodic=(True, True, True))
mesh.write_vtkhdf("gyroid.vtkhdf")
```

## Released features

- Gyroid, Schwarz P and Schwarz D, plus custom pointwise fields.
- Bands, spatial grading and blends, with automatic field derivatives.
- Optional matching periodic boundaries, including rigid mapped end transforms.
- Named intersections, constraint surfaces and feature edges.
- Quality optimization, sampled geometry checks, cancellation and resource limits.
- NumPy outputs and VTKHDF volume and boundary export.

[Install and generate a mesh](getting-started.md), then try the
[Python examples](examples.md). Read [quality and limits](quality.md) before
using a mesh in a simulation. Wheels include the native implementation and
expression compiler. Meshers is MIT licensed.
