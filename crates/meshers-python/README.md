# Meshers for Python

Generate tetrahedral meshes from bounded implicit solids. Use NumPy arrays directly for analysis, or install the `io` extra to export VTKHDF.

```sh
python -m pip install "meshers[io]"
```

```python
import meshers

mesh = meshers.generate("gyroid", band=(-0.5, 0.5), cells=24, periodic=(True, True, True))
points, tetrahedra = mesh.points, mesh.tetrahedra
mesh.write_vtkhdf("gyroid.vtkhdf")
```

The published 0.2.1 wheels use CPython 3.10+ and h5py for optional export. Wheels built from `main` use CPython 3.11+ and VTK for export. Source builds require Rust 1.95+.

- [Fields and grading](https://kmarchais.github.io/meshers/fields/) — built-in TPMS or a pointwise function with compiled gradients.
- [Named intersections](https://kmarchais.github.io/meshers/intersections/) — label cut faces and preserve periodic correspondence.
- [Quality and limits](https://kmarchais.github.io/meshers/quality/) — inspect orientation, retained material and sampled boundary error.
- [Python API](https://kmarchais.github.io/meshers/api/python/) — options, arrays, exceptions and cancellation.

Arrays own their storage. Editing the points or connectivity can invalidate diagnostics and boundary metadata.
