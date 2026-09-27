# Direct TPMS surfaces, experimental

The `experimental-surfaces` Cargo feature builds a surface-only triangle
extractor and exposes `meshers.generate_surface()` to Python. The function
returns points, triangles, wall/cap labels, and generation time. It never builds
tetrahedra. This is an experimental source build, not part of meshers 0.1.0.

```python
surface = meshers.generate_surface(
    field,
    bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
    cells=24,
    band=(-0.25, 0.25),
)
```

The current extractor accepts a scalar band in a box with one cell count for
all axes. It can smooth, improve connectivity, and polish triangle shapes.
`periodic=(True, ...)` raises `NotImplementedError`: those operations do not yet
apply the same changes to opposite faces. The extractor itself now makes the
same cap triangulation choice on opposite faces when the field is periodic.
This was verified on `split_p` at 24 cells per axis, including matching node
positions and all cap triangles on all three axes. After 10 smoothing passes,
the opposite faces no longer match. A polished surface is closed but is not a
periodic surface mesh.

The raw `split_p` surface had 27,500 points, 55,152 triangles, no open edges,
and minimum triangle angle about 1.3 degrees. With 10 smoothing passes, 12
connectivity-improvement rounds, and 40 polishing passes, it had 27,476 points,
55,104 triangles, no open edges, and minimum angle about 9.0 degrees. These
measurements do not establish simulation accuracy or a general quality floor.

The separate microgen experiment recovers an accepted `split_p` *volume* mesh
using meshers 0.1.0 with a larger background and intersection repair. This
surface API does not yet replace microgen's periodic surface mesh generator.

On Windows, all 18 core library tests and both experimental Python tests passed.
The broader meshers integration suite was not completed for this source branch.

Run the core check from the meshers repository, then install the feature build
into an active Python environment to run the Python check:

```sh
cargo test -p meshers-core --features experimental-surfaces triangles::tests
maturin develop --manifest-path crates/meshers-python/Cargo.toml --features experimental-surfaces
python -m pytest crates/meshers-python/tests/test_surface_experimental.py
```
