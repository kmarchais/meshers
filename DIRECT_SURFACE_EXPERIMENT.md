# Direct TPMS surfaces, experimental

The `experimental-surfaces` Cargo feature builds a surface-only triangle
extractor and exposes `meshers.generate_surface()` to Python. It returns
points, triangles, wall/cap labels, and generation diagnostics without building
tetrahedra. This is an experimental source build, not part of meshers 0.1.0.

```python
surface = meshers.generate_surface(
    field,
    bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
    cells=24,
    band=(-0.25, 0.25),
    periodic=(True, True, True),
)
```

The extractor accepts a scalar band in a box with one cell count for all axes.
For periodic fields, it moves opposite-face vertex copies together during
shape polishing and checks that cap triangles still match. The periodic
default uses 10 polishing passes and no smoothing or connectivity edits.
Explicit nonzero smoothing or improvement rounds are rejected because those
operations do not yet apply paired changes to opposite faces. Nonperiodic
requests keep the 10 smoothing, 12 improvement, and 40 polishing defaults.

Periodic results must have matching cap triangles and a minimum angle of at
least 5 degrees. A bad grid alignment retries up to four nearby cell counts;
`diagnostics["background_cells"]` gives the count actually used, and
`diagnostics["seconds"]` includes the retries. Zero polishing passes still
check cap matching but allow lower quality for raw inspection. If no candidate
succeeds, generation raises an error. The input field must be periodic on each
requested axis.

On the optimized local Windows build, raw `split_p` at 24 cells had 27,500
points, 55,152 triangles, no open edges, matching caps, and a 1.3-degree
minimum angle in 0.13 seconds. Ten periodic polishing passes raised the
minimum angle to 8.17 degrees in 1.33 seconds while keeping the caps matched
and the surface closed. Forty passes took 6.19 seconds and only reached 8.30
degrees. At 32 cells, the minimum angle fell to 3.45 degrees: refinement alone
does not cure slivers where the sheet meets a box edge at a shallow angle. The
two TPMS wall patches at 24 cells have minimum angles above 20 degrees; the
worst triangles are on the clipped cell caps.

A grid-aligned gyroid at 24 cells was degenerate and had unmatched caps. The
periodic retry selected 26 cells, generating matched caps with a 17.68-degree
minimum angle in 0.86 seconds for the accepted attempt. These measurements do
not establish simulation accuracy or a general quality floor.

The separate microgen experiment recovers an accepted `split_p` *volume* mesh
using meshers 0.1.0 with a larger background and intersection repair. This
feature-gated surface API has not yet been wired into microgen or released.
It cannot replace all of microgen's surface generation until the quality limit
at grazing cell-edge intersections and unsupported surface features are
addressed.

For identical repeated TPMS cells, `tile_periodic` avoids optimizing every
copy. The measurements and its limits are in `TPMS_SCALING.md`.

On Windows, all 18 core library tests and four experimental Python tests passed.
The broader meshers integration suite was not completed for this source branch.

Run the core check from the meshers repository, then install the feature build
into an active Python environment to run the Python check:

```sh
cargo test -p meshers-core --features experimental-surfaces triangles::tests
maturin develop --manifest-path crates/meshers-python/Cargo.toml --features experimental-surfaces
python -m pytest crates/meshers-python/tests/test_surface_experimental.py
```
