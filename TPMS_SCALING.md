# Scaling repeated TPMS meshes

`meshers.tile_periodic()` copies one validated periodic cell, welds shared nodes,
and removes internal clipping caps. It accepts a direct `Mesh` or the
feature-gated `SurfaceMesh`:

```python
unit = meshers.generate_surface(
    field,
    bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
    cells=26,
    band=(-0.25, 0.25),
    periodic=(True, True, True),
)
large = meshers.tile_periodic(
    unit,
    bounds=(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
    repeats=(3, 3, 3),
)
```

Use this only when the field, band, and unit-cell geometry repeat exactly.
Grading, changing thickness, and other nonrepeating cases still need direct
generation. The helper checks matching source caps and joins seam nodes. It
preserves a volume cell's minimum MMG quality and sampled surface error because
it does not alter any tetrahedron. It also rebuilds periodic pairs on the outer
box. It rejects mapped and intersection meshes whose extra constraint data it
cannot yet carry across cells.

One optimized Windows build, compiled fields, one worker, split-P band
`(-0.25, 0.25)`, 16 background cells per unit length:

| Mesh | Cells on each side | Elements | Direct whole box | Tile one cell | Tiled peak working set |
| --- | ---: | ---: | ---: | ---: | ---: |
| Surface, 10 polish passes | 1 | 25,008 triangles | 0.62 s | unit source | 40 MB |
| Surface, 10 polish passes | 2 | 191,616 triangles | 5.17 s | 0.68 s | 61 MB |
| Surface, 10 polish passes | 3 | 637,200 triangles | 21.25 s | 0.79 s | 130 MB |
| Surface, 10 polish passes | 4 | 1,499,136 triangles | not measured | 1.40 s | 291 MB |
| Volume, 4 optimize passes | 1 | 23,638 tetrahedra | 1.99 s | unit source | 43 MB |
| Volume, 4 optimize passes | 2 | about 190,000 tetrahedra | 20.49 s | 2.09 s | 70 MB |
| Volume, 4 optimize passes | 3 | about 630,000 tetrahedra | 76.19 s | 2.18 s | 148 MB |
| Volume, 4 optimize passes | 4 | 1,512,832 tetrahedra | not measured | 2.53 s | 304 MB |

Tiled and direct volume meshes can have different element counts: direct
whole-box optimization changes the interior seams, while tiling repeats the
validated cell exactly. The 16-cell split-P volume has minimum MMG quality
0.081, so these timings alone do not prove sufficient quality. At 32 cells per
unit length, the split-P unit cell reached 0.140 minimum MMG quality, had no
tetrahedra below 0.1, and had 0.00218 sampled surface error. Tiling that cell
to `(2, 2, 2)` produced 1,041,440 tetrahedra in 10.73 seconds with the same
quality and error. A gyroid surface at 26 cells per unit length tiled to
`(3, 3, 3)` produced 1,036,368 triangles in 1.26 seconds with a 17.68-degree
minimum angle.

Without optimization, direct extraction scaled nearly with output size in
these tests: split-P surface generation took 0.06 seconds for 25,008 triangles
and 1.50 seconds for 637,200; volume generation took 0.08 seconds for 23,634
tetrahedra and 2.44 seconds for 619,326. Optimization dominates larger direct
runs. Active-neighborhood bookkeeping reduced direct surface time at three
cells per side from 22.13 to 21.25 seconds without changing measured quality.
That remains too slow to call the direct optimized method efficient at every
size. More work is needed for large nonrepeating or graded domains.

The benchmark is `crates/meshers-python/examples/tpms_scaling.py`. Each row
above is one fresh-process run, not a median. Peak memory is process working
set on Windows. The listed times include generation and tiling but exclude
field compilation and process startup.

## One-axis density grading

A sheet gyroid with thickness `0.6 + 0.1*x` is nonperiodic along x, but still
periodic along y and z. The two isosurfaces are
`gyroid(2*pi*x, 2*pi*y, 2*pi*z) = +/- thickness/2`. Meshing the entire box
repeats the expensive optimization everywhere. Meshing a graded x slab that is
one unit wide in y and z, then tiling it only along y and z, preserves the
grade and gives matching lateral seams. A grade that also depends on y or z
cannot use that tiling direction.

One Windows release build, 16 grid points per unit length, four volume
optimization passes, one fresh process per measurement:

| Domain | Microgen legacy volume | Microgen to meshers volume | Direct meshers volume | Meshers graded slab plus lateral tiling |
| --- | ---: | ---: | ---: | ---: |
| 2 x 2 x 2 | 0.88 s, 68k mixed cells | 4.64 s, 133k tetrahedra | 4.20 s, 133k tetrahedra | 2.35 s, 129k tetrahedra |
| 3 x 3 x 3 | 0.87 s, 231k mixed cells | 13.57 s, 458k tetrahedra | 12.73 s, 458k tetrahedra | 3.87 s, 442k tetrahedra |

The direct meshers volume column uses the same pair of implicit constraints
as the microgen adapter. The slab uses one normalized band field so it can be
tiled by the current API. Its minimum MMG quality was 0.172 at two units and
0.174 at three; the full direct mesh had 0.122 and 0.138 respectively. The
slab's sampled surface error was 0.0043 and 0.0037, below the requested 0.01.
These are different tetrahedralizations of the same geometric sheet, so the
times are not a controlled measure of the optimizer alone. In particular,
the legacy microgen volume has mixed cell types and no MMG-quality gate.
The MMG executable was unavailable, so an end-to-end microgen-plus-MMG timing
was not measured.

When thickness is `0.6 + (0.1/3)*(x+y+z)`, no axis can be tiled. At two units,
legacy microgen took 0.91 seconds for 68k mixed cells; the meshers path took
4.61 seconds for 133k tetrahedra with minimum MMG quality 0.113. At three
units, legacy microgen took 0.98 seconds for 234k mixed cells; microgen's
meshers adapter took 14.15 seconds for 458k tetrahedra with minimum quality
0.098. Calling meshers directly took 13.95 seconds for the identical mesh.
The three-unit output would fail the adapter's default 0.1 quality gate; the
benchmark disabled that gate to measure its raw generation time. These results
do not support a general speedup for fully graded domains.

The direct periodic surface path has a separate scaling issue. At two units,
10 polish passes did not reach its 5-degree angle requirement; 40 passes
reached 13.16 degrees but took 23.03 seconds for 116k triangles. At three
units, even 40 passes and five background-resolution attempts did not meet
the 5-degree requirement. Microgen's
legacy surface took 1.21 seconds for 57k triangles, but its minimum angle was
0.0014 degrees. At three units its 190k-triangle surface took 1.30 seconds
and had a 0.00053-degree minimum angle. A fast, high-quality direct surface
for a large graded domain is not yet established.

Run `crates/meshers-python/examples/graded_comparison.py` for these cases.
The timings include geometry setup and generation, but exclude module import
and process startup. Each is one run, not a median; both output counts and
surface angles should be considered alongside time.
