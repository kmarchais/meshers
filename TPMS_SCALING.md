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
