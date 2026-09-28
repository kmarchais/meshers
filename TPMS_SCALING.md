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
The complete microgen plus MMG workflow is compared below.

When thickness is `0.6 + (0.1/3)*(x+y+z)`, no axis can be tiled. At two units,
legacy microgen took 0.91 seconds for 68k mixed cells; the meshers path took
4.61 seconds for 133k tetrahedra with minimum MMG quality 0.113. At three
units, legacy microgen took 0.98 seconds for 234k mixed cells. With its 0.1
quality gate enabled, direct meshers took 15.05 seconds for 458k tetrahedra
with minimum quality 0.117 and 0.00524 sampled surface error. The microgen
meshers adapter produced the same mesh in 14.55 seconds. The optimizer used
five passes, one more than its default four, because the fourth pass reached
only 0.098. This is a controlled improvement over rejecting the mesh or
rerunning the entire generator. Unremeshed VTK output remains much faster;
the complete workflow is the relevant comparison for FEM-ready tetrahedra.

For an MMG 5.8 comparison, `mmgpy` 0.16.2 provided the native library.
The PyPI package named [`mmg`](https://pypi.org/project/mmg/) is a Markdown
localization tool, not the mesh remesher.
The benchmark triangulates microgen's mixed cells, marks its boundary faces
as required, and performs one MMG remeshing pass through `mmgpy`. All required
faces were retained in the one-unit output. Four alternating fresh-process
runs of the one-unit, three-axis-graded case gave median generation times,
excluding imports, of 0.350 seconds for direct meshers and 1.766 seconds
for microgen plus `mmgpy`. Meshers returned 15,831 tetrahedra with minimum
MMG quality 0.120 and first-percentile quality 0.200. `mmgpy` returned
13,941 tetrahedra with minimum quality 0.0104 and first-percentile quality
0.110. Peak process working sets were about 44 MB and 497 MB respectively.
This is a one-pass `mmgpy` comparison, not a reproduction of microgen's
existing two-pass command-line workflow; the element counts and MMG size
controls differ. The `mmgpy` output's sampled surface error has not yet been
measured against the implicit field, so this is not a matched-accuracy result.
The two-unit `mmgpy` run exceeded a 120-second wall-time limit without
returning a mesh; its final quality and runtime remain unknown.

Microgen's existing `remesh_keeping_boundaries_for_fem` path was also run
unchanged with the MMG3D 5.8.0 executable bundled in `mmgpy` 0.16.2. The
executable lives in `mmgpy/bin/Release`, and its DLL directory `mmgpy/bin`
must be on `PATH`. The function makes two MMG3D command-line passes, first
with `-nofem`, then with `-ls -nr`. In three alternating fresh-process runs,
the one-unit case had median generation times, excluding imports, of 0.392
seconds for direct meshers and 3.607 seconds for the complete microgen plus
MMG workflow. The MMG output contained 9,109 tetrahedra with minimum quality
0.0137 and first-percentile quality 0.0971. A Windows temporary-file fix in
the microgen experiment branch was needed to let meshio reopen its input
file. The workflow results do not have matched element counts or measured
MMG surface error. At two units, MMG 5.8 failed in the first `-nofem` pass
after about 31 seconds. It reported that the neighborhood of one edge had
too many elements and could not complete Delaunay adaptation. Microgen's five
identical retries all failed, so there is no two-unit output to compare.

An earlier test used an MMG3D 5.7.0 executable supplied by `pymmg` 1.0.0.
Four one-unit runs had a median of 3.43 seconds, **including about 0.6
seconds of microgen import time**. Its 9,160-tet output had minimum quality
0.0123 and first-percentile quality 0.0989. The two-unit MMG 5.7 run was
stopped after 235 seconds of MMG CPU work without producing an output mesh.
This older test is retained only as historical context. Larger matched-quality
measurements remain necessary before claiming a general speedup.

The surface path previously spent most of its time improving an already
extracted mesh: raw extraction for the two-unit graded gyroid took 0.43
seconds, while the old default of 10 smoothing iterations, 12 topology
rounds, and 40 polishing passes took 23.03 seconds and still left a 13.16°
minimum angle. Four topology rounds without smoothing or polishing are now
the nonperiodic default. At matched grid counts, the resulting graded
gyroid took 1.96 seconds for two units (103k triangles, 15.34° minimum
angle) and 6.73 seconds for three units (349k triangles, 14.25°). A two-unit
graded split-P took 4.08 seconds (167k triangles, 7.28°). These single-run
times exclude imports. The three-unit gyroid produced 3.4 times as many
triangles as the two-unit case and took 3.4 times as long in this pair of
runs. The corresponding direct gyroid volumes took about
4.61 and 15.07 seconds. Microgen's raw VTK gyroid surface took 1.21 seconds
for two units and 1.30 seconds for three, but its minimum angles were
0.0014° and 0.00053°. Thus the improved surface path is faster than direct
volume meshing and avoids VTK's nearly degenerate triangles, though raw VTK
extraction remains faster. Larger graded cases and end-to-end FEA are still
unmeasured.

A two-unit gyroid graded only along x retains periodic y/z cap meshes. With
the shorter eight-pass default for partially periodic surfaces, it took 2.80
seconds, had a 6.24° minimum angle, and matched both cap nodes and triangles
on the two periodic axes. The corresponding direct volume took 3.89 seconds
with minimum MMG quality 0.122. This surface path is still considerably
slower than the fully nonperiodic four-round path; paired topology edits or
a faster periodic optimizer are needed for a larger speed margin.

Run `crates/meshers-python/examples/graded_comparison.py` for these cases.
Its `runtime_seconds` includes geometry setup and generation but excludes
module imports and process startup. `total_seconds` includes imports but not
process startup. Except for the stated four-run medians, each measurement
above is one run. Compare output counts and quality alongside time.

## Quality-led workflow comparison

These cases compare the native meshers path with microgen's VTK surface path
and its unchanged two-pass MMG 5.8 volume path. A regular gyroid is the simple
TPMS case; split-P and three-axis grading are the harder cases. Surface
quality is triangle angle in degrees. Volume quality is the MMG tetrahedron
shape metric on `[0, 1]`; larger values are better. The first percentile
shows whether poor quality affects more than one element. Periodicity was
checked by matching both nodes and cap triangles on each pair of box faces,
not just by requesting periodic output. All times exclude module imports.

For a one-cell gyroid with constant thickness `0.6` and 16 grid points per
axis, both paths returned fully periodic meshes:

| Output | Workflow | Time | Worst quality | First-percentile quality |
| --- | --- | ---: | ---: | ---: |
| Surface | microgen VTK | 0.035 s | 0.287° | 2.51° |
| Surface | meshers direct | 0.258 s | 15.05° | 28.74° |
| Volume | microgen + two-pass MMG | 3.421 s | 0.0225 | 0.1046 |
| Volume | meshers direct | 0.896 s | 0.1440 | 0.3888 |

Volume times are medians of three alternating fresh-process runs. Surface
times are single runs, with the meshers result rerun after the surface
optimizer change. The meshers volume passed a 0.1 minimum-quality gate
and had 0.0034 sampled surface error. Microgen's MMG output had fewer
tetrahedra, but its surface error against the implicit field has not been
measured, so this is not a matched-accuracy comparison. At thickness `0.5`,
the same microgen+MMG workflow failed its first pass on the uniform gyroid
because MMG reported a zero-quality element; meshers produced a periodic
volume above the 0.1 gate.

For a one-cell split-P sheet of thickness `0.5`, both surface outputs were
periodic. At 16 grid points, microgen VTK took 0.041 seconds with a 0.259°
worst angle and 2.72° first-percentile angle. With the revised surface
optimizer, meshers took 0.449 seconds with a 7.65° worst angle and 26.57°
first-percentile angle. The meshers volume at
16 grid points missed the 0.1 quality gate (minimum 0.068), but at 32
background intervals it produced a fully periodic volume in 10.9 seconds,
with minimum quality 0.140, first-percentile quality 0.413, and 0.00218
sampled surface error. At the 16-point input resolution, microgen's first
MMG pass ran for more than 120 CPU seconds without returning an output; it
was stopped. These different resolutions and the missing MMG result do not
support a numeric split-P speed ratio.

Three-axis grading makes the field nonperiodic, so no workflow should force
periodic boundary constraints in that case. Direct meshers produced a
`(3, 3, 3)` graded gyroid volume in 15.07 seconds with minimum quality 0.117,
first-percentile quality 0.230, and 0.00524 sampled surface error. On a
`(2, 2, 2)` graded case, microgen's MMG 5.8 first pass failed with a Delaunay
adaptation error after about 31 seconds and five identical retries. The
unremeshed microgen VTK volume is faster, but it has mixed cell types and no
comparable quality gate.

With a grade only along x, both workflows retained matching nodes and cap
triangles along the still-periodic y and z axes in a one-cell test. Direct
meshers took 0.409 seconds with minimum tetrahedron quality 0.132; microgen
plus MMG took 3.398 seconds with minimum quality 0.00044. The MMG output's
first-percentile quality was 0.0875, versus 0.210 for meshers. Microgen's
wrapper was called with `periodic=False` because its option applies to all
axes, but the required boundary triangles happened to retain the two lateral
periodic cap meshes in this run.

For large **surface-only** grading, meshers now meets a 10° minimum-angle
gate for the tested two- and three-unit gyroids with its shorter default
improvement path. The two-unit graded split-P exceeds 5°. Microgen's raw
VTK surfaces remain faster, but their worst angles are nearly zero. Results
at larger sizes and matched FEA accuracy remain open.
