# Adaptive and quadratic implicit meshing: results

Branch `claude/adaptive-implicit`. All numbers are element counts, quality, geometry
and FEA error; none is a timing claim. Timings wait for an idle `gpu` queue. Every
geometry metric is independent of the mesher's own diagnostics (`scorecard.py`):

- mesh to surface: Newton foot-point distance at 13 points per wall triangle;
- surface to mesh: exact point-to-triangle distance (VTK) from 2^17 projected samples;
- volume: analytic, or four scrambled Sobol sets of 2^23 points (uncertainty < 0.002%);
- quality: triangle shape 4√3A/Σl², minimum angle; MMG tet quality, minimum dihedral.

Cases: sphere r=0.35, ellipsoid (0.42, 0.25, 0.12), Gyroid sheet with raw band ±0.55
(the periodic-fea fixture), all in [-0.5, 0.5]^3, float64.

## Surfaces: curvature-adaptive remeshing (`remesh.rs`)

Target edge length from the chord sag of the exact surface, with gradation 0.3.
Split, collapse, flip and tangential relaxation, with exact projection onto walls,
box planes and their intersection curves. Matched Hausdorff distance H, against VTK
Flying Edges plus box clipping (the PR #2 control):

| Case | Flying Edges | Adaptive | Triangles |
|---|---|---|---|
| Sphere, H≈4e-4 | 10,616 tri, min angle 0.6° | 6,460 tri, min angle 33° | -39% |
| Ellipsoid, H≈4.5e-4 | 21,480 tri, min angle 0.2° | 7,234 tri, min angle 28° | -66% |
| Gyroid, H≈9e-4 | 44,938 tri, min angle 0.0°, ΔV +0.23% | 25,576 tri, min angle 18.7°, ΔV +0.10% | -43% |

With `balanced=True` the wall vertices move outward by half the local centroid sag,
so triangles deviate to both sides. Sizing then allows about twice the sag and the
measured deviation is still checked against the tolerance. Same independent H:

| Case | Adaptive | Balanced | Flying Edges at similar H |
|---|---:|---:|---:|
| Sphere, H≈4.2e-4 | 6,570 | 3,446 | 10,616 |
| Ellipsoid, H≈4.5e-4 | 7,194 | 4,130 | 21,480 |
| Gyroid, H≈9.8e-4 | 25,576 | 15,920 | 44,938 |

Balanced vertices are off the surface, so this mode is for printing and display;
curved quadratic FEA meshes keep vertices on the surface.

### Fast mode (`fast=True`): fewer triangles at Flying Edges speed

For printing and display, where triangle shape does not matter: collapse the 24-cell
extraction towards the size field (two passes, shape guard 0.05 instead of 0.3, no
flips or smoothing), then split only where the error exceeds the tolerance.
Supporting changes, which also apply to the quality modes:

- the error check per triangle is the exact maximum of the quadratic through its
  vertex and mid-edge samples (it was the largest sample, which read up to 2.7% low
  against the independent reverse distance on irregular triangles);
- balanced mode repairs locally: after balancing, faces still above the tolerance
  are split and balanced again (it used to tighten the tolerance everywhere, and
  in fast mode failed and fell back to unbalanced);
- lattice extraction solves all edge crossings in parallel and triangulates slabs
  in parallel with deterministic numbering; error scans, balancing and the final
  error check are parallel and incremental; manifold validation sorts instead of
  hashing.

Tolerance 1e-3, independent scorecard. Counts, H, ΔV and angles are deterministic.
Times: `timing_job.sh`, pueue job 1676 alone in the `gpu` queue on the OMEN (Intel
family 6 model 198, 24 threads), release build of 9dad555f frozen in `builds/`.
Compiled fields through Python, one warm-up, then the median of 5 runs interleaved
with all other configurations (`timing-threads1.json`, `timing-threads24.json`, with
every run). The pip VTK wheel uses the sequential SMP backend, so Flying Edges is
single-threaded in both columns; its two columns differ by run-to-run noise only.

| Case | Mesh | Triangles | H | ΔV | Min angle | 1 thread | 24 threads |
|---|---|---:|---:|---:|---:|---:|---:|
| Gyroid | Flying Edges 40 | 31,034 | 1.20e-3 | +0.34% | 0.0° | 0.107 s | 0.122 s |
| | Flying Edges 48 | 44,938 | 9.3e-4 | +0.23% | 0.0° | 0.175 s | 0.179 s |
| | fast | 26,780 | 9.9e-4 | +0.14% | 1.1° | 0.119 s | 0.086 s |
| | fast, balanced | 19,616 | 1.0e-3 | +0.08% | 1.8° | 0.125 s | 0.078 s |
| | balanced (quality) | 15,802 | 1.0e-3 | +0.05% | 18.0° | 0.226 s | 0.206 s |
| Schwarz P | Flying Edges 36 | 19,360 | 1.16e-3 | +0.26% | | 0.061 s | 0.072 s |
| | Flying Edges 40 | 24,064 | 9.5e-4 | +0.20% | 0.0° | 0.081 s | 0.091 s |
| | fast | 16,872 | 9.8e-4 | +0.19% | 0.0° | 0.076 s | 0.055 s |
| | fast, balanced | 11,880 | 1.0e-3 | +0.12% | 0.0° | 0.073 s | 0.048 s |
| | balanced (quality) | 8,946 | 1.0e-3 | +0.09% | 18.9° | 0.146 s | 0.116 s |
| Sphere | Flying Edges 32 | 4,808 | 1.01e-3 | -0.47% | | 0.016 s | 0.016 s |
| | fast, balanced | 1,954 | 9.9e-4 | -0.26% | 8.8° | 0.018 s | 0.015 s |
| Ellipsoid | Flying Edges 64 | 9,464 | 1.08e-3 | -0.41% | | 0.047 s | 0.049 s |
| | Flying Edges 80 | 14,856 | 6.9e-4 | -0.26% | | 0.077 s | 0.089 s |
| | fast, balanced | 2,950 | 1.0e-3 | -0.34% | 0.1° | 0.016 s | 0.013 s |

- Against the Flying Edges grid that meets H ≤ 1e-3 (48, 40, 32, 80), fast balanced
  has 2.0-5.0x fewer triangles and a smaller volume error except on the ellipsoid
  against grid 80. Speed: gyroid 1.4x faster on one thread and 2.3x on 24; Schwarz P
  1.1x and 1.9x; ellipsoid 4.7x and 6.8x (against grid 64, which misses 1e-3:
  2.9x and 3.8x); sphere 1.2x slower on one thread, parity on 24.
- The quality modes stay slower than Flying Edges at one thread (gyroid balanced
  0.226 s against 0.175 s) with 2.8x fewer triangles and 18° minimum angle.
- Volume block, same job (gyroid, tetrahedra): lattice-32 1.73 s / 0.58 s,
  adapt 1.81 s / 1.78 s (the tet adaptation barely uses the threads), MMG3D on the
  non-periodic seed 4.77 s / 4.65 s.
- Balancing offsets vertices by half the centroid sag, which minimizes H; the volume
  bias left is about a third (for a spherical cap the mean sag is 3/4 of the centroid
  sag). An offset of 3/4 would cancel the volume error at a larger H.
- Limits: the hard-min lattice (creases) fails the 1e-3 check in every adaptive mode,
  as before; Flying Edges 48 reaches only 1.26e-3 there, without a check. Fast mode
  keeps triangles valid (positive shape) but has slivers.

## Volumes: size-driven tet adaptation (`tet_adapt.rs`)

MMG-style split, collapse, edge removal, 2-3 swaps and relocation, with boundary
vertices projected onto the exact field. Against MMG3D (mmgpy 0.17.1) adapting the
same 24-cell seed with the same targets (hausd = tolerance, hmax 0.1, hgrad 1.3,
wall and box patches tagged):

| Case, tolerance | Adapted: tets / H / q_min / q_p01 | MMG3D: tets / H / q_min / q_p01 |
|---|---|---|
| Sphere 1e-3 | 13,689 / 7.7e-4 / 0.51 / 0.60 | 17,154 / 8.7e-4 / 0.42 / 0.52 |
| Ellipsoid 1e-3 | 16,059 / 9.6e-4 / 0.52 / 0.62 | 15,980 / 1.05e-3 / 0.36 / 0.51 |
| Gyroid 2e-3 | 42,441 / 1.5e-3 / 0.27 / 0.59 | 52,427 / 1.8e-3 / 0.14 / 0.47 |
| Gyroid 1e-3 | 74,662 / 1.0e-3 / 0.18 / 0.58 | 83,801 / 3.8e-2 / 0.15 / 0.46 |

MMG only sees the seed's facets, so its geometry cannot beat the seed (1.6e-3 for the
24-cell Gyroid); at 1e-3 it produced a 0.038 defect. Periodic cells are adapted on the
torus (copies merged, minimal-image geometry, box planes as internal constraints);
opposite faces match exactly. Periodic Gyroid at 1e-3: 74,938 tets, q_min 0.309,
against 171,679 tets, q_min 0.274 for the 32-cell lattice at the same H.

## FEA: homogenized stiffness error per DOF

Periodic Gyroid sheet, E=2.5, ν=0.25, the audited fedoo homogenization of the
periodic-fea experiment (six loads, single-thread CG). Error is the worst relative
strain-energy difference over all macroscopic strains against the quadratic
extrapolation of the archived native64/native80 linear tensors (observed order 2.1).

| Mesh | Order | Nodes | Elements | Error |
|---|---|---:|---:|---:|
| Lattice 24 | linear | 17,647 | 72,757 | 3.71% |
| Adapted 1e-3 | linear | 18,669 | 74,938 | 3.50% |
| Lattice 32 | linear | 39,717 | 171,679 | 2.12% |
| Lattice 8 | quadratic, straight | 6,914 | 3,257 | 4.24% |
| Lattice 12 | quadratic, straight | 19,931 | 10,832 | 1.67% |
| Lattice 8 | quadratic, curved | 6,914 | 3,257 | 1.31% |
| Lattice 12 | quadratic, curved | 19,931 | 10,832 | 0.29% |
| Adapted 4e-3 | quadratic, curved | 32,602 | 18,600 | 0.16% |
| Lattice 16 | quadratic, curved | 45,527 | 25,598 | 0.10% |
| Adapted 2e-3 | quadratic, curved | 70,024 | 42,768 | 0.04% |
| Lattice 20 | quadratic, curved | 79,848 | 46,881 | 0.04% |

- For linear tets in this thin sheet, adaptivity does not change error per DOF
  (3.50% vs 3.56% interpolated from the lattice): through-thickness resolution rules.
- Curved quadratic tets give 13-22x lower error at equal node count. A 1% target needs
  about 8k quadratic nodes and ~4k elements, against roughly 100k linear nodes and
  ~500k elements. Observed quadratic order is 3.7, close to the expected h^4.
- Curving the wall mid-edge nodes onto the exact field is most of the gain (3-6x over
  straight quadratic elements). Elevating a VTK or MMG mesh without the field gives
  only the straight result.
- The reference is cross-checked: extrapolating the quadratic 16- and 20-cell
  tensors with p=4 agrees with the linear 64/80 extrapolation within 0.007-0.012%,
  so quadratic errors are resolved down to ~0.04%. Quadratic lattice 20 (79,848
  nodes, 46,881 elements) is at 0.037%.

Balancing linear wall vertices (half the centroid sag outward, periodic copies
sharing one offset) cuts the volume bias about 2.5x but lowers the stiffness error
only ~5%: lattice 24 3.71% -> 3.51%, lattice 32 2.12% -> 2.01%, adapted 1e-3
3.50% -> 3.33%. Linear discretization, not geometry bias, dominates.

### Thick solid: Schwarz P network (`fea_solid.py`)

cos 2πx + cos 2πy + cos 2πz < 0, about 50% solid with thick ligaments. Reference:
p=4 extrapolation of the quadratic 16- and 20-cell tensors (they differ from it by
0.03% and 0.01%).

| Mesh | Order | Nodes | Elements | Error |
|---|---|---:|---:|---:|
| Lattice 16 | linear | 6,240 | 27,492 | 3.06% |
| Lattice 24 | linear | 19,004 | 90,792 | 1.35% |
| Lattice 32 | linear | 41,956 | 210,576 | 0.76% |
| Adapted 4e-3, hmax 0.2 | linear | 2,364 | 8,889 | 5.20% |
| Adapted 2e-3, hmax 0.2 | linear | 3,888 | 15,825 | 3.58% |
| Adapted 1e-3, hmax 0.2 | linear | 8,451 | 34,918 | 2.58% |
| Lattice 8 | quadratic | 6,864 | 3,852 | 0.26% |
| Lattice 12 | quadratic | 19,628 | 11,976 | 0.074% |
| Adapted 4e-3, hmax 0.2 | quadratic | 15,221 | 8,889 | 0.065% |
| Adapted 2e-3, hmax 0.2 | quadratic | 25,890 | 15,825 | 0.031% |

Geometric adaptivity lowers the linear error by about 20% at equal nodes on coarse
meshes and loses slightly at 8.5k nodes: FEA error follows stress resolution in the
ligament necks, which curvature sizing does not target. An error-driven size field
would be the next step for linear meshes. Curved quadratic elements again win by an
order of magnitude: 0.26% with 6.9k nodes and 3.9k elements against 0.76% with 42k
linear nodes and 211k elements.

Exact (quadrature) solid volume error of the same meshes:

| Grid | Linear | Curved quadratic |
|---:|---:|---:|
| 8 | +2.25% | +0.025% |
| 12 | +0.96% | -0.0007% |
| 16 | +0.47% | +0.0000% |

## Limits

- One FEA fixture (Gyroid sheet, one band). Thicker solids, lattices and spinodoids
  are untested; for thick solids interior grading should also reduce linear DOFs.
- Quadratic validity is checked by sampling the Jacobian at 25 points per element,
  not by Bernstein bounds; failing elements fall back to straight edges.
- No timing claims. Single-run wall clocks on a loaded machine suggest the adaptive
  surface is near Flying Edges' cost and the tet adapter near the lattice pipeline's;
  both are single-threaded.

## Issues #3 and #4 (8 cells, tolerance 0.01, minimum quality 0.1)

`issues.py`, full log in `issues.txt`. Independent two-sided distance H and volume
error of accepted meshes:

| Fixture | Lattice, 0 / 4 / 12 passes | Adapter (8-cell seed) |
|---|---|---|
| Gyroid | accepted all; H 1.0e-2, q_min 0.15 | H 9.9e-3, ΔV +1.5%, q_min 0.33 |
| Schwarz P | rejected at 0 passes | H 8.1e-3, ΔV +1.4%, q_min 0.39 |
| Neovius | rejected (2.9e-2) | H 9.8e-3, ΔV +1.8%, q_min 0.28 |
| Graded Gyroid | rejected at 0 passes (quality) | H 9.8e-3, ΔV +1.4%, q_min 0.34 |
| Thin Gyroid | rejected at 0 passes (quality) | H 9.0e-3, ΔV +1.8%, q_min 0.32 |
| Implicit BCC lattice | rejected (4.2e-2) | rejected (3.5e-2) |

Issue #3's non-monotonic Gyroid is fixed on this branch by PR #5's geometry guard.
Two adapter changes made the hard cases pass: projection falls back to bracketing
and bisection, and a wall vertex whose exact position would invert a neighbour is
inserted part of the way and moved onto the surface once refinement makes room.
That move keeps the split quality floor, so geometry is never bought with slivers.

The capsule lattice is a hard `min` of eight distances, so its surface has concave
creases. Triangles that straddle a crease keep a sag proportional to their size,
and the vertices that would resolve it cannot settle without slivers. The lattice
pipeline first passes at 24 cells. Capturing creases needs feature-aware
projection onto the curve where two `min` branches are equal; the compiled-field
JIT knows those branches. Spinodoid is untested here (needs Microgen).

### Curve chords and the spinodoid

The spinodoid (`spinodoid.py`, run in the PR #2 environment for Microgen's
`Spinodoid`) exposed a measurement gap. On a wall face, a chord that lies in a box
plane follows the wall/plane curve. The first-order check used the full field
gradient, which measures the distance to the unclipped wall; where the wall meets
the box at a shallow angle that wall passes close outside the box, while the curve
inside the plane is far. The 16-cell adapted spinodoid passed its own check (7.3e-3)
but its independent surface-to-mesh distance was 1.23e-2, at a curve on y = 0.5.

The adaptive paths now measure such chords with the in-plane gradient, size curve
edges by the curve's in-plane turning, and the surface remesher samples every edge of
every wall face (curve edges with a cap neighbour were sometimes never sampled).
The adapted spinodoid is now rejected at 16 cells (2.0e-2) instead of passing with a
hidden error. The lattice pipeline keeps its released check: with the in-plane
measure, three of its existing tests exceed their tolerance by 8% to 6x. Whether to
tighten the released check is a tolerance decision left open.

### Lattice with smooth unions

`smooth_lattice.py` rebuilds the issue #4 lattice with Microgen's cubic polynomial
smooth union (`implicit_ops.smooth_union`). Microgen's strut lattices currently go
CAD -> STEP -> gmsh; its F-rep path already offers smooth unions. A sharp concave
crease is also a stress singularity, so FEA would not converge there anyway.

| Blend k | Cells | Lattice pipeline | Adapter |
|---|---:|---|---|
| 0.08 | 12 | 5,696 tets, H 5.8e-3, ΔV -3.3%, q_min 0.32 | 6,225 tets, H 1.0e-2, ΔV -3.1%, q_min 0.48 |
| 0.08 | 16 | 24,480 tets, H 7.4e-3, ΔV -0.5%, q_min 0.20 | 8,129 tets, H 6.9e-3, ΔV -2.0%, q_min 0.34 |
| 0.03 | 8-16 | rejected | rejected (closest 1.09e-2 at 16) |

A blend much smaller than the cell size behaves like a crease. When the exact
position of a wall vertex would flatten a neighbour, the adapter now splits the
longest edge facing it in each such element and retries, which made the smooth
fixtures slightly better (Neovius q_min 0.28 -> 0.34). Equalizing sag to the
tolerance raises the mean deviation, hence the larger volume error of the adapted
16-cell mesh at equal H; balanced placement would cancel it for linear meshes.

## Real Microgen defaults on Linux G3, 2026-10-17

`MICROGEN_RESULTS.md` records the full default VTK-to-MMG workflow on the
Gyroid sheet, Schwarz P solid, fixed spinodoid and smooth-union lattice.
This replaces no earlier frozen results and does not relabel the historical
MMG-on-a-Meshers-seed measurements. Both Microgen paths now have real G3
execution evidence on all four fixtures, with unsupported FEA routes marked.

At about 20k nodes, curved Meshers has 0.2951% Gyroid and 0.0738% Schwarz-solid
full mixed-load stiffness disagreement, versus 3.7152% and 1.3458% for the
tested linear meshes. Their geometry passes the unchanged common gates,
but curved CG workflows take longer at these unmatched accuracies.
Every default Microgen control fails geometry or quality admission.

The negative findings matter. Meshers rejects periodic spinodoids at the
quality floor and rejects the smooth lattice's raw-field periodicity,
although Microgen completes diagnostic lattice FEA. Neither fixture has a
justified fine stiffness reference. Same-domain order controls separate
approximation-order changes from curved geometric-mapping changes without
claiming certified continuum error. Native meshing code remains unchanged.

## Public direct-fedoo Tet10 integration, 2026-10-17

`TET10_RESULTS.md` records the successful G3 integration through `Mesh.to_fedoo()`.
Four six-load cases pass: solid Tet4, solid Tet10, curved Gyroid and curved Schwarz
solid. Independent strain-controlled stress integration reproduces the saved
specimen tensors within 1.2e-11 relative mixed-load energy disagreement. Periodic
mid-edge displacement jumps remain below 6e-17. This validates conversion and
periodic FE integration, not local stress accuracy. Microgen's public converters
remain unsupported for Tet10.
