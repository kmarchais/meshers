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
