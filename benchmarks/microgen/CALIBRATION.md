# Calibration toward Microgen adoption

These measurements guide the next experiment. They are not a claim that Meshers
is ready to replace Microgen's meshing pipeline. The common benchmark measures
worst elastic energy error over a fixed quadratic displacement space on each
produced polyhedral domain. Engineering and homogenization accuracy still need
a converged reference on the same physical domain.

## Admission policy

The standard policy requires sampled distances in both directions below
`0.002 L`, relative volume error plus empirical reference uncertainty below 1%,
and minimum MMG tetrahedron quality 0.1 for FEA. Every measured repeat must pass.
FEA targets are 1%, 2% and 5%. Surface and geometry runs retain the geometry
policy and do not impose the FEA quality gate.

The Sobol volume uncertainty is an empirical consistency estimate, not a
confidence interval. Boundary samples and floating-point intersection checks
are not global certificates. Passing the volume test alone is insufficient.

## Initial Gyroid calibration

The runs below used one thread, one repeat and no warmup. Their generation times
are calibration samples, not stable performance claims. The local output
folders retain exact source snapshots, dependency hashes, meshes and all 24 FEA
solutions. Compact evidence files in `calibration/` retain the metrics and
provenance without mesh previews.

| Workflow | Grid | Tetrahedra | Generation | Worst energy error | Admission |
| --- | ---: | ---: | ---: | ---: | --- |
| Meshers balanced, wheel 0.1.0 | 28 | 65,386 | 2.853 s | 5.108% | Misses 5% FEA |
| Meshers balanced, wheel 0.1.0 | 36 | 130,438 | 6.010 s | 4.065% | Passes 5% |
| MMG with cleaned tetrahedral VTK seed | 36 | 166,937 | 15.764 s | 3.449% | Passes 5% |
| Meshers balanced, wheel 0.1.0 | 54 | 389,200 | 13.623 s | 2.821% | Passes 5%, misses 2% |
| MMG with cleaned tetrahedral VTK seed | 54 | 367,193 | 39.832 s | 2.583% | Passes 5%, misses 2% |

The original Microgen VTK path produced overlapping tetrahedra in the saved
Gyroid grid-16 case. A minimized regression exercises that actual defect.
Triangulating the background before clipping removes that overlap, but leaves
zero-volume cells. The explicit cleaned control removes those cells and merges
coordinates at its declared tolerance. Its grid-28 MMG output still misses the
geometry policy. No failed original pipeline is silently replaced.

AMG and direct LU were compared on the exact saved 130,438-tetrahedron mesh.
The worst energy error agrees to `3e-17` absolute; all displacement solutions
agree below `1e-6` relative. The AMG solve took 37.581 s versus 258.773 s for LU.
Patch and equilibrium checks remain mandatory. This makes finer-grid
calibration affordable without changing the FEA problem.

## Fast surfaces

The public Microgen surface workflow and a separate Flying Edges control are
measured independently from volume meshing. Meshers 0.1.0 still builds the full
volume and extracts its boundary; that full cost enters the surface comparison.

A direct Flying Edges prototype passed the strict geometry policy for Gyroid at
grid 36 and for the implicit lattice at grid 54. The prototype's first cap
triangulation left open edges. Its explicit cleanup now conforms existing
collinear cap edges without moving vertices or filling holes. Gyroid grid 54
still has open edges and remains a failure. Spinodoid grid 54 has accurate
volume but a sampled boundary distance of `0.002437 L`, so it remains excluded.
Common-runner repeats are needed before quoting a speed ratio.

This is the strongest reason to add a native surface-only API. Comparing full
native volume generation against a direct VTK surface algorithm would hide the
cost Microgen users actually pay for printing and visualization.

## Coarse implicit failures

Tracing issue #3 reproduced its reported Gyroid failures exactly. Quality
optimization moved already-admitted surface edges outside the geometry limit.
Refinement then left the offending midpoint unchanged. The correction rejects
quality proposals that worsen an incident implicit face beyond the larger of
its existing sampled error and the requested tolerance. It preserves the final
quality and geometry gates.

Uniform, graded and thin Gyroid now succeed at both four and twelve passes in
the native reproduction at grid 8, geometry tolerance `0.01 L`, minimum quality
0.1 and the original 250,000-tetrahedron budget. Schwarz P and Neovius still
reject. Issue #3 must stay open, and the lattice and Spinodoid in issue #4 need
further work. The common geometry validator must admit a returned mesh before
it enters benchmark comparisons.

## Next measurements

Measure short optimization budgets at grid 32, then choose the next FEA grid
from those results. Measure the 2% target rather than extrapolating it from
coarse meshes. Repeat selected qualifying configurations with warmups before
making timing claims. Keep surface timing separate and include all failures in
the report.
