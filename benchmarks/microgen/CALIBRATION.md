# Calibration toward Microgen adoption

These measurements guide the next experiment. They are not a claim that Meshers
is ready to replace Microgen's meshing pipeline. The common benchmark measures
worst elastic energy error over a fixed quadratic displacement space on each
produced polyhedral domain. Engineering and homogenization accuracy still need
a converged reference on the same physical domain.

The adoption baselines are Microgen's current VTK path and its public
`remesh_keeping_boundaries_for_fem` CLI pipeline, named `mmg-microgen`.
Direct `mmgpy` binding measurements are additional algorithm controls. They do
not measure Microgen's current remesher or imply support for that Python API.

## Current Microgen MMG pipeline

At grid 32, the public two-pass remesher returns 75,257 tetrahedra in 8.276 s.
The worst energy error is 4.213%, but the mesh remains excluded. Its sampled
forward boundary error is `0.006797 L`, reverse distance is `0.005353 L`,
and 266 tetrahedra fall below quality 0.1, with minimum 0.01023. Relative volume
error plus empirical uncertainty is 0.755%. Passing the FEA error test and
volume test does not compensate for these boundary and shape failures.

This is one calibration measurement with no warmup, retained in
`calibration/gyroid-microgen-mmg-32.json`. It used the installed command alias
that launches the bundled native MMG executable. Subsequent runs call that
native executable directly, record its exact path and hash, and can use an
existing installation through `--mmg-bin-directory`. Microgen's two-pass
function, parameters, required-boundary treatment and file exchange remain
unchanged.

The grid-54 run returns 332,322 tetrahedra in 30.401 s with 2.511% energy error.
It also remains excluded. Forward boundary distance is `0.004463 L`, reverse
distance `0.003591 L`, and 745 tetrahedra fall below quality 0.1, with minimum
0.004806. Volume error plus empirical uncertainty is 0.381%.
`calibration/gyroid-microgen-mmg-54.json` retains the full metrics and native CLI
hash. Both resolutions used `hausd=0.002 L`, `hsiz=L/cells`, `hgrad=1.3` and
Microgen's required-boundary policy. A smaller MMG approximation tolerance is a
useful next control, since its tolerance concerns the reconstructed boundary
rather than the analytic field. These failures provide no qualified speed ratio.

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

The early binding runs below requested cap protection but supplied no boundary
triangle entities to MMG. Their actual mode was free surface adaptation. The
independent geometry checks still apply to their returned meshes, but their
times must not be described as measurements of the corrected cap-protected
pipeline. The boundary-entity correction and fresh controls are described below.

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

At grid 32, the short-budget trial passes all strict admission gates with each
profile. All three have 96,060 tetrahedra. One pass takes 1.187 s with 4.574%
energy error, two take 2.031 s with 4.562%, and four take 4.537 s with 4.551%.
One pass is the measured cheapest candidate for the 5% target in this trial.
The small FEA change does not justify spending four passes for this case and
target. This does not establish that one pass suffices for other fields.

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
The completed common-runner trial uses one warmup and three measured repeats
per configuration, with all 72 measured rows and 24 warmups retained in
`calibration/surface-repeated.json`. Every repeat of the configurations below
passes the same strict geometry policy.

| Case | Workflow | Grid | Triangles | Median generation |
| --- | --- | ---: | ---: | ---: |
| Gyroid | Microgen public VTK surface | 54 | 72,512 | 0.1944 s |
| Gyroid | Direct Flying Edges control | 36 | 25,872 | 0.07543 s |
| Gyroid | Meshers zero passes, full volume then boundary | 24 | 20,944 | 0.4010 s |
| Gyroid | Meshers one pass, full volume then boundary | 24 | 20,944 | 0.8760 s |
| BCC lattice | Direct Flying Edges control | 54 | 34,904 | 0.1825 s |

At this geometry target, the fastest tested Meshers surface configuration costs
2.06 times the current Microgen VTK path and 5.32 times the direct Flying Edges
control for Gyroid. These are observed target-cost comparisons on this machine;
their triangle counts differ. One optimization pass costs more without improving
the zero-pass Gyroid's boundary distances. This supports the zero-pass profile
for this printing/visualization case.

The public VTK Gyroid surface fails topology at grids 24 and 36 but passes at
54. The direct control passes at 36 and fails topology at 54. Keep these
nonmonotonic results. The lattice public surface leaves the specimen box at
every tested grid. The original Meshers lattice field uses NumPy callback
fallback; neither zero nor one pass returns a mesh meeting `0.002 L` at grids
24, 36 or 54. Removing the FEA shape gate does not resolve that geometry failure.

This is the strongest reason to add a native surface-only API. Comparing full
native volume generation against a direct VTK surface algorithm would hide the
cost Microgen users actually pay for printing and visualization.

## Traceable BCC field control

The vector implementation of the BCC field uses `broadcast_arrays`, array
construction and reductions that the released native scalar tracer does not
support. A scalar-coordinate expression computes the same capsule distances
and also broadcasts over NumPy arrays. On 10,027 random and named points at
three scales, its NumPy values agree exactly with the original vector expression;
compiled values differ by at most `2.22e-16`, and compiled gradients are finite.
Regression tests retain the original expression and check gradients against
finite differences away from the capsule junctions.

The benchmark now uses that shared traceable expression for all workflows.
Meshers can use its public compiler and automatic derivatives instead of NumPy
callback fallback. Its preparation cost remains included. A fresh pilot must
measure the geometry and timing effect before any claim of improvement. The
earlier surface trial remains evidence for its original callback implementation.
The compiler-only pilot at grids 24 and 36 completed four attempts in 23 s,
but all four still reject `0.002 L`; its evidence is retained in
`bcc-traceable-band.json`. Faster rejection is not a usable-mesh speed result.

An optional `--meshers-solid-region negative` control removes the inactive
lower band wall for BCC and Spinodoid while keeping their signed fields,
material sets and independent admission policy fixed. It changes only the
native region choice; sheet bands and the legacy default remain unchanged.
The region is recorded even when native generation rejects the input.
At grid 24, zero and one pass still reject both fields under the strict policy.
Their evidence is retained in `implicit-negative-24.json`. The single-negative
control does not resolve these failures, and the legacy default stays unchanged.

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

The common-runner verification used native commit
`66b800cf80e791788b980c01a8438f0e7f8ec8c6`, proposed in
[PR #5](https://github.com/kmarchais/meshers/pull/5), without replacing the locked
wheel baseline. Its binary hash is retained in `native-geometry-correction.json`.
Both uniform Gyroid outputs also pass the independent smoke geometry policy.
The graded and thin outputs remain geometry warnings because reverse distances
range from `0.011138 L` to `0.011397 L`. Fixing the native acceptance regression
does not make those coarse meshes eligible for comparison. Native formatting,
lint and all 49 release tests passed; the real Python suite passed 102 tests
with one skip after correcting the scalar evaluator path.

The useful-resolution check of the same native commit passes all geometry,
quality and 5% FEA gates at Gyroid grid 32 with one and four passes. Both have
96,060 tetrahedra; minimum qualities are 0.1788 and 0.1824, and worst energy
errors are 4.574% and 4.551%. The single generation samples are 1.341 s and
4.680 s. Both use the compiled evaluator with zero Python callbacks. These
remain calibration times; the evidence is `native-guard-fine-32.json`.

## Measured 2% target and validator correction

One-pass Meshers at grid 80 passes the strict policy with 1,190,946 tetrahedra,
minimum quality 0.1828 and 1.942% worst energy error. Generation took 11.211 s.
The forward distance is `0.0001349 L`, reverse distance `0.00008669 L`, and
volume error 0.0382% before empirical reference uncertainty. This remains a
single calibration timing.

The first MMG grid-80 output was rejected by an absolute `1e-15 L³` determinant
cutoff. Inspection found 96 small, positive-volume cells with minimum MMG shape
quality 0.1350. Those cells are valid. The validator now tests collapse relative
to each cell's longest edge cubed, with a `1e-12` threshold, while retaining all
topology, overlap, geometry and FEA quality gates. A saved MMG cell is tested at
three sizes, with a collapsed version required to fail at each size. The
original false rejection is retained in `gyroid-two-percent-initial.json` and
must not be interpreted as an MMG defect. The corrected common-runner MMG
measurement is retained in `gyroid-two-percent-mmg.json`. It reaches 1.770%
energy error in 1,139,944 tetrahedra, but 739 cells fall below quality 0.1 and
some wall normals disagree with the analytic field. It remains a geometry
warning and cannot qualify for the FEA target. Its 109.088 s generation time
is retained without a speed claim.

The optional MMG control distinguishes wall patches and physical box planes with
references, allowing cap interiors to be retriangulated. It retains the same
cleaned seed and all independent admission gates. A corrected trial needs a
coarse proof before another expensive fine-grid run. The current Microgen CLI
baseline and repeated surface measurements take priority.

The boundary-entity regression exposed an adapter error in the early binding
controls. `mmgpy`'s volume-only conversion has no triangles until MMG analyzes
the mesh. Calling the required-triangle API or assigning references to that
empty list changes nothing. The corrected adapter extracts actual boundary
triangles, maps their point IDs back to the volume, and supplies a mixed
tetrahedron/triangle grid so MMG allocates both entity types. It verifies the
supplied count and records actual cap counts. The cube integration tests now
require nonempty cap constraints and all six distinct references in the tagged
control. An initial grid-32 tagged trial was stopped when this check failed;
its local output remains marked interrupted, and it provides no tagged result.

## Next measurements

Measure Microgen's current CLI pipeline at finer grids, then repeat selected
qualifying configurations with warmups before making timing claims. Keep
surface timing separate and include all failures in the report. Extend the
strict FEA calibration beyond uniform Gyroid, especially the coarse failures
in issues #3 and #4. The measured 2% native run used a two-million-tetrahedron
budget, a 60-million validation-pair cap and a 16 GiB worker RSS ceiling;
these larger limits remain explicit rather than becoming silent defaults.
