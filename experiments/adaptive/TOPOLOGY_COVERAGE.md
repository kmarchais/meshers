# Failure coverage beyond the four study specimens

Follow-up: [automatic polyhedral preservation](POLYHEDRAL_RESULTS.md) fixes the
default nonadaptive scalar route. The adaptive fused-octahedron failure below
remains recorded and uncorrected; it is not reclassified by that narrower change.

Kevin prioritizes failure coverage before optimization. The four original
specimens remain the comparison study. These additional analytic controls
exercise a handle, an internal cavity, disconnected solids and sharp corners.
They do not introduce new Microgen/MMG or stiffness comparisons.

The native Mac suite passes 143 tests with one optional wide-constraints skip.
After adding the sharp-corner regression, 11 focused native/fixture/recipe
controls pass. The new production test passes Ruff checks and formatting with
the repository's locked Ruff 0.16.7. No native library changed or rebuilt.

## Curved surfaces need separate volume admission

The torus has major radius 0.25 and minor radius 0.08. With the initial
0.002 adaptation and geometry tolerance, it passes distance 0.001965775,
minimum quality 0.443241 and the expected one-component, genus-one boundary.
Its volume deficit is 1.405454%, so it fails the unchanged 1% volume gate.
Independent signed boundary and cell volumes agree within 3.47e-18, and
wall-vertex field residuals are below 2.60e-14. The deficit comes from the
chordal approximation rather than missed topology or inconsistent volume sums.

One candidate tolerance is selected from `0.01 * exact_volume / exact_area`.
For this torus it is 0.0004. This scale is a refinement choice, not a certified
error bound. The resulting linear mesh has 76,214 tets, distance 0.000399645,
minimum quality 0.398901 and volume deficit 0.230699%. It passes every gate.
The material definition is unchanged. No surface offset or sizing sweep is used.

The spherical shell has inner radius 0.18 and outer radius 0.32. The initial
0.002 tolerance passes with 10,275 tets, distance 0.001443754, minimum quality
0.436980 and volume deficit 0.334739%. Both expected boundary components and
their combined Euler characteristic 4 are preserved.

Two separated spheres have radius 0.15 and centers at x = -0.23 and +0.23.
At tolerance 0.002, distance and topology pass but volume deficit is 1.813767%.
The same volume/area rule selects one tolerance 0.0005. It passes with 22,597
tets, distance 0.000494024, minimum quality 0.420379 and volume deficit 0.464510%.
Both sphere components remain separate, with combined Euler characteristic 4.

## Named constraints preserve sharp corners

The octahedron is `abs(x) + abs(y) + abs(z) <= 0.35`. Its fused scalar-field
route returns 61,302 tets with minimum quality 0.288827 and volume deficit
0.027431%. The native forward estimate is 0.001871043. Independent forward
distance is 0.001897237, but reverse distance reaches 0.002418269. The mesh
therefore fails the unchanged two-sided distance gate. A forward estimate
does not qualify the exact corners.

The same domain can be expressed as eight named linear half-spaces. The
existing `generate_intersection` route preserves their feature edges and the
four-plane corners. At the same seed24, quality0.1 and distance0.002 settings,
it admits 13,488 tets with minimum quality 0.317878. Sampled distance is
8.41e-17 and relative volume error is 2.22e-16. The production regression
checks all six corner coordinates, four active constraints at each corner,
assigned face planes, positive cell volumes and exact total volume.
See [the public example](../../docs/site/intersections.md).

This is an explicit supported representation, not an automatic scalar-field
feature fix. The fused-field reverse-distance failure remains recorded. Do
not replace true crease constraints with averaged tie derivatives, round the
corners, or infer arbitrary nonsmooth-field support from this named-plane case.

## Evidence and replay

Jobs 123 and 124 preserve the initial torus attempts. Job 123 failed in the
new auditor because NumPy returns shaped inverse indices; flattening the
indices corrected the auditor. That attempt is not a native meshing failure.
Job 124 records the actual torus volume rejection. Jobs 125–130 preserve
the admitted and rejected geometry controls described above.

The first eight-case replay, job 131, caught a missing historical test-source
snapshot after the recipe tests changed. Nine missing copies were restored
from commit `34d7540`, each matching the exact hash already recorded by its
run. Measured results and their manifests are unchanged. The restoration
receipt is
[coverage-source-restoration-20261024.json](../../handoff/macbook/coverage-source-restoration-20261024.json).
New benchmark and analytic runs now archive every Python input. Jobs 133–134
exercise both paths and their saved copies match all recorded source hashes.

Job 132 passes a fresh saved-array replay of all eight cases: the four analytic
controls and the original Gyroid, Schwarz P, spinodoid and smooth lattice.
It checks native/input hashes, fixed specimen/sample hashes, counts, quality,
positive volumes, oriented facets, independent geometry and periodicity where
requested. The analytic replay uses a second triangle sampling grid, exact
volume formulas and a boundary-volume identity. It is sampled evidence, not
a global surface-distance or continuum certificate. Results are in
[eight-case-coverage-20261024.json](../../handoff/macbook/eight-case-coverage-20261024.json).

Re-run `verify_topology_coverage.py --output <new-file>` in the recorded native
environment, through the shared bounded queue. The original meshes and fixed
sample archives must be present. Individual new probes use
`topology_coverage.py --run <existing-stage-directory> --case <name>`; the
admitted torus and disjoint recipes add `--tolerance .0004` and
`--tolerance .0005` respectively, and the octahedron adds `--named`.

No new comparative speedup, FE reference or full-goal completion is claimed.
Actual Gyroid/spinodoid admission and lattice reference convergence remain
unresolved. The next native feature investigation has a concrete red control:
the fused octahedron's reverse-distance miss, with an exact named-constraint
control on the identical domain. Preserve all four study cases and the
existing geometry, quality, periodicity and resource gates during that work.
