# Automatic polyhedral feature preservation, 24 October 2026

Meshers now recovers exact half-space constraints from supported compiled convex
piecewise-affine fields and uses its existing constrained cutter automatically.
This fixes the default nonadaptive Tet4 route for the fused octahedron without
requiring the caller to rewrite its eight planes.

The original four study specimens, definitions and admission gates are unchanged.
This is a failure-coverage improvement, not a new Microgen/MMG timing or
mixed-load stiffness result. The full comparison objective remains open.

## Diagnosis and change

The octahedron is `(abs(x) + abs(y) + abs(z) - 0.35) / sqrt(3) <= 0`.
Removing adaptation alone still rejects it. With both adaptation and optimization
disabled, the scalar cut returns 14,448 tets, minimum quality 0.141068 and
two-sided distance 0.002932698. Its largest witnesses lie on crease edges,
rather than the six tips. The initial piecewise-linear scalar cut loses face
identities. Optimization is not the sole cause.

The compiler can represent this value exactly as the maximum of eight affine
forms. The new bounded graph analysis retains those forms through affine
operations, positive constant scaling, sums, affine absolute values and maximum.
Nonlinear/nonconvex operations, zero-gradient branches and expansions beyond
the build's constraint capacity fall back. It does not infer arbitrary creases.
No derivative averaging, geometric smoothing or surface offset is introduced.

The specialization applies to compiled negative-side, nonperiodic, nonadaptive
Tet4 meshes with one worker, no interior refinement or separate gradient
override, and at most 20 optimization passes. Other options retain their scalar
path. Explicit gradient overrides keep their requested policy. The diagnostic
`polyhedral_constraints` identifies the specialization. See the
[public example](../../docs/site/intersections.md).

## Before and after on the same input

At seed24, tolerance 0.002, quality floor 0.1 and cap 250,000, the old default
nonadaptive route rejects with native sampled error 0.003240370. The preserved
old package reproduces this exact failure in job 146. This is an intentionally
failed regression, not an interrupted job and not a job to restart.

The candidate takes the same scalar function and admits 13,488 tets and 3,451
nodes. Independent sampled two-sided distance is 8.28e-17, minimum quality is
0.317878 and relative volume error is zero at reported precision. Oriented
facets, positive linear Jacobians and the expected boundary topology pass.
The external distance, volume, quality and resource limits are unchanged.

A rotated, translated octahedron also passes the production regression. Its
raw field is checked at surface vertices and facet centroids, along with exact
volume and assigned plane equations. Compiler controls compare reconstructed
values at fixed random points and reject nonconvex/nonlinear graphs. Fallback
modes, explicit gradients, cancellation, hard budgets and strict quality
rejection have end-to-end controls.

The final compatibility check caught integer-option regressions in the initial
dispatch. The corrected route retains the original float-type and zero-budget
exceptions and normalizes valid NumPy integers for the constrained backend.
An old/new package replay records those outcomes in
[polyhedral-option-compatibility-20261024.json](../../handoff/macbook/polyhedral-option-compatibility-20261024.json).

## Four-study regression and provenance

Fresh jobs 141–144 regenerate the admitted linear Gyroid, Schwarz P, spinodoid
and lattice controls. All four pass every existing admission gate. Against
their saved pre-change controls, Gyroid, Schwarz and spinodoid arrays are
bitwise identical. Lattice changes node/cell ordering only: coordinates,
tetrahedra and boundary triangles are exactly identical after canonical
relabeling. The comparison does not introduce a numerical tolerance.

Job 145 first required raw-array equality and failed on that lattice permutation.
Job 147 proves exact canonical equality and replays
the preserved historical eight-case data with original library hashes. That
audit takes 13.76 seconds at 1.263 GiB. Results are in
[polyhedral-change-audit-20261024.json](../../handoff/macbook/polyhedral-change-audit-20261024.json).

The baseline native SHA256 begins `8565629a`; the candidate begins `6ce9e397`.
The old library and Python package were copied locally before rebuilding.
Historical verifiers accept an explicit preserved library only when its hash
matches the original measured input. Frozen manifests and results are not
rewritten. New generation manifests also hash and archive the Meshers Python
facade, expression compiler and constrained wrapper, since dispatch is part
of the measured implementation.

Rust 1.95.0 and `experimental-surfaces` produce a native arm64 build. Clippy
with warnings denied, Rust formatting and Python lint pass. The final full
suite passes 171 tests with one optional wide-constraints skip. Two strengthened
raw-field geometry checks pass separately after their assertions were added.
No native Rust code changed after the checked build. Candidate source hashes
and native architecture are recorded in
[polyhedral-native-candidate-20261024.json](../../handoff/macbook/polyhedral-native-candidate-20261024.json).

## Remaining work

The adaptive fused-field octahedron remains an unsupported feature-preservation
case. The new dispatch deliberately keeps adaptation, periodic requests and
quadratic output on their existing paths. Extending adaptive feature handling
must carry actual constraint identities through cutting, collapse, relocation
and projection; another tolerance sweep or averaged subgradient is not a fix.

Keep all four original study cases during that work. Tested actual Gyroid and
spinodoid routes still lack admission, lattice reference convergence remains
unresolved, and no full mixed-load reference or comparative speedup is promoted.
