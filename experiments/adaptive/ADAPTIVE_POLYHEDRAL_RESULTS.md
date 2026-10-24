# Adaptive convex polyhedral coverage, 24 October 2026

Supported compiled convex polyhedra now retain exact affine planes during
nonperiodic adaptive Tet4 generation. The seed uses the existing constrained
cutter. Adaptive splitting, collapse and projection use active plane identities;
face and interior coarsening and spatial size callbacks remain active. Existing
crease vertices remain fixed except for checked collapse of numerical copies
within 64 machine epsilons times the box span, with identical active planes.
The final quality, distance and resource gates are unchanged.

## Failure and correction

The old package rejects the adaptive scalar octahedron. The first candidate
admits the axis-aligned specimen, but the full suite exposes a panic on the
rotated and translated specimen. This failure is retained in job 155.
A reconstructed Python seed has no aliases in job 156; it does not reproduce
all arithmetic of the actual Rust seed and cannot clear that seed.

Instrumentation in jobs 157–160 finds distinct native seed vertices about
6e-17 apart. The adapter's coordinate buckets merge them even without periodic
axes, creating a tetrahedron with a repeated vertex ID before any topology
operation. Nonperiodic initialization now preserves seed identities. Job 162
then rejects one sliver at q 1.6025e-15 instead of panicking. Numerical copies
of the same corner may now collapse through the ordinary link, boundary and
quality checks. Job 164 passes the rotated production regression, including
exact physical field checks at surface vertices and facet centroids and exact
volume. No length guard merely suppresses the original panic. Temporary
instrumentation is removed.

## Admission and regression evidence

Final native candidate SHA256 begins `55ae9ef2`, with Rust 1.95.0 and
experimental-surfaces. Job 173 admits the scalar adaptive octahedron:

| Metric | Measured value |
|---|---:|
| Tetrahedra | 2,565 |
| Nodes | 933 |
| Minimum MMG quality | 0.30043984 |
| Two-sided sampled distance | 8.9687e-17 |
| Relative volume error | 0 at reported precision |
| Boundary components / Euler characteristic | 1 / 2 |

Positive determinants, coherent facets, boundary-volume agreement and both
250000 admission caps pass. Adaptation performs 2518 collapses, 662 swaps and
3583 vertex moves. Spatial-size tests require actual refinement, bounded
callback batches, exact planes and q0.1; exceptions and hard refinement-budget
rejection also pass. Adaptive output does not yet export the named plane masks
available from the nonadaptive constrained route.

Jobs 169–172 regenerate all four original fixed study specimens. Their Tet4
counts remain 72757 Gyroid, 90792 Schwarz P, 72528 spinodoid and 218239 lattice.
Every original geometry, periodicity, quality and resource gate passes.
Job 174 independently replays both these meshes and the four preceding native
controls, resolving the old binary only by its original hash. Coordinates,
cells and boundary facets are exactly equal after canonical relabeling, with
no numerical tolerance. Fixture samples and specimen hashes are identical.
The same audit independently replays the adaptive octahedron. See
[the audit](../../handoff/macbook/adaptive-polyhedral-change-audit-20261024.json).

The full suite passes 176 tests with one optional wide-constraints skip. Native
core tests pass 49/49; Clippy with warnings denied, Rust formatting, changed
Python lint and typing pass. Repository-wide lint still reports 973 existing
findings, all in untouched handoff scripts and none in changed files.
A final facade-only correction converts nonadaptive bounds to an equal Python
list to match its declared Sequence API. The final suite is repeated after
that correction; the native binary and all four periodic control paths are
unchanged. Candidate source hashes are in
[the receipt](../../handoff/macbook/adaptive-polyhedral-native-candidate-20261024.json).

## Scope and next failure

This extends convex, compiled, nonperiodic linear meshes. It does not establish
arbitrary nonconvex CSG, nonlinear crease, periodic feature or quadratic feature
support. Next choose one bounded failing nonconvex sharp-feature control with
an independently known physical domain; first distinguish compiler rejection,
topology discovery and feature-preservation failures. Do not turn this into
another sizing sweep. Keep the four original specimens as regression controls.

No FE solve or reference promotion occurred. No new actual Microgen+MMG timing
ratio is claimed. Admitted actual Gyroid/spinodoid comparators and convergence
of all 24 lattice stiffness models remain open. The full objective is active
and unachieved.
