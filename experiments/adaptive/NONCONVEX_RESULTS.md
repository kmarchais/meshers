# Nonconvex sharp-feature investigation, 24 October 2026

Meshers now recognizes bounded min/max combinations of convex affine fields.
A shared background partition constructs a union without overlapping volume
or separately meshed interfaces. Adaptive projection uses union membership,
while actual boundary facets determine physical feature identities. Unsupported
expressions retain their scalar route. This extends coverage; it does not prove
that arbitrary inputs or all option combinations admit a quality mesh.

The control is an L prism. Its cross-section is the square [-0.35,0.35]^2
minus the upper-right square (0.05,0.35]^2, extruded over z in [-0.3,0.3].
Exact volume is 0.24. Its independently triangulated boundary has twelve
vertices, twenty triangles, one component and Euler characteristic two.
Rotation by 0.37 radians about z and translation by (0.02,-0.03,0.01) provide
a second placement of the same physical specimen. No smoothing changes its
concave corner or edges.

## Original failure and independent preflight

Job 176 rejects the original adaptive scalar route at internal sampled error
0.004322177 against 0.002. Diagnostic exports retain the original external
admission limits. The first export attempt uses an invalid zero internal
tolerance and is rejected; it is retained in job 177. The corrected diagnostic
uses an explicitly loose internal gate only to obtain geometry, while the
independent gates remain distance 0.002, absolute volume 1%, q 0.1, positive
determinants, coherent facets, exact topology and 250000 tets/nodes.

The raw scalar seed in job 178 has true two-sided distance 0.018332398.
Adaptation in job 179 improves distance to 0.003302337 but still fails it and
has qmin 0.02421324. These are rejected controls, not admitted competitors.

The shared-cut seed preflight in job 180 uses 72738 tets and 16633 nodes.
Independent replay in 181 admits it with qmin 0.11379553, distance7.23e-14
and relative volume error about -1.53e-14. Plane masks select two overlapping
convex regions; their overlap appears only once in the final volume.
This is a geometry/resource preflight, not a production timing comparison.

## Corrections and retained failures

The first production candidate passes axis-aligned and default-snap tests,
but the complete suite in 191 rejects both rotated unsnapped variants.
The diagnostic in 192 finds 218 artificial feature vertices on boundary patches.
Arrangement planes belonging to an inactive union branch must not become
physical curves. Boundary-facet support now defines the retained constraints,
both in seed repair and adaptation. This reduces the adaptive failure from
23 bad cells to 12 in 194, without clearing the quality gate.

The shared smoother also required all fields to be negative. That is correct
for an intersection but rejects valid union branches. After checking membership
in any complete region, the rotated adaptive route passes with zero snapping.
The remaining nonadaptive one-round unsnapped recipe still rejects in 196 at
qmin 0.01783753 with 113 cells below 0.1. The failed native image beginning
fa81963c and its source hashes are preserved. Its original failure is not
silently reclassified.

A single two-round control in 197 passes at qmin 0.13148723 with 63925 tets.
Trying relocation before collapse in a single round worsens the result to
qmin 0.01065045 with 166 bad cells in 199, so that ordering change is reverted.
The final implementation keeps the original requested operation counts.
Default snapping uses one round. The rotated, unsnapped, nonadaptive recovery
explicitly requests two; it does not establish admission for the one-round
recipe. No sizing sweep or threshold relaxation is used.

## Validation scope

The final full Python suite in 223 passes 190 tests with one optional skip. Controls
cover both placements, both snapping choices and adaptation on/off, with the
two-round recovery stated above. Metadata is compared against actual field
values at fixed random points. Budget, quality and cancellation rejection pass.
The complete native suite in 202 passes 91 tests, including shared-volume,
physical-facet and union projection controls. Final all-target Clippy checks
cover both the standard core and optional wide binding build. Two lint findings
were corrected without changing behavior; the native image was rebuilt and
the measured controls repeated. A callback carrying a private metadata
attribute is also checked: only actual CompiledField metadata may select the
union route. Uncompiled callbacks retain their requested geometry.

The default build supports at most eight distinct planes; the optional wide
build supports 32. Both cap the representation at 32 convex regions. Periodic,
quadratic, multithreaded and interior-refinement feature routes remain on their
existing path. Nonlinear leaves and arithmetic outside a nonconvex min/max
combination are unsupported by the new metadata. Union output does not yet
export named constraint masks, including on the nonadaptive route. Seed snapping fallbacks are propagated into
native diagnostics rather than reported as zero unconditionally.

Original Gyroid, Schwarz P, spinodoid and ordered smooth-union lattice specimens
remain fixed. No FE solve, lattice model removal, reference promotion or new
actual Microgen+MMG timing ratio occurs. General nonlinear/periodic feature
coverage, admitted actual Gyroid/spinodoid comparators and convergence across
all 24 lattice stiffness models remain open. The full objective is active and
unachieved.

## Admitted production controls

| Placement and route | Tets | Nodes | qmin | Two-sided distance |
|---|---:|---:|---:|---:|
| Axis aligned, adaptive, default snap 0.2, one round | 7354 | 1965 | 0.26939647 | 7.85e-17 |
| Rotated/translated, nonadaptive, snap 0, two rounds | 63925 | 14518 | 0.13148723 | 1.86e-16 |

Both pass independent volume, facet, positive-Jacobian and topology checks.
The rotated single-round failure remains rejected. These controls do not
provide an actual Microgen+MMG competitor or a speed ratio.

Final fresh controls 216–219 admit all four original specimens. Audit 222 replays
both these meshes and the previous controls with the preserved 55ae9ef2 image.
Coordinates, tetrahedra and boundary facets are exactly equal after canonical
relabeling; fixture hashes and configuration dictionaries are identical. No
comparison tolerance is introduced. The final image begins 83247d4d, after
all-target core and all-feature binding Clippy checks. Its predecessors,
including failed fa81963c and pre-lint 317e3daf, are preserved by hash. The earlier
verifier source was restored from its committed bytes without editing frozen
results or manifests. See [final replay](../../handoff/macbook/nonconvex-change-final-audit-20261024.json)
and [candidate receipt](../../handoff/macbook/nonconvex-native-candidate-20261024.json).

Final lint-clean native unit/intersection run 224 passes 62 tests. Full native
run 202 passed 91 before the lint-only syntax corrections. The final Python
run 223 passes 190 with one optional skip. Historical failures remain archived.
