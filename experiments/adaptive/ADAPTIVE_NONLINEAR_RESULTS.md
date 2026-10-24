# Adaptive nonlinear intersection failure recovery

The unchanged bored-prism adaptive recipe now passes. It retains the cylinder
and planar-cap identities through actual splitting and collapse, then projects
vertices onto their incident constraints together. This follows Kevin's priority:
fix rejected geometries before optimizing admitted cases.

## Matched failure control

Both runs use cells24, four adaptation iterations, one requested optimization
round, snap0.2, target distance0.002, minimum quality0.1 and cap250000. The
physical box and offset cylindrical bore, exact volume, independent surface
queries and admission gates are unchanged.

| Evidence | Before | Candidate |
|---|---:|---:|
| Status | rejected | admitted |
| Sampled distance | internal rejection0.003393814 | independent two-sided0.001977286 |
| Minimum MMG quality | unavailable on rejected output | 0.168559352 |
| Tetrahedra / nodes | unavailable on rejected output | 12955 / 3004 |
| Absolute relative volume error | unavailable on rejected output | 0.106744% |

The old result is `runs/macbook-bored-prism-adaptive-postcontrol-v1`; the first
candidate is `runs/macbook-bored-prism-adaptive-retained-v1`. The candidate also
passes positive Tet4 Jacobians, coherent boundary facets, used-node accounting,
one boundary component, Euler characteristic0 and boundary/cell volume balance.
These are sampled geometry checks, not a global intersection certificate.

## Implementation and limits

Only eligible compiled outer-maximum intersections use the new path. The seed
is made by the existing simultaneous constrained cutter. Actual boundary facet
labels propagate through topology changes; vertex labels come from incident
facets. Joint normalized Newton projection preserves the active surfaces, with
rank-aware handling of redundant locked-coordinate equations. Each incident
field supplies its own curvature; the jump of a maximum-field normal at a
crease is not interpreted as smooth-patch curvature. Spatial sizing, budget,
callback errors and final geometry/quality rejection remain active.

Periodic, quadratic, custom-gradient, multithreaded and interior-refinement
requests retain their existing routes. This does not establish arbitrary
nonlinear Boolean or periodic feature handling. It does not certify stress or
continuum accuracy. No actual Microgen+MMG timing ratio or FE reference is
promoted. All four original specimens, admitted actual Gyroid/spinodoid
comparators, and all24 lattice stiffness convergence remain in scope.

## Verification checkpoint

Python/fixture run250 passes185 tests with one optional skip. Run249 preserves
a test assertion mistake: the budget was correctly rejected by the seed cutter,
but the test expected a different error phrase. The corrected assertion accepts
the actual budget message; the rejection behavior was not weakened.
Core and binding Clippy245–246 and final current-source all-target/all-feature
checks pass. Native suite251 passes93 tests, including scaled nonlinear joint
projection with redundant locks and the existing large twisted-sheet regression.
Its complete debug run takes490 seconds and peaks at1.343GiB; that is test-suite
resource evidence, not a mesh workflow timing. Final build252 reproduces the
same native hash dc707190517fecbe8b8242593c5304a800fd265323335275fe8807085b9e6268.

Fresh controls253–256 admit all four original specimens. Independent replay258
checks their actual saved arrays against the preserved b2814b0e baseline:
coordinates, tetrahedra and boundary facets are exactly equal after relabeling,
and sample/specimen hashes are identical. Counts remain72757,90792,72528 and
218239 tets for Gyroid, Schwarz P, spinodoid and lattice respectively. The same
replay independently checks the bore and the affected hollow shell, including
geometry, topology, positive volumes, quality, resource caps and facet coherence.
See [the audit](../../handoff/macbook/adaptive-nonlinear-change-audit-20261025.json)
and [native receipt](../../handoff/macbook/adaptive-nonlinear-native-candidate-20261025.json).

The shell remains admitted with9469 tets /2207 nodes, qmin0.302376,
distance0.001862454 and volume deficit0.540544%. Its previous admitted scalar
route had10275 tets, qmin0.436980,distance0.001443754 and volume deficit0.334739%.
This is a metric tradeoff inside the gates; it is not an all-metric improvement.
Failure coverage takes priority over optimizing these admitted recipes.

All jobs through258 are terminal. The broader objective remains unachieved.
