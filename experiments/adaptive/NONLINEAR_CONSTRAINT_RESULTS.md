# Retained nonlinear intersections, 24 and 25 October 2026

The bored-prism control has a square section [-0.35,0.35]^2, caps at z=±0.3,
and a cylindrical bore of radius 0.12 centered at (0.03,-0.02). Its exact volume
is 0.7²×0.6−π0.12²×0.6. The boundary is connected with Euler characteristic 0.
An independent distance function minimizes distance to four finite rectangular
sides, two perforated caps and the finite cylinder. It handles void points and
exterior corners, rather than using the composite scalar as a distance oracle.
Deterministic exact surface samples include both circular rims. Ten fixture
controls validate surface membership, exact volume and analytic probe distances.

## Failure and controlled reference

Job 225 rejects scalar adaptive 24 at sampled error 0.003393814. The unchecked
raw seed 226 has reverse error 0.02006282 at a bore/cap rim, more than 0.21 away
from outer edges. This localizes an actual nonlinear rim miss. Diagnostic
exports use a loose internal gate only to save geometry; all independent
admission thresholds remain unchanged.

Seven named constraints describe exactly the same negative sublevel set.
The named 24 control 227 still rejects at error 0.003269072. Its saved diagnostic
228 has qmin 0.14558756 and volume error 0.10197%, but fails distance on the
cylinder near a cap. Outer planes and rim identities alone do not resolve
ordinary curvature error at that grid.

Before another run, a recorded quadratic chord-error estimate chooses cells 34:
ceil(24×sqrt(0.003269072/0.002)×1.1). It predicts distance 0.00162888 and about
152428 tets, below the 250000 cap. This is one curvature-derived control, not
a resolution sweep. Its one-round quality miss 0.09784093 in 229 remains a
failure. A second explicitly requested quality round in 230 admits 145924 tets,
30384 nodes, qmin 0.12139989, distance 0.001058141 and volume error 0.0502226%.
No gate is relaxed. At identical cells 34, snap 0.2 and two nonadaptive rounds,
the original fused field 231 still rejects at error 0.004960761.

## Implementation and evidence

The compiler retains the actual leaves of eligible outer maximum expressions,
with their own analytic gradient outputs. Native validation requires the roots
to match the parent field, so metadata cannot substitute another geometry.
Evaluators are compiled lazily, cached, and own their executable memory after
the parent is released. The parent evaluator is compiled once. Existing affine
and affine-union routes retain precedence; custom gradient policies are not
replaced. Retained field preparation and plane preparation are charged in the
public preparation diagnostics and complete generation call.

Eligible nonperiodic, nonadaptive Tet4 calls reuse the existing constrained
cutter. Affine leaves expand deterministically to their plane constraints.
The original scalar maximum remains the physical domain; no specimen, bore
radius, corner or field order is smoothed or altered. The outer-max identity
is checked at fixed random points, and retained gradient, ownership, caching
and invalid-root controls exercise the native interface.

The repaired scalar route 237 admits with the same reported counts, quality,
distance and volume as the explicit named control. It reports seven retained
constraints and exports the rim feature identities. The complete Python suite
236 passes 198 tests with one optional skip. Final all-target/all-feature binding
Clippy 238 passes; changed Python lint, formatting and typing pass.

## Limits and remaining work

This step covers a bounded outer maximum of analytic fields. It does not
decompose arbitrary nonlinear Boolean graphs or silently replace adaptive
requests with a nonadaptive mesh. Periodic, quadratic, multithreaded and
interior-refinement requests keep their existing routes. Adaptive nonlinear
constraint identities must still be carried through splitting, collapse and
projection; the original adaptive failure is not promoted by this result.

The original Gyroid, Schwarz P, spinodoid and ordered smooth-union lattice
remain in scope. No FE solve, reference promotion, model removal, local-stress
certification or new Microgen+MMG timing ratio occurs. Admitted actual Gyroid
and spinodoid comparators and all 24 lattice stiffness convergence remain open.
The complete objective is active and unachieved.

## Final replay and next step

Audit 243 proves exact canonical equality of the repaired scalar bore mesh and
the admitted named-control mesh, using identical cells, rounds, snapping,
quality and distance settings. It also replays all four original specimens
against their preserved 83247d4d controls. Coordinates, tetrahedra and boundary
facets are exactly equal after relabeling. Specimen/sample hashes remain fixed.
See [the audit](../../handoff/macbook/nonlinear-change-audit-20261024.json) and
[the native receipt](../../handoff/macbook/nonlinear-native-candidate-20261024.json).

The negative adaptive postcontrol 244 still rejects at exactly 0.003393814.
The implementation honors adapt=True; it does not substitute the admitted
nonadaptive recipe. The next bounded change must propagate actual retained
constraint identities through adaptive splits, collapse and joint projection.
Do not infer nonlinear feature identities from approximate nodal residuals or
repeat a sizing sweep. Keep the original four-study gates and the remaining
actual-comparator and lattice-reference problems in scope.
