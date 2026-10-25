# Public periodic affine normal-jump recovery, 25 October 2026

The original public tapered-prism recipe now returns an admitted mesh.
The fixed field is `max(abs(x) + 0.3*abs(z) - 0.35, abs(y) - 0.2)` on
`[-0.5,0.5]^3`, periodic in z. Its exact volume is 0.22. The wall-normal
jump at the seam remains 33.4 degrees. No smoothing or specimen change was
used.

| Public recipe | Previous image | Trace-aware construction |
| --- | --- | --- |
| Result | Reject sampled error 0.0106425143 | Admit |
| Tetrahedra / nodes | No returned mesh | 3,811 / 844 |
| Minimum quality | No returned mesh | 0.309002647 |
| Sampled two-sided distance | No returned mesh | 1.24e-16 |
| Relative volume error | No returned mesh | 2.22e-16 |
| Periodic node pairs / cap triangles per end | No returned mesh | 45 / 66 |

Both runs use cells 24, four adaptation rounds, one polish pass, min/max
length 0.002/0.1, gradation 0.3, quality guard 0.3 and one numerical thread.
Admission remains distance 0.002, absolute volume 1%, minimum q 0.1,
positive Jacobians, coherent oriented facets, expected boundary topology,
matching periodic nodes/caps and 250,000 node/tet limits. The new run is
`macbook-periodic-affine-traces-prism-v1`, job 315. Earlier rejected scalar
and mapped runs remain rejected historical evidence.

## Construction change

The complete max-of-affine field can have equal end traces while individual
planes exchange roles. Requiring each plane to be invariant excluded this
case. The shared native certificate now matches the complete multiset of
normalized restricted affine equations at each enabled pair of ends. It
requires exact matches in the represented coefficients and rejects incomplete
or mismatched multisets. It is conservative and is not a claim that all equal
max-field traces have the same representation.

For a certified noninvariant plane set, the mesher cuts the physical planes
with zero independent cleanup passes, verifies bijective end nodes and
matching cap triangles, then runs coupled periodic adaptation and the
requested final polish. The seed checks reuse the final periodic admission
function. The original raw-field and sizing periodic checks and all final
quality/geometry gates remain. Individual named intersection fields retain
their existing periodic trace requirement.

The feature-mobility correction committed as 7c69c54 is also required. Its
private saved-seed result has 3,770 tets; it is separate from this 3,811-tet
public generation result. The public constructor cuts normalized planes,
so its seed need not be bit-identical to that earlier research input.

## Controls

Public regression 311 fails before the change; 314 passes both the public
recovery and the named nonperiodic-field rejection check. Native controls 312
pass complete-trace multiplicity/mismatch checks and rejection of mismatched
seed nodes or opposite cap diagonals. Core and binding Clippy pass with
warnings denied. Production Python type checking and changed Python lint
pass. An exploratory type check including the pre-existing polyhedral test
file reports 53 diagnostics in existing dynamic option dictionaries and
missing compiled-property annotations; it is not recorded as a passing check.

This recovers one previously rejected public recipe. It does not establish
universal meshing coverage, Microgen timing dominance or stiffness accuracy.
All four fixed study cases and all 24 lattice model tensors remain in scope.
No FE solve or reference promotion is part of this experiment.

Full native suite 316 passes 99 tests with one explicit research skip.
Python suite 318 passes 208 tests with one optional skip. Independent audit
317 confirms the public recovery against the frozen rejected recipe at the
same options, admission gates and fixture sample hash.

Fresh Gyroid, Schwarz, spinodoid and lattice controls 319–322 all admit.
The periodic column control 323 remains admitted. Independent replay 325
proves all four canonical study meshes and the periodic column exactly match
the preceding image. Audit 324 stopped at the historical failed spinodoid CLI
launch directory. That failure remains; 325 explicitly selects its valid v2
replacement. It did not require any mesh regeneration or setting change.

Final arm64 image SHA256 is
`02d1c61873b8a2ac34e721fb9d47aaea5471f6dc7f3d325452124cae6d8faef6`.
The candidate, recovery, study and local artifact receipts use the prefix
`handoff/macbook/periodic-affine-traces-` and date `20261025`.
The source receipt for audit 317 resolves its helper hash to a hash-identical
archived copy after the later audit-directory CLI option was added. Frozen
audit output is unchanged.

## Rigid orientation control

Job 326 rotates the same prism 0.37 radians about z and translates it
by (0.02,-0.03,0), preserving the box ends, exact volume and unclipped walls.
At unchanged settings it admits 3,879 tets / 865 nodes, qmin 0.304040718,
distance 1.89e-16 and volume error at roundoff, with 46 node pairs and
68 matching cap triangles per end. Independent replay 327 passes and
13 fixture controls pass. This adds one orientation control to the
public recovery, with no universal coverage or timing claim. Native image
02d1c618 and the four fixed study specimens are unchanged.
