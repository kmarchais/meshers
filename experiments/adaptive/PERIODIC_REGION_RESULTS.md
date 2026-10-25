# Periodic re-entrant affine-region investigation, 25 October 2026

The next bounded control adds a re-entrant notch to the recovered periodic
prism. Its fixed field on `[-0.5,0.5]^3` is
`max(abs(x)+0.3*abs(z)-0.35, abs(y)-0.2, min(x-0.05,y-0.05))`, periodic in z.
The exact volume is `0.22 - 0.15*(0.275-0.05) = 0.18625`.
This auxiliary control does not replace any of the four study specimens.

Fourteen fixture controls pass. The explicit 18-vertex, 32-triangle target
has coherent closed boundary topology, exact oriented volume and six paired
nodes / four matching triangles on each periodic cap. Dense fixed samples
lie on its exact clipped boundary and raw scalar traces agree at the ends.

## Public failure

Job 328, `macbook-periodic-notched-wedge-scalar-v1`, rejects sampled geometry
error 0.0106362575 against tolerance 0.002 and returns no mesh. It uses
cells 24, four adaptation rounds, one polish pass, min/max length 0.002/0.1,
gradation 0.3, quality guard 0.3 and one numerical thread. External admission
remains distance 0.002, absolute volume 1%, q >= 0.1, positive Jacobians,
coherent facets, expected topology, matching periodic caps/nodes and 250,000
node/tet limits. No failure is reclassified or gate relaxed.

The compiled graph has eight planes and two region masks, 127 and 191.
The existing periodic convex-affine certificate does not admit a nonconvex
union of intersections. The scalar route remains the public fallback.

## Feasibility and physical seed

Job 329 independently audits a manufactured center-star mesh with 19 nodes
and 32 tets. It passes all geometry/quality/periodic controls with qmin
0.111167684 and exact volume. This is feasibility evidence only. It is
not a Meshers generator result, timing, FE reference or stress certificate.

Job 330 constructs one nonperiodic physical affine-region seed with zero
cleanup passes and a provisional generation quality floor of zero. It has
12,264 nodes and 54,534 tets, qmin 0.006212403 and roundoff volume error.
Its 160 end nodes and 252 cap triangles per end match exactly. The plane
trace permutation maps the complete region masks back to the same set,
and the represented end-plane trace coefficients agree exactly. Peak process
tree RSS is 0.402 GiB. The seed is explicitly unqualified and declares no
native periodic metadata; q >= 0.1 remains the final admission requirement.

The next bounded experiment is private coupled adaptation of these exact
seed bits with the actual region geometry and final gate. Its result must
be read before changing public eligibility. Native image 02d1c618 and all
four study meshes remain unchanged by the diagnostic harness and test code.
All 24 lattice models and the original stiffness/geometry/quality/periodic
gates remain in scope. No FE solve or reference promotion is part of this
investigation.

Private replay 331 passes the actual native final and independent geometry
gates at 3,896 tets / 898 nodes, qmin 0.301404653, roundoff geometry and
zero nominal volume error. It has 54 paired nodes and 82 matching cap
triangles per end. Its field and gradient reproduce the compiled expression
at 1,028 probes. Independent audit 332 verifies all three controls, preserved
public failure, exact seed bits and complete region trace multisets. The
private result remains explicitly unqualified as a public workflow. The next
step is a conservative complete-region end-trace certificate, not a sizing
change. Receipts are `handoff/macbook/periodic-region-investigation-audit-20261025.json`
and `handoff/macbook/periodic-region-local-artifacts-20261025.json`.

## Public recipe recovered

The conservative region certificate now compares the full multiset of
conjunction traces at opposite ends. It preserves plane membership within
each region and region multiplicities. Matching the unrestricted plane list
alone is insufficient. Signed zero is canonicalized; nonfinite equations,
invalid masks and unmatched branch multisets reject. Periodic regions must
use coupled adaptation. Nonadaptive periodic region generation is rejected.
Named individual-field periodic admission is unchanged.

Red public regression 333 fails on the preceding image. Focused native
controls 334 pass, including wrong branch grouping and multiplicity, invalid
masks, trace overflow and the mandatory adaptation gate. Build 335 installs
arm64 image `7bad4a454648f1bed28636929230411e68c9ed5c494664f134010809fcc5e59f`.
Public regression 336 passes. The full fixture record 337 admits 3,896 tets /
898 nodes, qmin 0.301404653, roundoff distance and zero nominal volume error,
with 54 paired nodes and 82 matching cap triangles per end.

Independent replay 339 verifies the public result at exactly the failed
recipe's options, gates and fixed sample hash. Its canonical mesh arrays are
identical to the private saved-seed recovery. Historical rejection 328 remains
unchanged. This is an actual public generation result, distinct from the
manufactured feasibility control and private adaptation.

Full native suite 338 passes 101 tests with one explicit research skip.
Python suite 340 passes 211 tests with one optional skip. Production Python
type checks, changed Python lint, and core/binding Clippy with warnings denied
pass. No FE reference, stiffness accuracy or comparative workflow speed
claim follows from this public failure recovery.

Fresh fixed controls 341–344 admit all four study cases. Column control 345
also admits. Independent audit 346 proves all four canonical study meshes
and the column exactly unchanged from native 02d1c618. These are geometry
controls, not refreshed stiffness or comparative timing results. The next
main-study control is the unchanged recorded lattice26 curved Tet10 recipe,
starting with geometry and resources. All 24 reference models remain retained.
