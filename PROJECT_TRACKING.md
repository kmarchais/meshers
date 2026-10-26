# Meshers project tracking

## Resumed investigation, 26 October 2026

Best completed mesh-workflow score is now **1.215175763%** across all 24
retained lattice models, after the curved quality pass described below. It
still fails 1%.

The preceding longest-edge closure prototype reaches
**1.216071201% worst disagreement across all 24 retained lattice models** at
250,000 nodes and 162,365 tets. This improves on the earlier directional
prototype's 1.608556732%, but still fails the unchanged 1% target. Six fresh
loads, independent reaction/energy/periodic and positive 64/125-point audits
pass. Exact preflight arrays, native binary and first operator are verified.
The full diagnostic workflow takes 502.321 seconds and peaks at 7.9395 GiB
under its explicitly declared 10 GiB ceiling. Neither earlier 8 GiB failure
is reclassified.

Composed-map verification of the actual parent fields passes, with a
Galerkin balance of 4.2921e-10. Only one of 8,697 remaining requested edge
orbits is quality-blocked. This replaces the earlier quality obstruction
with the unchanged node-cap obstruction; it does not satisfy the whole metric.

Against admitted modified actual Microgen+MMG, it has 162,365 versus 165,722
tets at the same 250,000 nodes, sampled distance 0.000278662 versus
0.001701239, absolute volume error 0.03124% versus 0.49053%, and minimum
quality 0.206873 versus 0.100001. However, median quality is worse, 0.664257
versus 0.796052, as is the first percentile, 0.383975 versus 0.429975.
Stiffness disagreement remains worse than actual's 0.994670612%. This is an
experimental Python route and an explicit tradeoff, not an all-metric win.

Residual screen 561 passes positive 64/125-point checks and finds a
0.102410865% possible energy correction, above the unchanged 0.05% solve
screen. Its hypothetical 1.116215063% all-model disagreement still fails 1%
and is not an equilibrated result. The largest-gain direction is 48.18% in
curve-vertex cells, 28.78% in other wall cells and 19.43% in neither group.
The remaining signal is less concentrated at the curves than in the earlier
dense-mesh investigation.

Consistent coupling preflight 562 reproduces the parent operator, passes
15/125-point Gram and gain checks within 8.83e-5, and completes four-level
AMG setup at 7.6076 GiB under the declared 10 GiB ceiling. Job 563 completes
six fresh zero-macro auxiliary responses on the unchanged physical domain.
Complete analytic-bubble/P2 field audit 564 passes positive 15/64/125-point
integration, fresh coarse/modal equilibrium, reactions and periodicity.
Galerkin balance agrees within 3.800e-11. The enriched diagnostic gives
**1.091221736%** worst disagreement across all 24 models, still failing 1%.
It is not a native mesh candidate or workflow timing result. Job 565 screens
the residual of this complete enriched field against the unchanged 0.05%
signal threshold. It finds only **0.015002691%** additional trial gain,
below 0.05%, so this residual family receives no further auxiliary solve.
Its hypothetical 1.076536305% disagreement is not a new equilibrium tensor
or error bound. No new density choice or reference is promoted. The complete
reference investigation is recorded in
`handoff/macbook/lattice-closure-reference-investigation-20261026.json`.
Reviewed evidence is in
`handoff/macbook/lattice-longest-edge-assessment-20261026.json`.

Job 543 was not restarted. Independent field audit 544 passes reactions,
Hill–Mandel energy, periodic jumps and positive 64/125-point integration.
Quadrature sensitivity is 4.1774e-8. Assessment 545 retains all 24 unchanged
lattice model tensors and gives **1.608556732% worst disagreement**, improving
on compact native's 2.412306220%, dense native's 1.935110817% and Uniform37's
2.085389679%. It still fails the unchanged 1% target. The prototype has 68,476
tets and 115,618 nodes. It includes experimental Python refinement and is
not a native-only Meshers performance result. No comparative timing or
reference promotion follows.

Read-only audit 546 composes all 15 saved parent maps and prolongs the actual
parent's six equilibrium fields for verification only. The freshly generated
seed arrays exactly match the parent. Mapping and field reproduction errors
are at most 4.441e-16; periodicity passes. The Galerkin energy balance agrees
within 1.3812e-9. The improvement is a consistent fixed-domain FE correction.
No parent solution entered mesh generation.

The same audit identifies 2,710 quality-blocked requested edge orbits out of
2,767. Of 2,706 blocking tets, 723 have no legal single-edge midpoint split.
A separate local full-edge template check, job 547, shows that full red splits
are not a sufficient repair. Across the whole final mesh, 1,011 tets have no
legal single-edge split; only 330 admit one of the three local red templates,
while 681 still fail quality. These counts cover different populations from
the 723 blocking tets. No conforming red/green mesh, curved child admission,
periodic closure or budget admission is inferred from this local check.

The next change must address how greedy directional splits create poor
shapes, rather than add another nearby sizing choice. The quality floor,
geometry gates and both 250,000 admission caps remain unchanged. Schwarz's
earlier win and the Gyroid/spinodoid actual-comparator gaps remain as recorded.
Evidence is retained in
`handoff/macbook/lattice-directional-assessment-20261026.json`.

### Quality-aware continuation preflight

The next policy keeps the same curve-distance law, ratio 5, physical domain,
60-pass resource guard and admission caps. It accepts a requested short-edge
split when both children pass q 0.1 and each child either satisfies the sampled
metric or retains a legal locally-longest split. Other requests use longer-edge
closure, whose leaves receive the same check. This is a one-step local guard,
not a guarantee that all later global refinement will be possible. No FE
fields enter selection.

A geometry-only regression extracted from level 3 of the earlier greedy run
has an admitted child at q 0.114833 that cannot split any of its six edges.
Frozen baseline 566 fails this regression. Green 567 passes 20 controls;
expanded controls 568 pass all 22, including periodic macro jumps, exact maps,
safe short-edge acceptance and already-completed metric handling. Source and
fixture hashes are frozen in
`handoff/macbook/lattice-quality-closure-selection-20261026.json`.

Geometry-only job **569** completes 60 levels at **224,358 nodes and 140,781
tets**. Every level passes geometry and periodic gates; maximum volume drift
is 4.663e-15. Final minimum quality is 0.100020510 and distance remains
0.000278662. It stops at the pass guard with two orbits still accepted in the
last pass. The 2,194 unresolved metric edges, maximum ratio 16.224, are not
reclassified as complete refinement. The node cap is not reached.

Exact operator/AMG probe 570 passes at 6.1954 GiB. Job **571** regenerates this
same frozen recipe and computes six fresh loads with the verified exec
transition and declared 10 GiB ceiling. Its combined deadline is explicitly
1,800 seconds after measured geometry admission takes 711.95 seconds. This
extends only the resource deadline, not an admission or accuracy gate. The
recipe is frozen in
`handoff/macbook/lattice-quality-closure-solve-selection-20261026.json` and
preflight evidence in
`handoff/macbook/lattice-quality-closure-preflight-evidence-20261026.json`.
The fresh workflow completes six loads at the exact preflight arrays and
operator hash. Field/positive-quadrature audit 572 passes, followed by all-24
assessment 573 at **1.315624132%** worst disagreement, still failing 1% and
worse than mandatory-longest-edge's 1.216071201%. Median quality is 0.577905,
also worse than the earlier prototype's 0.664257 and actual's 0.796052.
The full diagnostic workflow takes 934.233 seconds. No performance win is
claimed.

Nested audit 574 passes with Galerkin balance 2.7953e-10. None of the 1,385
blocking tets lacks every legal single-edge split, so the real dead-end
regression is addressed. However, 1,407 of 1,732 remaining requested edge
orbits remain quality-blocked. Fewer nodes/tets is a recorded tradeoff, not
promotion. This closure-policy variant receives no further density or nearby
policy retry. The best completed mesh-workflow score remains 1.216071201%
from the mandatory-longest-edge route. The next diagnostic concerns existing
edge-removal topology operations after refinement, with unchanged curved
cavity boundaries, rather than another sizing-law choice.
Evidence is recorded in
`handoff/macbook/lattice-quality-closure-assessment-20261026.json`.

### Curved edge-removal feasibility

Read-only probe 575 mirrors existing native valence-3/4 edge-removal patterns
on the mandatory-longest-edge prototype. Curved bipyramid controls pass.
It finds 10,876 locally feasible variants on 9,306 edges after checking
corner-quality improvement, exact cavity Tri6 boundary traces, volume and
positive Jacobians at the admission points plus a positive 125-point rule.
These variants overlap. No global mesh is changed or admitted by that probe.

Global preflight 576 then admits one disjoint pass with **8,508 operations**,
176 valence-3 and 8,332 valence-4. Counts decrease from 250,000 nodes and
162,365 tets to **249,824 nodes and 162,189 tets**. Curved boundary coordinates
and all boundary Tri6 associations are identical. Global geometry, volume,
coherent facets, point uniqueness, positive sampled Jacobians and periodic
counts pass. Every recorded corner-quality quantile improves or stays equal.
Minimum rises from 0.206873 to **0.222689**, first percentile from 0.383975 to
0.413661, and median from 0.664257 to **0.685789**. Median remains below
actual's 0.796052. Fresh fields and accuracy are assessed below.

Exact operator/AMG resource check 577 passes at **6.5114 GiB** under its
declared 10 GiB ceiling. This topology change requires fresh fields and is
not a nested-field correction or a timing result. Local and global evidence
are recorded in
`handoff/macbook/lattice-quadratic-edge-removal-feasibility-20261026.json` and
`handoff/macbook/lattice-quadratic-quality-preflight-evidence-20261026.json`.

The curved pass is integrated into the experimental fresh-generation driver.
Controls 578 pass all 23 checks, including no-op behavior and rejection of
overlapping operations. Live compact control 579 passes; its saved pre-pass
arrays reproduce the original native baseline exactly. A nonzero quality
pass explicitly marks the old approximation space as not preserved and is
rejected by the nested-field audit. The all-model assessor independently
checks the saved pre-pass boundary and volume instead.

Job **580** regenerates the mandatory-longest-edge source, applies the same
one-pass policy and solves six cold loads with the verified exec transition,
10 GiB and a 1,200-second whole-worker deadline. No saved mesh or FE field is
replayed in that live workflow. It completes at exact preflight arrays and operator hash. Independent field
audit 581 passes reactions, energy and periodic jumps, with positive125
quadrature sensitivity 4.0481e-8. Assessment 582 gives **1.215175763%** worst
disagreement across all 24 unchanged models. This tiny accuracy change from
1.216071201% accompanies the quality/count improvement, but still fails 1%.
The fresh diagnostic workflow takes 525.337 seconds and peaks at 8.0744 GiB
under its declared 10 GiB ceiling. Approximation spaces change, so no nested
field or Galerkin-transfer claim is made. No reference is promoted and no
comparative timing victory follows. Evidence is in
`handoff/macbook/lattice-quadratic-quality-assessment-20261026.json`.
The recipe and unchanged admission gates are frozen
in `handoff/macbook/lattice-quadratic-quality-solve-selection-20261026.json`.

### Longest-edge closure and memory handoff

The experimental refinement now propagates short-edge requests to longer
incident edges before bisection. Eighteen periodic, directionality and nesting
controls pass. The initial control attempt used two invalid test fixtures;
the corrected short edge and supported four-cell cube pass without changing
the closure algorithm. Both attempts remain recorded.

Geometry preflight 551 reaches **250,000 nodes and 162,365 tets** in 34 levels.
All levels pass geometry, coherent facets, periodicity and sampled curved
Jacobians. Quality rejects no orbit. Final minimum quality is **0.206873413**,
versus 0.100030317 for the earlier greedy prototype. Maximum volume drift is
3.775e-15 and map/control-field errors stay below 3.886e-16. It stops at the
unchanged node cap, not at the 60-pass guard. Its 10,077 remaining requested
edges make clear that the metric is still incomplete. No stiffness result is
yet available for this mesh.

Fresh operator/AMG setup 552 passes at 7.2485 GiB, but the combined live
workflow 553 hits the unchanged 8 GiB watchdog before completing load one,
with 8.595 GiB observed RSS. That attempt remains a resource failure. Exact
regenerated arrays match the preflight. The retained-generation-memory
hypothesis is tested by replacing the worker process image at the solver
boundary with `os.execv`, keeping the same PID, exact mesh hash and original
workflow clock. Imports and reload remain charged to the workflow.

Full compact control 554 completes six loads through this handoff. The
runnable verifier reproduces the baseline mesh, fields, stiffness tensor and
operator hash exactly. Boundary RSS falls from 574,603,264 to 215,105,536
bytes. This validates the transition, not the large candidate's memory use.
Job **555** completes the first fresh load with the exact expected operator
hash, then hits the same 8 GiB limit at 8.092239 GiB observed RSS on the next
load. This second attempt remains failed. The process handoff reduces
boundary RSS from 1,413,283,840 to 221,904,896 bytes, but full multi-load memory
still exceeds setup-only feasibility.

A recorded resource decision selects **10 GiB**, about 24% above that measured
failure, while the 24 GiB host has 13.9 GiB available. This uses Kevin's prior
authorization for a conservative measured Mac ceiling around 12 GiB. Default
workers remain at 8 GiB. Larger Mac workers require an explicit argument,
matching manifest/watchdog records and a matching declared assessment plan.
Control 556 verifies the 10 GiB configuration and result metadata. Job **557**
then completes the same frozen geometry and all six cold loads under this
declared ceiling and the unchanged 1,200-second bound. Jobs 558, 559 and 560
complete its independent field audit, all-24-model assessment and nested
energy audit, as reported at the top of this entry. Both earlier 8 GiB
failures stay failed. No new native-only feature, reference or comparative
performance winner is claimed.

Selection receipts, the failed resource attempt and compact control evidence
are recorded in `handoff/macbook/lattice-longest-edge-preflight-evidence-20261026.json`
and `handoff/macbook/lattice-exec-memory-protocol-20261026.json`.

## Daily pause summary, 25 October 2026

**Paused at Kevin's request.** Today improved Meshers correctness and produced
better geometry and quality on some routes. We have **not** demonstrated that
Meshers beats actual Microgen+MMG on all four cases. Lattice stiffness accuracy
remains the main unresolved native weakness; the Gyroid and spinodoid actual
comparator gaps also remain.

### Confirmed improvements

- The cap repair improved the curved spinodoid's
  sampled distance from 0.001333 to 0.000304 and minimum quality from 0.269
  to 0.301, with slightly fewer elements. Six fresh loads and independent
  audits give 0.183380673% worst retained-model disagreement.
- The strict linear spinodoid failure is repaired. Distance improves from
  0.004834255 to 0.001079175 at 131,011 tets. Minimum quality is 0.100653,
  above the unchanged 0.1 gate but below the old result's 0.152884.
- Commit `0fc144c` fixes an existing public option: scalar adaptation had
  silently ignored `interior_refine_nodes`. The rebuilt arm64 package passes
  105 core tests and 183 Python tests, with one ignored core test and one
  skipped Python test. Both feature Clippy checks pass. All four default
  point, tetrahedron and surface arrays reproduce exactly.

### Latest completed lattice comparison

Sizes were chosen before evaluating each new candidate's stiffness. All
figures below retain the fixed specimen, all 24 lattice model tensors and
the original geometry, quality, curved-Jacobian and periodic gates.

| Route | Tets | Nodes | Worst stiffness disagreement | Result |
|---|---:|---:|---:|---|
| Current compact native default | 41,340 | 65,721 | 2.412306220% | Fails 1% |
| Native adaptation at comparable count | 164,180 | 235,258 | 1.935110817% | Fails 1% |
| Uniform37 native alternative | 157,521 | 231,026 | 2.085389679% | Fails 1% |
| Admitted modified actual Microgen+MMG route | 165,722 | 250,000 | 0.994670612% | Passes 1% |

Both new native alternatives completed six fresh loads and independent field
and positive-quadrature audits. Dense adaptation has better minimum quality
and geometry than actual, but worse median quality, 0.680 versus 0.796.
Uniform37 improves median quality to 0.930 and sampled distance to 0.0001783,
but does not close accuracy. Its extreme upper quality tail is slightly worse
than actual. Neither is an all-metric win. No new paired timing ratio is
claimed. These two size choices are closed without nearby parameter retries.
Their evidence and tools are committed and pushed in `32389d2`.

### What the failure investigation established

The dense native domain supports a separately audited enriched FE result at
1.516977609% worst disagreement. It still fails 1%, and is not a native mesh
workflow. A further trial offers only 0.044029171% gain, below its existing
0.05% solve screen, so no additional auxiliary solve follows. The earlier
fixed-domain h5 result remains at 1.325163420%, also failing. No reference is
promoted.

A signed geometry audit finds re-entrant periodic wall junctions at all 278
examined single-plane points, with openings about 232–289 degrees. This is
geometry evidence, not a continuum or local-stress certificate. One fixed
curve-distance grading pilot exceeds the unchanged intermediate tet cap,
250,002 versus 250,000, and receives no solve.

A separate MMG tensor control on a saved Meshers seed is explicitly ineligible
as either native Meshers performance or an actual Microgen comparator.
Adapter errors and their corrections are preserved. After boundary roundoff
is measured and exact source coordinates restored, the control fails sampled
curved-Jacobian positivity. It receives no solve. It does not establish that
an arbitrary directional metric is feasible, especially because MMG can
modify metrics at required boundaries.

### Experimental Meshers refinement and stopping point

The new prototype combines a fresh native seed with Python exact-map,
periodic edge-star bisection. Its ranking uses geometry only, with no FE
fields or MMG calls. Each orbit is checked against physical quality and both
250,000 admission caps before allocation. A regression exposed midpoint-only
metric sampling that missed fine curve-endpoint targets. Endpoint-aware
sampling fixes it; all 16 directionality and nesting controls pass.

The corrected preflight reaches **115,618 nodes and 68,476 tets**. All 15
levels pass geometry, quality, sampled curved Jacobians and periodic checks.
Volume changes by at most 1.554e-15 relative, and map/control-field
reproduction errors stay below 3.331e-16. Minimum quality is 0.100030317.
The metric is still incomplete: 3,511 edges exceed its split threshold and
quality checks block further requested splits. Operator/AMG setup passes at
3.927 GiB peak RSS. This is an experimental route, not a shipped native
feature or a performance winner.

Local queue job **543** finished successfully before stopping. Its output is
[the directional workflow](/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-directional-restriction-amg-v1/result.json).
All six fresh load fields are complete and finite. Regenerated point, Tet10
and Tri6 arrays exactly match the admitted preflight. The reduced operator
hash also matches the resource preflight. Generation took
124.251 seconds, the six-load solve 142.622 seconds,
and peak process-tree RSS was 4.026 GiB. These are single
diagnostic timings, not a comparative speedup.

**The new prototype's stiffness accuracy is not yet assessed.** Independent
complete-field/positive-quadrature audits and the all-24-model assessment
remain for the next working session. No new reference or native winner is
promoted. No job from this session remains running or queued.

On resumption, read job 543's `worker.json`, `result.json` and complete fields
first; do not restart this completed job. Audit those fields, retain the exact
preflight-array and operator identity checks, then assess all 24 retained models using
the frozen endpoint-aware recipe. Read the outcome before selecting another
change. Do not resume closed size sweeps or overwrite frozen results.

Schwarz's earlier repeated win against an explicitly modified actual route
stands. Gyroid and spinodoid still lack admitted actual comparators in the
tested routes. All four cases remain in scope. The 1% stiffness target,
distance 0.002, absolute volume 1%, q 0.1, positive sampled curved Jacobians,
coherent facets, periodic gates and 250,000-node/tet admission caps stay fixed.
All earlier active or queued descriptions below are historical.

## Matched lattice count study still fails accuracy, 25 October 2026

The adaptive scalar path now honors its existing interior-node request.
Commit `0fc144c` passes 105 core tests, 183 Python tests and both feature
Clippy checks. One core test is ignored and one Python test is skipped.
All four default point, tetrahedron and surface arrays reproduce exactly.

A count-only selection produces 164,180 Tet10 elements and 235,258 nodes,
versus the admitted actual route's 165,722 and 250,000. All geometry, quality
and periodic gates pass. Six fresh loads and independent field/positive-rule
audits pass, but worst disagreement across all 24 retained lattice tensors
is **1.935110817%**, above 1%. Compact native remains at 2.412306220%.
Minimum quality improves on actual, but median quality is 0.680 versus 0.796.
No new reference, native winner or comparative timing claim follows.
Global maximum-length adjustments stop, and earlier failures stay rejected.

One existing uniform-seed alternative is selected from recorded quality and
count scaling before its stiffness evaluation. Its 37-cell preflight passes
at 157,521 tets and 231,026 nodes, minimum quality 0.339372, median 0.929516
and sampled distance 0.000178307. The upper quality tail remains slightly
worse than actual. Exact operator/AMG memory setup and six fresh loads pass.
Independent field audits and positive-rule checks also pass. Worst disagreement
is **2.085389679%**, still failing 1% and worse than the dense adaptive result.
This distinct topology choice is closed without a nearby resolution retry.
Both count-selected native alternatives remain rejected on accuracy. The source-link error in the initial selection receipt is preserved
and corrected without changing its counts or 37-cell choice.

Schwarz's repeated modified-actual-route win stands. Gyroid and spinodoid
still lack admitted actual comparators in the tested routes. Lattice reference
convergence remains unresolved. All four fixed cases, all 24 lattice models,
original admission gates and 250,000-node/tet caps remain in scope.
Earlier active and queued descriptions below are historical.

## Native cap repair passes four-case geometry checks, 25 October 2026

The fixed linear spinodoid now passes strict independent admission. Distance
improves from 0.004834255 to 0.001079175 at 30,336 nodes and 131,011 tets.
Minimum quality is 0.100652964, above the unchanged 0.1 gate but lower than
the old 0.152883712. This fixes a failure; it is not all-metric dominance.
The exact source passes 104 native tests, Clippy, 182 Python tests and all
four curved geometry/quality/periodic gates. Curved spinodoid improves distance
from 0.001333 to 0.000304 and minimum quality from 0.269 to 0.301, with slightly
fewer elements. Curved lattice also improves geometry and quality at nearly
the same count. Both changed meshes need fresh fields before accuracy claims.
Lattice preflight 476, six fresh loads 477 and independent audit 478 pass.
Worst disagreement across all 24 retained tensors is 2.412306220%, still
failing 1%. Fresh spinodoid loads and independent audit also pass; its six
distinct archived model tensors give 0.183380673% worst disagreement. Every
entry in all four archived spinodoid reports is checked without refitting.
Gyroid and Schwarz arrays are exactly unchanged. No new reference or
four-case performance winner is promoted.
See [native cap repair evidence](experiments/adaptive/NATIVE_CAP_REPAIR_RESULTS.md). Earlier checkpoints below are historical.


## Native periodic normals corrected, compact lattice accuracy still fails, 25 October 2026

Jobs 379 through 403 reproduce and fix physical periodic cap gradients and
cached edge normal images. The rebuilt arm64 native library passes 101 core
tests and 182 Python tests, with one skip in each suite. All four fixed-case
geometry, quality, and periodic gates pass. Gyroid, Schwarz, and spinodoid
keep their meshes within 4e-14. The unchanged lattice26 recipe drops from
123,694 to 41,131 tets and from 189,180 to 65,382 nodes. Minimum corner quality
improves from 0.20735 to 0.36161, and sampled distance improves slightly.

The smaller quadratic space still fails all-model stiffness at 2.427923243%.
The independently audited full cubic space on the same physical domain also
fails at 1.647040948%, correctly normalized against each of all 24 models.
Neither result promotes a reference, native cubic API, or matched-accuracy win.

Three independent FE controls verify exact variational gains of periodic
edge-star splits. Fresh coupled first-level solve 407 improves worst disagreement
to 1.962342587% at 116,569 nodes and 73,294 tets. Full second-level preflight 408
passes but is not solved, because most selected stars carry very little gain.
A fixed 90% bulk rule in gain-per-added-node order retains 618 stars. Preflight
409 and six fresh loads in 410 pass all independent field, positive quadrature,
periodic, and Galerkin checks. Only 4,224 added nodes reduce worst all-model
disagreement to 1.723104796% at 120,793 nodes and 76,027 tets, still failing 1%.
All regions remain eligible. No reference tensor enters allocation. Third-level
preflight 411 and coupled solve 412 pass all checks at 129,269 nodes and 81,751
tets. Worst all-model disagreement improves to 1.549508124%, still failing 1%.
Galerkin balance is within 1.600e-9. Fourth-level preflight 413 and fresh
coupled job 414 pass at 144,106 nodes and 91,996 tets, with minimum quality
0.100152291. Worst disagreement improves to 1.415999636%, still failing 1%.
Its Galerkin balance is within 1.601e-9. Fifth-level preflight 415 is stopped
by the unchanged 8 GiB watchdog during sparse setup. No solve starts. The exact
saved mesh's fresh-process operator and AMG setup pass in 416 at 5.040 GiB,
showing that retained analysis memory inflated the combined peak. Full replay
417 saves the trial checkpoint and replaces the process before setup. All
replay gates pass at 4.590 GiB, with exact 415 mesh-array equality and the same
operator hash as independent 416. The 168,615-node, 109,059-tet trial has a
0.065557630% largest gain and minimum quality 0.100032836. Failed 415 remains
failed. Coupled job 418 completes six fresh loads and independent field,
positive quadrature, Galerkin, and all-model audits. Worst disagreement improves
to 1.325163420%, still failing 1%. Galerkin balance is within 1.006e-9. A
sixth-level trial 419 preserves geometry and parent fields at 206,426 nodes
and 135,587 tets, but its largest gain is 0.039442352%, below the unchanged
0.05% signal. No resource setup or solve starts. Close this h sequence without
convergence or native/reference promotion. Native investigation isolates a
periodic cap-star obstruction in the fixed linear spinodoid route. Local repair
controls pass, but whole-field candidates still fail. Quarter-edge insertion
reduces native sampled error from 0.004610 to 0.003089, still above 0.002.
No repaired native candidate is promoted. The curved spinodoid remains admitted.
All four cases, the 0.1 quality gate, and all 24 lattice models remain in scope.
See [physical edge normal investigation](experiments/adaptive/PERIODIC_EDGE_NORMAL_RESULTS.md) and [native cap repair investigation](experiments/adaptive/NATIVE_CAP_REPAIR_RESULTS.md).
Earlier next-step descriptions below are historical. All four cases remain
in scope, and the full objective is incomplete.


## Measured edge refinement improves lattice error but still fails, 25 October 2026

Jobs 369–374 replay all 24 reference formulas, verify exact parent provenance,
and evaluate one residual-ranked periodic edge-star Tet10 refinement on the
unchanged current native domain. At 249,998 nodes / 163,802 tets, geometry,
quality, periodicity and the 8 GiB resource limit pass. Six fresh loads and
independent field/quadrature/Galerkin checks pass. Worst all-model disagreement
improves from 1.283822914% to 1.216217882%, still failing 1%. No production
default, native candidate, reference or timing victory is promoted. Do not
repeat this h recipe automatically. See [edge residual investigation](experiments/adaptive/EDGE_RESIDUAL_RESULTS.md).
All four cases and all 24 models remain in scope. Earlier next-step descriptions
below are historical; a separately reproduced physical-side gradient control
must pass before any native evaluation change is retained.


## Fixed native and actual wall geometry inspected, 25 October 2026

Three existing surface-projection controls pass in job 366, including exact
planar area/distance and explicit rejection reporting for feet outside the box.
Read-only jobs 367 and 368 compare the saved native lattice26 and admitted
actual lattice meshes against the same ordered implicit wall. The frozen
specimen JSON and sample NPZ hashes match between routes. No mesh coordinate,
connectivity, sample or material field is regenerated or changed.

Positive 9/16-point sampling initially changes native RMS by about 10%, so
job 368 adds 36/64 points per wall face. Native full-wall RMS settles at
1.579072e-5, changing 0.00222% from 36 to 64 points. All 987,776 native wall
probes project within the box. The curve-incident wall geometry is especially
close to the analytic wall, although this does not certify local stress.

Actual wall faces without curve vertices have RMS 4.306095e-4, stable across
the denser rules, with zero failed projections in that group. The two routes'
excluded curve strips differ, so these group values are not a comparison on
identical wall patches. Actual curve-incident probes include 4,659 projections
outside the box at 64 points per face. Their raw RMS and signed integrals are
saved with explicit failures and must not be treated as fully valid bounded
wall-projection metrics. This diagnostic does not replace the original
admission method or silently reclassify the admitted actual competitor.

The native mesh already has better original sampled distance, smaller absolute
volume error, higher minimum corner quality and fewer nodes/tets than this
actual route. It still fails the retained-model stiffness target. These new
geometric observations do not prove the cause of that disagreement, justify
fitting a geometry offset, or certify a converged reference. Worker 368 takes
13.35 seconds at 0.302 GiB, solely a saved-domain diagnostic, not new workflow
timing. See [fixed-domain diagnostic](handoff/macbook/lattice-fixed-domain-comparison-20261025.json).

Next audit the original lattice reference extrapolation assumptions and their
fixed-domain evidence before choosing another mesh allocation change. Retain
all 24 model tensors and all four cases. The residual screen remains below
its threshold for another solve, and no native/reference victory is declared.
Earlier selected-next descriptions below are historical.


## Remaining hierarchical residual screen closes another solve, 25 October 2026

Four reconstruction/Schur controls pass in job 364. They check the enriched
parent gradient and scattered residual against independent quadratic/bubble
energy decomposition. Job 365 then screens six new residual directions of
the complete audited field, including the analytic bubble contribution.
The geometry, original Tet10 space, cubic degree and all physical nodes/tets
remain unchanged. This is a fixed-parent trial, not another equilibrium solve.

Largest additional trial gain is 0.017533229%, below the existing 0.05% signal
screen. Positive 64/125-point gain matrices differ by 9.703e-6 relative, well
within 1%. The trial still fails all 24-model admission at 1.119427778% worst
disagreement. Worker time is 75.50 seconds at 1.238 GiB. No further auxiliary
solve is selected from this screen, and it does not bound the full remaining
FE error or establish convergence. No reference or native candidate is promoted.

The last independently equilibrated enriched result remains 1.136813109%;
the current native workflow remains 1.283822914%. Preserve that distinction.
The next investigation should compare geometric error on the fixed native
and admitted actual domains before another allocation or solve. Keep all 24
models and all four cases. No fitted geometry offset, new sizing sweep or
change to the 1% stiffness, distance, quality, periodic or count gates follows.
See [remaining residual trial](handoff/macbook/lattice26-hierarchical-remaining-trial-20261025.json).
Earlier selected-next descriptions below are historical.


## Complete hierarchical fields audited; lattice target still fails, 25 October 2026

Six fresh zero-macro auxiliary responses finish in job 361, using the exact
current lattice26 sparse parent operator. The 6x6 Schur correction preserves
the curved geometry, all original Tet10 functions and the 189,180 nodes /
123,694 tets. Its six extra scalar mode amplitudes do not turn it into a native
Tet10 mesh workflow. Nodal displacements alone omit its analytic cubic bubbles.
Worker time is 123.82 seconds at 5.148 GiB under the unchanged 8 GiB watchdog.

Independent job 362 reconstructs the complete fields with positive 15, 64 and
125 point integration and freshly assembles all coarse equilibrium rows.
The operator SHA reproduces exactly. Maximum coarse residual is 2.057e-10,
auxiliary residual 1.775e-11, modal residual 6.730e-13 and reaction/energy
relative difference 3.650e-12. Galerkin correction-energy balance differs by
2.169e-11 relative. Positive 125-point fixed-field energy differs from the
solved tensor by at most 1.185e-7 across mixed loads. Complete cubic traces
and zero-macro response periodic jumps pass. Audit takes 58.66 seconds at
4.820 GiB. No additional equilibrium solve occurs during the audit.

The largest stiffness correction is 0.147582534%. Worst disagreement across
all 24 unchanged models drops from 1.283822914% for the native parent to
1.136813109% for the enriched field. This still fails 1%. No reference,
native candidate or comparative timing claim is promoted. The admitted
modified actual Microgen route remains ahead on this lattice accuracy gate.

Saved-array replay 363 checks local positive semidefiniteness, global gain
balance and rule agreement. Curve-vertex incident cells carry 54.81% of the
largest-direction correction in 6.10% of cells, down from the unrelaxed trial's
60.18%. Wall cells without curve vertices carry 25.43%; cells touching neither
wall nor cap carry 17.32%. Groups overlap. This maps the resolved correction,
not the remaining error, and does not justify another curve-only sizing sweep
or copying scores to a different domain. Before another native allocation
change, distinguish unresolved approximation error from geometric error with
a same-domain residual check of these complete equilibrated fields.

See [complete-field audit](handoff/macbook/lattice26-hierarchical-coupled-audit-20261025.json) and
[localization replay](handoff/macbook/lattice26-hierarchical-localization-20261025.json). Schwarz's repeatable modified-actual-route
wins stand. Gyroid and spinodoid still lack admitted actual comparators in
tested routes. All four specimens, all 24 lattice models, fixed geometry,
quality, periodic gates and 250,000-node/tet caps remain in scope.
Earlier selected-next descriptions below are historical.


## Current-domain hierarchical coupling preflight passes, 25 October 2026

Eight virtual-work/force-balance/Schur controls pass. Job 359 then verifies
consistent 15-point mode coupling against positive 125-point integration:
mode Gram and gain changes are 6.655e-5 and 4.023e-5 relative, within 1%.
The actual sparse parent has 545,136 free rows and 44,320,304 stored entries,
534,024,196 CSR bytes. Six coupling columns need 26,166,528 reduced bytes.
All saved-field residuals pass, and the exact parent operator SHA reproduces
the six-load run. Four-level AMG setup passes at 5.018 GiB with CG intercepted
before entry. No dense global mode rows are passed to AMG or matrix saved.

Parser failure 357 and key-comparison failure 358 remain frozen. The latter
compared the (shape,SHA) tuple with a SHA string; comparing its SHA member
passes without changing the exact identity or physics gates. No numerical
operator discrepancy was established. Worker 359 takes 43.46 seconds.
See [coupling preflight](handoff/macbook/lattice26-hierarchical-coupled-preflight-20261025.json).
Next permit six fresh auxiliary zero-macro response solves and a 6x6 Schur
correction, followed by independent complete-field and positive-quadrature
audits. Geometry, original Tet10 space and node/tet caps remain fixed. No
reference, native candidate or timing victory is promoted; all four cases
and all 24 model tensors remain in scope.
Earlier selected-next descriptions below are historical.


## Lattice26 same-domain residual trial resolves a correction, 25 October 2026

After three passing analytic controls, job 355 builds six global trial modes
from shared cubic edge/face bubbles on the immutable current lattice26 mesh.
All original Tet10 nodal functions and geometry maps remain available; this
trial holds the original solved coefficients fixed and optimizes only six
extra scalar mode amplitudes per load. It is not a new equilibrium solve,
full cubic nodal discretization, native mesh candidate or workflow timing.
The 411,670 support features and 7,410,060 fixed coefficient entries are
reported explicitly. No new geometric nodes/tets are created or admitted.

The largest fixed-domain energy gain is 0.117081127%, above the existing
0.05% signal screen. Positive 64/125-point gain matrices differ by 8.564e-6
relative, within the unchanged 1% gain-precision gate. The complete cubic
trace check covers 494,776 shared/periodic face records, using the six exact
P2 zeros and four complementary points; maximum jump is zero. About 60.18%
of the largest-direction correction energy lies in the 7,548 curve-vertex
incident cells, 6.10% of the mesh. This localizes this trial correction only.

Worst disagreement remains 1.167137138% over all 24 unchanged models, failing
1%. Worker time is 68.12 seconds at 1.091 GiB under the unchanged 8 GiB cap.
The failed first control 352 exposed an installed Tet4 derivative-list type;
its source is preserved, conversion is corrected, and controls 353–354 pass.
See [hierarchical trial](handoff/macbook/lattice26-hierarchical-trial-20261025.json).
Next check consistent 15-point coupling, actual reduced assembly and AMG
setup feasibility before any fresh auxiliary response solve. Do not feed
six dense global mode rows to AMG. No reference or native workflow is promoted;
all four cases and all 24 models remain in scope.
Earlier selected-next descriptions below are historical.


## Fresh lattice26 fields pass; all-model target still fails, 25 October 2026

The unchanged recorded curved Tet10 recipe on native 7bad4a45 admits 189,180
nodes / 123,694 tets, versus 199,972 / 128,300 historically. Minimum corner
quality rises from 0.149329 to 0.207354 and q_p01 from 0.293528 to 0.522231.
Sampled distance is slightly worse, 0.000343162 versus 0.000323640, within
0.002. Volume error is also nominally slightly worse; all original geometry,
positive Jacobian, facet and periodic gates pass. No all-metric win follows.

Fresh Mac workflow 348 completes all six AMG loads in 151.245 seconds including
generation/workflow checks, at 5.380 GiB RSS under the unchanged 8 GiB cap.
Independent audit 349 passes reactions, energy and periodic jumps; positive
125-point fixed-field reintegration shifts worst mixed-load energy by only
4.10e-8, below 1e-4. Assessment 351 preserves all 24 model tensors and reports
worst disagreement 1.283822914%, versus previous 1.268698522%. The 1% target
still fails. No reference is promoted, and no comparative timing ratio is
claimed. The existing admitted modified actual lattice route remains a case
where Meshers is behind on the retained-model target.

Assessment 350's raw-label identity requirement failed on numbering alone.
Its source/failure is preserved; 351 proves exact coordinates and all local
Tet10/Tri6 associations after node/whole-row canonicalization, without sorting
local shape nodes or changing any scientific gate. The pilot now uses the
active venv, records actual host/dependency/thread metadata, enforces both
250,000-node/tet caps and fails the queue stage on incomplete requested loads.
See [current lattice assessment](handoff/macbook/lattice26-current-native-assessment-20261025.json).
Next use these fully audited fields for one same-domain hierarchical residual
trial, with no sizing sweep, geometry change, new reference or native workflow
claim. Keep all four cases and all 24 models in scope; the full goal is open.
Earlier selected-next descriptions below are historical.


## Public re-entrant periodic recipe recovered, 25 October 2026

Complete region-conjunction end traces now receive a conservative native
certificate before physical seed construction and the existing periodic gate.
Red 333 fails; green 336 and independent public replay 339 pass. The unchanged
public notched recipe admits 3,896 tets / 898 nodes, qmin 0.301405 and roundoff
geometry/volume, with 54 paired nodes and 82 matching cap triangles per end.
The public arrays exactly match the earlier private recovery. Historical
failure 328 is unchanged. Full native suite 338 passes 101 tests / one explicit
research skip; Python 340 passes 211 / one optional skip. Independent audit
346 proves all four fresh fixed study meshes and the periodic column exactly
unchanged from native 02d1c618. Final arm64 image is 7bad4a45.
See [periodic region results](experiments/adaptive/PERIODIC_REGION_RESULTS.md).
No FE solve, reference promotion or comparative timing victory follows.
Next replay the recorded lattice26 curved Tet10 recipe as one geometry and
resource preflight on the new native image, before considering fresh loads.
Keep its settings, 250,000-node/tet caps, fixed specimen and all 24 models.
The full objective remains open. Earlier selected-next descriptions below
are historical.


## Re-entrant periodic failure has a controlled private recovery, 25 October 2026

The new auxiliary notched prism fails the public scalar recipe 328 at error
0.0106363 against 0.002. Its fixed volume is 0.18625 and its compiled graph has
eight planes / two region masks. Manufactured 329 proves feasibility at qmin
0.111168. Physical seed 330 has 54,534 tets / 12,264 nodes and exact geometry,
160 matching end nodes and 252 cap triangles per end, but qmin 0.0062124 and
is explicitly unqualified. Private coupled adaptation 331 of the exact seed
bits passes native final and independent geometry gates at 3,896 tets / 898
nodes, qmin 0.301405, roundoff distance/volume and 54 paired nodes / 82 matching
cap triangles. Independent investigation replay 332 passes. This does not
recover the public route or promote a reference. Native 02d1c618 and all four
fixed study meshes are unchanged. See [periodic region results](experiments/adaptive/PERIODIC_REGION_RESULTS.md).
Next certify complete region end traces before considering public dispatch.
All four study cases and all 24 lattice tensors remain in scope.
Earlier selected-next descriptions below are historical.


## Periodic prism orientation control passed, 25 October 2026

Job 326 admits the recovered prism rotated 0.37 radians about z and translated
(0.02,-0.03,0), without clipping the side walls or changing its volume/end
traces. It uses 3,879 tets / 865 nodes, qmin 0.304041, roundoff distance/volume,
46 node pairs and 68 matching cap triangles per end. Independent replay 327
passes at identical options/gates. Thirteen fixture controls pass. Native
02d1c618 and all four fixed study specimens are unchanged. No FE/reference
or speed claim follows. Next run one periodic re-entrant notched-prism
control to test nonconvex branch coverage; no sizing sweep is queued.
Earlier selected-next descriptions below are historical.


## Public periodic affine normal-jump recipe recovered, 25 October 2026

The unchanged public tapered-prism recipe now admits 3,811 tets / 844 nodes,
qmin 0.309003 and roundoff geometry/volume error, with 45 periodic node pairs
and 66 matching cap triangles per end. Complete normalized affine end-trace
multisets must match before physical cutting; the actual seed must then pass
the existing node/cap gate before coupled torus adaptation. Named fields keep
their individual trace requirement. Red 311 rejects; green 314 and independent
replay 317 pass. No old failure was reclassified. Full suite 316 passes 99 native
tests / one explicit research skip; Python 318 passes 208 / one optional skip.
Fresh controls 319–323 and independent replay 325 prove all four study meshes
and the periodic column exactly unchanged from native 99817e40. Audit 324
selected the earlier failed spinodoid CLI directory; its failure is retained,
and 325 explicitly selects the valid v2 predecessor. No FE solve, reference
promotion or comparative timing claim follows.
See [periodic affine trace results](experiments/adaptive/PERIODIC_AFFINE_TRACE_RESULTS.md).
Next check one rotated version of the recovered prism for orientation
robustness before optimization. Original four specimens and all 24 lattice
model tensors remain fixed and in scope. The full objective remains open.
Earlier selected-next descriptions below are historical.


## Verified affine curve mobility, 25 October 2026

The adapter incorrectly froze vertices on two affine planes. Rank-aware
classification now permits their valid line motion and holds endpoints fixed
by an independent third constraint. The exact saved periodic-prism seed
recovers from ten cells below q = 0.1 to none, qmin 0.304656, 3770 tets / 838 nodes,
with geometry, positive Jacobians, periodic caps and coherent facets admitted.
This is private-core research; the public prism route is still rejected.
Full native suite 297 passes 97 tests / one explicit research skip; Python 301
passes 206 / one optional skip. Independent study audit 308 proves all four
fresh fixed study meshes exactly unchanged against native 30506b13. The
periodic affine column remains admitted, now 2614 tets / 604 nodes, qmin 0.321550
and roundoff distance. No timing, FE accuracy or reference claim follows.
See [affine curve results](experiments/adaptive/AFFINE_CURVE_RESULTS.md).
The next failure is the public construction of a verified periodic seed with
permuting end traces. Keep the individual named-field periodic gate intact.
All four cases and all 24 lattice models remain in scope; the full objective
is unachieved. Earlier selected-next descriptions below are historical.


## Verified periodic image-query correction, 25 October 2026

The red/green normal-jump regression proves local out-of-cell error queries
needed the prescribed periodic image. The correction changes no field, gradient
policy or admission gate. Native suite 283 passes 95 tests / one explicit
research skip; Python suite 286 passes 206 tests / one optional skip. Fresh
controls 287–290 admit all four cases, independently replayed in 292. Gyroid,
Schwarz and spinodoid are exactly unchanged. Lattice uses 213,800 tets / 46,131
nodes instead of 218,239 / 48,077, with unchanged distance 0.001768701 and
improved reported quality quantiles. Nominal volume error is slightly worse
(0.429176% versus 0.429123%), within the same 1% gate. No timing ratio or FE
accuracy claim follows. The original periodic affine control remains exact.
See [periodic distance results](experiments/adaptive/PERIODIC_DISTANCE_RESULTS.md). Saved-seed replay 293 removes 111
false error splits, but ten bad cells remain and the public tapered-prism
recipe still rejects. All ten touch affine curves; inspect their actual
constraint freedom before another native change. Original four cases and
all 24 lattice models remain in scope; the full objective is unachieved.
Earlier selected-next descriptions below are historical.


## Periodic normal-jump failure isolated to a feasible input,25 October2026

The invariant affine fix is committed/pushed as6c843b6. The next exact tapered
prism has matching periodic scalar traces and a33.4degree wall-normal jump.
Raw adaptive272 rejects error0.0106425; mapped273 rejects0.00234375; zero-round
mapped274 rejects q0.045969. None is reclassified. A separate manufactured
20tet mesh passes all gates at q0.171170,proving geometric/quality feasibility
without claiming a native workaround. Independent replay276 passes;12fixture
controls pass. No production source,native image,FE result or reference changed.
See [periodic normal-jump investigation](experiments/adaptive/PERIODIC_SEAM_RESULTS.md). Next inspect one physical
six-plane seed and actual periodic end correspondence before a native change.
Keep the mapped recipe closed and preserve all four cases and24 lattice models.
Earlier selected-next descriptions below are historical.


## Explicit resume and verified periodic affine recovery, 25 October 2026

Kevin resumed from the pause checkpoint. Native test265 had finished; it was
not restarted. The periodic column now admits 2,663 tets / 655 nodes, qmin 0.302251
and distance 6.21e-17 at unchanged settings. Independent replay 271 reproduces
the rejected scalar corner witness and exact named-domain control, and proves
all four freshly generated study meshes exactly unchanged after relabeling.
Python/fixture/tools suite 264 passes 205 tests / 1 optional skip; full native
release suite 265 passes 94 tests. Final build 266 reproduces the measured image.
See [periodic affine results](experiments/adaptive/PERIODIC_POLYHEDRAL_RESULTS.md). No new FE/reference/timing claim follows.
The supported periodic planes must each be invariant along requested axes;
non-affine normal jumps and branch permutations remain open. Next diagnose
one actual periodic-seam feature mechanism before any further native change.
All four cases and all 24 lattice models remain in scope. Earlier pause and
selected-next descriptions below are historical.



## 25 October continuation — adaptive nonlinear failure recovered

The unchanged cells24 adaptive bore recipe now admits12955 tets and3004 nodes,
qmin0.168559,distance0.001977286,volume error0.106744%. Actual surface labels
survive splits/collapse and joint projection. No option or admission gate was
relaxed. Python/fixture suite250 passes185 tests / one optional skip; full
all-feature native suite251 passes93 tests. Core and binding Clippy pass.
All four fresh study controls253–256 admit and independent replay258 proves
exact canonical equality with preserved b2814b0e meshes and fixed fixtures.
The hollow-shell regression257 remains admitted, with fewer elements but worse
quality/distance/volume error than its prior admitted result; no dominance claim.
See [adaptive nonlinear results](experiments/adaptive/ADAPTIVE_NONLINEAR_RESULTS.md)
and [the audit](handoff/macbook/adaptive-nonlinear-change-audit-20261025.json).
No FE solve, reference promotion, model removal or new actual Microgen timing
ratio. General nonlinear Boolean/periodic feature handling, admitted actual
Gyroid/spinodoid comparators and all24 lattice convergence remain open. The full
goal stays active. Earlier selected-next descriptions below are historical.


Eligible nonadaptive Tet4 calls now retain analytic outer-max branches.
The bored-prism scalar route that failed at matched settings now admits with
145924 tets, qmin 0.121400, distance 0.00105814 and volume error 0.05022%.
Its mesh is exactly equal after relabeling to the admitted seven-field named
control. All four original study meshes also replay exactly against their
preserved baseline. Full suite 236 passes 198 tests / 1 skip; final all-feature
binding Clippy 238, Python lint and typing pass. See
[nonlinear constraint results](experiments/adaptive/NONLINEAR_CONSTRAINT_RESULTS.md).
Adaptive nonlinear postcontrol 244 still rejects at 0.003393814; no adaptive
option is bypassed. Next carry actual field identities through adaptive
splitting, collapse and projection. No FE solve, reference promotion, model
removal or new actual Microgen timing ratio. Actual Gyroid/spinodoid admission
and all 24 lattice convergence remain open. The full goal stays active.


Meshers now supports bounded nonconvex min/max combinations of convex affine
fields for nonperiodic Tet4. Shared cuts preserve one conforming union; actual
boundary facets define constraints, and smoothing accepts any complete region.
The L-prism default adaptive control admits with 7354 tets, qmin 0.269396 and
near-zero geometric error. Rotated adaptive controls pass; the rotated,
unsnapped nonadaptive recipe needs two explicit optimization rounds. Its
one-round quality failure remains rejected, and the failed reordering trial
is reverted. See [nonconvex results](experiments/adaptive/NONCONVEX_RESULTS.md).
Final native 83247d4d controls admit all four original specimens and exactly
match prior meshes after canonical relabeling, with identical fixture hashes
and configurations. Full Python suite 223 passes 190 tests / 1 skip; full native
suite 202 passes 91 tests; final lint-clean native check 224 passes 62 unit/intersection tests.
All-target core and all-feature binding Clippy pass. No FE solve, reference
promotion or new actual Microgen timing ratio. The full goal stays active.
Next investigate one nonlinear wall/cap intersection failure relevant to the
lattice boundary limits, keeping all four original controls. General nonlinear
and periodic feature handling, admitted actual Gyroid/spinodoid comparators
and all 24 lattice reference convergence remain open.


Adaptive convex polyhedral Tet4 generation now preserves exact plane constraints.
The rotated control exposed a nonperiodic seed-welding panic; preserving vertex
identities and checked collapse of numerical corner copies fixes it. The final
adaptive octahedron admits with 2565 tets, qmin 0.300440 and distance 8.97e-17.
All four fresh original study controls admit and exactly match their previous
meshes after canonical relabeling. Final suite 176 pass / 1 skip, core 49 pass and
Clippy pass. See [adaptive polyhedral results](experiments/adaptive/ADAPTIVE_POLYHEDRAL_RESULTS.md).
General nonconvex/nonlinear/periodic feature handling, actual Gyroid/spinodoid
comparators and lattice reference convergence remain open. No FE solve,
reference promotion or new timing ratio; the full goal stays active.


Meshers now preserves supported convex polyhedral features automatically in
default nonperiodic, nonadaptive Tet4 generation. The compiler recovers exact
affine half-spaces and reuses the constrained cutter. The formerly rejected
scalar octahedron passes with 13,488 tets, qmin0.317878 and distance8.28e-17;
rotated/translated and option-gating regressions pass. Final suite171pass1skip.
Fresh controls for all four original specimens pass, with exact geometric
mesh equality to their baselines; lattice differs only in labels/order.
See [polyhedral results](experiments/adaptive/POLYHEDRAL_RESULTS.md).
Adaptive/nonlinear/periodic crease handling remains open. The original actual
comparator and lattice reference gaps persist; the complete goal is unachieved.

The failure-first continuation now has an eight-case saved-artifact replay:
all four fixed specimens plus torus, hollow shell, separated spheres and a
sharp octahedron. New curved controls required independent volume admission;
one analytic volume/area-derived tolerance each resolves the torus and sphere
deficits. The fused octahedron still misses reverse distance, while the exact
named-plane route passes and has a production corner regression. No native
algorithm change or broader comparative win is claimed. See
[topology coverage](experiments/adaptive/TOPOLOGY_COVERAGE.md).
Next investigate the remaining fused-field feature miss using its exact-domain
named control; preserve all four original cases. Nine missing historical
source copies were restored by recorded hash, and new runs archive Python
inputs. The complete objective remains open.

Current priority is failure coverage, at Kevin's request, before optimization.
The resumed session recovered two omitted historical Tet4 recipes. Lattice
and spinodoid each pass three fresh native Mac generations, including two
required-admission runs and independent replay. All four fixed specimens
therefore have verified Tet4 routes. No native code or gate changed.
See [Tet4 coverage correction](experiments/adaptive/TET4_COVERAGE.md).
This corrects the missing-native-admission assessment below. General input
coverage, admitted actual Gyroid/spinodoid comparators and lattice reference
convergence remain unresolved. The full goal is unachieved.

The mesh-only comparison is complete across all four fixed cases. See
[Mesh comparison results](/Users/kmarchais/dev/meshers/experiments/adaptive/MESH_COMPARISON_RESULTS.md).
Matched Tet10 shows 45.31 times faster generation on Schwarz with 4.97 times
fewer tetrahedra, and 10.90 times faster on lattice uniform 25 with 3.25 times
fewer. Lattice p95 quality is slightly lower despite better minimum, low-tail
and median quality. Tet4 remains the primary comparison and has trade-offs;
native lattice/spinodoid admission and admitted actual Gyroid/spinodoid controls
remain gaps. All 138 saved meshes were replayed. No FE reference or full-goal
completion is promoted. Next follow the mesh-only quality/fidelity evidence,
rather than the historical FE investigation sequence.

The complete parent25 correction map passes both positive rules and independent
saved-array replay. Wall cells without curve vertices carry74.60% of its
normalized energy; curve-vertex cells16.37%. Stop curve-only allocation and
keep closed cap-length/energy-sizing recipes closed. Next review native periodic
normal/feature handling against existing branch evidence before selecting a
controlled regression. The audited field still fails all24 accuracy at1.711485%.
No reference, native candidate or timing ratio is promoted. All four cases and
the original unachieved objective remain in scope; work resumed explicitly.


Resumed at the user's explicit instruction. Parent25's six fresh coupled
responses pass complete-field, positive-rule and freshly assembled equilibrium
audits, but all24-model disagreement remains1.711485%, failing1%. No reference
or native candidate is promoted. Next localize the resolved correction before
choosing a native boundary change. Schwarz repeats stand; Gyroid/spinodoid
actual-comparator gaps and lattice convergence remain unresolved. All four
cases, original gates and the unachieved objective stay in scope.


Paused at the user's request after job90 finished successfully. No successor
is queued. The parent25 coupling/rule preflight passes; reduced assembly,
AMG setup and auxiliary solves have not run. Resume only on explicit user
instruction, beginning with the pause section in `experiments/adaptive/PLAN.md`
and `handoff/macbook/pause-20261024.json`. The four-case goal is unachieved.


Parent25's fixed-domain six-mode h/p trial reduces energy by 0.374103%, but
still misses all-24-model accuracy at 1.859718%. Both positive rules agree;
no reference or native candidate is promoted. Most trial gain lies near the
wall/caps outside curve-vertex cells, so stop the curve-only allocation.
Cell batching matches the old integrator and makes this diagnostic practical.
Next inspect consistent-coupling and memory feasibility before choosing any
six-mode re-equilibration solve. Schwarz repeats stand; Gyroid/spinodoid actual
comparator gaps and lattice convergence remain unresolved. The full four-case
goal is active and unachieved.


Preferred-only star refinement improves the audited native26 control to
1.183637% across all 24 models, still failing 1%. Default refinement and all
gates are unchanged. A coarser parent25 budget preflight stops at round13
because two Newton projection feet leave the box. Valid constrained target
feet are demonstrated on the same probes, but the mesh remains rejected.
Next validate an opt-in constrained-projector fallback, preserving historical
results and thresholds, then re-audit before any refinement or FE. The fresh
modified actual lattice route still passes at 0.994671%; no lattice timing
win follows. Schwarz repeats stand; Gyroid and spinodoid gaps remain. The
full four-case goal is active and unachieved.

One fresh Mac actual lattice workflow now passes all 24 models at 0.994671%.
The route declares tuned VTK51/two-pass MMG, a one-vertex fixed-boundary repair,
manual straight Tet10/direct fedoo and star refinement. Every added stage is
charged: 314.03 seconds including diagnostics. Independent fields and positive
quadrature pass; this is one diagnostic workflow, without a timing ratio.
The native star-refined control still fails at 1.230455%. Next inspect and
preflight preferred-star-only refinement under the unchanged native26 domain,
250000-node/tet caps and all scientific gates. Schwarz repeats stand; Gyroid
and spinodoid comparator gaps remain. The four-case goal is unachieved.

A star-refined actual Microgen-derived saved-mesh lattice control now passes
all 24 models at 0.988498%. The native star-refined control still fails at
1.230455%. Both pass full field, periodic and positive-quadrature audits.
The passing actual route includes explicit extra refinement; its frozen G3
seed is not fresh Mac generation or a workflow timing. Prepare that fresh
actual generation next, reading geometry before any solve. Parent23 reference
convergence remains unresolved. Repeated Mac Schwarz wins stand; Gyroid and
spinodoid comparator gaps remain. The full four-case goal is unachieved.

Ranked nested lattice correction now passes complete-field, positive-rule,
periodic and freshly reassembled equilibrium audits. All twelve fresh
auxiliary responses fit at 5.769 GiB. The largest fixed-domain stiffness
decrease is 0.093128%; worst disagreement across all 24 models improves to
1.741320%, still failing 1%. No reference or native candidate is promoted.
Inspect remaining correction and geometry evidence before another native
change. Gyroid and spinodoid comparator gaps remain; repeated Mac Schwarz
wins stand. All four cases remain active and the full goal is unachieved.

Correction-energy localization now distinguishes complete curve-vertex
neighborhoods from literal curve edges. One measured-score-per-node nested
preflight selects parents carrying 88.8088% of the resolved six-load score
and passes geometry at 249,999 nodes. Actual reduced assembly fits at
5.394 GiB, without a solve. Independent saved-array replay passes. Next is
coupling/AMG resource validation; no new accuracy or reference is claimed.
The last coupled 24-model disagreement remains 1.826184%.

The coupled parent23 reference investigation now passes complete-field,
positive-quadrature, periodic, pin and freshly reassembled equilibrium audits.
Its six fresh auxiliary responses fit the unchanged 8 GiB resource ceiling.
Worst disagreement remains 1.826184% across all 24 models, so no reference or
native candidate is promoted. Next is correction-energy localization before
any new mesher change. The current Mac checkpoint records all four cases.

The latest lattice reference diagnostic proves that saved nested and cubic
fields share the original parent23 physical domain, then measures a resolved
0.2593816% trial-energy gain from their six periodic difference fields.
Independent constrained-basis replay passes, but the trial still misses the
24-model target at 1.899837%. No reference or native candidate is promoted.
Next is one quadrature/resource preflight for the coupled space, before any
new equilibrium solve. See the current Mac checkpoint for all four cases.

Five fresh Mac Schwarz pairs now pass the supplementary 0.1% full-tensor
target, with one excluded warmup and alternating arm order. Median charged
workflow time is 131.85894 seconds for modified actual Microgen and 19.14923
seconds for Meshers. All five paired ratios favor Meshers, from 6.7632 to
6.9203. The explicit one-pass route, charged SciPy repair and manual straight
Tet10 remain declared. Independent replay verifies all twelve mesh/field sets
and 96 model comparisons. The other three fixtures remain unresolved.

One fresh Mac Schwarz pair now passes the supplementary 0.1% full-tensor
comparison. The actual Microgen route includes explicit one-pass overrides,
a charged one-vertex SciPy repair, and manual straight Tet10. Workflow times
are 132.63564 versus 19.27508 seconds, a single-pair 6.8812 ratio. Independent
field/tensor/artifact replay and positive-rule reintegration pass. This is
one pair; repeated Mac timings are next. No four-case victory is claimed.

Mac jobs 24–28 locate an independent spinodoid wall-geometry failure and
reject simple wall-node projection. The first fresh Schwarz pair rejects
its Microgen parent at quality 0.0956791; no FE or timing claim follows.
A one-vertex fixed-boundary repair preflight then passes all geometry gates.
Job 29 tests this explicitly added, charged repair in a fresh paired workflow.
The current Mac checkpoint records the evidence and remaining four-case gaps.

Mac jobs 21–23 complete the selected two-orbit spinodoid cap preflight.
The initial control failure is preserved; corrected derivative and known-repair
controls pass. Eight starts miss cap quality 0.1, reaching about 0.00504.
Independent saved-coordinate replay confirms the rejection. No move, MMG,
FE, admitted comparator, or reference follows. Stop this search recipe and
localize the separate spinodoid distance failure before more cap repair.
See the latest Mac checkpoint for the active continuation.

Latest Mac comparison work is in
[fixed-spinodoid extraction and cap feasibility](experiments/adaptive/MICROGEN_RESULTS.md).
The exact frozen coefficients now survive cross-platform FFT roundoff.
Paired preclipping removes 282 spinodoid source face conflicts, but the
coherent seed still fails geometry and retained-cap feasibility. A checked
cap-only obstruction rules out one fixed-neighbor periodic vertex move;
another free cap vertex participates, so coupled motion remains untested.
No MMG or FE follows. Twenty-nine focused controls and expanded provenance
verification pass. The next selected preflight is one closed active cap patch.

Native MacBook continuation is recorded in the
[24 October Mac checkpoint](handoff/macbook/CHECKPOINT-20261024.md).
Native setup, 54 focused controls, twelve solver-control loads, the second
nested preflight, and reduced assembly have passed. Level2 has 218,135 nodes
and 141,144 tets. Its six fresh loads and independent audits pass, with a
0.328272% maximum additional FE correction and 2.135117% worst retained-model
disagreement. No reference is promoted. The longest-edge lattice26 refinement
also passes audits but misses at 1.256515%, despite reaching 250,000 nodes.
That recipe is stopped. Verified Apple Accelerate single-thread controls
now pass 55 tests and twelve loads; earlier Mac numerical timings are excluded
from one-thread comparison baselines. Curve-priority refinement also passes
six fresh loads and all audits, but misses at 1.255625%. Stop these one-stage
bisection recipes. Independent checkpoint verification checks all 3,509
transferred files, 150 input hashes, and 72 full model comparisons. Return
to missing admitted actual-comparator coverage before another lattice recipe.

Updated 24 October 2026. Covers the recorded project history on
`handoff/g3-microgen-20261017` through commit `9a20151`, plus completed
audit 473 and nested-refinement experiments 475–476 below. The history starts with the initial
release candidate on 19 September and includes the resumed 20 October geometry,
reference, and fixed-field sensitivity studies.
This is a consolidated index and status ledger, not a replacement for frozen
experiment records. Historical results on other machines remain historical.
Work in separate support repositories is included only where this checkout
records its integration or evidence.

## Reading this tracker

- Sections 1–3 cover the product foundation, solver verification, and early
  adaptive meshing work before the G3 migration.
- Sections 4–10 cover the fixed comparison protocol, implemented changes,
  experiments, failed approaches, and reference uncertainty.
- Section 11 is the current four-fixture scoreboard.
- Sections 12–14 index commits, validation evidence, and preserved artifacts.
- Section 15 records unfinished work and the selected next experiment.

Completed means that the cited implementation or measurement exists. It does
not mean that the four-fixture objective has been achieved. Planned and draft
work has no measured result until its bounded run and independent checks finish.

Cap-topology preflights 345–349 find local alternatives for all 24 obstructing
Gyroid caps, with two stars needing larger cavities. Jobs 350–354 then complete
a rank-guarded replacement using full one-rings. The new seed preserves its
geometry and periodic caps, but still fails actual volume quality. A final
admitted MMG competitor has not yet been established by this preparation.

## Current conclusion

Kevin requested primary continuation on the MacBook Pro on 24 October 2026,
with Sol 6.1 medium. Read the [MacBook handoff](handoff/macbook/HANDOFF.md)
for current state and selected continuation. G3's last job, 476, has finished.
The objective remains unachieved. Work resumed on 23 October after the earlier
pause. Job 473 finishes the saved lattice26 field and quadrature audits.
They pass, while worst retained-model disagreement remains 1.2686985%.
Jobs 475–476 establish a controlled nested FE correction on a fixed curved
domain; reference convergence remains unresolved.

Meshers has a Rust/Python implicit mesher, adaptive surface and volume paths,
curved Tet10 output, direct fedoo conversion, and independently checked
elasticity experiments. Several correctness and geometry-admission defects
have been fixed and regression-tested.

The current research objective is **not achieved**: demonstrate that Meshers
beats the actual Microgen+MMG workflow at matched accuracy on Gyroid sheet,
Schwarz P solid, a fixed spinodoid, and a fixed smooth-union lattice.

- Schwarz has one demonstrated, repeated matched-target workflow advantage
  against an explicitly modified actual one-pass Microgen route.
  Five new interleaved stronger Tet10 pairs also pass the supplementary 0.1%
  target with common AMG. Their paired workflow ratios are 7.3858–7.5401.
  This remains an explicit one-pass/manual-Tet10 comparison, not defaults.
- Gyroid and spinodoid have admitted Meshers candidates, but no admitted actual
  Microgen/MMG comparator in the tested routes.
- Lattice stiffness accuracy remains unresolved across all retained reference
  models. Neither the best tested Meshers candidate nor the strongest tested
  Microgen-derived Tet10 control passes the fixed 1% target against every model.
- The latest eight-round curve-local sizing pilot fails that target at
  1.439306%. Its results and audits
  are preserved rather than promoted as an improvement.
- Three neighboring cubic reference levels fail the declared common-error
  extrapolation screen. No new reference replaces the retained model set.
- The quartic curve-edge diagnostic resolves a trial-energy reduction of only
  0.648 ppm. Allowing cubic curve-corner relaxation raises the incremental
  quartic signal to only 0.676219 ppm against the 500 ppm screen. This does not
  justify a native sizing change or establish a globally equilibrated tensor.
- Computing stopped for the requested 19 October shutdown. Work resumed on
  20 October after checking the checkout, memory, and idle shared queue.
  The eight-round curve control improves size realization and single-run cost,
  but its audited fresh full-tensor result still fails the stiffness target.

Latest continuation record: [20 October checkpoint](handoff/g3/CHECKPOINT-20261020.md).
The [19 October shutdown record](handoff/g3/SHUTDOWN-20261019.md) remains historical.

### Work status at this update

| State | Work |
|---|---|
| Implemented and validated | Direct fedoo Tet4/Tet10 conversion, optional spatial sizing, periodic projection and quality fixes, partial curved settling, exact-tie derivative correction |
| Measured with limited scope | Repeated Schwarz advantage against the labeled one-pass Microgen variant; single-run recovery and geometry improvements |
| Fresh stronger-comparator pilot | Common-AMG Schwarz Tet10 pair, 409.87 versus 55.10 seconds; both pass supplementary 0.1%, independent fields and positive-rule audit pass |
| Completed stronger-comparator repetitions | Job 333, five measured pairs after warmup; both arms pass 0.1%, medians 409.28 versus 55.19 seconds; independent verifier passes |
| Reporting fix validated | Explicit standard Job completion retains values; red/green success/failure tests and queued checks 334–335 pass; old warning remains frozen |
| Retired Gyroid repair source | The one-pass FEM source has 146 same-orientation shared faces, already present before repair; no further repair or FE on it |
| Necessary validity guard | Future common admission also requires coherent oriented facets; numerical thresholds are unchanged and all frozen rows remain unchanged |
| Comparison integrity checked | All 24 archived Schwarz arm meshes match coherent audited controls; six manual Tet10 lifts preserve parent corners; earlier timing conclusions remain supported |
| Completed negative studies | Default Microgen admissions, unsuccessful sizing and repair recipes, volume-offset rationale, three-level cubic extrapolation, edge-only quartic signal |
| Completed negative diagnostic | Coupled cubic curve-corner relaxation reproduces original edge gains exactly but fails the useful-signal screen; polynomial-mode route closed |
| Still open | Four-fixture matched-accuracy comparison, lattice reference uncertainty, admitted actual MMG controls for Gyroid and spinodoid |
| Completed extraction preflight | Gyroid31 five- and six-tet preclipping remove all 288 oriented-face conflicts, but each introduces 24 cap faces with quality ceiling below 0.1; no MMG or FE follows |
| Completed cap preparation | Twenty safe cones and four full rings replace all 24 blocked Gyroid caps; points, wall, volume, incidence, and periodicity pass; actual volume quality still requires remeshing |
| Actual prepared-seed MMG replays | Jobs 355–357 hit intermediate caps; 359 fits and passes geometry/periodicity except quality 0.047711; 104 bad cells need diagnosis; no admitted comparator or FE solve |
| Incomplete valid-source repair | Job 366 reduces 104 bad tets to 16 with 64 fixed-boundary moves; quality 0.051264 still rejects; 368 finds a single-vertex obstruction that does not yet exclude coupled moves |
| Completed Gyroid obstruction checks | Boundary-only witness 370 excludes interior-only repair at fixed connectivity; eight local flips fail, complete existing-vertex star solve rejects, and one proposed Steiner point lies outside the cavity |
| Completed lattice feasibility check | Job 377 finds a 0.9019313% necessary lower bound and a trial matrix within 0.9084715% of all 24 models; the 1% target is algebraically feasible, not a validated physical reference |
| Last completed workflow | Job 379, fresh curved lattice26, six AMG loads, 199,972 nodes and 128,300 tets, 736.51 seconds; preliminary worst disagreement 1.2686985%, independent field/quadrature audit pending |
| Historical pause | Work paused on 20 October; resumed on 23 October; user requests MacBook primary continuation on 24 October |
| Latest independent audit | Job 473 passes lattice26 fields and positive quadrature; its 1.2686985% retained-model miss remains |
| Latest controlled refinement | Jobs 475–476 preserve parent23 geometry, pass controls and audits, and measure a 0.499513% nested FE correction; no reference promoted |

### Final work on 20 October, jobs 366–379

The earlier sections retain today's extraction, cap preparation, Schwarz
repetitions, reference sensitivity, and rejected enrichment work. The final
sequence adds the following results without changing native code, specimen
definitions, references, or admission thresholds.

- Fixed-boundary repair 366 accepts 64 moves and reduces 104 bad Gyroid tets
  to 16. Quality 0.05126434 still fails. A second sweep accepts no moves.
- Full-star check 368 finds a fixed-neighbor obstruction. Boundary-only check
  370 strengthens it to constraints unaffected by other interior vertices.
  Its tangent upper bound is -0.000976627527263, independently reproduced
  with extended-precision arithmetic. This blocks this fixed-connectivity
  repair, not arbitrary remeshing.
- All eight 2-to-3 cavities in 372 fail. Complete-star trial 374 uses 13
  existing vertices, 20 original tets, and 477 quality-filtered candidates.
  The bounded binary solve reports infeasible. Trial 375 rejects one analytic
  cap-optimal Steiner point because it is outside the original star.
  No rejected source proceeds to FE. Larger cavities and other inserted-point
  constructions remain untested.
- Lattice feasibility check 377 retains all 24 stiffness models. The necessary
  pairwise lower bound is 0.9019312871%. A harmonic trial matrix has worst
  disagreement 0.9084715371%. This proves algebraic feasibility of the fixed
  1% condition, not mesh realizability or continuum accuracy. Six queued
  controls pass.
- Geometry preflight 378 admits the current-native lattice26 configuration.
  Fresh workflow 379 regenerates it and completes all six engineering-strain
  loads using curved Tet10, 15-point quadrature, and cached elastic-block AMG.

Job 379 finishes successfully at 21:41 Europe/Paris. Its final record is
`experiments/adaptive/runs/g3-lattice-adapt26-tol0028-amg-solve-v1/result.json`.
The worker exits zero with no limit failure.

| Final lattice26 measurement | Result |
|---|---:|
| Nodes / tetrahedra | 199,972 / 128,300 |
| Independent sampled two-sided maximum distance | 0.0003236401 |
| Relative volume error against frozen estimate | -0.03139036% |
| Minimum corner MMG quality | 0.14932885 |
| Minimum sampled curved Jacobian | 2.5104482e-9 |
| Periodic nodes, cap triangles, oriented facets | Pass |
| Generation / six-load solve | 9.03 / 716.94 seconds |
| Charged total workflow / bounded worker wall | 736.51 / 738.82 seconds |
| Worker peak RSS | 5,781,733,376 bytes, below 8 GiB |
| CG iterations per load | 108, 108, 108, 109, 108, 110 |
| Largest reported relative residual | 9.41e-13 |
| Preliminary worst disagreement across all 24 models | 1.2686985%, fails 1% |

The last comparison is a lightweight read-only generalized-eigenvalue check
of the saved symmetric tensor against each unchanged model, using
`max(abs(eigvalsh(C_candidate - C_model, C_model)))`. Its worst model is
index 0 of `g3-lattice-reference-sensitivity-v2.json`. Three models exceed
1%. It is not an independent displacement, reaction, Hill-Mandel, or
higher-positive-quadrature audit. The run's `qualified=true` means geometry
admission, while its final status remains `solved-needs-reference`.
Saved meshes and all six fields stay local and ignored.

This candidate does not resolve lattice accuracy. Its binary differs from
the archived lattice24 binary, and the old solve used a different solver
configuration. Neither a controlled refinement gain nor a matched-accuracy
speedup follows from their raw timings. No four-fixture victory is claimed.

The corner-coupling extension includes five algebraic tests, manufactured
periodic controls, and the frozen-domain lattice comparison. Job 330 measures
594 displacement DOFs on 1,408 tets in 263.45 seconds at 0.6863 GiB. Precision
passes, but useful signal fails. No native API or benchmark reference changes.

### Resumed work, 23–24 October, jobs 473–476

Job 473 independently audits the lattice26 fields. Reaction error is
1.0611e-11, Hill-Mandel error 4.1723e-12, and periodic jumps zero. Positive
125-point reintegration shifts energy by 4.0667e-8; final quadrature change
is 9.9686e-9. The 117.41-second audit passes at 0.7454 GiB on G3. These
checks do not explain the candidate's 1.2686985% miss.

The smallest archived admitted parent has 43,844 tets, so uniform eight-child
refinement exceeds the fixed cap. Job 475 uses conforming periodic edge-star
bisection with exact parent quadratic map restriction and six-field
prolongation. Ten controls and independent admission pass. Parent23 grows
from 68,095 nodes and 43,844 tets to 123,514 nodes and 76,947 tets.
Map/field errors are 3.33e-16, volume change 1.40e-14, sampled distance
0.0006492808, and minimum quality 0.10385108. Worker is 14.82 seconds.

Job 476 completes six fresh refined-space loads in 245.83 seconds, worker
247.63 seconds at 3,304,448,000 bytes peak RSS. Independent fields and positive
quadrature pass. All generalized parent-to-refined corrections are negative,
largest magnitude 0.499513%. The full Galerkin energy balance agrees within
6.6756e-8 relative. Worst disagreement across all 24 models is 2.4480846%,
still above 1%. This coarser fixed-domain diagnostic does not replace the
latest best native candidate or establish continuum convergence.

Records under `experiments/adaptive/runs` are
`g3-lattice-adapt26-field-audit-v1`, `g3-lattice-nested23-preflight-v1`, and
`g3-lattice-nested23-level1-v1`. Maintained code is `nested_refinement.py`,
`nested_preflight.py`, and `nested_solve.py`. Selected Mac continuation is
one second-level nested preflight, read before any further solve.

## 1. Product foundation and release preparation

### Initial Rust/Python implementation

The initial release-candidate lineage provides:

- CPU float64 tetrahedral meshing of bounded implicit solids and TPMS bands.
- Built-in Gyroid, Schwarz P, and Schwarz D fields, custom pointwise fields,
  callable grading, thickness functions, and blends.
- Automatically compiled supported NumPy expressions, automatic derivatives,
  cached native evaluation, and bounded NumPy callbacks for unsupported graphs.
- Owned point, tetrahedron, boundary, and periodic-pair arrays in Python.
- Named multi-constraint intersections, retained surface identities, feature
  edges, coordinate maps, and explicit rigid periodic transforms.
- Quality optimization, minimum-quality admission, topology validation,
  boundary-link checks, snapping rollback, cancellation, and preserved callback
  exceptions.
- VTKHDF export, including boundary and periodic metadata. Quadratic triangle
  and tetrahedron export was added during the adaptive work.

The core is independent of Python, CUDA, OpenCascade, and HDF5. Its code forbids
unsafe Rust. Python bindings and the native field compiler are separate layers.
Some intersection and surface features have separate build-feature limits.
The constrained path does not imply general multithreaded support for every API.

Release preparation added wheel and source-distribution checks, platform/Python
test matrices, installed-wheel tests, executable examples, API documentation,
coverage workflows, performance comparisons, and PyPI Trusted Publishing
configuration. Workflow configuration alone does not prove a particular remote
run succeeded or that every branch change has been published.

Evidence: [changelog](CHANGELOG.md), [project README](README.md),
[core README](crates/meshers-core/README.md),
[Python README](crates/meshers-python/README.md),
[release checks](docs/site/release.md), [testing](docs/site/testing.md), and
the workflows in `.github/workflows`.

### Early integration and performance fixes

The 3 October work reduced scalar-evaluation cost, removed unnecessary gradient
work and mutex locking, reduced shared signal-counter overhead, added bulk
compiled point evaluation without Python callbacks, and directly lowered scalar
min/max and squares into JIT instructions.

It also normalized optimizer penalties across physical scales, protected sampled
interior components during snapping, bounded interruption-test workloads,
combined optional surface and wide-constraint APIs, allowed preview surfaces to
skip distance sampling, and identified the integration API independently of the
wheel version. The default constraint width remains eight; the optional wide
build supports 32.

Evidence: commits `e645847` through `2ca2b35` and
[integration build notes](docs/implicit-microgen-integration.md). The matching
Microgen port and older research worktree are separate work, not proof that the
current default Microgen pipeline uses Meshers.

## 2. Solver and periodic integration verification

### Manufactured elasticity solutions

An independent P1 elasticity assembly tested a cube, cylindrical Gyroid sheet,
and toroidal Gyroid sheet at three resolutions each. Affine patches checked
assembly and boundary handling. The manufactured displacement study observed
approximately second-order displacement and first-order energy convergence.
All physical tetrahedral volumes were positive and explicit solver residual
checks passed.

Fedoo 0.8.4 repeated the same nine mesh solves. Error norms agreed with fresh
independent assemblies within 5e-13 relative. This validates assembly, loading,
and boundary integration, but both studies share the SciPy CG backend.

These historical measurements were on an Apple M4, not G3 timing baselines.
They do not establish periodic or nonlinear correctness, continuum geometry
accuracy, or parity with Microgen/MMG.

Evidence: [solver verification](docs/site/solver-verification.md),
`demos/solver/convergence.py`, `demos/solver/fedoo_convergence.py`,
`docs/solver-convergence-results.json`, and `docs/fedoo-convergence-results.json`.

### Periodic homogenization

Fedoo 1.0.0b2 was checked on a solid cell and four Gyroid refinements using six
engineering-strain loads and the automatic stress-controlled route. All 60
systems converged. Checks included opposite nodes and cap triangles, periodic
displacement jumps, force balance, mixed-load Hill-Mandel energy, averaged
stress versus reactions, symmetry, positive stiffness eigenvalues, and the
solid-material limit. Loading-path disagreements were below 4e-12 relative.

This established translational periodic integration, not rotational periodic
homogenization or continuum certification. Full-cell volume, including void,
normalizes effective stiffness.

Evidence: [periodic homogenization](docs/site/periodic-homogenization.md),
`demos/solver/fedoo_homogenization.py`, and
`docs/fedoo-homogenization-results.json`.

## 3. Adaptive and quadratic meshing, 11–15 October

| Workstream | Implemented or measured | Remaining qualification |
|---|---|---|
| Geometry-preserving optimization | Preserve implicit tolerance during quality moves; prototype quadratic interpolation and nested interior-refinement objectives | Prototype objectives are not general accuracy guarantees |
| Adaptive surfaces | Curvature sizing, split/collapse/flip, tangential relaxation, exact wall and box-curve projection | Thin features and sharp creases remain difficult |
| Adaptive tetrahedra | Size-driven split/collapse, edge removal, swaps, relocation, curvature/thickness targets | Adaptation is largely serial |
| Periodic adaptation | Torus/minimal-image treatment with internal box-plane constraints | Per-specimen periodic admission is still required |
| Curved Tet10 | Exact wall mid-edge projection with sampled Jacobian checks and validity fallbacks | No proof of positive Jacobians everywhere |
| Balanced surfaces | Offset vertices to reduce two-sided chord error and triangle count | Display/printing mode, not on-surface curved FEA geometry |
| Fast surfaces | Collapse-focused coarsening, error-triggered splits, exact quadratic triangle-error estimate, local balancing repair | Slivers are allowed; sphere one-thread speed remains worse |
| Extraction and diagnostics | Parallel crossings/slabs, deterministic numbering, incremental error scans, sorted manifold validation, phase timings | Coarser initial extraction remains open |

Historical independent surface scorecards measured forward and reverse distance,
volume, shape, and angles. At similar sampled distance, adaptive surfaces used
39% fewer sphere triangles, 66% fewer ellipsoid triangles, and 43% fewer Gyroid
triangles than the tested Flying Edges controls. Balancing reduced counts further.

The controlled OMEN timing study used a warm-up and five interleaved repetitions
at one and 24 threads. Fast balanced output used about 2–5 times fewer triangles
than the compared Flying Edges grids. It was faster on Gyroid, Schwarz P, and
ellipsoid, but about 1.2 times slower on sphere at one thread. Tet adaptation
barely scaled with threads. These are surface-route results, not G3 workflow
comparisons.

The early MMG volume comparison remeshed a **Meshers seed**. It was useful as an
adapter comparison but was not the actual default Microgen pipeline. Early FEA
experiments also showed much lower stiffness error per node for curved Tet10,
while balanced linear and curvature-only sizing did not establish a general
linear-FEA improvement.

Evidence: [adaptive results](experiments/adaptive/RESULTS.md), historical sections
of [adaptive plan](experiments/adaptive/PLAN.md), `scorecard.py`, `timing.py`,
`timing_job.sh`, `timing-threads1.json`, and `timing-threads24.json`, all under
`experiments/adaptive`.

## 4. G3 migration and comparison rules

On 17 October, work moved from OMEN to the dedicated Linux G3. Commit `2e35c7e`
preserves the handoff. Linux dependencies and native code were rebuilt rather
than copying Windows binaries or virtual environments. Small historical tensor
records were restored; large original OMEN meshes and fields were not copied.

The sibling support checkouts are:

- `/home/kmarchais/dev/periodic-microgen` at
  `2b1ea80b36d4ba03842a7d9087638abf3868fe3e`.
- `/home/kmarchais/dev/microgen-benchmark` at
  `602dc900cd028d205b36e29dcdb5b4a9cf4f581a`.

### Fixed specimens

All four occupy `[-0.5, 0.5]^3` with isotropic `E=2.5`, `nu=0.25` and six
engineering-strain loads, including shear.

| Specimen | Definition retained throughout the comparison |
|---|---|
| Gyroid sheet | Raw field band `[-0.55, 0.55]`; Microgen offset 1.1 |
| Schwarz P solid | Negative sublevel of `cos(2*pi*x)+cos(2*pi*y)+cos(2*pi*z)`; not a Schwarz sheet |
| Spinodoid | Fixed FFT16 realization, seed 20261010, `k0=4*pi`, bandwidth 3, retained-mode threshold 0.05, translated into the comparison box |
| Smooth lattice | Eight radius-0.12 origin-to-corner capsules, sequential cubic smooth union with `k=0.08` and frozen corner order |

Changing VTK sampling does not change the spinodoid realization. Reordering the
nonassociative smooth unions or adding neighboring struts would change the
lattice specimen and is not an admissible shortcut.

### Admission and execution

- Independently sampled two-sided boundary distance at most 0.002.
- Absolute relative volume error at most 1%.
- Minimum MMG corner quality at least 0.1, with separate positive sampled curved
  Jacobians for quadratic meshes.
- Matching periodic nodes and cap triangles, with no unresolved projections
  promoted as reliable complete distance maxima.
- Since 20 October, newly measured outputs also require coherent oriented
  facet incidence and matching extracted/saved boundaries. Numerical thresholds
  and frozen flags stay unchanged; all prior Schwarz arm replays pass.
- Full mixed-load stiffness disagreement at most 1% against justified retained
  references. The supplementary Schwarz 0.1% target is separate.
- Heavy work through the shared `gpu` pueue group, `Job` progress, and
  `qrun --check-start`, one bounded experiment read before selecting its successor.
- About 8 GiB per worker, one numerical thread, two build threads, explicit
  timeouts, and staged mesh/solve caps. The final tet cap stays 250,000. A later
  measured-safe solve stage allows 250,000 nodes without reclassifying earlier
  stage-cap failures. Separate seed-only preflights do not relax final caps.
- Same-host timing only. Live workflows charge coarse solves, size construction,
  meshing, audits, and final solves when used. Saved-mesh replay is labeled.

The stiffness score is the maximum absolute generalized eigenvalue of
`C_candidate - C_reference` relative to `C_reference`, not six diagonal errors.
Distance, Jacobian, reference-model, and QMC checks are empirical, not continuum
or local-stress certificates.

Evidence: [original G3 handoff](handoff/g3/HANDOFF.md),
[actual workflow definition](experiments/adaptive/MICROGEN.md), and
[decision trail](.audit/g3-microgen.tsv).

## 5. Actual Microgen workflow and stronger controls

### Defaults measured, not substituted

The first G3 comparison used Microgen 2.0.0b1's actual volumetric VTK extraction
and `remesh_keeping_boundaries_for_fem`, backed by mmgpy 0.17.1. Actual calls and
source hashes were recorded. The default MMG route first uses `nofem=True`, then
`ls=True` and `nr=True`, with unset sizing options and Microgen's retry policy.

Every default Microgen/MMG specimen fails at least one unchanged admission
gate. Original Meshers Gyroid linear24 and curved12, and Schwarz linear24 and
curved12, pass geometry. Initial Meshers spinodoid routes fail quality, while
raw smooth-lattice input fails periodic scalar compatibility. Later work fixes
or explicitly transforms some of those Meshers routes; it does not erase the
original failures.

| Initial admitted Meshers route | Nodes | Mixed-load disagreement | G3 total workflow |
|---|---:|---:|---:|
| Gyroid linear24 | 17,647 | 3.7152% | 17.48 s |
| Gyroid curved12 | 19,931 | 0.2951% | 47.08 s |
| Schwarz linear24 | 19,004 | 1.3458% | 15.21 s |
| Schwarz curved12 | 19,628 | 0.0738% | 38.07 s |

These initial disagreements use the initial declared references, not the later
conservative multi-model envelopes. Faster rejected default MMG rows are not
matched-accuracy victories for either workflow.

### Why tuning and seed inspection were needed

An analytic retained-cap quality ceiling was derived and self-checked:
`q_tet <= 4*sqrt(3)*A/S`, where `A` is a fixed triangular cap's area and `S` its
sum of squared edge lengths. Caps below 0.1 cannot support an admitted adjacent
tet while those faces remain fixed. This obstructs particular Gyroid, spinodoid,
and lattice seeds, not every possible Microgen mesh.

Tuned Schwarz defaults exposed first-pass MMG failures and geometry loss.
Intermediate meshes and actual MMG output were audited. Seed scans identified
cap-feasible alternatives. One explicitly modified one-pass recipe eventually
qualified. Gyroid first-pass FEM, interior-quality, skewness, and required-entity
grading controls still failed quality. Fixed-realization spinodoid seed grids
16–41 remained obstructed in the tested screens; grid41 used a separate bounded
seed-only inspection. Tuned lattice VTK51/two-pass MMG eventually passed geometry.

Unsupported routes, failed retries, timeouts, and stage-cap stops remain records,
not omitted samples. Microgen's public Tet10 converters remain unsupported.

Evidence: [Microgen results](experiments/adaptive/MICROGEN_RESULTS.md),
`cap_quality_bound.py`, `g3-cap-quality-bound.json`, `microgen_seed_screen.py`,
`microgen_tuned.py`, and the frozen `g3-microgen-*` records under
`experiments/adaptive/runs`.

### Fixed-boundary Gyroid interior diagnostic

The existing boundary-face audit already rules out cap and wall-face quality
ceilings on two tested outputs. A new single-star feasibility check reaches
quality 0.1 on all 26 incident tets while preserving boundary and signed
volume. The two-pass parent would lift to 325,911 Tet10 nodes, over the stage
cap, whereas the earlier FEM43 one-pass parent would lift to 153,010 nodes.

One bounded repair of that smaller saved parent accepts 108 interior moves
and reduces bad tets from 222 to 49. Minimum quality remains 0.0079758, and
the second sweep accepts no moves. Geometry and periodic caps remain unchanged.
The explicitly added SciPy postprocess stays unqualified and receives no FEA.
It is not Microgen defaults or a Meshers-seeded MMG comparison. These failed
searches motivated the feasibility and source-consistency audits described
next. Detailed records are in `MICROGEN_RESULTS.md`.

The next single-star diagnostic finds a negative weighted concave-tangent
bound at quality 0.1. Extended-precision replay confirms it. A topology
preflight then stops at a same-side source cavity. The subsequent lineage
audit finds 684 same-orientation shared faces in the VTK43 seed and 146 in
the one-pass FEM output and partial repair. Positive determinants and matching
boundary arrays did not catch them. Retire that source rather than apply
more moves or bypass the cavity guard.

The two-pass Gyroid and all representative admitted/native controls pass this
necessary incidence check. Every archived Schwarz timing arm and manual lift
is checked separately or by exact array identity, so those timing conclusions
are unaffected. The common future admission now fails closed on incoherent
facets. Legacy metric replays and the old 96-comparison, 205-hash verification
still pass. Next compare one fixed-resolution extraction-order variant before
another MMG recipe. This does not certify absence of curved or nonlocal
intersections and does not supply an admitted Gyroid comparator.

### Latest source-validity study ledger

These checks explain why work moved upstream from interior relocation to
Microgen extraction. They do not establish an admitted new competitor.

| Job | Work and result | Evidence under `experiments/adaptive/runs` |
|---|---|---|
| 334–335 | Single-star feasibility succeeds, but whole-mesh repair leaves 49 bad tets and fails quality; later source audit retires the repair input | `g3-gyroid-firstpass-fem43-interior-relocation-v1` and the corresponding single-star records |
| 336 | A weighted concave-tangent bound rejects the exact 0.1 quality target for one fixed-neighbor star; extended-precision replay agrees | `g3-gyroid-convex-star-feasibility-v1`, `g3-gyroid-convex-star-witness-audit-v1.json` |
| 337 | The proposed two-to-three cavity preflight stops on same-side source neighbors; no flips occur | `g3-gyroid-two-to-three-preflight-v1` |
| 338 | Twelve-mesh lineage/control audit finds 684 oriented-face conflicts in the Gyroid VTK43 seed and 146 in the one-pass FEM output; representative admitted controls have none | `g3-facet-incidence-lineage-controls-v1` |
| 339 | All 24 archived Schwarz timing-arm meshes match coherent controls; six manual Tet10 lifts preserve their parent corners | `g3-schwarz-admitted-facet-replay-v1.json` |
| 340 | New admission rejects incoherent facets while reproducing legacy geometry metrics; the previous 96 tensor comparisons and 205 hashes still verify | `g3-facet-admission-check-v1` |
| 341 | The relocation entry point rejects the actual invalid source before any move, new mesh, or FE solve | `g3-gyroid-invalid-source-rejection-v1` |
| 342 | Six-tet preclipping preserves all sample hashes and removes face conflicts, but introduces 24 obstructing caps and fails distance | `g3-gyroid-preclip31-v1` |
| 343–344 | Exact-point-order guard catches a five-tet conversion mismatch; corrected retry preserves samples, passes incidence and distance, but retains the same 24 cap obstructions | `g3-gyroid-preclip31-five-v1`, `g3-gyroid-preclip31-five-v2` |
| 345–346 | First attempt fails JSON serialization; corrected retry finds passing planar flips for all 24 blocked cap faces | `g3-gyroid-cap-flip-preflight-v1`, `g3-gyroid-cap-flip-preflight-v2` |
| 347–348 | Complete two- and three-tet star checks find 22 disjoint valid volume plans; no partial replacement occurs | `g3-gyroid-cap-cavity-preflight-v1`, `g3-gyroid-cap-cavity-preflight-v2` |
| 349 | Each unresolved six-vertex cavity has 13 tet candidates and 8,191 nonempty subsets, with no matching signed boundary | `g3-gyroid-cap-cavity-enumeration-v1` |
| 350–351 | Single-neighbor expansions fail, but complete nine-vertex rings solve both unresolved stars with independent integer-face and volume replay | `g3-gyroid-cap-cavity-expansion-v1`, `g3-gyroid-cap-cavity-fullring-v1` |
| 352–354 | First complete copy exposes near-flat cones; guarded rebuilding produces a globally consistent, periodic seed with no obstructing caps, but quality still below 0.1 | `g3-gyroid-cap-replacement-v1`, `g3-gyroid-cap-cavity-rankguard-v1`, `g3-gyroid-cap-replacement-v2` |
| 355–357 | Actual MMG on the explicitly prepared seed produces over-cap outputs; the fitting first stage has distance 0.00203065 and quality 0.033234, so count alone cannot qualify it | `g3-gyroid-prepared-two-pass-mmg-v1`, `g3-gyroid-prepared-two-pass-mmg-hausd0015-v1`, `g3-gyroid-prepared-two-pass-mmg-hausd00175-v1` |
| 359, 361 | Fitting two-pass output passes distance, volume, incidence, and periodic caps but fails quality 0.047711; localization finds 104 bad cells with no fixed-face ceiling obstruction | `g3-gyroid-prepared-two-pass-mmg-second0015-v1`, `g3-gyroid-prepared-mmg-quality-localization-v1` |
| 364 | All 104 bad cells have free interior vertices; none is all-boundary; this supports a bounded repair trial but does not prove feasibility | `g3-gyroid-prepared-mmg-quality-localization-v2` |
| 366, 368 | Repair accepts 64 moves but leaves 16 bad tets; the worst fixed-neighbor star has a negative checked bound, independently replayed, with three other interior vertices in its support | `g3-gyroid-prepared-mmg-interior-repair-v1`, `g3-gyroid-prepared-repair-star-feasibility-v1`, `g3-gyroid-prepared-repair-star-witness-audit-v1.json` |

The source defect predates MMG in this saved lineage. The recorded evidence
does not yet isolate which clipping or tetrahedralization operation introduces
it. Coherent corner-face incidence is necessary, but it does not exclude all
nonlocal intersections or certify curved elements everywhere.

## 6. Direct fedoo Tet10 integration

`Mesh.to_fedoo()` was added for Tet4 and Tet10, preserving curved coordinates
and ordering in independent arrays. It does not transfer boundary conditions,
surface tags, or periodic transforms. It avoids, rather than fixes, the tested
Microgen `SingleMesh.from_pyvista` and `BoxMesh.from_pyvista` Tet10 rejection.

The integration study checked two solid patches and admitted curved Gyroid and
Schwarz meshes, with 24 prescribed-strain solves. Independent stresses, energy,
reactions, periodic jumps, and saved stress-controlled tensors agree within
about 1.2e-11 relative. Mid-edge jumps are below 6e-17. VTK/fedoo reference
coordinates required an explicit mapping; the initially incorrect test mapping
was corrected in a fresh run rather than changing conversion to fit the test.

Stronger competitors were then constructed by manual straight-Tet10 elevation
of genuine admitted Microgen meshes, preserving their planar domains. These
are explicitly manual elevation plus direct fedoo, not default Microgen Tet10
support and not Meshers-seeded MMG.

Evidence: [Tet10 results](experiments/adaptive/TET10_RESULTS.md),
`experiments/adaptive/tet10_integration.py`, and the public conversion tests.

## 7. Stress recovery and spatial sizing

An optional periodic `size_field` input was added to adaptive generation.
It accepts positive physical lengths, validates finite values and callback
errors, checks opposite periodic traces, and combines the requested sizes with
existing curvature, thickness, bounds, and gradation. Default generation is
preserved when no size field is supplied.

| Study | Measured outcome | Decision |
|---|---|---|
| Tet4-to-straight-Tet10 recovery diagnostic | Six-load compliance-error matrices satisfy fixed-domain enrichment identities; useful Schwarz localization | Justified one live linear-sizing pilot, not stress certification |
| Schwarz linear recovery sizing | At common 1% target, 25,854 nodes and 31.76 s versus uniform 35,393 and 37.77 s | Useful single-run improvement within Meshers; retain optional API |
| Gyroid linear recovery sizing | Coarse trial adds cost; tighter trial reduces nodes but worsens disagreement | No matched-accuracy gain established |
| Lattice recovery ranking/enrichment | Some localization, but the later predeclared incremental recovery screen fails | Do not promote recovery from ranking alone |
| Initial live lattice sizing trials | No demonstrated workflow benefit; some runs hit resource caps | Stop blind sizing-strength tuning |
| All-load non-recovery energy screen | Top 20% captures 94.02% correction trace and at least 92.24% in every mixed load | Justified one separately labeled energy-sizing study |
| Charged energy-sizing workflow | 204,956 nodes, 812.17 s, 5.573 GiB, 1.352008% worst retained-model disagreement | Fails 1%; worse than previous best; stop this recipe |

The last pilot includes all twelve fresh coarse and final loads. Independent
field and positive-rule audits pass; fixed-field energy shifts by only 2.49e-7
relative. Passing localization does not imply an efficient remeshed workflow.
An earlier stage-cap rejection remains rejected even when a later fresh run
uses the previously validated larger solve stage.

Evidence: [stress-sizing results](experiments/adaptive/STRESS_RESULTS.md),
`stress_recovery.py`, `sizing_pilot.py`, `p2_recovery.py`,
`lattice_energy_screen.py`, and `energy_sizing_pilot.py` under
`experiments/adaptive`.

## 8. Periodic geometry and native correctness fixes

Implemented changes include:

- Keeping adaptive projection inside the periodic cell partition rather than
  crossing a retained internal box plane.
- Evaluating the fixed spinodoid with its exact compiled spectral field instead
  of changing its realization or accepting an unsafe coarse projection.
- Using absolute box-plane tolerances in independent audits and revalidating
  frozen inputs without rewriting their reports.
- Retaining valid partial curvature at lattice junctions and separately checking
  sampled curved Jacobians.
- Honoring the requested tetrahedral quality floor during refinement.
- Preserving wall projection targets for initially off-field vertices.
- Enforcing tet budgets inside topology operations and before audit/solve.
- Allowing quality-preserving partial settling of pending wall vertices before
  attempting local subdivision.
- Correcting exact min/max tie derivatives to average operand derivatives,
  including overflow-safe half weighting. Primitive values and explicit user
  gradient callables are unchanged.

The smooth lattice's nonassociative scalar construction has a mismatch in void
across periodic faces. An explicitly recorded positive void clamp at 0.02 keeps
the material and zero boundary unchanged and admits tested routes. A symmetric
algebraic representation was also studied. Fully raw unclamped input still
fails periodic preflight. This is not general raw-field periodic support.

Partial settling admits a previously rejected curved24 spinodoid. A coarser
sizing tolerance then produces a 73,951-node admitted candidate in a fresh
364.07-second workflow, with 0.10333–0.18439% disagreement across six retained
models. Archived-binary replay confirms that the same recipe previously failed
distance. This improves Meshers' admitted route, not its missing Microgen
comparison coverage.

Exact-tie gradient regressions were red before the compiler correction. The
full Python suite afterward passed 133 tests with one optional wide-constraints
skip. All four geometry controls remained admitted. Actual lattice gradient
error dropped to 1.81e-10. A raw-clamped before/after control improved sampled
distance from 0.00039027 to 0.00026671 with 516 fewer nodes. Both already passed
geometry, and no new six-load tensor was measured for that control.

Rejected diagnostics include unsafe coarse spinodoid snapping, collapse,
collective projection, classic VTK clipping, quadratic damping shortcuts,
ineffective coarse lattice repair objectives, selective midpoint restoration,
and affine corner postprocessing. Do not revive them without new evidence.

Evidence: [periodic results](experiments/adaptive/PERIODIC_RESULTS.md),
native and Python regressions, and archived before/after builds.

## 9. References, solver memory, and independent audits

### Fine references and sensitivity

Gyroid and Schwarz initially used sensitivity-checked curved16/20 extrapolation,
with archived linear comparisons where available. Spinodoid later gained a
field-audited fine reference and order/length sensitivity checks.

Lattice uniform and adaptive refinement families disagree enough that close
adjacent-level agreement is insufficient. All 24 empirical lattice models are
retained. A candidate may agree with the finest actual tensor yet fail their
combined 1% envelope. No convenient family is selected after seeing the answer.

Independent field audits compare integrated stresses, reactions, mixed-load
energy, and periodic jumps. Positive 64-point and 125-point reintegration checks
have separate affine and monomial controls. They evaluate fixed fields, not
higher-quadrature equilibrium re-solves. Reduced quadrature as a shortcut was
tested and rejected where it failed the declared sensitivity checks.

### Reference-only cubic enrichment

A research-only Tet20 registration was built and checked on analytic patches.
It enriches the FE displacement on the same quadratic physical geometry,
separating one FE correction from a change in geometric domain. It is not a
public Meshers Tet20 product feature or a continuum reference by itself.

The lattice assembly fit at 5.893 GiB, but initial diagonal and translation-only
block solver attempts failed bounded continuation tests. The actual reduced
system was captured and inspected. Elastic translation/rotation candidates,
node-major block AMG, explicit residual tests, and integrated control solves
justified a new preset.

Further failures exposed redundant symmetry allocations, reaction-access
copies, and retained completed reduced matrices. Complete matrix hashes guarded
symmetry reuse, vector-first reactions avoided unnecessary copies, and explicit
completed-cache release bounded memory. Reused load prefixes were rechecked
against the current equations; they were not charged as fresh workflow solves.

The complete six-load cubic result has 222,331 nodes and 44,976 elements, peaks
at 6.341 GiB, and passes stress/energy, periodic, residual, and positive-rule
checks. It retains status `solved-needs-reference`. Its disagreement with the
24 models, normalized by cubic24, spans about 0.154–1.7505%. This reverse
diagnostic is distinct from candidate-to-model admission and does not settle
lattice accuracy.

On 20 October, a separate cubic geometry probe projects 30,972 higher-order
wall nodes while preserving topology, corner coordinates, and the six saved
displacement coefficient arrays. Both domains pass unchanged geometry gates;
sampled distance improves from 0.00054354 to 0.00013101. Positive 125-point
fixed-field energy changes by 0.0105922%, below the predeclared 0.05% signal
screen. Quadrature sensitivity is 9.965e-7 relative. No projected-domain
equilibrium solve is selected. This does not bound total geometry error or
create a new stiffness reference. Three added interpolation/admission controls
pass. This led to a finer seed25 geometry and node-budget preflight, not
another curve-metric strength trial.

That preflight passes at 75,113 quadratic nodes and 48,696 elements. Exact
shared-edge and face counting predicts 241,555 cubic nodes and matches the
constructed mesh. Assembly-only peak RSS is 5.818 GiB under the same 8 GiB
ceiling. A new budget check prevents cubic construction above 250,000 nodes.
A fresh six-load parent25 control is selected before any higher-order
equilibrium solve. Resource feasibility alone is not reference convergence.

The fresh parent25 control completes in 128.64 seconds and passes independent
field checks. The first cubic25 solve then fails before completing a load
because full-matrix symmetry subtraction needs a 1.47 GiB temporary beyond
the 8 GiB address-space cap. A blocked complete Frobenius audit preserves the
same 1e-12 threshold and matrix-hash reuse gates. Ten focused tests and twelve
integrated control loads pass. A fresh cubic25 retry is selected with the
same resource and first-load limits; no finer-reference accuracy result is
claimed before its completion and audits.

The retry completes all six fresh cubic25 loads at 6.753 GiB peak RSS, with
1,071.18 seconds of worker time. Its own-rule, periodic, higher-rule, and
fresh-process audits pass. Cubic24 and cubic25 differ by 0.248177%, leaving
reference convergence unresolved. Against cubic25, the curve pilot differs
by 0.483742%, the finest adaptive Meshers tensor by 0.651535%, and the strongest
Microgen-derived Tet10 tensor by 0.917416%. Do not choose this single reference
to override the unchanged 24-model admission failures. Cubic25 itself fails
their model-normalized envelope at 1.920819%.

A third neighboring parent23 passes geometry at 68,095 quadratic nodes and
fits the existing cubic stage at an exact 218,508 nodes. One fresh parent
control completes in 119.43 seconds and passes independent field checks.
Actual cubic23 assembly fits at 5.014 GiB. The subsequent fresh higher-order
solve preserves the unchanged caps. These additional nonnested-mesh data
expose reference uncertainty rather than remove it.

Before reading the third cubic tensor, a full-matrix consistency analysis is
defined for a common leading-error model. It accounts for unequal spacings,
checks correction direction and proportionality, tests an order range of
0.25 to 6, and separates correction signal from empirical quadrature shifts.
Synthetic controls precede dataset analysis. Both existing length proxies
remain visible, and the exact retained model tensors are checked. This is a
diagnostic protocol, not a new benchmark admission rule.

The third reference completes all six fresh loads at 5.976 GiB and passes
independent field and quadrature audits. The subsequent three-level analysis
passes seven synthetic controls but rejects both length proxies. Corrections
reverse with cosine -0.959707, and the common-matrix residual is 28.1%.
The exact retained 24-model set and field provenance are verified. Do not
discard the middle level or fit only the closer odd levels. No extrapolated
tensor is promoted. Next measure unresolved quartic curve-edge modes on a
single frozen domain, with coupled overlap terms and independent trial-energy
and periodic-trace controls. That diagnostic subsequently runs and is not an
equilibrium reference or a native API feature.

Three unit tests and manufactured periodic controls pass. The real lattice
uses 102 edge orbits, 306 mode DOFs, and 332 support tets. Its initial 64/125
precision test fails at 6.03%. A single higher-rule 125/216 check passes at
0.017934%, but the trial-energy reduction is only 0.648 ppm, versus the
predeclared 500 ppm signal. No native sizing change is justified. The small
mode-only witness does not bound global coarse-space re-equilibration or
exclude other FE errors. The joint curve-corner relaxation check subsequently
passes its Schur-complement identity and reproduces the original edge gains
exactly. Its incremental signal is still only 0.676219 ppm; higher-rule change
is 0.018333%. Stop this polynomial-mode route without escalating degree.

Evidence: [cubic results](experiments/adaptive/CUBIC_RESULTS.md),
`cubic_reference.py`, `cubic_assessment.py`, `amg_reference.py`,
`capture_reference_system.py`, and `block_amg_probe.py` under
`experiments/adaptive`.

### Volume reference audited before correction

A proposed lattice volume offset was investigated rather than added. Four
nested Sobol sequences were extended from 2^18 to 2^24 points each, totaling
67,108,864 classifications, with exact reproduction of the original prefixes.
The refined mean is 0.2482374460, replicate spread 0.0000137091, and last nested
mean shift 0.0000007004. Both predeclared precision screens pass.

The finest Meshers volume error becomes approximately -0.001898%, with zero
inside the empirical replicate range. About 94% of its previously apparent
deficit came from the old reference estimate. A geometry offset is unjustified.
The Microgen-derived control remains about -0.4582% by volume, inside admission.
Old frozen estimates and stiffness models are not overwritten.

Evidence: [volume audit](experiments/adaptive/VOLUME_RESULTS.md) and
`experiments/adaptive/volume_reference_audit.py`.

## 10. Latest periodic-curve study

The fixed lattice's periodic extension has measured wall-normal jumps of
48.19–109.47 degrees, with matching scalar traces. Full same-domain FE correction
is concentrated near actual wall/box intersection curves: 3.1306% of cells carry
75.0159% of whitened correction trace and at least 70.5900% in every mixed load.
Cap cells away from those curves carry little of this measured correction.
This is distinct from the corrected compiler tie derivative and does not prove
a continuum stress singularity.

One geometric metric was tested, using a fresh coarse mesh and exact projected
curve samples, without coarse FE or oracle fields:

`target = min(0.1, 0.0025 + 0.8 * periodic_distance_to_curve_samples)`.

The matching native gradation is explicitly tuned to 0.8. Both the same-gradation
geometry control and the curve candidate pass geometry admission. The fresh
six-load curve workflow has:

| Measurement | Result |
|---|---:|
| Nodes / Tet10 elements | 118,475 / 68,415 |
| Sampled two-sided distance | 0.00032671 |
| Minimum corner quality | 0.16959 |
| Charged total workflow | 342.19 s |
| Peak RSS | 3.032 GiB |
| Worst disagreement across all 24 models | 1.456109% |
| Positive 125-point fixed-field energy shift | 3.8591e-7 relative |

The pilot fails 1% and is worse than the previous best Meshers result. Independent
field and quadrature audits pass, so their measured errors do not explain the
miss. All six signed energy differences from the finest saved tensor are
positive.

The saved mesh has 2,822 curve edges, median chord 0.00260035. Of these, 288
remain above the native 4/3 split ratio after four rounds; the worst ratio to
the requested endpoint target is 2.85078. This does not prove those splits are
admissible. It identifies a concrete realization question before strengthening
the size field or paying for another solve.

Evidence: [curve study](experiments/adaptive/CREASE_RESULTS.md),
`periodic_crease_audit.py`, `curve_sizing_pilot.py`, and frozen
`g3-lattice-curve-*` records under `experiments/adaptive/runs`. Last job 279
completed the positive-rule audit. No new production/native feature was added
by this geometric pilot.

On 20 October, job 305 tests eight rather than four adaptation rounds while
preserving the curve metric and all gates. It passes geometry at 111,605 nodes
and minimum corner quality 0.30065. Oversized curve edges fall from 288 to 32,
and the worst requested-size ratio from 2.85078 to 1.59475. This establishes
a size-realization gain, not stiffness accuracy. One fresh charged six-load
workflow then completes in 274.06 seconds at 2.851 GiB peak RSS. Its independent
field audit passes, but worst disagreement across all 24 unchanged models is
1.439306%. All eight fine-family models reproduce exactly and both other
families remain retained. Positive 125-point fixed-field energy shifts by only
5.0073e-7 relative. More adaptation rounds alone do not close the gap.

## 11. Current comparison scoreboard

Percentages below are retained-model disagreement, not certified continuum
errors. Single-run and saved-mesh measurements are not repeated workflow ratios.

| Fixture | Best relevant Meshers evidence | Actual Microgen/MMG evidence | Current verdict |
|---|---|---|---|
| Schwarz, 1% target, G3 | Curved8, 6,864 nodes, 0.25870% primary disagreement, 7.7916 s median | Explicit VTK30 one-pass variant, 32,070 Tet4 nodes, 0.92806%, 40.5036 s median | Five interleaved pairs pass; measured speedup 5.1155–5.2908 times, only for this modified comparator |
| Schwarz, supplementary 0.1%, G3 | Curved16, 42,928 nodes, 0.031469–0.093778%, 55.19 s median common-AMG workflow | Actual VTK30 one-pass plus manual straight Tet10, 207,200 nodes, 0.053985–0.062716%, 409.28 s median | Five measured pairs pass; total ratios 7.3858–7.5401, phase ratios 7.2806–7.4331; explicitly modified comparator |
| Schwarz, supplementary 0.1%, Mac | Fresh curved16, 42,928 nodes, 0.031469% to 0.093778%, 19.14923 s median | Actual VTK30 one-pass plus charged interior repair and manual straight Tet10, 207,099 nodes, 0.053480% to 0.062162%, 131.85894 s median | Five fresh interleaved pairs pass; workflow ratios 6.7632 to 6.9203; explicitly modified comparator |
| Gyroid, 1% | Curved9, 11,857 nodes, 0.72921–0.95990%, 20.56 s single run | Tested default and tuned routes remain unqualified | Meshers passes tested target; comparative win not established |
| Spinodoid, 1% | Partial-settle adaptive24 at tolerance 0.003, 73,951 nodes, 0.10333–0.18439%, 364.07 s | Tested fixed-realization seeds 16–41 have retained-cap obstructions | Meshers has an admitted route; comparative win not established |
| Lattice, 1% | Latest audited native result 1.2686985%, 199,972 nodes; eight-round curve pilot 1.439306%, 111,605 nodes | Admitted tuned two-pass parent; strongest manual straight Tet10 control, 238,411 nodes, 1.004354% | Neither passes all 24 models at 1%; no matched-accuracy win |

The strongest lattice Microgen-derived Tet10 control uses validated explicit AMG.
Its six-load solve takes 530.80 s, with 555.78 s for that saved-mesh study and
4.378 GiB peak RSS. Default-CG timeout remains recorded. These numbers do not
include a fresh VTK/MMG workflow and must not be used as such.

## 12. Commit milestones

This groups the committed lineage by substantive work. Formatting, test,
documentation, and failure-recording commits between milestones remain in Git.
Use `git log --reverse --format='%h %ad %s' --date=short` for the complete ledger.

| Date | Commits | Milestone |
|---|---|---|
| 19 Sep | `3e2cbf4`, `d1324b7` | Initial v0.1 release candidate and publication preparation |
| 3 Oct | `e645847` through `2ca2b35` | Scalar/batch/JIT performance, scale and component fixes, optional integration APIs |
| 11–14 Oct | `66b800c` through `dfd0e33` | Geometry guards, optimizer prototypes, adaptive surfaces/tets, torus adaptation |
| 15 Oct | `09598fe` through `80cf40c` | Phase timings, curved elements, export, crease robustness, balanced/fast surfaces, controlled timings |
| 17 Oct | `2e35c7e`, `004ed9c` | G3 migration and actual Microgen defaults study |
| 17 Oct | `31d9cf8`, `4e6e137`, `de64e18` | Direct fedoo conversion, cap quality bound, tuned MMG failures |
| 17 Oct | `3bcca84` through `5f7d74c` | Spatial sizing, intermediate MMG audits, cap screens, Schwarz recovery gain and Gyroid limitation |
| 18 Oct | `d5da490` through `64e6681` | Periodic projection, exact spinodoid and lattice representations, audit tolerance, finer references, partial curvature |
| 18 Oct | `bff2c47`, `af4e25a` | Cap-feasible actual comparator and repeated Schwarz matched-target advantage |
| 18 Oct | `4a24c1e` through `01590ae` | Quality-floor and projection-target fixes, lattice field/reference audits |
| 18 Oct | `115e4a7` through `26cbdfe` | Cheaper Gyroid, admitted lattice comparator, negative correction/repair controls, stronger Tet10 and stricter Schwarz controls |
| 18 Oct | `06558a2` through `dad565e` | Recovery clarification, enforced tet budgets, unsafe coarse spinodoid diagnoses and rejected repairs |
| 18 Oct | `94a14b7` through `14a946c` | Cubic feasibility, isolated systems, bounded AMG, hash reuse, safe reactions/cache lifetime, completed six-load reference |
| 18 Oct | `f8d7534` through `300b5af` | Stronger admitted lattice parent, retained timeout, Tet10 AMG control, error localization and shutdown checkpoint |
| 19 Oct | `090a71c`, `56e1b59`, `f28fa01` | Resumed screens, partial-settle fix, smaller admitted spinodoid, seed41 bound |
| 19 Oct | `5a91f82`, `1ac72d9` | Negative charged energy-sizing study and higher-precision volume reference |
| 19 Oct | `234e497`, `347a010` | Exact-tie gradient correction, negative curve pilot, final audits and shutdown checkpoint |
| 20 Oct | `0ebaa71` | Eight-round curve-size realization improves, but full-tensor admission still fails; consolidated project tracker |
| 20 Oct | `4587cd6` | Fixed-field cubic geometry sensitivity and bounded complete symmetry verification |
| 20 Oct | `b928d1c` | Fresh finer cubic25 reference and independent audits; retained reference uncertainty |
| 20 Oct | `b186d39` | Fresh cubic23 reference, rejected three-level extrapolation, unchanged 24-model set |
| 20 Oct | `06e4f27` | Controlled quartic edge diagnostic, rejected initial precision, resolved weak higher-rule signal |
| 20 Oct | `ca388d2` | Validated corner-coupling diagnostic and weak incremental signal; polynomial-mode route closed |
| 20 Oct | `9f9ee96` | Fresh stronger Schwarz Tet10 pilot with common AMG and independent field/quadrature audits |
| 20 Oct | `2f0363f`, `5789b90` | Repeated exact-recipe Schwarz Tet10 advantage; independent mixed-load, timing and provenance verification; late reporting warning retained |
| 20 Oct | `c6ee7db` | Progress-transition fix, exact archived-source verification, feasible Gyroid star and incomplete whole-mesh interior repair |
| 20 Oct | `3d37ccc`, jobs 336–341 | Checked local obstruction, failed cavity preflight, upstream Gyroid incidence defect, all-arm comparison replay, and fail-closed admission/source guards |
| 20 Oct | `051adb8` | Paired fixed-field preclipping and five-tet conversion controls; coherent seeds still have obstructing caps |
| 20 Oct | `e99fa92` | Local cap flips, complete-star checks, and bounded exhaustive cavity obstruction; no partial replacement |

## 13. Verification and evidence register

Historical test counts belong to their specific source/build snapshots. They
are not claims that a fresh full suite ran while writing this document.

| Checkpoint | Recorded validation |
|---|---|
| Spatial-size API | Linux rebuild, 20 core and 35 focused Python tests; constant-size default preservation, periodic sizes, invalid values and exceptions |
| Direct fedoo API | Focused conversion tests plus 24 integration solves and independent interpolation/energy checks |
| Reference admission | Regression cases keep borderline 0.1% failures rejected rather than rounding them into passes |
| Partial settling | Red native regression before fix; 25 core and 37 focused Python tests after fix; four admitted geometry controls and archived-binary replay |
| Exact-tie derivative | Two red smooth-composition regressions; full Python suite 133 passed, one optional skip; four geometry controls |
| Curve-local pilot | Fresh six-load solve; independent field and positive 64/125-point audits; source/input hashes and bounded worker records |
| Cubic solver memory | Ten focused blocked-symmetry tests and twelve integrated control loads; fresh six-load cubic25 and cubic23 solves and independent audits |
| Three-level reference screen | Seven synthetic controls, source-field provenance, exact retained-model equality, and both declared length proxies checked; scientific model rejected |
| Quartic edge diagnostic | Three unit tests and manufactured periodic controls; direct trial-energy identity, overlap coupling, pin gauge, and periodic traces checked; higher-rule precision passes but useful signal fails |
| Corner-coupling diagnostic | Five algebraic tests, manufactured periodic controls, exact edge-only regression, Schur identity, and 125/216-point precision pass; useful signal fails |
| Live Schwarz Tet10 pilot | 33 controls, fresh VTK/MMG and native meshes, twelve fresh loads, independent field audits, all eight reference models, and native positive-rule sensitivity pass |
| Repeated stronger Schwarz comparison | Five measured pairs and warmup, 38 controls per pair, 72 fresh loads, independent fields; separate verifier checks 96 mixed-load comparisons, timing statistics and 205 hashes |
| Progress and Gyroid repair | Two reporting regressions are red then green; queued single-star and complete-repair controls pass; 18 focused checks pass; repair remains scientifically rejected |
| Source consistency and admission | Six convex controls, extended-precision witness replay, twelve-mesh lineage audit, 24 timing-arm and six-lift replay; admission/repair guards are red then green; 69 focused tests pass |
| Extraction order | Both paired runs verify identical point/scalar hashes; seven current extraction/progress controls and sixteen combined focused tests pass; both altered seeds fail fixed-cap feasibility |
| Rank-guarded cap preparation | Binary controls match exhaustive results; actual near-flat regression and fifteen queued controls pass; fresh 24-plan copy passes independent incidence, periodicity, wall/point identity, volume, and cap-ceiling checks; 103 input/source hash records verified across eight completed runs |

Canonical detailed reports are under `experiments/adaptive`:

| Report | Subject |
|---|---|
| `PLAN.md` | Requested sequence, dated decisions, open work |
| `RESULTS.md` | Cross-study summary and historical adaptive measurements |
| `MICROGEN.md`, `MICROGEN_RESULTS.md` | Exact defaults, actual calls, seed/tuning controls, admissions and timings |
| `TET10_RESULTS.md` | Public direct-fedoo integration and stronger manually elevated competitors |
| `STRESS_RESULTS.md` | Recovery diagnostics, spatial sizing, negative live energy pilots |
| `PERIODIC_RESULTS.md` | Geometry defects, representations, native fixes, rejected repairs |
| `CUBIC_RESULTS.md` | Reference-only enrichment, solver/memory failures and validated completion |
| `VOLUME_RESULTS.md` | Independent nested QMC audit and rejected volume-offset rationale |
| `CREASE_RESULTS.md` | Physical periodic normal jumps, full correction localization, negative curve pilot |

The append-only decision ledger is `.audit/g3-microgen.tsv`. Frozen run
directories retain selected JSON, driver snapshots, and logs. Manifests record
versions, hashes, limits, and inputs. Do not edit historical outputs to match
new conclusions or remove failures from a timing sample.

## 14. Environment and artifact preservation

The G3 has 14 GiB RAM. The tested environment uses Python 3.12, Rust 1.95.0,
Microgen 2.0.0b1, mmgpy 0.17.1, fedoo 1.0.0b2, and PyAMG 5.3.0. Exact per-run
requirements and hashes are authoritative. User-local GLU/OpenGL libraries
satisfy import dependencies without a privileged system installation.

The native binary at the shutdown checkpoint has SHA256:

`908e326211792142f83f854ad2696cd879d927bd81deacf38eb6a0cb63bc4dc1`.

The exact-tie correction changes the Python graph generator, not that binary.
Older binaries are archived locally for controlled before/after comparisons.

Large NPZ meshes, displacements, Gram arrays, curve samples, local packages,
and native binaries remain ignored on G3. Git contains selected evidence, not
a complete binary backup. Some historical large inputs remain only on OMEN.
Small archived tensors cannot support claims of saved-field revalidation.

Use the runtime library path recorded in
[the G3 checkpoint](handoff/g3/CHECKPOINT-20261020.md). Do not store private
machine addresses or credentials in project tracking. GitHub-backed operations
use the verified `kmarchais` account.

## 15. Open work and restart order

Primary work moves to MacBook Pro at the user's request. Read
`handoff/macbook/HANDOFF.md`. G3 jobs 473, 475, and 476 are terminal, and the
lattice26 audits are complete. The selected continuation is a second-level
nested geometry/resource preflight from the solved level1 fields. No G3
successor is queued. Earlier restart instructions below are historical.

The following paragraphs preserve the earlier Gyroid investigation sequence.
Its formerly selected successors have finished, as summarized above.

The latest completed experiment is `experiments/adaptive/preclip_probe.py`,
with `experiments/adaptive/test_preclip_probe.py`. It compares paired Gyroid
resolution31 extraction using identical background coordinates, scalar
samples, and offset 1.1. Normal clipping has 288 oriented-face conflicts.
Pretriangulation removes those conflicts for both five- and six-tet background
decompositions, but each introduces 24 retained-cap obstructions at minimum
quality ceiling 0.08710052. The five-tet variant passes sampled distance at
0.00163287. Neither qualifies, and no MMG or FE solve follows.

Cap-flip and complete-star preflights then identify two stars needing larger
cavities. Complete-ring binary solves supply those replacements. Independent
quality checking catches two additional near-flat cones before MMG. A relative
determinant guard rebuilds all 24 plans as 20 safe cones and four full rings.
The resulting fresh seed has no obstructing caps, preserves points, walls,
volume, incidence, and periodic triangles, and passes sampled distance. Actual
minimum quality remains 0.00109506, below final admission 0.1.

The successor at that checkpoint was one audited, explicitly modified actual Microgen/MMG
geometry replay on this prepared seed. It is not defaults or a fresh complete
workflow timing. All intermediate outputs retain the 250,000-tet cap and
8 GiB ceiling. No FE follows before final geometry and periodic admission.

The first three replays, jobs 355–357, reject at intermediate tet caps.
Replay 359 subsequently fits at 133,826 tets and passes geometry/periodicity
except quality. Repair 366 reduces 104 bad tets to 16 with 64 moves, but still
fails quality. No FE follows. The subsequent boundary-only witness rules out
coupled interior moves at this connectivity. Bounded local topology and one
Steiner-point trials also reject, without proving general impossibility.
Read the latest checkpoint for completed task identifiers and pause state.
The previous 41-test focused suite and seven queued repair controls pass at
their recorded source snapshots.

| Priority | Next action | Completion criterion |
|---|---|---|
| On resume | Review completed lattice26 fields and decide whether to audit them | Independent checks before any stiffness promotion; preliminary error already fails 1%, no automatic solve or parameter sweep |
| Completed negative | Boundary-only constraints and bounded Gyroid topology trials | Witness independently replayed; fixed-connectivity relocation blocked, tested local remeshing rejected; no admitted comparator |
| Retired source | One-pass FEM43 interior repair | Invalid oriented source; current code rejects it before relocation; do not retry sweeps or bypass cavity guards |
| Completed prerequisite | Retain progress values across worker-stage transitions | Red/green success/failure regressions and queued checks pass; shared helper unchanged |
| Completed comparison | Stronger Schwarz Tet10 repetitions | Five interleaved pairs pass, both cost definitions favor Meshers; no further Schwarz tightening selected |
| Closed diagnostic | Curve-corner relaxation coupled to the same quartic modes | Completed; incremental signal remains small, so no further polynomial-mode trial is selected |
| 2 | Resolve lattice FE versus geometry/reference uncertainty | A justified reference envelope and a candidate passing every retained applicable model, without post-hoc model deletion |
| 3 | Establish admitted actual Microgen/MMG comparators for Gyroid and spinodoid, or clearly bound unsupported tested routes | Genuine geometry-qualified controls, not MMG on Meshers seeds or altered specimens |
| 4 | Run fresh repeated matched-target workflows where both competitors qualify | Same-host interleaved repetitions, charged preparation and solves, all scientific gates met |
| 5 | Add true branch-aware crease and box-curve handling only when evidence justifies a native change | Regression plus independent geometry and full-tensor gain; averaged subgradients alone are not crease constraints |
| 6 | Reconsider volume correction only with resolved bias; separately assess parallel relocation | Measured need and benefit, not fitting an uncertain volume estimate; no current justification for a lattice offset |
| 7 | Reduce fast-surface initial extraction density and improve one-thread sphere cost | Matched sampled geometry with fewer triangles and measured same-host speed improvement |

The user's original order was actual workflow comparison, Tet10 integration,
linear stress-recovery sizing, volume correction/parallel relocation,
branch-aware creases/spinodoid curves, then coarser extraction. The first three
have substantial implemented and measured work, but comparison coverage remains
incomplete. Later work was selected from measured blockers, not treated as a
license to batch speculative changes.

General CAD/BRep replacement, universal geometry support, certified Hausdorff
bounds, globally proven curved Jacobian positivity, nonlinear mechanics, and
continuum/local-stress certification remain outside the demonstrated results.
