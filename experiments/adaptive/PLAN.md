# Adaptive implicit meshing: plan and status

## A modified actual lattice control passes, 24 October 2026

One saved actual Microgen-derived lattice control now passes all 24 retained
models at 0.988497945% worst disagreement. The new native curve-star control
still fails at 1.230455074%. These are fresh Mac solves of saved meshes, not
live regeneration or a matched-accuracy timing comparison. No lattice speed
win, new reference, or continuum/local-stress certificate is claimed.

Jobs 58–59 localize the complete ranked-space increment between the two
previously audited coupled solutions. The own-rule Galerkin identity error is
7.64e-14, and positive 125-to-216-point local-trace change is 8.08e-8.
Curve-vertex incident cells carry 66.8117% of the correction in 2.4994% of
cells. Literal-curve-edge cells carry 23.3435% in 0.4557% of cells. Groups
overlap. The top 1% of cells carry 73.2238%. This is resolved discrete
correction, not an estimate of remaining error.

The first diagnostic finishes all three FE rules, then fails its wall
projection assertion. Its source, arrays and log remain frozen. The second
run reuses those hash-identical completed arrays and reports the projection
failure explicitly. Eleven Newton projections on six faces leave the box,
by up to 6.76e-5. Surface-rule stability passes, but complete trimmed-wall
geometry interpretation does not. Neither the original geometry admission
nor its samples are silently changed. No co-location conclusion follows.

The optional `--prefer-box-curve-stars` preflight mode prioritizes all edges
incident on actual wall/cap curve vertices, including radial edges, before
the existing descending-length fallback. Complete periodic orbits, disjoint
stars, exact parent maps, q0.1, positive sampled Jacobians and 250,000-node/tet
limits remain unchanged. It is mutually exclusive with literal-edge priority;
the default selection is unchanged. A radial-edge selection control and all
18 queued controls pass. No size parameter sweep or oracle field is used.

| Saved-mesh control | Nodes | Tets | Worst error across all 24 models | 1% target |
|---|---:|---:|---:|---|
| Original live native lattice26 | 199,972 | 128,300 | 1.2686985% | Fail |
| Earlier literal-curve-edge nested native | 249,999 | 158,967 | 1.2556246% | Fail |
| New curve-vertex-star nested native | 250,000 | 159,226 | 1.2304551% | Fail |
| Original tuned actual Microgen straight Tet10 | 238,411 | 158,252 | 1.0043541% | Fail |
| New star-refined actual Microgen-derived control | 249,998 | 165,772 | 0.9884979% | Pass |

Jobs 60–61 retain the original native lattice26 domain and parent field to
3.34e-16. They bisect 78.1472% of its curve-vertex incident cells. Distance
stays 0.0003236401; quality is 0.1023041 with positive sampled Jacobians.
All six stiffness changes are negative, largest 0.03868354%. Galerkin balance
is 3.71e-10 and positive-rule energy change is 9.64e-9. Six fresh loads take
212.13 seconds at 5.758 GiB. Stop this allocation recipe without promotion.

The actual control's initial worst signed error is +1.0043541%; its most
negative error across all models is -0.8937705%. That motivates one bounded
refinement test, not rounding the original failure into admission. Jobs
62–63 retain its affine domain exactly. Distance stays 0.0016215853, volume
change is zero, quality is 0.1015540, and sampled Jacobians are positive.
They bisect 73.6880% of curve-vertex incident cells. Six fresh loads take
232.18 seconds at 6.198 GiB. Galerkin balance is 3.86e-14, the largest
stiffness decrease is 0.09635235%, and positive-rule energy change is
1.64e-14. Full field and periodic audits pass. The worst final error is
negative against sensitivity model 7. More FE refinement is not automatically
helpful to this empirical target.

The passing control retains actual tuned VTK51/two-pass MMG provenance,
manual straight Tet10 elevation/direct fedoo, and an explicitly added
periodic star-refinement step. It is not the default Microgen pipeline or
MMG on a Meshers seed. All 24 original models and original failures remain.
The full parent23 reference investigation remains unresolved at 1.7413203%.

Next prepare fresh Mac actual VTK51/two-pass MMG generation with this declared
added refinement and charge every stage. Read geometry and counts before a
new solve; the frozen G3 seed does not prove fresh Mac MMG reproduction.
The native control still needs improvement before a matched-target pair.
Schwarz retains the five repeated Mac workflow wins. Gyroid and spinodoid
still lack admitted actual comparators in tested routes. All four cases and
the unchanged objective remain active.

Records are `runs/macbook-lattice-hp-ranked-localization-v1`,
`runs/macbook-lattice-hp-ranked-localization-v2`,
`runs/macbook-lattice-adapt26-curve-star-preflight-v1`,
`runs/macbook-lattice-adapt26-curve-star-v1`,
`runs/macbook-lattice-actual-microgen-curve-star-preflight-v1` and
`runs/macbook-lattice-actual-microgen-curve-star-v1`.
Jobs through 63 are terminal. Earlier selected-next descriptions are historical.

## Ranked nested coupled correction, 24 October 2026

Mac jobs 53–56 complete the measured-score nested-space investigation.
The 249,999-node, 162,599-tet mesh retains the original parent23 curved
physical domain and the previous approximation space. Root-map geometry
reproduction is within 1.67e-16. Fifteen versus positive 216-point integration
changes the mode Gram matrix by 0.020572% and the fixed-lift trial gain by
0.025114%, below the declared 1% precision gate.

The prolonged quadratic base is not equilibrated on the refined space.
Its six reduced residuals are 19.29–31.69% of their internal-force norms.
The generalized correction therefore includes six fresh lift-residual
responses as well as six fresh mode responses. A dense full augmented-system
control verifies the condensation and rejects dropping the lift residual.
Combined AMG setup fits at 5.584 GiB. All twelve responses finish in 467.17
seconds at 5.769 GiB, under the unchanged 8 GiB ceiling and explicit timeouts.
These are reference-research auxiliary responses, not twelve macro workflows.

The complete modal field audit takes 174.96 seconds at 5.580 GiB. Own-rule
energy agrees with the predicted tensor within 3.96e-14, and stress-energy
balance is 2.46e-14 relative. Positive 125-to-216-point energy change is
4.37e-9. Fresh operator replay checks all twelve auxiliary responses below
5.56e-12, all six complete equilibrium states below 1.66e-12, and mode
stationarity within 2.56e-12. Periodic jumps are below 4.45e-16 and the pin
is exactly zero. All six signed stiffness changes from the preceding coupled
space are negative; the largest decrease is 0.09312825%.

Worst disagreement across every retained model improves from 1.82618412%
to 1.74132028%. It still fails 1%. No model is deleted, reference promoted,
or native candidate promoted. This controlled fixed-domain result is not a
live Meshers regeneration, same-host competitor timing, or continuum/local
stress certificate. The best live lattice candidate remains the audited
adapt26 result at 1.2686985%. The reference problem remains open.

Read the remaining correction and geometry evidence before selecting a
native change. Do not repeat the exhausted literal-curve-edge sizing recipe
or assume this parent23 investigation resolves geometry bias. Gyroid and
spinodoid still lack admitted actual comparators in tested routes. The five
fresh Mac Schwarz pairs retain their explicitly modified actual-route win.
All four cases and their original admission gates remain in scope.

Records are `runs/macbook-lattice-hp-ranked-coupling-v1`,
`runs/macbook-lattice-hp-ranked-setup-v1`,
`runs/macbook-lattice-hp-ranked-solve-v1` and
`runs/macbook-lattice-hp-ranked-audit-v1`. Jobs through 56 are terminal.
Earlier selected-next descriptions are historical.

## Measured-score lattice allocation preflight, 24 October 2026

Job 49 localizes the complete coupled correction, using all six loads and
two declared mixed directions. Every element tensor contributes to the full
Galerkin identity; the 15-point identity error is 1.15e-13 relative to the
baseline. Positive 125-to-216-point local-vector changes are 1.066e-5 for
the six-load normalized trace, 1.364e-5 for the largest-correction direction,
and 1.019e-5 for the worst retained-model direction. All pass the 1% precision
gate. The worker takes 155.91 seconds at 0.622 GiB.

| Positive 216-point localization | Element fraction | Six-load correction fraction |
|---|---:|---:|
| Incident on a curve vertex | 1.98379% | 59.10487% |
| Containing a curve edge | 0.43006% | 17.30880% |
| Top 1% ranked by correction | about 1% | 54.1449% |
| Top 5% ranked by correction | about 5% | 87.4500% |

Geometry groups overlap; these fractions are not a partition. The trace is
normalized by the complete baseline stiffness matrix. Localization also
retains the worst correction and worst model mixed-strain witnesses. This
is resolved discrete correction, not certified remaining continuum error.
It does not justify repeating the earlier literal-curve-edge-only recipe.

The nested-refinement research helper now optionally ranks complete periodic
stars by measured element-score sum per added-node upper bound. Its default
ordering and every quality, overlap, count, geometry and periodic guard stay
unchanged. A disjoint-star budget regression proves that the measured score
can select a shorter high-score star instead of the longest edge. All eleven
nested-refinement controls pass. No product Meshers API is added.

Job 50 independently replays both positive-rule localization arrays and all
geometry-group shares, then uses the six-load normalized trace for one ranked
preflight. It bisects 21,455 parent cells through 4,723 accepted periodic
orbits. The selected 15.2008% of parents carry 88.80881% of the measured
score. Final counts are 249,999 nodes and 162,599 tets, within unchanged caps.
Quality is 0.10003495, sampled Jacobian minimum is 2.47990e-7, and distance
remains 0.0006492808. Volume changes by 6.67e-16 relative. Geometry and the
quadratic base reproduce within 4.45e-16. Periodicity, coherent facets and
the original pin node 4847 are preserved. The unchanged cubic function span
and included old quadratic space retain the old augmented space.

Job 51 captures the actual reduced operator without solving. It has 740,964
mesh unknowns and 58,642,444 nonzeros, with 706,673,188 CSR bytes. Peak RSS is
5.394 GiB under the unchanged 8 GiB watchdog. This is mesh-operator assembly,
not proof that coupling columns and AMG together fit. Job 52 independently
replays all child nodal geometry and field values, score coverage, sampled
geometry and periodic base fields from saved arrays. Maximum nodal replay
error is 1.12e-16. All new meshes and arrays remain saved, with hashes.

Next is one coupled-force and AMG resource preflight on this ranked space.
The saved prolongation is the old coupled field's quadratic base, not a new
Tet10 equilibrium solution. Account for its six-load residual through a
general coupled correction or fresh macro solves; do not reuse the earlier
formula that assumes an equilibrated base. Retain the complete modal field,
fixed cubic span, original pin, 8 GiB and 250,000-node/tet limits. Read those
checks before choosing any solve. No new accuracy result or reference is
promoted; the last audited coupled disagreement remains 1.826184% across
all 24 models. Gyroid and spinodoid actual-comparator gaps remain open, and
the repeated Mac Schwarz result stands. The full four-case goal is active.


## Audited coupled lattice correction, 24 October 2026

Jobs 41 to 43 establish a consistent quadrature and resource path before any
coupled solve. The 15-point mode Gram matrix differs from the positive
216-point result by 0.0478437%; its trial gain differs by 0.0572437%, both
below the fixed 1% precision gate. Assembled mode forces satisfy direct
virtual work within 5.03e-16 and annihilate rigid translations. Full and
reduced coupling storage is 31,411,440 and 31,009,392 bytes.

Reduced assembly has 646,029 unknowns and 51,030,214 nonzeros. The six global
scalar modes give 646,035 augmented unknowns, with the same 218,135 physical
nodes and 141,144 tets. The sparse mesh operator stays separate from those
six columns. Assembly peaks at 6.044 GiB; setup of the validated four-level
block AMG peaks at 5.359 GiB and stops before CG. The 8 GiB ceiling and all
250,000-node/tet admission caps stay unchanged.

Job 44 then solves six fresh auxiliary mode responses with a 240-second
limit per response and a 1,800-second worker timeout. It uses saved baseline
fields as explicit reference inputs; it is not six fresh macro workflows.
Every baseline load is checked against the assembled operator. The positive
six-by-six Schur system gives a 0.3254513% own-rule energy reduction from the
nested baseline. The worker takes 177.49 seconds at 5.991 GiB. No dense
augmented graph is passed to AMG, and no solved prefix is replayed.

The complete displacement representation combines quadratic response fields
with the frozen cubic-minus-nested modes. It is saved as `hp-fields.npz`;
ordinary nodal Tet10 interpolation would lose its cubic part. Job 45 checks
that representation independently, including the constant pin gauge:

| Check | Result |
|---|---:|
| Own-rule energy versus Schur prediction | 1.23e-13 maximum mixed-load difference |
| Own-rule stress/energy relative difference | 2.14e-12 |
| Positive 216-point energy shift | 3.003e-6 |
| Positive 125-to-216-point energy difference | 6.001e-9 |
| Freshly reassembled auxiliary reduced residual | at most 1.075e-11 |
| Freshly reassembled complete-field reduced residual | at most 2.038e-10 |
| Mode-equilibrium relative residual | 1.594e-12 |
| Six-load periodic jump | at most 4.45e-16 |
| Pin displacement | 0 |

The full-field audit takes 139.39 seconds at 5.720 GiB. Every sampled
Jacobian remains positive. All 24 original model tensors are retained.
Worst disagreement is 1.82618412%, so the 1% target still fails. This is an
audited equilibrium in the stated augmented space on the frozen approximate
parent23 domain. It is not a promoted reference, a best native Meshers
candidate, a live workflow timing, or continuum/local-stress certification.

Job 46 records a provenance-verifier failure: its old source lookup tried
a commit predating the h/p helper. Job 47 verifies the exact committed
source bytes. Job 48 generalizes that lookup to exact same-path Python blobs
on the current branch and records the verifier hash and repository head.
A regression checks historical lookup after HEAD changes and rejects wrong
bytes. Data hashes never use this fallback. No solve or field audit is rerun.

Next, localize the resolved complete-field correction energy, using all six
loads and the full Galerkin energy identity. Check positive-rule precision
and concentration on actual geometry groups before selecting any new
refinement or Meshers change. This is resolved correction, not a certified
estimate of remaining continuum error. Do not revive curve-only sorting,
raise caps, tune an offset or add another polynomial degree from this result.
Gyroid and spinodoid comparator coverage remain open; the verified repeated
Mac Schwarz comparison stands. The full four-fixture goal remains active.


## Shared-domain lattice field combination, 24 October 2026

Job 37 proves that the saved second nested Tet10 level and cubic23 fields
use the same original curved parent23 physical domain. Both solve-to-mesh
provenance links and eight transferred NPZ hashes pass. All ten child Tet10
nodes reproduce the composed quadratic map, and all twenty cubic nodes
reproduce the original quadratic geometry, within 1.67e-16. Reference-volume
fractions sum to one for every original parent. The two prolongation checks
cover all 76,947 and 141,144 children, reproducing fields within 3.34e-16.
Periodic caps/nodes and all six affine field jumps pass.

The nested mesh has 218,135 nodes and 141,144 tets; the cubic mesh has
218,508 nodes and 43,844 tets. Neither mesh, physical geometry, field array,
fixture, quality threshold or count limit changes. The original cubic mesh
is read from its hash-identical prepared directory, not regenerated.

Job 38 integrates all 144 cross terms between the twelve saved load fields.
It optimizes six global periodic cubic-minus-nested difference fields. Source
coefficient arrays remain immutable; mode weights vary, without independent
nodal re-equilibration. Two controls check exact mixed polynomial integrals
under nonidentity child maps and cancellation of manufactured opposing errors.

| Positive rule | Maximum relative trial-energy reduction | Direct energy-identity error |
|---|---:|---:|
| 125 points | 0.259381777% | 1.31e-16 |
| 216 points | 0.259381631% | 1.98e-16 |

The gain-matrix rule change is 6.71442e-7 relative, below the existing 1%
precision gate. The gain exceeds the existing 0.05% signal screen. Own-field
energy differences from solved tensors stay below 7.69e-6; minimum sampled
Jacobian is 4.94816e-7. The worker takes 90.98 seconds at 0.472 GiB, under
its 600-second timeout and unchanged 8 GiB watchdog. These are diagnostic
costs, not live meshing-workflow timings.

Job 39 independently solves the 18-by-18 constrained-basis system instead
of reusing the six-mode formula. Both rule results reproduce, with saddle
residual below 1.44e-16. Every original lattice model is retained. Trial
worst disagreement is 1.89983697% across all 24, still failing 1%. This is
not a new equilibrium tensor, native candidate, reference or continuum bound.

The resolved gain justifies one coupled-space feasibility preflight. First
check a consistently assembled quadrature rule against the positive-rule
screen, then measure coupling storage and reduced assembly under the existing
8 GiB and 250,000-node/tet caps. Record the six extra scalar mode unknowns
explicitly. Do not apply AMG blindly to globally dense mode rows. No full
coupled solve follows before those checks pass. This uses already solved
spaces, with no higher degree, offset or speculative sizing sweep.
Gyroid and spinodoid comparator coverage remain open; the repeated Mac
Schwarz result stands. The full four-fixture objective remains active.


## Verified repeated Mac Schwarz comparison, 24 October 2026

Job 33 completes one excluded warmup and five fresh interleaved pairs, with
alternating arm order. Every pair passes all eight declared full mixed-load
Schwarz models at the supplementary 0.1% target. Native errors span
0.031469% to 0.093778%; Microgen errors span 0.053480% to 0.062162%.
All geometry, quality, sampled Jacobian, coherent-facet and periodic gates pass.

| Charged cost | Modified actual Microgen median | Meshers median | Five paired ratios |
|---|---:|---:|---:|
| Science-audited workflow | 131.85894 s | 19.14923 s | 6.7632 to 6.9203 |
| Generation, repair/lift and solve | 119.40758 s | 17.48343 s | 6.6949 to 6.8605 |

Workflow ranges are 130.16436 to 132.30789 seconds for Microgen and 19.11868
to 19.54234 seconds for Meshers. Meshers is faster in every measured pair
under both cost definitions. No pair is discarded, and warmup is excluded.

The comparator is actual VTK30/MMG with explicit one-pass overrides, a
charged fixed-boundary SciPy interior repair, and manual straight Tet10/direct
fedoo. It is not Microgen's default or public Tet10 conversion. Meshers uses
fresh nonadaptive curved16. Both use the same block-AMG preset and one
numerical thread. Stage limits stay at 250,000 nodes/tets and 8 GiB; the
whole study takes 980.56 seconds within its one-hour deadline. No previous
mesh, solved prefix or preconditioner is replayed.

All six pairs pass 42 controls. Five verifier unit controls also pass.
Job 34 stops before geometry replay because the verifier omitted support-path
initialization; its exact failed driver and log remain preserved. Job 35
corrects initialization and independently verifies 262 hash records, all
96 full-matrix model comparisons, every timing statistic, and all twelve
saved mesh/field sets. Geometry and periodic field jumps are recomputed
from actual arrays. Minimum quality is 0.100001001 for Microgen and 0.4499953
for Meshers. Maximum reaction/energy error is 7.60e-12; periodic jumps are
at most 5.56e-17. The timing study is not rerun after the verifier failure.

This establishes the repeated Mac Schwarz result for the declared modified
route. It does not resolve the other three fixtures or certify continuum
or local-stress accuracy. Keep all 24 lattice models and all fixed specimens.
The next reference investigation uses the already solved nested and cubic
fields on their original parent23 domain; first prove the physical mappings
and field provenance agree before computing a joint cross-energy diagnostic.
No new solve, polynomial degree, geometry offset or sizing sweep is selected.
Earlier selected-next descriptions below are historical.


## First verified Mac Schwarz pair, 24 October 2026

Job 29 completes one fresh same-host pair after the separately preserved
quality rejection. The actual Microgen route is VTK30 plus one MMG call
with `nofem`, hausd 0.0005, hmax 0.1 and hgrad 1.3, explicit omission of the
second call, one fixed-boundary SciPy interior move, and manual straight
Tet10/direct fedoo. The repair is charged. It is not a Microgen default.
Meshers freshly generates curved16 with one optimization pass and no adaptation.
Both use the same explicit six-load block-AMG preset, one numerical thread,
250,000-node/tet limits and the 8 GiB process-tree RSS watchdog.

Both pass geometry, quality, positive Jacobians, coherent facets and periodic
caps/nodes. Full mixed-load error against all eight declared Schwarz16/20
models is 0.053480–0.062162% for Microgen and 0.031469–0.093778% for Meshers.
Both pass the supplementary 0.1% target. Independent reactions and energy
agree within 7.60e-12, with periodic jumps at most 5.56e-17.

Charged science-audited workflow time is 132.63564 seconds for Microgen and
19.27508 seconds for Meshers, a single-pair ratio of 6.8812. The generation,
lift/repair and solve subtotal is 119.88952 versus 17.60300 seconds, ratio
6.8107. Peak RSS is 4.910 GiB for the Microgen solve and 2.244 GiB for Meshers.
This is one pair with no excluded warmup, not a repeated timing claim.

Job 30 independently replays saved artifacts, full model errors, charged
workflow totals, the earlier geometry witnesses and unchanged repaired
boundary/connectivity. Job 31 reintegrates the fresh native six-load fields
with positive 64- and 125-point rules. Maximum energy shift from the solve
is 2.88454e-9; successive-rule shift is 5.80694e-10. These are fixed-field
checks, not new equilibrium solves or continuum certification.

The selected next experiment is one bounded repeated Mac study: use this
verified pilot and quadrature record, exclude one warmup per arm, then run
five fresh interleaved pairs with alternating arm order and the same explicit
charged repair. Port the repeated runner's historical G3 evidence lookup
without rewriting frozen records. Every pair must pass before its successor.
Gyroid and spinodoid still lack admitted actual comparators in tested routes;
lattice remains at 1.2556246% against all 24 retained models. No new reference
or four-fixture victory is promoted.


## Spinodoid geometry and fresh Mac comparison, 24 October 2026

Job 24 replays the fixed spinodoid samples against the coherent sampling36
seed. All 113 reverse samples above 0.002 are wall samples; none is a cap
sample. The worst 32 nearest faces are wall triangles, including witnesses
well inside the box. Reverse maximum is 0.0044944787. Closest-point replay
agrees with the existing VTK distance evaluator and a manufactured plane
control. Forward Newton failures remain separate diagnostics.

Job 25 tests hypothetical damped wall-node projection on the exact frozen
field. All 22,602 periodic representatives converge below 1e-9 residual
distance, with box coordinates locked and periodic copies consistent.
Reverse maximum improves to 0.0028468542, but 27 samples still fail and 146
tets invert. No projected mesh is applied or admitted. Stop this simple
projection recipe as well as the previous two-orbit cap search.

Job 26 attempts the recorded Schwarz pair afresh on the Mac. The Microgen
parent passes distance, volume, facets and periodicity, but quality 0.0956791
rejects it before FE. Historical G3 wins remain historical. Job 27 localizes
one failing tet, with one free interior vertex and no fixed-face obstruction.
Job 28 reuses the existing complete-star relocation control: one move of
0.000748602 gives global quality 0.100001001 while boundary coordinates,
connectivity and volume stay unchanged. All geometry gates pass. This is
saved-mesh feasibility, not a live workflow or a timing win.

Job 29 selects one fresh paired workflow with that explicit SciPy repair
charged to Microgen, followed by manual straight Tet10 and direct fedoo.
Keep the existing supplementary 0.1% eight-model mixed-load test, all original
admission gates, and independent six-load field audits. Read its actual
result before deciding on repetitions or another experiment. Gyroid and
spinodoid comparator coverage and lattice reference convergence remain open.


## Independent two-orbit cap search, 24 October 2026

Mac job 21 stops at its manufactured fixture: the optimized cube does not
retain the assumed quarter-grid node. Its exact test and failed log are
preserved. Job 22 selects the actual neighboring periodic orbit, passes both
analytic-derivative and known-repair controls, then tests four independent
planar unknowns on the spinodoid. All eight deterministic starts miss the
unchanged 0.1 cap gate. Direct replay finds minimum cap quality about
0.00503593, positive incident tetrahedra, and exactly matching periodic
translations across all candidates. The patch has 24 cap triangles and 51
incident tets. No move, MMG call, or FE solve follows.

This stops this two-orbit search recipe; nonlinear search failure does not
prove infeasibility. Job 23 independently replays every returned candidate
from saved coordinates, all transferred hashes, run inputs, and all 72 lattice
comparisons. No reference or comparator is promoted. Before expanding cap
motion again, localize the coherent spinodoid seed's independent reverse
distance failure (0.00449448 versus 0.002) and unresolved forward projections.
Use that evidence to decide whether the extraction geometry is worth repairing;
do not spend another cap sweep on a source that still fails geometry.
All four cases remain in scope, and fresh same-host workflows are still needed.


## Spinodoid coverage and next cap patch, 24 October 2026

Mac jobs 14–20 preserve the exact FFT16/52-mode realization and compare one
sampling36 extraction pair. Five-tet preclipping removes all 282 oriented
face conflicts, but geometry and 164 retained-cap obstructions still reject.
Only 58 blocked faces have a passing single flip. Polygon-corner bounds find
no acute-corner obstruction, and every blocked face has a cap-interior vertex.

One worst periodic cap-interior star nevertheless has a checked fixed-neighbor
obstruction. An independent 80-digit replay at quality 0.1 passes and identifies
another free cap vertex, 13894, in the active support. Stop single-vertex cap
sweeps. Select one bounded preflight giving both active free-vertex orbits
independent planar motion, with outer cap curves fixed and source volume
orientation guarded. The active boundary edge is 4.04e-6, compared with a
0.02778 free-vertex edge that rigid translation would preserve. Test four
unknowns with analytic-derivative and known-feasible controls. These joint
constraints are nonlinear; failed search is not an impossibility proof.
Read direct candidate checks before mesh edits, MMG, or FE. Distance still fails.

Job 17's unsupported two-cell manufactured control is preserved as failure;
the corrected controls and final 29-test suite pass. Full results and source
provenance are in `MICROGEN_RESULTS.md` and the current Mac checkpoint.
No reference, admitted comparator, or workflow win is promoted.

## Curve-priority diagnostic also misses, 24 October 2026

Mac jobs 11–13 pass sixteen controls, geometry and memory preflights, six
fresh loads, independent fields, positive quadrature, and full Galerkin
assessment. Curve-priority refinement fits 249,999 nodes and 158,967 tets,
but changes worst stiffness by only 0.0132094%. Worst 24-model disagreement
is 1.2556246%, still failing 1%. Stop these one-stage bisection recipes.
No native feature, reference, or workflow advantage is promoted.

Independent checkpoint verification passes all 3,509 transferred file hashes,
150 input hashes with exact historical source resolution, and 72 full model
comparisons. All jobs are terminal and no successor is queued. The full goal
remains active. Return to missing admitted actual Gyroid and spinodoid
comparators, reading preserved coherent-seed and cap evidence before choosing
one new preflight. Do not revive a rejected recipe without new evidence.
Read `handoff/macbook/CHECKPOINT-20261024.md` for the verified Mac environment,
threading correction, full completed sequence, and four-case limits.

## Negative bulk refinement and corrected Mac threading, 24 October 2026

Mac jobs 6–8 finish the fixed-domain lattice26 longest-edge test. It fits
exactly 250,000 nodes and passes every geometry, field, positive-rule, and
Galerkin check. Its largest additional FE correction is only 0.0122671%,
and worst 24-model disagreement remains 1.2565145%. Stop that recipe.

The Mac NumPy wheel uses Apple Accelerate. The initial OpenBLAS/OMP settings
did not explicitly select its single-thread API. Preserve those results,
exclude their numerical timings from verified one-thread comparisons, and
use the corrected policy recorded in `handoff/macbook/threading-clarification-20261024.json`.
Jobs 9–10 pass 55 focused controls and twelve fresh numerical loads after
setting and querying native Accelerate mode 1.

Job 11 selects one curve-edge-first nested preflight on the same original
lattice26 domain. Earlier full FE correction is concentrated near wall/box
curves; the bulk test spends its cap for little gain. Change edge ordering
only, preserving every guard and default behavior. Read the actual preferred
orbit acceptance, counts, geometry, and field reproduction before any solve.
Do not strengthen a sizing recipe, relax caps, or infer a new reference.

## Native Mac results and selected candidate check, 24 October 2026

The native arm64 environment, shared serial queue, 54 focused tests, and
twelve fresh solver-control loads pass. Mac job 3 preflights the second
nested parent23 level at 218,135 nodes and 141,144 tets. Job 4 measures actual
reduced assembly at 4.464 GiB. Job 5 completes six fresh loads and independent
field, positive-rule, and Galerkin checks. The additional FE correction is
0.328272%, and the balance error is 3.78355e-8. Worst disagreement across all
24 models is 2.135117%, still failing 1%. No reference is promoted.

Select one fixed-domain nested preflight on the best saved native lattice26
candidate. Its worst retained-model error is positive. The controlled
softening seen on parent23 makes an FE-resolution check useful without
another sizing-strength sweep or a change to physical geometry. Mac job 6
uses the same periodic edge-star rule, fixed samples, 250,000-node/tet caps,
and 8 GiB RSS watchdog. Read actual admission and counts before any solve.
This is a saved-mesh diagnostic. Fresh same-host workflows remain necessary
for a timing claim. Read `handoff/macbook/CHECKPOINT-20261024.md`.

## Primary MacBook continuation, 2026-10-24

Kevin requests the primary session move to MacBook Pro with Sol 6.1 medium.
Read `handoff/macbook/HANDOFF.md`. G3 jobs 473, 475 and 476 are terminal.
The saved lattice26 audit passes but its 1.2686985% model score fails 1%.
One nested Tet10 refinement preserves geometry and old fields exactly;
its six-load correction is 0.499513% and independent energy balance passes.
This has not established an asymptotic reference. After native Mac setup,
verify transferred data and prepare one second-level nested geometry/resource
preflight before deciding on another solve. Older active/paused instructions
below describe their historical checkpoints.

## Current G3 continuation, 2026-10-18

Active branch: `handoff/g3-microgen-20261017`. The stopping objective is a
demonstrated matched-accuracy advantage over actual Microgen/MMG, not an
isolated meshing improvement. The four-fixture objective remains unmet.

- Schwarz has a repeated 5.1155 to 5.2908 times workflow advantage against an
  explicitly modified actual one-pass Microgen route. It is not the default
  two-pass pipeline. Both representative meshes remain below 1% across the
  conservative fine-window reference models.
  A stronger manual straight-Tet10 Microgen control also passes a supplementary
  0.1% stiffness target. Fresh Meshers curved16 passes that target with 42,928
  nodes and an 88.01-second single-run workflow. There is no new paired live
  timing claim against Tet10; a coarser genuine Microgen control also passes
  the target but remains a saved-mesh study.
- Gyroid curved9 passes admission and field checks with 11,857 nodes,
  0.72921% to 0.95990% reference-model disagreement and 20.56 seconds total.
  The same-helper curved12 control takes 39.17 seconds. These are single runs.
  Tested genuine Microgen variants remain quality-rejected and their unranked
  stiffness diagnostics also exceed 1%.
- Spinodoid's exact FFT16 field is preserved. All integer VTK sampling grids
  16 through 40 have analytic retained-cap quality obstructions; 41 exceeds
  the tet cap. Meshers has an empirical fine reference but no admitted
  Microgen comparator in the tested range.
- Lattice accuracy remains unproven across uniform and adaptive reference
  families. The genuine tuned two-pass Microgen grid51 candidate is admitted,
  but differs by 3.89090% from the finer adaptive tensor. A fresh current-native
  Meshers linear32 control differs by 3.95566%. Neither supports 1% accuracy.
  Recovery localization did not justify another sizing solve. A finer admitted
  adaptive level still reaches 1.26988% across retained reference models.
  Uniform volume correction worsens stiffness and is rejected. A coarse-seed
  projection-repair hypothesis also fails its target and is reverted. A strict
  coarse20 adaptive budget passes geometry but exceeds the current solve-node
  stage. Fine-mesh quadrature sensitivity passes and enables uniform34 within
  8 GiB, but that level still does not settle the retained reference models.
  Lattice Microgen straight-Tet10 routes remain capped or quality-rejected.
  Same-domain enrichment validates recovery localization, but two live sizing
  pilots show no benefit. A final scale trial exceeds the fixed tet cap and is
  excluded, including a mistakenly solved uniform arm. Explicit post-adaptation
  budget checks now precede audit and solve. Native adaptation also stops at
  the first topology operation crossing the tet budget, validated by replay.
  Stop blind sizing-scale tuning. An actual Gyroid first-pass FEM override
  still fails quality and remains unranked.
  Coarse spinodoid diagnostics reject direct snapping, local collapse and
  collective projection. The existing current 24-cell curved route also fails
  distance. Do not promote these repairs. Before another large lattice run,
  assess a reference-only cubic FE enrichment on an analytic solid patch to
  separate quadratic discretization error from geometric reference uncertainty.
  That patch now passes, and the admitted lattice's cubic assembly fits at
  5.893 GiB. Diagonal-CG cannot meet the first-load continuation gate and is
  stopped. Installed AMG passes analytic and archived tensor controls without
  showing a small-case speed gain. A same-geometry reference retry uses a hard
  first-load deadline and retains field and quadrature gates, but also times
  out. Its actual reduced system is captured without solving. Translation-only
  node-block AMG still misses the same deadline; a documented rotational
  candidate extension converges on that fixed system and passes twelve integrated
  control solves. All six lattice loads are now being retried with this validated
  preset under unchanged resource and accuracy gates. Its first load passes,
  but a redundant second symmetry allocation hits the address-space cap. Exact
  matrix hashing now guards reuse of the already verified symmetry property;
  control solves and twelve small regressions pass. The full reference retry
  also archives each field immediately. It retains two solved fields, then
  fedoo's reaction extraction hits a separate matrix-product allocation. A
  vector-first equivalent passes native reaction comparisons; verified prefix
  replay avoids solving those fields again. The resumed reference has recovered
  four columns, then repeated problem caches prevent the fifth allocation.
  Corrected delegate binding alone is insufficient. Explicitly releasing the
  completed reduced cache passes a real six-matrix lifetime regression and
  both six-load control cases. Another same-limit reference run resumes the
  four verified fields and still needs the final pair and audits.
  `CUBIC_RESULTS.md` records these stages. There is no complete new tensor,
  Tet20 product API or workflow win yet.

Every heavy stage uses the shared queue, an 8 GiB worker ceiling, conservative
threads, a timeout and staged mesh caps. Admission gates and frozen outputs
remain unchanged. `MICROGEN_RESULTS.md` and `PERIODIC_RESULTS.md` contain the
recipes, failures and reference limits. No continuum or local-stress accuracy
certification is claimed.

## Historical starting point

Goal (Kevin, 2026-10-14): surface and volume meshes of implicit geometries that are
faster than VTK/MMG, with better element quality, fewer elements, and better
geometry fidelity (volume and two-sided boundary distance).

Historical constraint: a long GPU job owns the `gpu` queue. Until it is idle, run only small
single-threaded checks. Element counts, quality and geometry error are deterministic
and can be measured now; timing waits for the queue.

Historical branch `claude/adaptive-implicit` = `implicit-microgen-ready-20261003` (surface API,
compiled fields) + `fea-accuracy-optimizer` (PR #5 geometry guard, quality floor).

## Ideas, in order

1. Scorecard: independent two-sided distance (Newton projection, VTK exact
   point-to-mesh), volume error, triangle/tet quality, counts. No timing.
2. Size field from the implicit field: sag ~ k l^2 / 8, so l = sqrt(8 tol / k_max),
   k_max from the Hessian of the normalized field. Clamp and grade (hgrad).
3. Surface: coarse extraction, then adaptive isotropic remeshing on the implicit
   surface (split long / collapse short / flip to valence 6 / tangential relax +
   exact projection), with box caps and cap curves as features.
4. Volume: same size field drives tet adaptation (split, collapse, swaps, smoothing)
   with exact boundary projection; interior grades coarser.
5. Optional quadratic (tet10) output with mid-edge nodes projected on the surface.

## Status

- [x] Integration branch builds; 69 core tests pass.
- [x] Scorecard: `scorecard.py` (no timing). Reference volumes cached.
- [x] Surface remesher `remesh.rs`, Python `generate_surface(adaptive=True)`.
  At matched Hausdorff vs Flying Edges: sphere -39% triangles, ellipsoid -66%,
  gyroid -39%; min angle 17-35 deg vs 0-3 deg. Not yet faster than FE.
- [x] Tet adapter `tet_adapt.rs`, Python `generate(adapt=True)`; same seed and
  targets as MMG3D (hausd, hmax 0.1, hgrad 1.3, tagged patches):
  sphere 1e-3: 13.7k tets q_min .51 vs MMG 17.2k q_min .42;
  ellipsoid 1e-3: 16.1k .52 vs 16.0k .36; gyroid 2e-3: 42k .27 vs 52k .14;
  gyroid 1e-3: 75k H 1.0e-3 vs MMG 84k H 3.8e-2 (MMG defect). MMG cannot beat
  its seed's facet error; we project to the exact field.
- [x] Periodic adaptation on the torus (merged copies, internal box planes).
- [x] FEA error per DOF: linear adaptivity ties the lattice; curved quadratic tets
  give 13-22x lower error at equal nodes. Native `generate(order=2)`.
- [x] Issues #3/#4: adapter passes all smooth fixtures at 8 cells; crease lattice open.
- [x] Balanced vertex placement for surfaces (~40% fewer triangles at equal H).
- [x] Curve chords measured in-plane on adaptive paths (found via the spinodoid).
- [x] Thick solid FEA: linear adaptivity +-20%; curved quadratic ~10x better per node.
- [ ] Crease capture for min/max fields (compiled-field branch information); smooth
  unions (k ~ strut radius) already work.
- [ ] Error-driven size field for linear meshes (stress in ligament necks).
- [ ] Spinodoid at 16 cells: remaining in-plane curve error 2.0e-2 does not refine.
- [~] Speed: phase timers; targeted edge removal, cached normals (surface 2.5x,
  volume 1.7x, single-run indications). Parallel relocation not done.
- [x] Fast surface mode (Kevin, 2026-10-15: no quality needed, beat VTK on speed with
  fewer triangles and good fidelity): collapse-only coarsening + error splits, exact
  quadratic error per triangle, local balance repair, parallel extraction. Indicative:
  parity with Flying Edges on 1 thread, 1.6-2.5x faster on 24, 2-3x fewer triangles.
- [ ] Fast mode, single thread: collapse (~26 ms) works on the 2x-dense Freudenthal
  extraction; a coarser or marching-cubes extraction would cut it.
- [ ] Box-edge wedges limit q_min (gyroid ~0.18-0.27) in non-periodic boxes.
- [x] Timing study, surface vs Flying Edges and volume vs lattice and MMG, 1 and 24
  threads (`timing_job.sh`, job 1676; RESULTS.md, fast mode).
- [ ] FEA wall time at a 1% target (`timing.py --fea`, not run yet).
- [ ] Tet adaptation does not scale with threads (gyroid 1.81 s on 1, 1.78 s on 24).

## G3 continuation, 2026-10-17

The real default Microgen VTK/MMG pilot and same-skeleton controls are in
`MICROGEN_RESULTS.md`; protocol and exact defaults are in `MICROGEN.md`.
This is a reviewable diagnostic result, not a completed four-fixture
matched-accuracy comparison. Gyroid and Schwarz-solid references are
supported by refinement checks. Meshers rejects the tested spinodoids and
the smooth lattice, and neither has an admitted fine reference here.
No native meshing code or frozen admission rule was changed. Do not infer
local stress accuracy or a tuned Microgen workflow speed ranking.

## G3 Tet10 integration, 2026-10-17

The public `Mesh.to_fedoo()` route passes Tet4 and Tet10 solid patches and
six-load curved Gyroid/Schwarz reproduction. Full periodic mid-edge pairing,
six-node cap matching and sampled curved Jacobians pass. See `TET10_RESULTS.md`.
Microgen's own `from_pyvista` converters still reject Tet10. Next is the requested
coarse stress-recovery sizing experiment for linear meshes.

## G3 stress sizing, 2026-10-17

`STRESS_RESULTS.md` records the fixed-domain recovery diagnostic and two actual
Schwarz remeshing pilots. The optional `size_field` API preserves existing sizing
limits and passes its native/Python checks. At target 0.04, recovery reduces nodes
by 32.7% and total workflow time by 29.2% against uniform sizing, including its
coarse solve. Both narrowly miss 1% stiffness disagreement. An admitted Microgen
baseline and repeated matched-accuracy timings remain open.

At nominal 0.038, both Schwarz candidates meet the 1% target. Recovery uses
27.0% fewer nodes and takes 15.9% less total time in a single run. The actual
Microgen first-pass audit identifies geometry degradation. An explicit tighter
first-pass override fixes distance but leaves impossible retained-cap triangles
on the 24-point VTK grid. A bounded seed-grid screen is next.

The six-grid screen found no feasible retained caps at resolutions 21, 22, 23,
25, 26 or 27. Gyroid recovery pilots at nominal 0.065 and 0.04 pass admission but
do not establish a matched-accuracy gain. Keep spatial sizing optional. Next
assess the requested volume-correction and relocation priorities against the
measured phase costs and volume bias, rather than assume the Schwarz result
generalizes. The four-fixture Microgen comparison is still incomplete.

## G3 periodic projection and spinodoid, 2026-10-18

`PERIODIC_RESULTS.md` records a reproduced periodic projection bug, its core
regression and the fixed specimen rerun. Projection now rejects a target outside
the current cell before wrapping. All 21 core and 35 focused Python tests pass.
The same spectral spinodoid can use a verified compiled expression of its 52
Microgen Fourier modes. At 32 cells, adaptation returns an admitted linear mesh
with 16,721 nodes, H=0.001864, volume error +0.04458% and minimum quality 0.2543.
The spinodoid stiffness reference and smooth-lattice admission remain open.

The admitted spinodoid now completes all six loads in a live 27.71-second
workflow. The smooth lattice also has an admitted linear route after a positive
void-value clamp that preserves its exact material domain and wall. At 32 cells
it has 48,761 nodes, H=0.001773, volume error -0.43152% and q_min=0.10073, and
completes six loads in 84.14 seconds. All four specimens therefore have admitted
Meshers routes with explicit input representations. Next obtain justified fine
stiffness references for spinodoid and lattice, then complete the matched-accuracy
comparison. Their current tensors alone do not establish stiffness accuracy.

The admitted curved spinodoid now has a complete tensor at 114,659 nodes.
Its two-sided distance is 0.000177 and its full workflow takes 490.86 seconds.
The admitted linear mesh differs by 5.3472% in the worst mixed load. A tighter
0.0015 adaptation target admits 147,894 nodes and is selected for the next
reference-sensitivity solve under the unchanged 8 GiB ceiling and an explicit
1,200-second reference stage. Lattice's coarse curved failure is localized to
mid-edge nodes straightened by the sampled Jacobian validity fallback.

The tighter spinodoid solve completes and passes independent saved-field
stress/energy/reaction checks. Its tensor differs by 0.020306% from the preceding
fine tensor. The empirical model envelope puts the curved candidate at
0.04643% to 0.12747% disagreement, below the 1% target, and the linear candidate
at 5.39348% to 5.47426%. Next test a targeted curved-validity fallback improvement
for the lattice, whose coarse curved geometry still fails. A qualified actual
Microgen competitor and repeated matched-accuracy timings remain outstanding.

The lattice diagnostics identify both excessive straightening after curved
validity failure and incorrect composite gradients at smooth-union ties.
Partial-curving corrections retain the unchanged sampled validity gate and
preserve periodic copies in a new regression. A verified algebraically symmetric
field fixes the gradient without changing compiler primitive conventions or the
physical specimen. The 24-cell nonadaptive Tet10 route is admitted with 69,025
nodes, H=0.0005435, volume error -0.03105% and q_min=0.3410. Its six-load workflow
takes 239.60 seconds. Next test a finer lattice reference candidate and audit its
saved fields. Stiffness accuracy and a Microgen timing advantage remain unproven.

Three admitted curved lattice levels now have complete, independently audited
tensors. Their successive worst mixed-load differences are 0.21527% and
0.30689%, so they do not establish a stable asymptotic rate. The conservative
finer-pair model envelope leaves the cheapest curved candidate at 0.86871% to
2.66674% and the linear candidate at 3.37702% to 5.23059%. Keep lattice accuracy
unproven. The 189,187-node reference peaks at 6.56 GiB under the hard 8 GiB
ceiling. Stop blind grid refinement and return to the genuine Microgen first-pass
diagnostic on the cap-feasible default Schwarz seed. A matched-accuracy workflow
victory remains unestablished.

## Schwarz matched-target result, 2026-10-18

A cap-feasible 30-point actual VTK seed and explicitly modified one-pass MMG
recipe now meet all admission gates and the 1% stiffness target. The default
two-pass pipeline still fails. Against that labeled comparator, an eight-cell
curved Meshers mesh passes at 0.25870% stiffness disagreement with 6,864 nodes,
versus 0.92806% and 32,070 nodes for Microgen. Five interleaved live pairs with
the same prescribed-strain CG helper give median total times of 7.7916 and
40.5036 seconds. Meshers is faster in every pair by factors 5.1155 to 5.2908.
Independent saved-field integration agrees within 1.32e-11 relative.
See `MICROGEN_RESULTS.md` and `runs/g3-schwarz-matched-v1`.

This is one supported matched-target result, not a four-fixture completion.
Gyroid and spinodoid still need admitted comparators; lattice still lacks a
stable stiffness reference. Its adaptive quality rejection blocks useful graded
reference meshes. Next reproduce that rejection in a regression test before
testing a single seed-quality hypothesis. Do not raise resource limits or remove
slow-order reference models to claim a global win.

## Adaptive quality-guard correction, 2026-10-18

An exact failing lattice regression and phase probes identify geometry-error
enforcement as the first quality loss. The main adaptation rounds have no tets
below 0.1, but enforcement passes guard zero and creates over 1,300 bad elements.
Respecting the relative guard alone is insufficient because repeated operations
can halve quality. Splitting and pending wall projection now also honor the
existing caller-required minimum. All 23 core and 36 Python tests pass, including
the exact regression. Temporary source probes and diagnostic export flags are
removed. See `PERIODIC_RESULTS.md` for rejected hypotheses and the decision trail.

The new lattice candidate passes quality at 0.11744 but still fails independent
distance at 0.006986 because a corner remains off the wall. No solve runs.
Spinodoid remains admitted, and the repeatedly timed nonadaptive Schwarz mesh
is byte-for-byte identical in its point and connectivity arrays. Next diagnose
the stalled corner projection while preserving both quality and geometry gates.
This does not resolve lattice stiffness uncertainty or complete the four cases.

The optimizer also misclassifies an off-field wall vertex as unconstrained when
it relies only on a near-zero field residual. A one-tetrahedron regression
reproduces that defect. Shared wall-topology-based projection selection passes
24 core and 36 Python tests and reduces the lattice distance failure from
0.006986 to 0.002284. Four passes reach 0.002245, still above 0.002, so no solve
runs. Schwarz arrays remain byte-identical and spinodoid admission metrics
remain unchanged. Continue the remaining wall-projection investigation without
calling this a lattice reference or a four-fixture win.

The broader wall-repair objective policy passes its controls but gives no
actual lattice improvement, so it is reverted. The verified 24-cell adaptive
seed instead admits 123,570 Tet10 nodes and completes six loads in 639.81
seconds total. Its tensor differs by 0.82695% from uniform32. Both saved-field
audits pass, and positive 64-point and 125-point reintegration changes energy
by at most 5.65e-7 relative. This does not resolve reference convergence.
Next tighten only adaptive tolerance from 0.004 to 0.003 with admission before
FEA, the existing 200,000-node and 1,800-second stage, and the unchanged hard
8 GiB ceiling. Read that tensor before choosing any further refinement.

The tighter 175,818-node candidate completes at 5.64 GiB peak RSS, and all five
saved-field audits pass. Its adaptive-family model envelope is 0.15983% to
0.94546%, but disagreement with earlier uniform-family models reaches 1.30614%.
Retain both families and keep lattice 1% accuracy unproven. The cheaper
adaptive candidate also fails its own worst model at 1.07075%. Before another
costly lattice solve, return to missing actual Microgen coverage by screening
Gyroid VTK seeds for retained-cap feasibility. Read the six bounded seed stages
before choosing one MMG successor. This does not replace the frozen defaults.

The Gyroid study finds two cap-feasible VTK seeds among seventeen tested.
Both actual tuned two-pass MMG routes remain quality-rejected. Separate,
explicitly unranked six-load diagnostics also miss the 1% stiffness target at
2.58544% and 2.29279%; quality repair alone is insufficient. Failed ordinary
and skewness third-pass controls are removed from the maintained runner and
preserved in frozen drivers. Do not claim a matched comparison from them.

The admitted curved10 Meshers Gyroid pilot reaches 0.45586% disagreement in
24.71 seconds, with 13,953 nodes and independently checked fields. Run a fresh
curved12 control with the same prescribed-strain helper and native build before
further coarsening. This distinguishes parameter savings from the older
stress-controlled timing. Lattice cross-family uncertainty and admitted
Gyroid, spinodoid and lattice Microgen comparators remain unresolved.

## Cubic lattice reference checkpoint, 2026-10-18

The reference-only cubic24 enrichment now completes all six loads under 8 GiB.
Fresh-process stress, energy, positive 64-point and 125-point integration and
periodic-jump checks pass. Explicit completed reduced-cache release fixes the
measured fedoo retention failure without changing equations or benchmark
defaults. The six-load tensor is not a live workflow timing or continuum bound.
See `CUBIC_RESULTS.md` and `runs/g3-lattice-cubic-assessment-v1`.

Retained-model disagreement reaches 1.75050% when normalized by the cubic
tensor. All 24 empirical models remain retained. Finest adaptive Tet10 remains unproven at the 1%
target across every model. The next geometry-only stage tests second-pass
hausd 0.001 on the actual VTK51 Microgen lattice seed, retaining first-pass
0.0015, to assess a stronger straight-Tet10 comparator inside 250,000 nodes.
Read its quality, distance and exact lift count before any solve. The existing
Schwarz matched-target result remains the only demonstrated workflow win.

The actual tuned lattice parent passes admission at second-pass hausd 0.001
and lifts to 238,411 straight Tet10 nodes. Default diagonal CG hits its unchanged
1,800-second timeout after five loads. A separately validated explicit AMG
retry completes six fresh solves and independent field checks within 4.378 GiB.
The worst retained-model disagreement is 1.004354%, versus 1.269882% for the
best tested Meshers candidate. Neither passes 1%; do not round or remove models.
No timing victory follows from saved solves under different solver policies.
See `TET10_RESULTS.md` and `runs/g3-lattice-tet10-retained-comparison-v1`.

The worst error direction is predominantly deviatoric normal strain. Next
localize same-domain Tet10-to-cubic enrichment energy in that direction, with
the Galerkin identity and independent quadrature checks, before selecting a
Meshers sizing or geometry change. No speculative native patch is queued.

The measured same-domain cap concentration justifies one analytic periodic
cap-size pilot. It passes geometry and six-load FE consistency at 159,147
nodes, but its worst retained-model error is 1.411382%. This is little better
than the smaller untargeted baseline and worse than the best tested Meshers
mesh. No further cap-length tuning is selected. A confirmed host suspend
invalidates that pilot's timing comparison, without changing its frozen tensor.

Stop study work for Kevin's requested computer shutdown after pushing a
reviewable checkpoint. Do not mark the four-case goal achieved. Read
`/home/kmarchais/dev/meshers/handoff/g3/CHECKPOINT-20261018.md` when Kevin resumes.

## Resumed bounded study, 2026-10-19

Kevin resumes the goal after the shutdown. The Tet10 recovery screen misses
its fixed incremental threshold, so do not promote it to live sizing. The
actual Gyroid required-entity grading variant worsens quality and is stopped.

Spinodoid saved-mesh probes reject selective curvature restoration and affine
corner displacement. They identify a narrower native opportunity: partial
pending-wall progress before quadratic elevation. A red regression confirms
the missing progress. The candidate retains the original quality floor and
branch, containment, periodicity, and Jacobian gates.

Job 230 passes 25 core and 37 focused Python tests and admits the formerly
rejected fixed curved24 spinodoid. Job 231 regenerates it for six prescribed
strain loads. Jobs 232 and 235 pass independent field and higher-rule checks.
All six unchanged reference models give disagreement 0.05430% to 0.13529%.
Job 234 retains all four geometry admissions, including bitwise unchanged
Schwarz and Gyroid controls. Keep the fix, but make no speed claim from the
slower 640.08-second single workflow. Test one 20-cell spinodoid geometry
before another solve. No admitted actual Microgen spinodoid comparator exists.

The 20-cell seed rejects, including its two-pass optimizer control. Stop seed
coarsening. Tolerance 0.003 on the fixed 24-cell seed admits 73,951 nodes and
passes the full tensor target under all six unchanged reference models, with
independent field and positive-rule checks. Total live workflow is 364.07
seconds. The pre-patch binary rejects the exact recipe, confirming that the
native fix enables its admission. Do not claim an interleaved speed benchmark.

Return to comparison coverage rather than continue spinodoid tuning. Actual
Microgen seed41 is cheap enough for a separate 300,000-tet seed-only preflight
under the unchanged 8 GiB ceiling. It still has 190 obstructing cap triangles.
Stop that route without MMG or FEA. No larger final mesh or solve gate is
introduced. Gyroid's missing admitted MMG control and lattice's retained-model
accuracy disagreement remain higher priorities than another spinodoid solve.

## Charged all-load energy pilot, 19 October 2026

Job 250 validates a cheaper non-recovery indicator using all six parent
fields, with no oracle strain in its input. Its top 20% captures 94.02% of
the full correction trace and at least 92.24% for every mixed load. The full
Galerkin matrix identity and basis invariance pass. This does not reverse
the failed incremental recovery screen.

The first live energy-sizing attempt hits the unchanged native tet budget.
One geometry-only floor adjustment passes at 204,956 nodes, above the initial
200,000-node stage. Job 253 then runs a fresh, fully charged workflow in the
already validated 250,000-node stage under the same 8 GiB ceiling. All twelve
coarse and final loads complete in a charged 812.17-second workflow. No old
cap rejection is reclassified.

Jobs 254 and 255 pass independent field and positive-rule checks. All 24
reference models remain unchanged. Worst disagreement is 1.352008%, worse
than the best prior Meshers result and above 1%. Stop energy-sizing tuning.
Read `STRESS_RESULTS.md` before selecting another stress-based change.

Before a volume-correction feature, tighten the independent lattice volume
estimate. Its current four scrambled Sobol estimates use only 2^18 points
each, with absolute spread 0.00015640. That exceeds the apparent Meshers
volume deficit, about 0.000078. A higher-precision, bounded volume audit must
leave old records and all admission gates unchanged. Do not offset geometry
to fit an uncertain volume mean or to manufacture stiffness agreement.

Job 258 completes that audit using four nested 2^24-point sequences. Frozen
prefixes reproduce exactly. New mean volume is 0.2482374460, spread is
0.0000137091, and the final mean-prefix shift is 0.0000007004. Both precision
screens pass. The finest Meshers volume error is -0.001898% at the new mean,
and its empirical replicate range includes zero. Do not add a volume offset.
Keep old records and every stiffness model unchanged. See `VOLUME_RESULTS.md`.

Parallel relocation cannot resolve the measured accuracy gap, and final
meshing costs only about nine seconds of the 812-second failed pilot. Before
changing curved geometry, reproduce the known raw smooth-union gradient
mismatch at tied branches with a focused regression. The compiler currently
chooses one operand's derivative at a min or max tie, whereas its absolute
value derivative is zero there. The existing symmetric-minimum workaround
is not default raw-field support. Test any correction separately, keeping
periodicity, quality, distance, and full stiffness gates fixed.

Jobs 261, 264, and 266 complete that compiler check. Both smooth-composition
gradient regressions fail before correction. Exact min/max ties now average
operand derivatives, consistently with zero derivative for abs(0), without
changing primitive values. The documented subgradient policy changes; an
explicit gradient callable can preserve another policy. This does not provide
true crease constraints.

The full Python suite passes 133 tests with one optional wide-constraints
skip. All four geometry controls stay admitted. The actual raw lattice gradient
matches centered differences within 1.81e-10. A raw-clamped before/after
comparison improves distance and removes 516 nodes, but both arms already
pass geometry. Fully raw scalar input still rejects at periodic axis 0.
There is no new raw-clamped six-load tensor or matched-accuracy timing win.

Return to the unresolved lattice stiffness/reference discrepancy rather than
claiming that compiler correctness solves it. Volume bias and fixed-field
quadrature sensitivity are too small in the measured candidates to justify
an offset or explain their model misses. Further reference work must separate
finite-element enrichment from geometric interpolation error. Preserve the
same-domain cubic result and all 24 empirical models; do not promote another
stress-sizing strength trial from localization alone.

## Curve-local pilot and shutdown checkpoint, 19 October 2026

Jobs 271, 274, 276, and 278 measure periodic-extension normal jumps, localize
the full fixed-domain FE enrichment near wall/box curves, admit a geometric
curve-local size field, and audit its fresh six-load tensor. The workflow
takes 342.19 seconds at 118,475 nodes but fails all-model admission with
1.456109% worst disagreement. The previous best remains 1.269882%. Preserve
all 24 models and do not claim that localization alone improves Meshers.

An endpoint metric audit finds 288 of 2,822 curve edges above the native
split ratio after four rounds. Before another solve after restart, inspect
split eligibility and use one geometry-only increased-round control with
the same metric and gates. Read `CREASE_RESULTS.md`. Do not start a fresh
long solve during the requested shutdown window. Preserve local fields and
record the final independent positive-quadrature audit before stopping.

## Eight-round realization follow-up, 20 October 2026

The resumed geometry control 305 changes only adaptation rounds from four
to eight. It admits 111,605 nodes, raises minimum corner quality to 0.30065,
and reduces oversized curve edges from 288 to 32. Fresh six-load workflow 306
then takes 274.06 seconds at 2.851 GiB. Independent field audit 307 passes,
but all 24 unchanged models still give 1.439306% worst disagreement. The eight
fine-family models reproduce exactly and their envelope itself reaches
1.301824%. More rounds are not a stiffness-target solution.

Do not select another round count or stronger curve metric from localization
alone. Complete positive-rule reintegration, then investigate geometric-domain
approximation separately from FE enrichment or justify a bulk-refinement
control with independent evidence. Preserve all references and resource gates.
Read `CREASE_RESULTS.md` and the new consolidated `PROJECT_TRACKING.md`.

Job 309 subsequently tests a separate cubic boundary mapping while retaining
the original displacement coefficients. Both domains pass the unchanged gates
and sampled distance improves fourfold. Fixed-field energy changes by only
0.0105922%, below the predeclared 0.05% signal screen; quadrature sensitivity
is 9.965e-7 relative. Do not call these projected fields equilibrated or launch
a projected-domain solve from that signal. Preserve all 24 models. Select one
finer nonadaptive seed25 geometry preflight, then count the exact cubic nodes
before considering an independent higher-order refinement solve within the
existing 250,000-node stage. See `CUBIC_RESULTS.md`.

The seed25 parent admits and its exact cubic count is 241,555. Assembly 311
fits at 5.818 GiB. Fresh parent loads 312 and independent audit 313 pass.
However cubic solve 314 fails in the full-matrix symmetry audit before its
first completed load. Blocked exact Frobenius verification preserves the
1e-12 gate and passes ten focused and twelve integrated control cases.
One fresh retry 316 uses the same 8 GiB ceiling and first-load deadline.
Read its actual result before any successor; no reference model changes.

Job 316 completes all six fresh cubic25 loads at 6.753 GiB. Own-rule checks,
higher-rule sensitivity, and fresh-process audit 317 pass. Cubic24-to25 change
is 0.248177%, and the relative ranking of the two native Tet10 candidates
depends on reference choice. Do not select the favorable reference or discard
any of the 24 models. Cubic25 itself misses their envelope at 1.920819%.

Job 318 admits a third neighboring parent23, whose exact cubic node count
is 218,508. One fresh parent23 six-load control is selected in the existing
8 GiB stage. Audit its fields and read resource feasibility before a cubic23
successor. Three nonnested meshes may still fail an asymptotic-convergence
test; do not assume that an additional point proves a continuum reference.

Parent23 six-load control 319 and independent audit 320 pass. Assembly 321
confirms 218,508 cubic nodes at 5.014 GiB. A fresh cubic23 reference is selected
as job 322, with no prefix replay and unchanged first-load and worker limits.
Verify that job's authoritative state before any retry. See the current
`handoff/g3/CHECKPOINT-20261020.md` for continuation and audit instructions.

Job 322 and independent audit 323 pass. The predefined three-level analysis
324 passes seven synthetic controls but rejects both length proxies: the
correction direction reverses, cosine -0.959707, and common-matrix residual
is 28.1%. Exact field provenance and all 24 retained tensors are verified.
Do not promote an observed order, fit only the odd levels, or run another
nonnested grid from this result alone.

Select a bounded fixed-domain residual witness on the saved cubic24 fields.
The proposed quartic edge modes share coefficients across periodic curve-edge
orbits and require their full overlap matrix, independent trial-energy checks,
and periodic trace controls. Read the declared protocol in `CUBIC_RESULTS.md`
before implementation. It remains a diagnostic, not an equilibrated stiffness
reference, a revived recovery-strength recipe, or a new native feature.

Jobs 325 and 326 pass the quartic derivative, affine, coupled-energy, gauge,
and periodic-trace controls. The lattice edge-only witness 327 fails its
64/125-point precision gate. One higher-rule follow-up 328 passes 125/216
precision but has only 0.648 ppm trial-energy reduction, far below the
predeclared 500 ppm signal. Preserve the failed initial route and all sources.

No native feature or stiffness reference follows. Select only one missing
coarse-response check: existing cubic curve-corner functions relaxed jointly
with the same edge modes. Separate coarse-only and incremental quartic gains,
verify the Schur complement and nested-space ordering, and reproduce the
original edge-only result. Keep the same higher-rule precision and signal
screens, memory and support caps. If still small, stop this polynomial route.

Jobs 329 and 330 complete that check. Five algebraic tests and manufactured
controls pass. The lattice incremental signal is 0.676219 ppm versus 500 ppm,
with 0.018333% higher-rule change. Original edge-only results reproduce
exactly and the Schur identity passes. Close this polynomial-mode route.
Do not select another degree, strength, or reference tensor from it. All
geometry and benchmark gates remain unchanged, and the four-fixture goal
is still active and unachieved. Detailed evidence is in `CUBIC_RESULTS.md`.

Next close one missing comparison measurement rather than revive lattice
degree or size-strength tuning. Run a single fresh Schwarz pair against the
stronger genuine Microgen one-pass plus manual straight-Tet10 route, using
the same validated explicit AMG preset on both arms. The declared protocol
is in `TET10_RESULTS.md`; `schwarz_tet10_live.py` charges the fresh VTK/MMG
parent instead of treating its saved mesh as free. Retain the supplementary
0.1% target, all eight Schwarz models, and all geometry/resource gates.
Read the pair and its independent field audits before choosing repetitions.

Job 331 completes that fresh pair. Both arms pass the supplementary 0.1%
target and all eight Schwarz models with the same explicit AMG. The total
single-pair ratio is 7.4382, or 7.3368 for generation/lift/solve only. Both
independent field audits pass. After reading it, positive-rule audit 332
passes with native energy shift 2.89e-9. Select interleaved repetitions of
these exact recipes before any repeated stronger-Tet10 timing claim. Keep
both cost definitions and the explicit one-pass/manual-elevation labels.
No additional lattice sizing or degree trial follows from this result.

The repeated stronger-Tet10 protocol is now implemented in
`schwarz_tet10_repeated.py`: one fresh warmup per arm, then five fresh pairs
with alternating arm order. Preserve the pilot recipes, common AMG, eight
reference models, both cost definitions, and all admission/resource gates.
Eight focused ordering/summary controls pass. Queue this as one bounded
experiment with a 3,600-second whole budget, then read its authoritative
completion before choosing a successor. Host clock discrepancy above five
seconds excludes timing. No lattice or missing-comparator conclusion changes
until new fixture-specific evidence exists.

Job 333 completes all five measured pairs and warmup within 2,923.13 seconds.
All geometry, field, clock and eight-model 0.1% gates pass. Median common-AMG
audited times are 409.28 versus 55.19 seconds; paired ratios are 7.3858–7.5401.
The independent verifier recomputes 96 mixed-load comparisons and every timing
statistic, checks 205 hashes, and passes. Both timing definitions favor Meshers
in every measured pair. Preserve explicit one-pass/manual-Tet10 labels.

Before another experiment, fix the recorded `no_chips` warning. A direct
standard-Job capture proves the end-of-part context clears values. No shared
helper changes or successor job were submitted. After that reporting fix,
return to missing admitted actual Gyroid/spinodoid comparators and lattice
reference uncertainty rather than tightening Schwarz again. All four-fixture
requirements and admission gates remain unchanged.

The progress-transition fix is now tested in `job_parts.py`, using public Job
methods and explicit completion values. Both success/failure capture tests
are red before the fix and green afterward. Queued checks 334 and 335 pass.
Frozen driver verification remains possible after source changes by exact
archive-hash matching; mesh and field changes are still rejected.

Do not repeat the existing Gyroid cap/wall-face audit. Job 334 demonstrates
one feasible complete-star relocation on the actual two-pass output. Exact
Tet10 node counting instead selects the smaller FEM43 precursor for one
whole-mesh repair trial. Job 335 reduces bad tets from 222 to 49 with 108
accepted moves. Worst quality remains 0.0079758 and the second sweep makes
no progress. Boundary, volume, geometry and periodic caps stay unchanged;
the route is still rejected, and no FE solve follows.

Stop that recipe. Next diagnose a single remaining worst star with a convex
slack feasibility problem and a checked supporting-plane bound. A failed
SLSQP search alone does not establish infeasibility. Read this diagnostic
before choosing any coupled relocation, topology repair, or MMG setting.

Job 336 and extended-precision witness replay pass. Three boundary-only
constraints obstruct quality 0.1 for the frozen worst star. The resulting
2-to-3 preflight 337 stops at an invalid oriented source cavity; no flips
are applied. Audit 338 finds the source issue already in the VTK43 seed,
with 684 same-orientation shared faces, and 146 in the one-pass FEM output.
The two-pass Gyroid and all representative admitted controls pass this
necessary incidence check. Retire the invalid one-pass repair route rather
than tune its optimizer or bypass the cavity guard.

All 24 archived Schwarz timing-arm meshes and six lifts pass replay 339.
The common geometry helper now rejects incoherent incidence uniformly for
future measurements, without changing numerical admission thresholds or
frozen flags. Check 340 reproduces legacy metrics and verifies all prior
96 mixed-load comparisons and 205 hashes using exact archived source.

Next queue one paired extraction preflight at resolution 31: normal scalar
clipping versus the same sampled background triangulated before the two
clips. Keep the analytic field, offset, sample coordinates and scalar data
identical, with explicit non-default labeling and unchanged 250,000 background
and material-tet limits. Read incidence, cap, periodic and geometry results
before choosing MMG or FEA. Do not start another resolution or sizing sweep.

Check 341 confirms the retired repair entrypoint rejects the actual invalid
source before any moves. Its seven controls include a red/green positive-
overlap regression. No new mesh or FE solve is produced. Shared geometry
admission also rejects this source, while every legacy scalar geometry metric
and the admitted Schwarz control reproduce. Do not reuse the invalid one-pass
route to tune relocation or cavity operators.

Jobs 342–344 complete that extraction-order check. Normal Gyroid31 clipping
has 288 face conflicts. Both five- and six-tet background variants eliminate
them with identical coordinate/scalar hashes, but each creates 24 retained
caps with analytic quality ceiling below 0.1. The five-tet arm passes distance
at 0.00163287. Neither seed qualifies. No MMG or FE solve runs.

First five-tet attempt 343 stops correctly on point-order mismatch; retry344
uses a checked bijective connectivity mapping and restores exact coordinates
before clipping. Seven queued controls and sixteen combined focused tests pass.
Original failed records and matching source snapshots remain preserved.

Next inspect local cap-triangulation feasibility on the coherent five-tet
seed, preserving periodic cap boundaries. A valid replacement must precede
volume remeshing. Do not send analytically blocked caps to MMG or start another
grid-resolution or sizing sweep. The four-fixture goal remains unachieved.

Jobs 345–349 test local cap flips before volume remeshing. All 24 obstructing
cap faces have passing planar alternatives. Complete-star checks find 22
valid disjoint cavity plans, leaving two three-tet stars. Their existing-apex
cones fail. Independent enumeration checks all 8,191 subsets of 13 positive
tet candidates on each six-vertex cavity and finds no matching signed boundary.
No partial mesh replacement, MMG, or FE runs. Preserve the failed JSON-ID
serialization attempt and exact driver/control snapshots as well as successes.

Next assess one bounded one-ring cavity expansion around those two stars,
keeping target cap polygons and periodic boundaries fixed. Do not repeat cones
or brute-force enumeration on the same six vertices. This is not a general
remeshing impossibility and supplies no admitted new comparator.

Jobs 350–354 complete a cap-feasible copy. All ten single-neighbor expansions
fail, but complete nine-vertex rings solve both missing stars. The first copy
exposes two near-flat cones in independent quality checking, so it is not used
for MMG. A numerical relative-determinant guard, with an actual-source
regression, rebuilds 20 safe cones and four full rings. The fresh copy preserves
points, wall triangles, volume, incidence, periodic caps, and distance while
removing every cap obstruction. Its actual quality remains below 0.1.

Select one geometry-only actual `remesh_keeping_boundaries_for_fem` replay on
that saved prepared seed. Label preclipping/cap retriangulation non-default,
record both actual MMG calls and intermediate geometry, and exclude replay
costs from end-to-end timing claims. First-pass sizing is explicitly overridden
to hausd 0.0005, hmax 0.06, hgrad 1.3 while retaining nofem. The actual second
level-set call remains. No FE before full final admission and periodic checks.

Actual prepared-seed MMG controls 355–357 hit intermediate tet caps. First-pass
hausd 0.0005 returns 649,230 tets and 0.0015 returns 259,134. At 0.00175 the
first pass fits 229,778 tets, but misses distance at 0.00203065 and quality at
0.033234. Its actual second pass at hausd 0.0005 returns 446,978 tets and stops.
No FE or admitted competitor follows. Select only second-pass hausd 0.0015
next, with all first-pass settings and admission gates unchanged. Saved-seed
replays remain explicitly excluded from complete-workflow timing comparisons.

Replay 359 fits 133,826 tets and passes distance, volume, incidence, and periodic
cap checks, but quality 0.047711 rejects. No FE runs. Localization 361 finds
104 bad cells, all touching cap and wall vertices, without any fixed-face
analytic obstruction. Next classify actual free interior vertices before
choosing relocation or another topology/geometry route. Do not infer that a
passing face ceiling guarantees a feasible complete-star move.

Classification 364 finds all 104 bad cells have free interior vertices and
none is all-boundary. Queue one fixed-boundary repair using the existing
guarded relocation algorithm with an explicit source argument. Keep two
sweeps, 400 attempts, complete-star target 0.100001, connectivity, boundary,
volume, and every final gate fixed. Seven focused controls pass. Read this
trial before another MMG or FE experiment. Applicability is not feasibility.

Repair 366 is terminal but incomplete. It accepts 64 moves, reducing bad tets
104 to 16, with minimum quality 0.05126434. The second sweep accepts none.
All boundary, volume, geometry, periodic, and incidence checks remain fixed.
No FE follows. Stop this recipe rather than repeating the same sweep.

Diagnostic 368 finds a checked fixed-neighbor obstruction at worst node19446,
with independent extended-precision replay. Three other interior vertices
participate in the active witness. Check the boundary-only subset of necessary
constraints before choosing coupled relocation or connectivity changes.

Check 370 finds a boundary-only obstruction independent of other interior
moves, with extended-precision replay. Preflight 372 rejects all eight 2-to-3
alternatives. Complete-star trial 374 reports infeasible for 477 quality-filtered
candidates on 13 vertices. Trial 375 rejects one analytic cap-optimal apex because
it lies outside the fixed star. No additional point, source edit, or FE follows.
Stop these particular recipes without claiming universal remeshing failure.

Before more expensive tuning, assessment 377 checks the retained lattice
target. Its pairwise lower bound is 0.9019313%, and a trial SPD matrix passes
all 24 models at 0.9084715%. The 1% rule is algebraically feasible. Keep every model,
the same admission metric, and the distinction between a mathematical trial
and a physical stiffness reference.

Preflight 378 selects current-native lattice26 under the fixed solve budget.
It passes at 199,972 nodes, 128,300 tets, distance 0.0003236401, and quality 0.149329.
Run one fresh regeneration and six-load 15-point curved-Tet10 AMG workflow before assessing
accuracy. The old lattice24binary differs; no pure-refinement or timing claim
comes from the rough count projection. No volume-offset, size-strength,
parallelization, or speculative feature change is selected.

## Pause checkpoint, 20 October 2026

The user requested a pause after the current job finishes. Job 379 is now
terminal success: six fresh lattice26 AMG loads, 199,972 nodes, 128,300 tets,
736.51 seconds charged workflow. Preliminary worst disagreement across all
24 unchanged models is 1.2686985%, so the 1% condition is still not met.
Independent saved-field and higher-quadrature checks remain pending.
No successor experiment or audit was queued. Resume only on explicit user
instruction, starting from the root `PROJECT_TRACKING.md` and the final
20 October checkpoint. Older selected-next instructions above are historical.
