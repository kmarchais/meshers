# Native MacBook continuation, 24 October 2026

This MacBook session owns primary Meshers continuation on
`handoff/g3-microgen-20261017`, starting at `0b24b564`.
The four-case objective is paused at the user's request and remains unachieved. G3 job 476 is complete
and has not been restarted. No delegates or cross-model reviews were used.

## Native environment

The M4 has 10 CPU cores and 24 GiB unified RAM. The native environment is
`/Users/kmarchais/dev/meshers/.local-work/macbook-runtime/venv`, using Python
3.12.9. Every requested G3 numerical package version installed from arm64
wheels. No platform version replacement was needed. Each new run records
its installed requirements. Native facts, binary hash, and the resource
probe are in [environment-20261024.json](environment-20261024.json).

Meshers was rebuilt with Rust 1.95.0, the repository lockfile, two build
workers, and `experimental-surfaces`. Mac job 0 passed in 56.26 seconds at
1.307 GiB peak process-tree RSS. Job 1 passed 54 focused tests. Job 2 passed
twelve fresh AMG loads on a cubic affine cube and the saved curved Schwarz
control. Tensor disagreements were below 3.45e-12 and 4.42e-11 respectively.
Periodic, reaction, own-rule, and higher-rule checks passed. These are
environment controls, not workflow comparisons.

The first Mac numerical runs did not explicitly control Apple Accelerate,
which is NumPy's BLAS backend here. A fresh process reports native BLAS mode
0, allowing automatic threading. Their numerical results and audits remain
valid, but numerical timings from jobs 1–8 are excluded from verified
one-thread comparison baselines. Frozen values and original manifests remain
unchanged. The correction is recorded in
[threading-clarification-20261024.json](threading-clarification-20261024.json).

The shared workflow now sets `VECLIB_MAXIMUM_THREADS=1` and calls the native
`BLASSetThreading(1)` API on its calling thread. `BLASGetThreading()` verifies
mode 1 before numerical work. Jobs 9–10 pass 55 focused controls and twelve
fresh cube/Schwarz loads with the corrected policy, without changing tensors.
New manifests and worker records include the actual policy. The Apple API
is thread local; these workers perform numerical work on their main thread.

Pueue 4.0.4 runs through its Homebrew user service. Its shared `gpu` group
has one worker. The restored workbench checkout and external environment
support `qrun --check-start`. A user launch agent serves the local Queue tab
at `http://127.0.0.1:7811/#queue`. Standard Job progress and results were
verified through actual queued stages. No old timers or unrelated jobs were
changed.

macOS rejects the former 8 GiB `RLIMIT_AS` call with
`ValueError: current limit exceeds maximum limit`. On macOS,
`microgen_workflow.bounded` now uses its existing process-tree RSS watchdog
without the Linux pre-exec address-space limit. Worker records declare this
policy and the 0.25-second polling interval. Small success, timeout, and
memory-stop controls pass. Summed RSS is conservative when processes share
pages; polling does not prevent every transient allocation. The cap remains
8 GiB. No 12 GiB protocol or mesh-cap change has been needed.

The transfer was rehashed independently: all 3,509 files, totaling
1,650,689,800 bytes, match. The new receipt is
`/Users/kmarchais/dev/meshers-handoff-20261024/macbook-primary-verification.json`.
Both support repository commits match the handoff. Frozen Linux manifests
and their hashes remain unchanged.

## Second nested level

Mac job 3 completes the selected geometry and field preflight from the
finished G3 level1 mesh. Its fourteen controls pass. The output is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-nested23-level2-preflight-v1`.

The preflight has 218,135 nodes and 141,144 Tet10 elements. It bisects
64,197 parent elements in 12,861 complete periodic edge orbits. Neither
fixed 250,000 cap is reached. The curved physical domain and the six old
displacement fields are retained, with maximum interpolation errors
2.78e-16 and 3.33e-16. Relative volume change is -2.52e-14.

Sampled distance is 0.0006492808, absolute relative volume error is
0.0310173%, minimum corner quality is 0.10010709, and minimum sampled
curved Jacobian is 4.9374e-7. Periodic nodes, cap triangles, and oriented
facets pass. Worker time is 87.72 seconds at 1.047 GiB peak RSS.

The old level1 solve lacks fixture files. The new `--fixture-source` option
reads their hash-checked copies from the original parent23 directory.
Only the new output receives copies. A regression rejects changed samples.
Future nested solves also copy their prepared fixed samples and record
their hashes, so the omission does not recur.

After reading the preflight, job 4 measures the actual first reduced system
without solving or saving large matrix archives. It has 646,029 free DOFs,
51,030,214 stored entries, and 614,946,688 CSR bytes. Complete three-component
node blocks, finite entries, and positive diagonals pass. Worker time is
9.61 seconds at 4.464 GiB peak RSS. Assembly fitting does not guarantee
hierarchy or six-load solve memory, so the same watchdog remains active.

## Selected solve and interpretation

After reading assembly, job 5 completes one fresh level2 six-load solve with
validated elastic-block AMG, 15-point curved Tet10 integration, and explicit
completed-cache release. The solve has a 1,200-second timeout. Independent
field and positive-quadrature audits and Galerkin assessment each have a
300-second timeout. Its immutable output directory is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-nested23-level2-v1`.

All six residuals pass. Solve worker time is 191.45 seconds at 4.508 GiB
peak RSS. Independent fields and positive quadrature pass. Every generalized
parent-to-level2 correction is negative, with maximum magnitude 0.328272%,
versus 0.499513% for the preceding level. The full Galerkin balance agrees
within 3.78355e-8. Worst retained-model disagreement falls to 2.135117%,
still above 1%. No reference is promoted. Nonuniform bisection and two
corrections do not establish an asymptotic order or a continuum bound.
This saved-mesh diagnostic supplies no live workflow timing or cross-host
speed ratio.

After reading that result, job 6 tests one nested preflight on the best saved
native lattice26 candidate, whose original worst model disagreement is
1.2686985%. Its worst-model generalized error is positive, so measured
same-domain FE softening is a relevant hypothesis. The edge-selection rule
is unchanged and uses lengths, periodicity, quality, and budgets. It does
not use oracle fields to tune a size function. The new preflight is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-adapt26-nested-preflight-v1`.
Job 6 passes at exactly 250,000 nodes and 160,291 tets. It accepts 7,971 edge
orbits, bisects 31,991 parents, and stops further orbits at the unchanged
node cap. Map and field errors are 3.33e-16. Volume change is -1.11e-15.
Distance 0.0003236401, quality 0.14932885, positive sampled Jacobians,
periodic caps, and oriented facets pass. Worker time is 11.27 seconds at
1.140 GiB peak RSS. The source's three independent field-audit input hashes
were rechecked.

After reading that result, job 7 measures actual reduced assembly at
726,645 DOFs, 57,890,122 stored entries, and 697,588,048 CSR bytes. It passes
in 10.85 seconds at 5.270 GiB peak RSS. Job 8 then completes six fresh loads and
the same independent audits, under unchanged timeouts and the 8 GiB watchdog.
Its directory is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-adapt26-nested-v1`.
All fields, positive-rule checks, and Galerkin assessment pass. Peak solve
RSS is 5.038 GiB. The largest FE correction is only -0.0122671%, with full
Galerkin balance error 6.07082e-10. Worst disagreement across all 24 models
is 1.2565145%, still failing 1%. Stop this longest-edge recipe. Spending the
remaining node budget there does not establish a useful candidate improvement.

After reading the negative result and passing corrected thread controls,
job 11 preflights one change to edge ordering. It prioritizes actual
wall/cap intersection edge orbits, then uses the existing length order.
Earlier same-domain correction localization places 75.02% of the correction
trace near curve vertices, whereas the longest-edge test resolves little.
There is no sizing-strength parameter or change to physical geometry.
Every periodic, quality, map-restriction, and field-nesting guard remains.
The default length-only path is unchanged. A short-edge control checks the
requested bisection and exact parent map; nonexistent preferred edges reject.

The new preflight directory is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-adapt26-curve-nested-preflight-v1`.
Job 11 passes sixteen controls. It accepts 544 of 549 preferred edge orbits
and fits 249,999 nodes and 158,967 tets. Map and old-field errors are below
3.34e-16. Relative volume change is -3.77e-15. Distance is unchanged at
0.0003236401, corner quality is 0.10271976, sampled curved Jacobians are
positive, and periodic caps and oriented facets pass. Worker time is 13.37
seconds at 1.167 GiB RSS, with the corrected numerical thread policy.

After reading that preflight, job 12 measures actual reduced assembly at
721,914 DOFs, 57,483,976 stored entries, and 692,695,372 CSR bytes. It passes
in 9.97 seconds at 6.230 GiB peak RSS. Job 13 then completes one fresh six-load
solve, independent field and positive-rule audits, and full Galerkin and
24-model assessment. The 8 GiB cap and timeouts remain unchanged. The directory is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-adapt26-curve-nested-v1`.
All six loads, independent fields, positive quadrature, and full Galerkin
assessment pass. Solve worker time is 216.33 seconds at 5.472 GiB RSS. The
largest additional FE correction is only -0.0132094%, and full energy balance
agrees within 4.82986e-10. Worst disagreement across all 24 models remains
1.2556246%, failing 1%. Prioritizing the curve edges gives little additional
gain over the 1.2565145% length-only result. Stop these one-stage bisection
recipes. No public native feature, stiffness reference, or workflow win is
promoted from either negative result.

## Initial verification and return to comparison coverage

[verify_checkpoint.py](verify_checkpoint.py) independently rechecks all
3,509 transferred files and 1,650,689,800 bytes, 150 Mac input hashes, and
72 full generalized-eigenvalue comparisons for the three completed nested
studies. All pass. Twenty-six changed source references resolve to exact
frozen driver or Git content, without relaxing any mesh or field hash.
The receipt is [verification-20261024.json](verification-20261024.json).

Mac jobs 0–13 are complete. No successor is queued at this checkpoint.
The active goal remains the full four-case objective. Return to the missing
admitted actual Microgen comparator coverage, beginning with the preserved
Gyroid and fixed-spinodoid failures in `MICROGEN_RESULTS.md` and the current
coherent seed/cap evidence. Select one bounded preflight from that evidence
before any MMG or solve. Do not repeat rejected interior moves, local cavities,
seed grids, sizing strengths, or polynomial degrees without a new reason.
Lattice reference uncertainty remains open; these negative refinements do
not justify selecting a favorable tensor or removing any retained model.

Schwarz retains its historical repeated wins against explicitly modified
actual Microgen routes. Gyroid and spinodoid still lack admitted actual
comparators in tested routes. Lattice reference convergence remains unresolved.
Fresh Mac comparisons still need same-host, admitted, matched-target workflows.
The fixed specimens, distance 0.002, volume 1%, quality 0.1, sampled positive
Jacobians, periodic gates, coherent facets, and node/tet caps are unchanged.

The append-only decision trail is
[macbook-continuation.tsv](../../.audit/macbook-continuation.tsv).

## Spinodoid comparator work, jobs 14–20

The subsequent primary continuation moves to missing actual-comparator
coverage. Read the new opening section of
`/Users/kmarchais/dev/meshers/experiments/adaptive/MICROGEN_RESULTS.md`
for full recipes, gates, failures, and interpretation.

The shared fixed-specimen loader now restores the exact transferred 52-mode
coefficients before spinodoid sampling. The raw Mac FFT differs in 33 entries
by up to 8.88e-16; the new regression is red before restoration and green
afterward. Mode indices, threshold, and the stored fixture remain unchanged.

Job 14 compares actual Microgen sampling36 extraction with five-tet
pretriangulation. Background coordinates and scalar hashes match. The
233,280-tet background fits the cap. Original extraction has 282 oriented
face conflicts. Pretriangulation removes them at 49,032 nodes and 180,277
tets, but leaves 164 analytically obstructing cap faces. Its finite reverse
sample maximum is 0.00449448 and 148 forward projections are unresolved.
Quality also fails. Both seeds remain rejected, with no MMG or FE.

Job 15 finds a passing individual flip for only 58 of 164 blocked faces.
Job 16 scans 2,388 unambiguous cap boundary corner instances and finds no
acute-corner obstruction. Every blocked face has at least one cap-interior
vertex. A planar motion diagnostic is therefore relevant, but these necessary
checks do not establish a valid cap or volume replacement.

Job 17 fails its manufactured control before production analysis, because
`cells=2` is unsupported. The failed log and exact original test are preserved.
Job 18 corrects the control to `cells=4`, passes derivatives, concavity,
periodicity, and known feasible/obstructed cases, and tests the worst cap face.
Its one free periodic orbit, nodes 13306/13320, cannot meet the cap target with
all neighbors fixed. No move is applied, and the source mesh is unchanged.

Job 19 independently verifies a cap-only supporting-plane witness with
80-digit Decimal arithmetic. At the unchanged 0.1 gate, its residual upper
bound is -0.00063315085191554. The second active cap triangle contains another
free cap vertex, 13894. This does not rule out coupled cap motion or changes
to connectivity. Do not run a whole-cap single-vertex sweep from this result.

Job 20 passes 29 focused controls. The expanded verifier independently
replays both saved seed incidence, cap bounds, and periodicity, verifies
210 input hashes including archived failed-control sources, rehashes all
3,509 transferred files, and reproduces all 72 earlier lattice comparisons.
The new receipt is
[verification-20261024-spinodoid.json](verification-20261024-spinodoid.json).
All jobs through 20 are terminal. Job 17 remains a recorded failure.

The selected next experiment is one bounded coupled planar-patch preflight:
give the two active free cap vertices independent planar motion and move
their periodic copies consistently, while outer curves stay fixed. The
active boundary edge is 4.03873e-6, versus a 0.0277778 edge between free
vertices. Rigid translation preserves that long edge and is too restricted.
The four-variable joint constraints are nonlinear; controls must verify
derivatives and known feasible cases. A failed search must remain inconclusive.
Read direct cap-quality, volume-positivity, and periodic checks before any
mesh edit or MMG. This tests motion identified by the checked witness.
The independent geometry failure remains, so cap feasibility alone cannot
admit this comparator. All four fixtures and the full goal remain in scope.

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

Records: `runs/macbook-spinodoid-cap-coupled-v1` (failed control),
`runs/macbook-spinodoid-cap-coupled-v2` (completed search), and
`runs/macbook-coupled-checkpoint-verification-v1`. Independent receipt:
[verification-20261024-coupled.json](verification-20261024-coupled.json).
Jobs through 23 are terminal. No successor is queued. Earlier selected-next
instructions above are historical; the geometry localization above is current.

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

The final provenance receipt is
[verification-20261024-schwarz-final.json](verification-20261024-schwarz-final.json).
It includes code-source dictionaries as well as artifact/input manifests.
Jobs through 32 are terminal. The failed job 26 and its unqualified parent
remain preserved. Earlier selected-next descriptions are historical.

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

Independent receipt:
`runs/macbook-schwarz-tet10-repeated-verification-v2/result.json`.
The complete checkpoint receipt is
[verification-20261024-schwarz-repeated.json](verification-20261024-schwarz-repeated.json).
Jobs through 36 are terminal. Job 34 remains a recorded verifier failure.
No successor is queued, and no reference is promoted. The full goal is active.

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

Records are `runs/macbook-lattice-hp-shared-domain-v1`,
`runs/macbook-lattice-hp-cross-energy-v1` and
`runs/macbook-lattice-hp-cross-assessment-v1`.
The provenance receipt is
[verification-20261024-hp.json](verification-20261024-hp.json).
All jobs through 40 are terminal. No successor is queued. Earlier selected-next
instructions above are historical. Keep the full goal active.

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

Records are `runs/macbook-lattice-hp-consistent-rule-v1`,
`runs/macbook-lattice-hp-coupled-assembly-v1`,
`runs/macbook-lattice-hp-amg-setup-v1`,
`runs/macbook-lattice-hp-coupled-solve-v1` and
`runs/macbook-lattice-hp-coupled-audit-v1`.
The final provenance and full-model replay receipt is
[verification-20261024-hp-coupled-sources.json](verification-20261024-hp-coupled-sources.json).
The earlier successful receipt and failed job 46 are preserved separately.
Jobs through 48 are terminal. No successor is queued. Earlier selected-next
instructions above are historical. Keep the full goal active.

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

Records are `runs/macbook-lattice-hp-correction-localization-v1`,
`runs/macbook-lattice-hp-ranked-preflight-v1` and
`runs/macbook-lattice-hp-ranked-assembly-v1`.
The independent receipt is
[verification-20261024-hp-ranked.json](verification-20261024-hp-ranked.json).
Missing measured helper snapshots for the localization were restored from
exact committed bytes matching their original input hashes. Original
results and hashes are unchanged. Jobs through 52 are terminal. No successor
is queued; earlier selected-next descriptions are historical.

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

Job 57 independently rechecks all 3,509 transferred files, 1,388 measured
input hashes and 120 full mixed-load model comparisons. The saved receipt is
[verification-20261024-hp-ranked-solve.json](verification-20261024-hp-ranked-solve.json).
The verifier distinguishes raw solved records from complete model audits;
it preserves strict source/data hashes and every retained model. Its exact
source snapshot is retained with job 57. All jobs through 57 are terminal.
No successor solve is queued. The four-case goal remains active.

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
stars, exact parent maps, quality 0.1, positive sampled Jacobians and 250,000-node/tet
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

Job 64 independently replays both curve-star preflights, local correction
arrays, wall projections and all new model comparisons. All 3,509 transfer
files, 1,699 measured input hashes and 168 full mixed-load comparisons pass.
Its receipt is [verification-20261024-curve-star.json](verification-20261024-curve-star.json).
The failed job 58 and unqualified geometry interpretation remain explicit.
The exact verifier source is retained with job 64. Jobs through 64 are
terminal; no successor is queued. The next selected action is fresh Mac
actual-route generation, with every added stage declared and charged.
The four-case goal stays active and unachieved.

## Fresh Mac actual lattice workflow passes, 24 October 2026

The declared actual lattice route now passes fresh Mac generation and all
24 retained stiffness models. Worst disagreement is 0.9946706124%, negative
against sensitivity model 7. This closes the frozen-seed reproduction gap
for one explicitly modified route. It does not establish a lattice Meshers
win: the native star-refined control still fails at 1.230455074%.

Job 65 regenerates the fixed smooth-union specimen through actual VTK51 and
two MMG calls. First pass uses nofem, hausd0.0015, hmax0.06 and hgrad1.3;
second pass uses hausd0.001, hmax0.06, hgrad1.3, ls and nr. The VTK seed has
51,309 nodes and 249,432 tets. First MMG has 24,936 nodes and 115,536 tets;
second has 33,953 nodes and 158,430 tets. Every intermediate count passes.
Full intermediate geometry diagnostics are saved and their cost is charged.

Fresh final quality is 0.0998809921, so the unrepaired Mac output is rejected.
Distance 0.0017012390, volume error -0.4905337%, positive Jacobians, coherent
facets and periodicity pass. Job 66 identifies one bad cell with one free
interior vertex and no fixed-face quality obstruction. No FE precedes repair.
Job 67 adds the existing guarded SLSQP interior relocation, explicitly enabled
only for this exact lattice recipe. Eight controls pass, including rejection
of a changed recipe or failed geometry. One vertex moves 1.43275e-5; complete
star quality reaches 0.100001001. Boundary coordinates and connectivity stay
exactly fixed, and signed volume changes by 3.70e-17. The initial failure is
retained rather than rounded into admission.

Job 68 counts straight Tet10 elevation before allocation: 238,729 nodes.
Shared manual elevation is factored from the existing control without changing
its midpoint construction or default behavior. A shared-edge/node-cap control
and twelve nested controls pass. The declared curve-vertex-star priority
then produces 250,000 nodes and 165,722 tets. Parent geometry and six affine
control fields reproduce exactly within recorded tolerances. These controls
are not equilibrium fields. No unnecessary parent FE solve is performed.
Distance and volume are preserved; q0.100001001, positive sampled Jacobians,
periodic nodes/caps and coherent facets pass under unchanged admission limits.

Job 69 performs six fresh 15-point loads with the same cached block-AMG preset,
new fields and a new hierarchy. Each load has a 240-second continuation gate,
the stage has a 1,200-second timeout, and the RSS watchdog stays at 8 GiB.
Solve plus final field archive takes 218.55733 seconds; peak worker RSS is
6.10088 GiB. Independent job 70 verifies reactions within 8.74e-12,
Hill-Mandel energy within 2.66e-12, and zero periodic jumps. Positive
64-to-125-point energy change is 1.79e-14. Job 71 retains every original
model and passes the unchanged 1% condition at 0.9946706124%.

| Charged fresh Mac stage | Seconds |
|---|---:|
| VTK/MMG generation and intermediate/final geometry audits | 15.290462 |
| Quality localization | 0.811610 |
| Explicit interior repair and admission | 3.422781 |
| Straight elevation, star refinement and admission | 12.487351 |
| Six-load solve and final field archive | 218.557328 |
| Independent field and positive-quadrature audit | 63.405356 |
| All-model assessment | 0.053772 |
| Total | 314.028660 |

Job 72 preserves this full cost receipt. The job 71 named-stage subtotal is
313.217050 seconds; the supplemental receipt also charges quality localization.
Timers retain diagnostics and serialization within their named boundaries.
Synthetic controls, untimed provenance checks, bookkeeping, queue and inspection
gaps are excluded. Intermediate audit overhead is deliberate, and
`timing_comparable` remains false. This is one fresh diagnostic workflow,
not production/repeated paired timing, the default Microgen pipeline, MMG on
a Meshers seed, or a continuum/local-stress certificate. No reference is promoted.

The next Meshers action is a bounded preferred-star-only refinement budget
preflight on the same native lattice26 domain. Its prior priority pass
bisects 6,934 curve-vertex parent cells, then its fallback bisects 23,992
other parents: 77.5787% of all 30,926 bisections are outside curve-vertex
cells. This is a cell-split census, not exact node-cost attribution. Inspect
node allocation and test repeated preferred-star refinement with bulk fallback
suppressed, preserving the parent space and every geometry, quality and
periodic guard. Use the unchanged 250,000-node/tet caps. No new sizing sweep,
model-fitting geometry offset, or solve before reading the preflight is selected.

Schwarz retains its five repeated Mac modified-actual-route workflow wins.
Gyroid and spinodoid still lack admitted actual comparators in tested routes.
Lattice native accuracy and reference convergence remain unresolved. All four
cases remain in scope, and the full objective stays active and unachieved.

Records are `runs/macbook-lattice-actual-vtk51-stage-audit-v1`,
`runs/macbook-lattice-actual-vtk51-quality-localization-v1`,
`runs/macbook-lattice-actual-vtk51-interior-repair-v1`,
`runs/macbook-lattice-actual-live-preflight-v1`,
`runs/macbook-lattice-actual-live-solve-v1`,
`runs/macbook-lattice-actual-live-field-audit-v1`,
`runs/macbook-lattice-actual-live-assessment-v1` and
`runs/macbook-lattice-actual-live-costs-v1`. Jobs through 72 are terminal.
Earlier selected-next descriptions are historical.

Job 73 independently replays all three fresh VTK/MMG stages, the shared
straight elevation, affine geometry/field prolongation, final admission,
all 24 fresh-model comparisons and every charged-cost component. All 3,509
transferred files, 1,947 measured input hashes and 192 full model comparisons
pass. The receipt is
[verification-20261024-fresh-lattice.json](verification-20261024-fresh-lattice.json).
Both exact verifier sources and their hashes are preserved. Jobs through 73
are terminal; no successor is queued. The next selected native budget
preflight does not change the 1% target, fixed geometry or admission caps.
The four-case goal remains active and unachieved.

## Preferred-only refinement improves the control but exposes a projection failure, 24 October 2026

The new native26 fixed-domain control improves worst disagreement across all
24 models to 1.183636899%, from 1.230455074% with a one-round curve-star
priority plus bulk fallback. It still fails 1%. No native candidate, reference
or workflow timing advantage is promoted. The fresh modified actual route
remains admitted at 0.9946706124%, with its separately recorded 314.03-second
diagnostic workflow cost.

The research refiner adds an optional `preferred_only` switch. It filters
complete periodic edge orbits and retains all quality, overlap, count,
geometry and periodic guards. Default and measured-score behavior stay
unchanged; preferred-only mode requires explicit edges and rejects mixed
score ranking. A disjoint-tet control proves that a remote cell remains
untouched despite spare budget, while default refinement still splits it.
All thirteen nested controls and nineteen queued controls pass.

Job 74 repeatedly refines only edges incident on actual wall/cap curve
vertices. Every accepted parent is in that neighborhood; no bulk fallback
is used. Composed maps retain the original native26 domain, old approximation
space and pin 50999. Original geometry and six displacement fields reproduce
within 4.45e-16. Distance stays 0.0003236401, volume changes by 8.88e-16,
quality is 0.10021975, and minimum sampled Jacobian is 1.52083e-9.

| Native26 round | Nodes | Tets | Curve parents bisected | Added nodes |
|---|---:|---:|---:|---:|
| 1 | 213,098 | 135,234 | 6,934 | 13,126 |
| 2 | 227,732 | 143,312 | 8,078 | 14,634 |
| 3 | 244,696 | 152,748 | 9,436 | 16,964 |
| 4 | 249,996 | 155,829 | 3,081 | 5,300 |
| 5 | 249,996 | 155,829 | 0 | 0 |

The implementation accepts no additional star in the remaining budget and
stops. Job 75 completes six fresh loads in 218.97 seconds at 5.832 GiB peak
RSS. Independent field, periodic and positive-rule audits pass. All six
stiffness changes from the original parent are negative, largest 0.08623610%.
Galerkin balance is 3.92e-11 and positive-rule energy change is 9.97e-9.
This supports the allocation change, but the target remains unmet.

The already audited nonadaptive parent25 offers comparable geometry accuracy
at 75,113 initial nodes, rather than native26's 199,972. Its distance is
0.0003819821 and initial quality is 0.322484. It is a distinct fixed domain;
no nested relation or monotone energy comparison to native26 is claimed.
One bounded preflight uses the identical preferred-only rule and 250,000
node/tet caps, without a new size parameter or seed sweep.

Job 76 reaches 108,850 nodes and 64,716 tets after eight rounds. All gates pass,
but the declared round limit stops work while stars are still being accepted.
This is not convergence. Job 77 continues that verified geometric checkpoint
for one more bounded eight-round block. Original source fields and composed
maps are retained; no equilibrium solve or solved prefix is replayed.
Rounds 9–12 pass. Round 13 stops at 161,494 nodes and 88,203 tets because two
forward projections fail. Distance, volume, quality, sampled Jacobians,
coherent facets and periodicity otherwise pass. The mesh remains rejected;
no FE is selected, and no earlier round is substituted to hide the failure.

Job 78 independently replays exactly those two forward samples on face 17087.
Unbounded Newton feet leave the box by 4.00488e-6 and 2.35737e-6. Box-constrained
SLSQP trials find valid zero-level targets at distances 2.68001e-5 and
2.41392e-5, with field residuals at roundoff. A trimmed-plane control proves
that a constrained target can exist when the unbounded normal foot leaves
the box. These feasible points provide distance upper bounds, not certified
closest points. The original rejection, samples and mesh remain unchanged.

Next validate an explicitly versioned, opt-in box-constrained projection
fallback on the same samples. Preserve the default historical auditor,
0.002 distance, 1% volume, q0.1, positive-Jacobian and periodic gates. Keep
failed samples in the calculation; a constrained candidate must satisfy the
same target level and box constraints. Re-audit the unchanged rejected mesh
before any further refinement or solve. Do not silently reclassify the old
failure. This is an observation-method investigation, not permission to
waive geometry admission or fit a mesh to empirical stiffness models.

Schwarz retains its repeated Mac wins against the declared modified actual
route. Gyroid and spinodoid comparator gaps remain open. Lattice native
accuracy and reference convergence remain unresolved. All four cases and
the unchanged objective remain active and unachieved.

Records are `runs/macbook-lattice-adapt26-preferred-only-preflight-v1`,
`runs/macbook-lattice-adapt26-preferred-only-v1`,
`runs/macbook-lattice-parent25-preferred-only-preflight-v1`,
`runs/macbook-lattice-parent25-preferred-only-preflight-v2` and
`runs/macbook-lattice-parent25-projection-failures-v1`.
Jobs through 78 are terminal. Earlier selected-next descriptions are historical.

Job 79 independently replays all saved preferred-only rounds, composed
original-domain maps, edge-neighborhood membership, field reproduction,
geometry/periodic results and the two feasible constrained target witnesses.
It preserves the round13 geometry rejection. All 3,509 transferred files,
2,039 measured input hashes and 216 full mixed-load comparisons pass.
The receipt is
[verification-20261024-preferred-only.json](verification-20261024-preferred-only.json).
All three exact verifier sources and their hashes are retained. Jobs through
79 are terminal; no successor is queued. No Meshers native candidate,
reference or lattice timing win is promoted. The four-case goal remains
active and unachieved.

## An opt-in constrained auditor admits the geometry, but parent25 still misses accuracy, 24 October 2026

The unchanged rejected parent25 mesh now passes a separately versioned
geometry audit. The subsequent cap-budget correction also passes geometry
and FE consistency, but worst disagreement across all 24 models is
2.214326639%, failing 1%. Stop the coarse-bulk, curve-only allocation recipe.
No native candidate, reference or timing advantage is promoted. Native26's
best saved control remains 1.183636899%; the fresh modified actual route
passes at 0.9946706124%.

`geometry(..., box_projection_fallback=True)` opts into
`newton-box-slsqp-v1`. The default historical auditor is unchanged. Only
failed unbounded Newton samples enter SLSQP; every attempted sample remains
in the result. A recovered target must meet the same finite level, normalized
field-residual tolerance of 1e-9, fixed box and original locked coordinates.
The forward distance is measured to that feasible point, a conservative
upper bound rather than a certified closest point. Distance 0.002, volume 1%,
quality 0.1, positive sampled Jacobians, coherent facets and periodic gates stay fixed.
Failed constrained trials remain failed. Ten controls cover a trimmed plane,
locked coordinates, an unreachable level, fixed corners and facet admission.

Job 80 re-audits the exact rejected round13 mesh, prolongation and sample
files. Their hashes are identical. The default still returns two failures;
the opt-in method recovers both. Counts, distance, volume, quality and sampled
Jacobians are unchanged. No refinement or FE occurs in this re-audit, and
job77's rejected result remains frozen. The corrected record explicitly
identifies the changed observation method before continuation is allowed.

Job 81 stops before refining because the checkpoint reader eagerly evaluates
a missing legacy-field fallback. Its exact source and failed log are retained.
The reader and independent replay are corrected; a regression covers a
record with `total_rounds` but no `completed_rounds`, legacy records, and
cumulative counts. Job 82 passes twenty queued controls and continues only
from the explicit same-mesh re-audit, using the same selection rule and caps.

| Continued round | Nodes | Tets | Curve parents bisected |
|---|---:|---:|---:|
| 14 | 177,268 | 95,236 | 7,033 |
| 15 | 196,416 | 103,617 | 8,381 |
| 16 | 218,704 | 113,308 | 9,691 |
| 17 | 244,530 | 124,797 | 11,489 |
| 18 | 249,998 | 127,549 | 2,752 |
| 19 | 249,998 | 127,549 | 0 |

Every accepted parent is curve-vertex incident; no bulk fallback is used.
The implementation finds no further accepted star in the remaining budget.
All fifty-one unbounded projection failures at the final mesh are recovered
by valid box-constrained targets. Distance stays 0.0003819821, quality is
0.1000034664 and minimum sampled Jacobian is 2.75613e-11. The original parent25
domain, approximation space and pin are retained through composed maps.
This domain is distinct from native26; their controls are not nested.

Job 83 completes six fresh loads in 162.95 seconds at 5.715 GiB peak RSS.
All reduced residuals are below 9.79e-13. Independent field, periodic and
positive-rule audits pass. The largest stiffness decrease from parent25 is
0.64684004%; all six signed changes are negative. Galerkin balance is
4.20e-8, and positive-rule energy change is 8.02e-8. These verified corrections
still do not meet the retained target. The shorter saved-mesh solve is not
a matched-accuracy or live timing win.

Next inspect the available cubic25 fields for one shared-domain bulk
approximation probe, before another sizing change. The prepared mesh is
`runs/g3-lattice-cubic25-assembly-memory-v1/mesh.npz`; all six fields are in
`runs/g3-lattice-cubic25-reference-blocked-symmetry-v2/fields.npz`, with audit
`runs/g3-lattice-cubic25-assessment-v1`. First verify exact parent25 embedding,
composed maps, periodic fields, unique-map count and resource feasibility.
Then decide on a bounded positive-rule cross-energy trial without a new
equilibrium solve. Keep every model and distinguish a trial subspace gain
from remaining-error, reference or live-candidate certification. Do not
assume that still deeper curve-only refinement fixes the unresolved bulk.

The fresh actual lattice control remains admitted with its explicit extra
repair/refinement and diagnostic cost of 314.03 seconds. Schwarz retains five repeated
Mac wins against its declared modified actual route. Gyroid and spinodoid
actual-comparator gaps remain. All four cases, fixed specimens, original
thresholds and the full unachieved objective stay in scope.

Records are `runs/macbook-lattice-parent25-box-projection-reaudit-v1`,
`runs/macbook-lattice-parent25-preferred-only-preflight-v3`,
`runs/macbook-lattice-parent25-preferred-only-preflight-v4` and
`runs/macbook-lattice-parent25-preferred-only-v1`. Jobs through 83 are terminal.
Earlier selected-next descriptions are historical.

Job 84 independently replays both projector versions, the same-mesh re-audit,
every recovered target in the continued rounds, original-domain maps,
periodic fields, all model errors and preserved failed workers. All 3,509
transferred files, 2,132 measured input hashes and 240 full mixed-load
comparisons pass. Its receipt is
[verification-20261024-box-projector.json](verification-20261024-box-projector.json).
The historical rejection and schema failure remain explicit. Exact verifier
sources and hashes are retained. Jobs through 84 are terminal; no successor
is queued. No native candidate, reference or lattice timing win is promoted.
The four-case goal remains active and unachieved.

## Parent25 saved h/p trial resolves a surface signal, 24 October 2026

The exact-domain parent25 trial passes positive quadrature checks, but worst
mixed-load disagreement across all 24 retained models is 1.859718171%, above
1%. No reference or native candidate is promoted. The six periodic cubic-minus-
quadratic modes reduce the fine quadratic field energy by at most 0.374103334%.
This is a saved-field trial with fixed coefficients outside those six modes,
not a new equilibrium solution, remaining-error bound or live workflow timing.

Preflight85 verifies transferred cubic data hashes, original parent embedding
within 1.67e-16, composed maps, complete partition, original field reproduction,
and periodic fields on both node sets. Fine quadratic counts remain 249,998
nodes and 127,549 tets; the cubic mesh has 241,555 nodes on 48,696 parent tets.
The physical domain is unchanged. Frozen G3 evidence is read through resolved
path prefixes; no source hashes or inputs are rewritten.

The grouped integrator finds 42,838 map patterns and projects 5,198 seconds for
the two positive rules. Bounded cell batches avoid those singleton map loops.
Three polynomial controls pass. On the actual fields, benchmark86 measures
30.224 seconds grouped and 2.206 seconds batched, with all 144 entries agreeing
within 7.50e-16. This 13.70x integrator ratio is a resource result, not a
meshing-workflow speedup. Four-point tensors are discarded for science claims.

Trial87 fails before integration because its provenance list names the wrong
`periodic.py` path. Its log, failed source and failure note are preserved.
Corrected trial88 finishes in 160.155 seconds at 0.8571 GiB. Positive 125- and
216-point rules change the trial tensor by 3.86e-9 and the gain tensor by
4.58e-7 relative. The high-rule quadratic/cubic diagonal blocks agree with
their saved own energies within 3.17e-7 and 2.68e-5 respectively, below the
unchanged 1e-4 cross-energy check.

The normalized six-load gain puts 20.6100% in curve-vertex cells, occupying
0.08627% of volume. Wall-vertex cells without a curve vertex carry 70.6453%;
cap-vertex cells without a curve vertex carry 36.7881%. These last two groups
overlap and must not be added. Cells touching neither wall nor cap carry
3.71554%. These fractions localize this trial's gain only. They support stopping
the coarse-bulk curve-only allocation, but do not certify the remaining error
or justify transferring a score between different physical mesh domains.

Next inspect a bounded resource preflight for re-equilibrating these six
periodic modes against the fine quadratic space, with consistent assembly
quadrature and full residual checks. Read that preflight before any solve.
Keep the 8 GiB ceiling and 250,000 physical node/tet caps. Do not repeat a
curve-only sizing sweep or treat the trial field as equilibrated. Any subsequent
native change needs its own geometry and full-tensor evidence.

The fresh explicitly modified actual lattice route remains at 0.994671%;
native26's best saved control remains at 1.183637%, and its live candidate at
1.268699%. No lattice matched-target timing advantage follows. Schwarz retains
five fresh paired wins against its declared modified actual route. Gyroid and
spinodoid still lack admitted actual comparators in tested routes. All four
fixed specimens, every admission gate and the original unachieved goal remain.

Records are `runs/macbook-lattice-parent25-hp-preflight-v1`,
`runs/macbook-lattice-parent25-hp-batched-benchmark-v1` and
`runs/macbook-lattice-parent25-hp-trial-v2`. Older next-action entries below
are historical.

Independent job89 passes in 156.434 seconds at 2.5626 GiB, checking all 3,509
transferred files, 2,261 input hashes and 288 full model comparisons. Its new
replay reconstructs mode stationarity, trial/gain congruence, per-cell gain
sums, group fractions and all 48 positive-rule model comparisons. On 64
actual cells it compares both full positive rules to the original grouped
integrator, including a batch size of 17. Failed trial87 remains in the
receipt. The receipt is
`handoff/macbook/verification-20261024-parent25-hp.json`. Jobs through89 are
terminal; no auxiliary response solve has been selected.

## User-requested pause after coupling preflight90, 24 October 2026

The user requested stopping after the current job. Job90 is terminal success,
78.101 seconds at 0.50165 GiB. No successor is queued. The goal is paused,
not achieved. Resume only on explicit user instruction.

The parent25 15-point coupling preflight passes. Its mode Gram and gain differ
from the positive216 rule by 0.0265102% and 0.0288212% respectively, below the
unchanged 1% consistency gate. Virtual-work error is 6.66e-16 and all six mode
Gram eigenvalues are positive. The saved coupling has 749,994 by six entries,
35,999,712 bytes. Four focused controls pass before submission; three cross-
energy controls also pass inside the worker. Saved coupling algebra, force
balance and exact input hashes were rechecked directly after completion.

This is a force/rule preflight only. No reduced matrix assembly, AMG setup,
auxiliary solve, reference promotion or native-candidate promotion follows.
A further full checkpoint replay was not queued after the stop request.
The prior full replay89 remains the last complete evidence audit. Receipt:
`handoff/macbook/pause-20261024.json`.

On explicit resume, first read job90's saved result and source metadata, then
inspect one bounded reduced-assembly memory preflight using the existing
`hp_coupled_preflight.py --parent25 --assembly` route. Read it before choosing
AMG setup or a solve. Do not restart completed jobs or run older next-action
instructions automatically. All four cases and every gate remain in scope.
Schwarz's five modified-actual-route paired wins stand; Gyroid/spinodoid actual
comparators and lattice native accuracy/reference convergence remain unresolved.
