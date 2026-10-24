# Native MacBook continuation, 24 October 2026

This MacBook session owns primary Meshers continuation on
`handoff/g3-microgen-20261017`, starting at `0b24b564`.
The four-case objective remains active and unachieved. G3 job 476 is complete
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
