# G3 continuation checkpoint, 20 October 2026

## Active checkpoint after classification 364

Job 364 is terminal success. All 104 bad cells in the fitting prepared-seed
MMG output have interior vertices, and none is all-boundary. The worst cell
has two free vertices. No fixed cap or wall face is analytically obstructing.
This permits a relocation trial without establishing that it will succeed.

Job 366 is queued behind another project's shared work. Observe pueue task
366 or qrun terminal session 21723. Do not restart it or touch other jobs.
Its run directory is `runs/g3-gyroid-prepared-mmg-interior-repair-v1` under
`experiments/adaptive`, created when it starts. There is no repair result at
this checkpoint. Read it before choosing a successor.

The reused `interior_relocation_trial.py` now takes `--source`, defaulting to
the retired original input. Source-incidence rejection remains intact. The
new explicitly accepted source is
`runs/g3-gyroid-prepared-two-pass-mmg-second0015-v1`, with two recorded actual
MMG calls. Seven focused controls pass. The queued limits are two sweeps,
400 attempts, 600-second worker, one thread, and 8 GiB. Every accepted move
must raise its complete incident star above 0.100001, reduce bad-cell count,
and preserve boundary, connectivity, and signed volume. Full admission and
periodicity are rerun before any FE work. No native change or default-pipeline
claim follows. Earlier sections below are historical.

## Active checkpoint after jobs 359 and 361

Both jobs are terminal success but the mesh remains scientifically rejected.
Actual two-pass MMG replay 359 returns 29,876 nodes and 133,826 tets. Its
sampled distance 0.0019272611, volume +0.72004462%, positive sampled Jacobians,
coherent oriented facets, and periodic caps pass. Quality 0.04771095 fails
0.1. No FE runs and no comparative timing claim follows from saved-seed work.

Localization 361 finds 104 bad cells, all touching cap and wall nodes. There
are 75 with a cap face, 26 with a wall face, and zero with analytically
obstructing cap or wall faces. The next classification checks free interior
vertices rather than assuming relocation applies.

Job 364 is queued behind other shared work. Observe pueue task 364 or qrun
terminal session 95223. Do not restart or cancel other jobs. Its output will
be `runs/g3-gyroid-prepared-mmg-quality-localization-v2` under
`experiments/adaptive`, created when it starts. No v2 result exists at this
checkpoint. Read it before selecting a successor. Current script has the
necessary source-incidence guard and standard worker-part progress retention.

Forty-one focused tests pass, plus Ruff and whitespace checks. Commit
`fb0f587` preserves the completed cap preparation and earlier MMG stage-cap
failures. No native code, references, specimen definitions, or admission
thresholds change. Earlier sections below are historical.

## Active queue checkpoint after job 357

Jobs 355–357 are terminal and scientifically rejected. Job 356 returns
259,134 tets at first-pass hausd 0.0015. Job 357 fits its first pass at
229,778 tets with hausd 0.00175, but distance is 0.0020306511 and minimum
quality 0.03323404. Its actual second pass at hausd 0.0005 returns 446,978
tets and stops at the unchanged intermediate cap. No FE runs.

Job 359 is queued behind another project's shared-queue build. Do not touch
that other job. Observe qrun terminal session 60619 or pueue task 359 rather
than restarting it. The selected replay keeps first-pass hausd 0.00175,
hmax 0.06, hgrad 1.3 and changes only second-pass hausd to 0.0015. Its run
directory is `runs/g3-gyroid-prepared-two-pass-mmg-second0015-v1` under
`experiments/adaptive`, which is created only when the job starts. No result
exists at this checkpoint. Read it before selecting any successor.

Pueue's Running and Queued states are JSON objects, like Done. To inspect
active tasks, exclude the Done key rather than excluding all object states.
Earlier abbreviated status filters would miss active work. No local heavy
computation was launched from those observations. Heavy work stays queued.

The prepared seed, all native inputs, reference tensors, and admission gates
remain unchanged. Saved-seed MMG times exclude preparation and do not support
a workflow comparison. Earlier sections below are historical.

## Latest complete cap preparation, jobs 350 through 355

Jobs 350–354 are terminal. Single-neighbor expansions fail for the two old
unresolved stars. Joint nine-vertex rings solve both. The first complete
copy, job 352, passes incidence and periodicity but exposes roundoff-scale
cells from two earlier cones. Do not use that copy as an MMG input.

A relative determinant guard rejects new tets at most 1e-12 of the original
cavity's summed determinant. An actual-source regression catches the issue.
Job 353 rebuilds 20 safe cones and four full rings. All four binary models
solve with 102 candidates, exact signed-face replay, and positive volume.
Job 354 independently applies and audits all 24 plans to a fresh copy:

- 30,271 points, 101,472 tets, original points and wall triangles unchanged.
- All oriented facets and periodic cap triangles pass.
- Sampled distance 0.0016328718, relative volume error +0.71897165%.
- Relative volume change is zero and cap ceiling minimum is 0.1053009971.
- Minimum replacement relative determinant is 0.0015171466.
- Actual minimum corner quality is still 0.0010950575, so it is not admitted.

Prepared source is `runs/g3-gyroid-cap-replacement-v2`. Its plan source is
`runs/g3-gyroid-cap-cavity-rankguard-v1`. No native code, benchmark references,
or specimen definitions change. Preparation uses an explicitly non-default
Microgen background triangulation and cap retriangulation, not a Meshers seed.

Job 355 invokes actual `remesh_keeping_boundaries_for_fem` on the saved prepared
seed, with recorded first-pass overrides `nofem=True`, hausd 0.0005, hmax 0.06,
hgrad 1.3. That first pass returns 649,230 tets and stops at the unchanged
250,000 intermediate cap. No second MMG call or FE runs. Its 58.02-second
saved-seed diagnostic is not a full-workflow timing or admission.

The selected next control changes only first-pass hausd to 0.0015, retaining
the actual public second call and hausd 0.0005 there. Inspect its live queue
and results before choosing any successor. All stages retain 8 GiB, one
thread, and 600-second MMG worker limits. Earlier sections below are historical.

## Latest cap-topology findings, jobs 345 through 349

All jobs are terminal. Job 345 is a preserved JSON serialization failure,
corrected with Python integer IDs and an actual roundtrip regression in 346.
The five-tet Gyroid seed remains unchanged. No MMG or FE call runs.

Cap-flip preflight 346 finds 24 passing local planar alternatives among 48
candidates. All 24 obstructing faces are covered. Patch boundaries and areas
are unchanged. Two-tet volume-star check 347 finds 20 passing cavities and
four three-tet stars. Extended check 348 finds 22 valid disjoint plans, but
two three-tet stars have no valid cone from their existing non-patch vertices.
No partial replacement is applied.

Independent finite-vertex enumeration 349 checks 13 nondegenerate oriented
tet candidates and 8,191 nonempty subsets on each of those six-vertex stars.
Neither has a target signed-face solution. Positive controls reproduce the
original cavity and a changed-boundary convex cavity. A missing boundary
face is rejected. This is finite-precision local evidence, not an obstruction
to all remeshing, larger cavities, or Steiner points.

Next assess bounded one-ring cavity expansion around the two unresolved
stars. Keep target cap polygons and periodic boundaries fixed. Do not repeat
the same cones or enumeration and do not send blocked caps to MMG. The source
mesh is in `runs/g3-gyroid-preclip31-five-v2/pretriangulated`. Latest results are
`runs/g3-gyroid-cap-flip-preflight-v2`,
`runs/g3-gyroid-cap-cavity-preflight-v2`, and
`runs/g3-gyroid-cap-cavity-enumeration-v1` under `experiments/adaptive`.
All workers stay below 0.294 GiB with one thread and 120-second stage limits.
No native code, references, specimen definitions, or admission gates change.
The four-fixture goal remains active and unachieved. Earlier sections below
are historical.

## Latest extraction-order findings, jobs 342 through 344

All three jobs are terminal. No MMG or FE solve runs in this phase. Normal
Gyroid31 scalar clipping produces 288 same-orientation shared faces.
Triangulating the background before clipping removes all such conflicts in
both tested decompositions, with identical original coordinate and scalar
hashes. Both pass periodic cap triangulations, but each introduces 24 fixed-cap
quality obstructions, minimum ceiling 0.08710052. Neither is admitted.

The six-tet version has 34,743 nodes, 123,344 tets, and sampled distance
0.00253049. The five-tet version has 30,271 nodes, 101,472 tets, and sampled
distance 0.00163287. Volume errors are +0.719315% and +0.718972%.
Five-tet background count is 135,000; six-tet count is 162,000. All material
and background limits remain 250,000, with 8 GiB and one thread per worker.

Job 343 is a preserved implementation failure at an exact-input-order guard,
not a scientific result for five-tet extraction. Microgen structured order
differs from VTK rectilinear order. Retry 344 explicitly maps connectivity
back to original point order and restores exact coordinates before clipping.
Seven controls and sixteen combined focused tests pass. Original driver and
test snapshots match their declared hashes. Both completed paired-run progress
audits pass. No native code, benchmark references, or specimen changes occur.

Do not send either blocked seed to MMG. Next inspect local cap-triangulation
feasibility on the coherent five-tet seed, preserving periodic cap boundaries,
before any volume-remeshing or FE experiment. No such cap trial has run.
Frozen records are `g3-gyroid-preclip31-v1`,
`g3-gyroid-preclip31-five-v1`, and `g3-gyroid-preclip31-five-v2` under
`experiments/adaptive/runs`. Earlier sections below are historical.

## Latest source-validity findings, jobs 336 through 341

All are terminal; preflight 337 fails at its correct source guard. No job is
queued. Do not restart observing sessions 54541, 89586, 56535, 39179, 10787
or 13041. Recheck the shared queue before work. The full four-fixture goal
remains active and unachieved. Native code and specimen definitions did not
change in this phase.

The worst remaining star at node 704 has a checked negative tangent bound at
quality 0.1. The normalized bound is -0.0017874161, independently reproduced
with 64-bit-mantissa arithmetic within 1.3e-16. Three supporting tets use only
boundary vertices opposite the moving point. The subsequent source audit is
essential to interpreting that algebraic obstruction.

The 2-to-3 preflight 337 refuses a non-oriented source cavity. Cells 1297 and
1298 share `[704, 9126, 9369]`, have positive determinants and same-side
opposite vertices. No flips are applied. Do not bypass the guard or treat the
failed run as a successful topology result.

Lineage audit 338 finds 684 same-orientation shared faces in the saved VTK43
Gyroid seed, before MMG, and 146 in the one-pass FEM stage, final output and
partial interior repair. All have positive corner determinants and matching
saved boundary arrays. The input is invalid under this necessary incidence
check. Retire that repair source. The two-pass Gyroid has zero conflicts but
still fails quality. The audit does not yet isolate which extraction step
introduces the seed conflicts.

Representative admitted controls pass. Replay 339 covers all 24 archived
Schwarz timing-arm meshes by exact identity to audited controls and checks
all six manual lifts preserve parent corners/surfaces. Earlier Schwarz timing
conclusions are unaffected. Necessary corner-incidence consistency is not
a certificate against curved or nonlocal intersections.

The common `microgen_workflow.geometry` admission now requires oriented facet
consistency, uniformly for future outputs. Numerical H, volume, quality,
positive-Jacobian, projection and periodic thresholds are unchanged. Check 340
passes sixteen controls, reproduces every legacy metric on bad-Gyroid and
good-Schwarz records, and retains/rejects the appropriate controls. Its old
study verifier still passes 96 comparisons and 205 hashes, using exact driver
archives and declared measurement Git ref `2f0363f` for thirteen source records.
Frozen flags, tensors and timings are not rewritten.

Check 341 verifies early source rejection in the retired interior-repair
entrypoint. Seven controls pass; status is `source-incidence-rejected`, with
no moves, new mesh or FEA. Archive snapshots added for earlier checker/control
versions match their already recorded SHA256 exactly. The original failed
preflight, partial repair and witness outputs remain intact.

The final focused suite passes 69 tests; targeted Ruff and whitespace checks
also pass. No native rebuild is needed for these Python research changes.

Next queue one paired extraction preflight at resolution 31: unmodified
Microgen scalar clipping versus the same sampled background triangulated
before its two scalar clips. Keep offset 1.1, analytic field, sample coordinates
and scalar arrays identical. Explicitly label the variant non-default. Bound
background and returned material tets at 250,000, and read incidence, cap
feasibility, periodicity and independent geometry before selecting MMG or FEA.
No such extraction trial has run yet. No new resolution or sizing sweep.

Detailed records are under `experiments/adaptive/runs`:

- `g3-gyroid-convex-star-feasibility-v1` and its witness audit v1/v2 JSON files.
- `g3-gyroid-two-to-three-preflight-v1`, the failed source-guard run.
- `g3-facet-incidence-lineage-controls-v1`.
- `g3-schwarz-admitted-facet-replay-v1.json`.
- `g3-facet-admission-check-v1`.
- `g3-gyroid-invalid-source-rejection-v1`.

Read `MICROGEN_RESULTS.md`, `PLAN.md`, `RESULTS.md`, `PROJECT_TRACKING.md` and
the append-only `.audit/g3-microgen.tsv`. Earlier sections below are historical.

## Latest completed diagnostics, jobs 334 and 335

Both are terminal success; scientific admission is not. No job is queued.
Do not restart observing sessions 79681 or 32411. Check the shared queue
before work. The full four-fixture goal remains active and unachieved.

The progress transition issue below is now fixed for selected controllers
with `job_parts.worker_part`, using public standard Job methods. Success and
failure capture tests are red before the fix and green afterward. Completion
retains elapsed time and the 8 GiB ceiling. Both queued checks pass, with six
live reports in job 335. No shared workbench helper or frozen output changes.

Existing audit 212 already rules out low cap/wall-face ceilings on the Gyroid
outputs. Job 334 instead tests one worst free interior vertex, node 30242,
on the original two-pass mesh. All 26 incident tets reach strict quality 0.1
with a 1e-9 positive guard, while boundary and signed volume stay unchanged.
Global quality still fails. That mesh would lift to 325,911 Tet10 nodes,
above the stage. The smaller one-pass FEM43 precursor lifts to 153,010.

Job 335 makes one explicitly added SciPy interior-relocation attempt on that
smaller saved genuine Microgen mesh. Two sweeps, 400-attempt cap, complete
star target 0.100001, no boundary or topology changes. Six controls pass.
It accepts 108 moves in sweep one and none in sweep two. Bad tets fall from
222 to 49, but minimum quality remains 0.0079758. All determinants stay
positive, volume change is zero, and original H=0.00151229 and periodic caps
reproduce. The route remains `repair-incomplete`, `qualified=false` and no
FEA runs. This is not a Microgen default or a Meshers seed. Do not call it
an admitted or matched-accuracy competitor.

Read `runs/g3-gyroid-interior-star-feasibility-v1` and
`runs/g3-gyroid-firstpass-fem43-interior-relocation-v1` under
`experiments/adaptive`. The maintained drivers are `interior_star_probe.py`
and `interior_relocation_trial.py`. Stop this sweep recipe. Next diagnose one
remaining worst star using an independent convex-slack feasibility problem
and checked supporting-plane bound. A failed SLSQP search alone is not an
infeasibility proof. Read that diagnostic before coupled or topology changes.

Schwarz's repeated comparison remains valid. Verification v2 still checks
96 comparisons and 205 hashes after source changes, using exact measured
driver archives for seven source records. Changed mesh/field data cannot use
that fallback. Read
`runs/g3-schwarz-tet10-live-amg-repeated-verification-v2.json`.
Eighteen focused tests pass. Native code and specimen definitions are unchanged.
The sections below retain earlier checkpoint reasoning as history.

## Latest completed experiment, job 333

Job 333 is terminal success. Do not restart its observing session 95526 or
completed jobs 331 and 332. No successor is queued. Recheck the authoritative
shared queue before new work. Its output is
`experiments/adaptive/runs/g3-schwarz-tet10-live-amg-repeated-v1`.

One fresh warmup per arm and all five alternating measured pairs pass 38
controls each, geometry, periodic, field and eight-model supplementary 0.1%
gates. All 72 loads are fresh. Common explicit AMG, recipes and admission
limits match the pilot; no saved mesh, load or previous-run hierarchy replay.

Median audited times are 409.2799 seconds for actual VTK30/one-pass MMG plus
manual straight Tet10/direct fedoo and 55.1914 seconds for native curved16.
Every paired workflow ratio favors Meshers, 7.3858–7.5401. Phase-only ratios
are 7.2806–7.4331. Microgen has the smaller model disagreement. This is a
common-target result, not identical measured errors or the default pipeline.

The experiment takes 2,923.13 seconds within its 3,600-second budget. Peaks
are 3.8143 and 1.4571 GiB. Clock discrepancy is below 8.1e-7 seconds. No
host interruption or resource gate fails. All measurement input hashes pass.

The new `verify_schwarz_tet10_repeated.py` independently recomputes 96 full
mixed-load comparisons and every timing statistic, verifies 205 hashes,
and reads every geometry/field gate. Its two controls include an off-diagonal
error invisible to diagonals. Read its output
`experiments/adaptive/runs/g3-schwarz-tet10-live-amg-repeated-verification-v1.json`.

One operational issue remains before another queued experiment. The initial
progress audit passes, but a late `no_chips` warning appears after 16 minutes.
A direct standard-Job capture proves `Part.__exit__()` ends the part with
empty values after populated worker reports. Fix stage-transition reporting
using standard Job APIs before any successor. Preserve the warning; no shared
workbench helper was changed. Scientific results and timings remain valid.
The warning capture is in
`experiments/adaptive/runs/g3-schwarz-tet10-progress-warning-v1.json`.

After that reporting fix, return to missing admitted actual Gyroid/spinodoid
comparators and lattice reference uncertainty. Do not tighten Schwarz again
from this result alone. The full four-fixture goal remains active and unmet.

## Latest live comparison, jobs 331 and 332

After closing the polynomial diagnostic, one fresh stronger Schwarz Tet10
pair now completes. Both arms use the same explicit AMG preset. Actual
Microgen VTK30/one-pass MMG at hausd .0005, hmax .1, hgrad 1.3 is followed by
manual straight elevation and direct fedoo, not public Microgen Tet10 support.
The competing native mesh is nonadaptive curved16. No saved meshes, solved
load prefixes, or previous-run preconditioners are replayed.

Both pass the supplementary 0.1% target across all eight Schwarz models.
Microgen has 207,200 nodes and 0.053985–0.062716% disagreement, native has
42,928 nodes and 0.031469–0.093778%. Both pass unchanged geometry and periodic
gates. Fresh field reaction errors are below 7.60e-12 and jumps 5.56e-17.

Science-audited total times are 409.87 versus 55.10 seconds, ratio 7.4382.
Generation/lift/solve subtotals are 388.25 versus 52.92 seconds, ratio 7.3368.
Checkpoint I/O within those phases remains included, despite the frozen v1
subtotal label saying serialization is excluded. Frozen values are unchanged;
the maintained runner corrects that wording. Microgen is more accurate against
the models. These are common-target, single-pair times without warmups, not
identical measured errors or a repeated/default/four-fixture speed claim.

After reading the pair, positive-rule audit 332 checks fresh native fields.
The 125-point energy shift is 2.89e-9 and 64-to-125 change 5.81e-10, below
1e-4. The Microgen four-point affine rule already agrees with fifteen-point
reintegration within 1.24e-11. Worker peaks are 3.8152 and 1.4316 GiB.

Read `runs/g3-schwarz-tet10-live-amg-pilot-v1` and
`runs/g3-schwarz-tet10-live-amg-pilot-quadrature-v1`, under
`experiments/adaptive`. The maintained runner is `schwarz_tet10_live.py`.
Its 33 combined controls pass. A separate artifact-hash record preserves the
pilot's generated meshes and fields without rewriting its frozen summary.

The selected next experiment is interleaved repetitions with these exact
recipes, common AMG and both timing definitions. Do not claim repeated
stronger-Tet10 performance before measuring it. No successor is queued here.
Lattice reference uncertainty and missing admitted actual MMG comparators
for Gyroid and spinodoid remain unresolved. The full goal remains active.
The earlier sections below are historical.

## Latest completed diagnostic, jobs 329 and 330

The corner-coupling check proposed later in this checkpoint is now complete.
Five algebraic tests and manufactured periodic controls pass. On the frozen
cubic24 lattice, 102 edge orbits and 96 cubic corner orbits give 594 DOFs on
1,408 tets. Edge-only gains reproduce exactly. The incremental quartic gain
after corner relaxation agrees with the Schur complement within 5.37e-16.

The incremental worst mixed-load gain is 0.676219 ppm, versus the unchanged
500 ppm screen. Its 125/216-point change is 0.018333%, passing precision.
Trial-energy identity error is 3.86e-13 and periodic trace mismatch 3.39e-21.
Worker time is 263.45 seconds, peak RSS 0.6863 GiB. The 300-second stage
and 8 GiB cap remain unchanged. All source hashes reproduce after the run.

Stop this polynomial-mode route. No higher degree, mode-strength tuning,
new equilibrium reference, or native feature is selected. All 24 models
remain retained and the four-fixture objective remains unachieved. No
successor is queued by this checkpoint. Recheck the shared queue before work.

Read `runs/g3-quartic-corner-controls-v1` and
`runs/g3-lattice-quartic-corner-witness-v1` under `experiments/adaptive`.
The remaining sections retain the earlier reasoning and proposal as history.

The four-fixture matched-accuracy goal is active and unachieved. The source
and documentation changes through `b186d39` are pushed. This checkpoint adds
the controlled quartic edge witness and its resolved weak signal.
No delegates, unrelated jobs, timers, or remote hosts are involved.

## Current execution state

Jobs 322 through 328 all completed. No Meshers successor is
queued at this checkpoint. Recheck the authoritative shared queue before
running anything. Do not restart cubic23 or treat its old observing terminal
session 71321 as a current worker.

The completed cubic23 job uses `cubic_reference.py`, source
`runs/g3-lattice-cubic-parent23-fields-v1`, prepared mesh
`runs/g3-lattice-cubic23-assembly-memory-v1`, and the validated elastic-block
AMG preset. No saved load prefix is replayed. All paths beginning with `runs`
are relative to `experiments/adaptive`.

The source parent is independently admitted, has six fresh solved loads, and
passes field audit 320. Actual cubic assembly 321 has 218,508 nodes and
43,844 elements and peaks at 5.014 GiB. The worker retains 8 GiB, the existing
250,000-node and 250,000-tet stages, a 240-second first-load deadline, and a
1,800-second total worker timeout. One numerical thread is used.

Its sampled geometry comes from the admitted quadratic parent through an
exact cubic embedding, not a projected cubic boundary. The next selected
question is now one coarse-corner coupling check on the same quartic curve-edge
witness. The edge-only controls and lattice measurements are complete; the
joint coarse-corner check has not been implemented or run.

## What changed the next action

The eight-round curve pilot improves size realization and single-run cost but
still fails the fixed model envelope at 1.439306%. Do not tune its strength
or round count again without new evidence.

A separate projected cubic boundary probe passes geometry and improves
sampled distance fourfold. With saved displacement coefficients fixed, energy
changes by only 0.0105922%, below the predeclared 0.05% signal screen. It is
not an equilibrium solution or a bound on total geometry error. No projected
domain solve was selected. Records are in
`runs/g3-lattice-cubic-geometry-sensitivity-v1`.

The actual cubic25 assembly fits, but first solve 314 fails before completing
a load because full symmetry subtraction requests a 1.47 GiB temporary beyond
the 8 GiB address-space ceiling. The AMG helper now checks the complete
relative Frobenius symmetry norm in row blocks. Its 1e-12 threshold and
hash-guarded reuse are unchanged. Ten focused tests and twelve integrated
control solves pass. The failed run stays preserved.

Fresh retry 316 completes all six cubic25 loads without replay. Worker time
is 1,071.18 seconds, peak RSS is 6.753 GiB, and own-rule error is 1.03e-11
relative. Periodic jumps are below 2.23e-16. Positive 125-point energy shift
is 2.6354e-5, below the unchanged 1e-4 sensitivity gate. Fresh-process audit
317 passes, including 64-to-125-point change 6.8313e-6. All fields remain local.

## Reference uncertainty is still real

Cubic24 and cubic25 differ by 0.248177%. Cubic25 is stiffer in every
generalized mode. Relative to cubic25, the finest native adaptive tensor
differs by 0.651535%, the eight-round curve tensor by 0.483742%, and the
strongest actual Microgen-derived Tet10 tensor by 0.917416%.

Do not use that single favorable reference to claim a timing victory. All
24 earlier models remain retained. Cubic25 itself has a model-normalized
envelope of 0.114924% to 1.920819%, so it does not qualify as a target-passing
replacement. Neither two nonnested levels nor a higher polynomial degree
proves continuum convergence. The third level is intended to test that
uncertainty, not guarantee that an extrapolation will pass.

Read `CUBIC_RESULTS.md`, `PLAN.md`, `RESULTS.md`, `PROJECT_TRACKING.md`, and
the append-only `.audit/g3-microgen.tsv`. Detailed new records are in:

- `runs/g3-lattice-cubic25-reference-blocked-symmetry-v2`.
- `runs/g3-lattice-cubic25-assessment-v1`.
- `runs/g3-lattice-cubic25-reference-envelope-v1`.
- `runs/g3-lattice-cubic-parent23-preflight-v1`.
- `runs/g3-lattice-cubic-parent23-fields-v1`.
- `runs/g3-lattice-cubic-parent23-field-audit-v1`.
- `runs/g3-lattice-cubic23-assembly-memory-v1`.

## Completed third-level result

Cubic23 completes all six fresh loads in a 1,018.27-second worker, peaking
at 5.976 GiB. Its own-rule error is 1.29e-11 and periodic jumps are below
2.23e-16. Positive 125-point energy shift is 7.6754e-6, and the minimum
parent-minus-cubic eigenvalue is +0.000137181. Fresh-process audit 323 passes,
including 64-to-125-point change 1.1212e-6. These three eight-model reports
remain retained:

- `runs/g3-lattice-curve-grading-reference-comparison-v1.json`.
- `runs/g3-lattice-reference-sensitivity-v2.json`.
- `runs/g3-lattice-microgen-reference-comparison-v1.json`.

Before reading the full cubic23 tensor, `cubic_refinement.py` defines a
common leading-matrix model, order range 0.25 to 6, at most 10% relative
matrix-fit residual, and signal at least ten times observed quadrature shifts.
Job 324 passes seven synthetic controls, verifies all field hashes and parent
provenance, and checks exact equality of the retained 24 model tensors.

Both length proxies reject the model. Correction-direction cosine is
-0.959707 and the best proportionality coefficient is negative. Matrix-fit
residual is 28.1%. Cubic24 is softer than both neighboring levels. Do not fit
only the closer odd pair, drop the middle level, promote an observed order,
or select a favorable reference for a timing claim. This rejects the stated
local model, not every possible convergence mechanism.

Records are in `runs/g3-lattice-cubic23-reference-v1`,
`runs/g3-lattice-cubic23-assessment-v1`, and
`runs/g3-lattice-cubic-three-level-screen-v1`.

## Completed edge-only diagnostic

Read the measured protocol and limits in `CUBIC_RESULTS.md`. The witness uses
the saved cubic24 fields in `runs/g3-lattice-cubic-reference-cache-release-v1`
and the prepared mesh in `runs/g3-lattice-cubic-assembly-memory-v1`. Freeze
their geometry, topology, and source fields.

The scalar mode on each periodic wall/box-curve edge orbit is
`16*lambda_i*lambda_j*(lambda_i-lambda_j)**2`. Give each orbit one shared
three-component displacement coefficient, supported on the complete incident
element star. Full overlap terms are included. Three unit tests and analytic
periodic controls pass, including nonzero quartic face traces, explicit
off-orbit coupling, and a nonzero pin-gauge test. No regularization is added.

The actual lattice has 102 edge orbits, 306 mode DOFs, and 332 support tets.
Job 327 passes energy and trace identities but fails the 64/125-point
precision gate at 6.0323% change. Preserve that rejected run. One higher-rule
125/216-point check, job 328, passes at 0.017934% change. Its worst mixed-load
trial-energy reduction is only 6.4791e-7 relative, or 0.648 ppm, versus the
500 ppm useful-signal screen. The direct energy identity agrees within
2.04e-13; continuity, periodic quartic traces, and the actual pin offset are
zero. Worker time is 26.57 seconds and peak RSS is 0.3404 GiB.

This resolved small witness does not justify a native sizing feature. It is
not an equilibrium tensor or a bound on all FE error: the original coarse
coefficients were fixed. Inputs, geometry, and all 24 models remain unchanged.
Records are in `runs/g3-quartic-edge-controls-v1`,
`runs/g3-quartic-edge-controls-v2`,
`runs/g3-lattice-quartic-edge-witness-v1`, and
`runs/g3-lattice-quartic-edge-witness-higher-v2`.

## Next single coarse-coupling check

Before closing this polynomial route, allow only the existing cubic corner
basis functions at periodic curve vertices to relax with the same quartic
edge modes. Include complete vertex stars and all cross terms. The earlier
Tet10-to-cubic audit covered 1,408 curve-vertex-incident tets versus 332 for
the edge subset; this motivates a test, not a claim about current cubic error.
Keep geometry and source fields fixed. No new mode degree or grid is selected.

Report corner-only, edge-only, and joint trial gains. Use the incremental
quartic gain after coarse relaxation for the signal, not the joint total.
Verify it by the Schur complement, and check nonnegative gain increments
under space enlargement. Reproduce the original edge-only Gram before the
comparison. Any nonzero source-pin mode value needs an explicit uniform
translation gauge, which changes no strains or periodic jumps.

Retain the validated 125/216-point pair, 1% precision gate, 0.05% signal
screen, 8 GiB ceiling, at most 512 scalar modes and 5,000 support tets.
Run controlled coupled-space checks before a lattice measurement. If the
increment remains small, stop this polynomial route rather than escalating
degree or strength. Trial energy is still not globally equilibrated or a new
homogenized reference. No such joint lattice check has run yet.

## Remaining objective and runtime

Schwarz still has the only narrow repeated matched-target workflow win.
Gyroid and spinodoid still lack admitted actual MMG comparators. Lattice
still lacks a justified target-passing comparison across the retained models.
Local stress and continuum certification remain outside the demonstrated scope.

Use the pinned `.venv`, Rust 1.95.0 if rebuilding, and
`LD_LIBRARY_PATH=/home/kmarchais/dev/meshers/.local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu`.
The native binary SHA256 remains
`908e326211792142f83f854ad2696cd879d927bd81deacf38eb6a0cb63bc4dc1`.
Large NPZ archives and native binaries are ignored local artifacts, not backed
up by the committed text. Heavy work remains in the shared `gpu` queue through
`qrun --check-start` and standard `Job` progress. Before GitHub operations,
switch to and verify `kmarchais`. Do not touch unrelated queue entries.
