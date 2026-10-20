# Adaptive implicit meshing: plan and status

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
