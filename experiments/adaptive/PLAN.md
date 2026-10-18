# Adaptive implicit meshing: plan and status

## Current G3 continuation, 2026-10-18

Active branch: `handoff/g3-microgen-20261017`. The stopping objective is a
demonstrated matched-accuracy advantage over actual Microgen/MMG, not an
isolated meshing improvement. The four-fixture objective remains unmet.

- Schwarz has a repeated 5.1155 to 5.2908 times workflow advantage against an
  explicitly modified actual one-pass Microgen route. It is not the default
  two-pass pipeline. Both representative meshes remain below 1% across the
  conservative fine-window reference models.
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
  stage. Test finer-mesh quadrature sensitivity before considering larger solves.

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
