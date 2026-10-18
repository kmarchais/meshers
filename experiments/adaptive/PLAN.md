# Adaptive implicit meshing: plan and status

Goal (Kevin, 2026-10-14): surface and volume meshes of implicit geometries that are
faster than VTK/MMG, with better element quality, fewer elements, and better
geometry fidelity (volume and two-sided boundary distance).

Constraint: a long GPU job owns the `gpu` queue. Until it is idle, run only small
single-threaded checks. Element counts, quality and geometry error are deterministic
and can be measured now; timing waits for the queue.

Branch `claude/adaptive-implicit` = `implicit-microgen-ready-20261003` (surface API,
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
