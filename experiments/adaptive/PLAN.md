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
