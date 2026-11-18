# Meshers tasks

Start with the first open task under **Next**. Check a task when its evidence is linked here. The [project review](docs/project-review-2026-11-17.md) gives the reasons for this order.

## Completed experiments

- [x] Check mapped-cap tolerance and final element quality after fitting. Add regression cases for a translated component and a folded Tet10. [PR #9](https://github.com/kmarchais/meshers/pull/9)
- [x] Freeze and replay uniform Gyroid, adaptive Gyroid and Schwarz P Tet4 on the released wheel. Test directional sizing on an analytic cap. [PR #9](https://github.com/kmarchais/meshers/pull/9)
- [x] Run the Tet4 solver replay and all 11 release fixtures on all five CPU configurations. All 45 replays (540 linear solves) and 165 meshing samples passed. [Dated measurements](demos/solver/release-wheel/configurations-20261118/README.md)
- [x] Audit exported walls and caps in both directions, material samples and periodic faces. All 15 exported meshes passed across the five configurations. [PR #10](https://github.com/kmarchais/meshers/pull/10)
- [x] Run the three-way Schwarz P study at a 1% stiffness target. Across five configurations, all 45 native runs solved above the target and all 90 Microgen meshes failed geometry admission before FEA. [PR #11](https://github.com/kmarchais/meshers/pull/11)

## Next

- [ ] **Locate the Microgen wall error.** At resolution 24, audit the VTK seed and each MMG pass. Record the first failing stage and worst location. Replay the historical single-pass repair as a control.
- [ ] **Bracket the native 1% stiffness target.** Vary adaptive size or tolerance, measure free DOFs and local stress, and compare with curved Tet10. Record any target that remains out of reach within the 250,000-node and element limits.
- [ ] **Compare cost at matched accuracy.** Once both routes qualify, repeat on the same host and include preparation, repair, assembly, solve and peak memory. Publish the cheapest qualifying setting for each route.
- [ ] **Test one optimization on an unseen specimen.** Choose it from the measured bottleneck, freeze the settings, then replay on a held-out phase, thickness or spinodoid seed before changing a default.

## After the comparison

- [ ] Cover thin or missing components, empty Tet10 output, interior folds in the runtime validator, damaged periodic pairs and stale boundary labels after export.
- [ ] Add an acceptance test for algorithm changes alongside the frozen release replay. Keep geometry, quality and FEA gates while allowing connectivity to change.
- [ ] Benchmark sensitivity to phase, seed and thickness. Measure thread scaling with the geometry algorithm held fixed. Include cold starts, field evaluation, export memory, failures and cancellations.
- [ ] Remesh an imported Tet4 model twice. Verify material volumes, boundary labels, periodic pairs and field transfer after refinement and coarsening.
- [ ] Check a nonperiodic graded solid under traction, then compare triangle and quad shell meshes at matched bending error.
- [ ] Mesh a trimmed CAD solid with a hole or fillet and verify feature and face labels. Compare swept Hex8, Tet4 and curved Tet10 on a beam or sector at matched FEA error.
