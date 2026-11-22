---
icon: lucide/list-checks
---

# Development benchmarks

Meshers development is benchmark-based. AI coding agents use these suites to add support for implicit geometries and improve algorithms. Compare a pinned base and candidate on the same solid, hardware and numerical settings.

## What counts as an improvement?

| Goal | Evidence to report |
| --- | --- |
| Shorter execution time | Measure generation, assembly and solve time separately. Report peak memory and compare repeated runs. |
| Better mesh quality | Check orientation and conformity. Report the minimum and lower tail of tetrahedron shape quality and the dihedral-angle range. |
| More accurate boundaries and volumes | Measure mesh-to-wall **and** wall-to-mesh distance, including caps. Compare signed volume by component or material with an independent reference. Mark sampled checks. |
| Better FEA convergence | Refine the same physical problem. Plot displacement, energy or stiffness error against mesh size or free degrees of freedom. |
| Fewer nodes and elements at comparable quality | Count both on meshes that meet the **same** geometry, shape and FEA error limits. Record rejections and adaptive fallbacks alongside successes. |

## Automated gates

| When and what runs | Cases | What the result establishes |
| --- | --- | --- |
| **PR performance** · Linux · [runner](https://github.com/kmarchais/meshers/blob/fc54032b/tools/performance_gate.py) | Seven 24-cell fixtures: compiled and callback gyroid, graded sheet, mapped intersection, Schwarz P, Fourier spinodoid and smooth lattice. One warm-up pair, five alternating measured pairs. | Generation time, peak memory, nodes, tets, shape distribution, sampled wall/cap error and volume change against base. Volume agreement with base alone cannot establish physical accuracy. [Limits and replay](performance.md#pull-request-gates) |
| **Release performance** · manually selected · [fixtures](https://github.com/kmarchais/meshers/blob/fc54032b/tools/performance_worker.py) | All seven PR fixtures plus adaptive gyroid, spinodoid and lattice, and a quadratic gyroid: **11 total**. | Adds curved-element volume and Jacobian checks. The Tet10 shape and wall-error metrics describe the linear parent. [Archived runs](benchmark-history.md) |
| **Independent geometry admission** · Python 3.11 on Linux, Windows and macOS · [CI recipes](https://github.com/kmarchais/meshers/blob/fc54032b/.github/workflows/ci.yml) | Gyroid, Schwarz P, Neovius, graded gyroid, thin gyroid and BCC lattice. Eight cells per axis, with 0, 4 and 12 optimization passes: **18 recipes per OS**. | Independent topology, normals, two-sided wall distance ≤ 0.01, minimum shape quality ≥ 0.1, volume error plus reference uncertainty ≤ 10%. The separate seven-solid study includes the seeded spinodoid. [Validator and replay](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/implicit-coverage/README.md) |
| **Candidate-wheel FEA** · same three platforms · [runner](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/schwarz_workflow.py) | One adaptive periodic Schwarz P Tet4 mesh, reopened from VTKHDF. Six strain and six stress loads. | Sampled wall/cap distance ≤ 0.002, volume within 1%, minimum shape quality ≥ 0.1 and worst relative energy discrepancy ≤ 1% against the frozen stiffness tensor. [Replay](performance.md#candidate-wheel-fea-gate) |

## Geometry and solver studies

| Study and replay | Coverage | Evidence |
| --- | --- | --- |
| [Seven-solid coverage and cap witness](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/implicit-coverage/README.md) | The six CI geometries plus a seeded, 52-mode spinodoid. Compares uniform and adaptive recipes; bounds the distance at a failing cap sample. | Spinodoid passes at 20 cells, misses the shape floor at 24, and the eight-cell target remains open. [Per-attempt records](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/implicit-coverage/followup-2026-11-18.json) |
| [Tet4 displacement and energy convergence](solver-verification.md) | Cube at 6/9/12 cells per axis; cylindrical and toroidal gyroids at 12/18/24. Independent P1 assembly and Fedoo, with affine patch controls. | Error versus maximum physical edge length on nine meshes. Final displacement orders 1.92–2.01 and energy orders 1.00–1.05. [Raw errors and rates](https://github.com/kmarchais/meshers/blob/fc54032b/docs/solver-convergence-results.json) |
| [Manufactured elasticity on implicit solids](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/implicit-coverage/README.md#refinement-at-the-2-manufactured-energy-target) | Uniform gyroid and adaptive BCC. Six affine patch loads plus the worst combination of 18 quadratic displacement fields after geometry admission. | Best recorded worst energy errors were 5.78% and 10.49%, respectively. Both miss the 2% target. [Refinement records](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/implicit-coverage/followup-2026-11-18.json) |
| [Periodic homogenization](periodic-homogenization.md) | Solid cube material control and periodic gyroid refinement. Six strain and six stress solves per mesh. | Effective stiffness, symmetry, positive stiffness and Hill–Mandel energy balance. Refinement stability measures tensor change; continuum error remains unknown. [60 solve records](https://github.com/kmarchais/meshers/blob/fc54032b/docs/fedoo-homogenization-results.json) |
| [Frozen 0.2.1 wheel replay](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/release-wheel/README.md) | Uniform and adaptive gyroid, adaptive Schwarz P. Three repeats per case on five hardware/OS configurations. | Exported volumes, node/element counts, solver residuals, free DOFs and load responses. [Two-sided geometry audit](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/release-wheel/geometry-admission.md) · [Configuration records](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/release-wheel/configurations-20261118/README.md) |

## Cost comparisons and profiling

| Question | Benchmark | Scope |
| --- | --- | --- |
| Which Microgen or Meshers route meets the application limits? | [Historical ten-case harness at its runnable revision](https://github.com/kmarchais/meshers/blob/602dc900/benchmarks/microgen/README.md) | Includes Diamond, Neovius, scaled and thin gyroids, lattice and spinodoid. The standard sweep uses 24/36/54/80 intervals, one warm-up and three repeats. The eight-cell smoke run checks coverage. |
| Which step first loses the Schwarz P wall or misses stiffness? | [Three-route refinement series](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/release-wheel/schwarz-series.md) · [Stage diagnostics and repair controls](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/release-wheel/schwarz-diagnostics-20261118/README.md) | Native Tet4, default and tuned Microgen/MMG. Audit the seed, intermediate passes, repair and exported mesh. A rejected route has no qualifying solve time. |
| How much does reaching the same 1% stiffness target cost? | [Schwarz P matched-cost trials and raw runs](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/release-wheel/schwarz-matched-cost-20261118/README.md) | Three alternating timed pairs per configuration. On M4: native **19.7 s, 133,730 tets, 30,070 nodes**; repaired **29.1 s, 154,002 tets, 32,099 nodes**. Both qualify against the pinned reference. The errors differ. |
| Which local operation consumes time? | [Adaptive gyroid volume](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/examples/adapt_gyroid.rs) · [Experimental triangle surface](https://github.com/kmarchais/meshers/blob/fc54032b/crates/meshers-core/examples/adapt_surface.rs) | Per-phase Rust timings for one gyroid. The surface example reports triangle counts and sampled error as well as extraction, split, collapse, flip and relaxation times. Use paired runs for a speed claim. |

### Suites on separate revisions

- [Periodic stress suite at `b721a33f`](https://github.com/kmarchais/meshers/blob/b721a33f/demos/solver/periodic-bench/README.md) covers 17 fields, 33 settings and two requested routes: **66 cases**. Score cap matching, wall distances, volume, dihedrals and shape; distinguish a returned mesh from one that passed on the requested route. [Observed failures](benchmark-gaps.md#periodicity-candidate)
- [Gate calibration at `fe3bbf31`](https://github.com/kmarchais/meshers/blob/fe3bbf31/benchmarks/gate-calibration-20261121/README.md) measures same-build noise and tests detection of injected time and memory regressions. It evaluates the gate itself.

## Benchmarks to add

These are proposed studies. The [roadmap](roadmap.md) sets development priorities. [Open targets](benchmark-gaps.md) records failures already demonstrated on pinned inputs.

| Priority | Gap in current evidence | Experiment |
| --- | --- | --- |
| 1 | A lost cavity can leave the material component count unchanged. One witness gives little evidence about nearby thicknesses or positions. | Sweep cavity radius, film thickness and gap width across grid offsets. Probe material and void, count retained features, compare analytic volume and both wall distances. Record a lost feature returned as valid separately from an explicit rejection. |
| 2 | One frozen spinodoid cannot establish accuracy or cost across realizations. The release timing case has no independent stiffness target. | Reserve spectra and phases outside the tuning set. Compare uniform and requested adaptive routes against the same geometry, shape and converged stiffness budgets. Report strict passes, fallbacks, time, nodes, tets and free DOFs. |
| 3 | The release suite checks Tet10 orientation, but its wall error measures the parent. | On a sphere and a mapped cap, sample the **curved** exported faces in both directions, check volume, then compare Tet4 and Tet10 at matched FEA error. Retain rejected or uncertified fits. |
| 4 | The triangle phase timer and historical surface comparisons leave direct quad extraction without a recurring, matched-accuracy performance gate. | Compare triangles and experimental quads on the same implicit wall and labelled crease. Require boundary accuracy, manifold topology and valid face orientation before comparing vertex/face counts, time and memory. |

### First controls for planned features

| Feature | Suggested first comparison |
| --- | --- |
| Solver-driven adaptive refinement | Compare error-driven and uniform refinement on a problem with a localized stress gradient. Measure error versus free DOFs and total solve–adapt time. Verify transfer of fields, loads and boundary labels. |
| CAD/BRep input | Mesh a plate with a cylindrical hole. Recover the analytic radius and volume, preserve the named bore surface, and reapply its load after export. |
| New volume cell families | Begin with an affine patch test on a hex mesh, then a curved geometry at the same physical error as Tet4/Tet10. For mixed meshes, check conformity and load transfer across cell-type transitions. |

<details markdown="block"><summary>Other FEA capabilities worth benchmarking</summary>

- **Material interfaces:** mesh two bonded phases with shared interface faces. Verify each phase’s volume and material ID, then recover the analytic stiffness of a layered specimen.
- **Directional sizing and hybrid cells:** resolve a thin layer through its thickness with prisms/wedges or swept hexes. Compare DOFs at equal error. Test pyramid transitions where quad and triangular faces meet.
- **Discrete inputs:** replay a labelled triangle surface or segmented voxel volume alongside its analytic source. Measure surface deviation and feature retention after remeshing.

</details>
