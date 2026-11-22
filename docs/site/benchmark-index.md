# Development benchmarks

| If you change… | Start with | What it measures |
| --- | --- | --- |
| Meshing speed, memory or element count | [PR gate and 11-case release suite](performance.md) | Five measured base/candidate pairs. Time, peak memory, nodes, tetrahedra, shape distribution, sampled boundary error and volume. The release suite adds adaptive cases and Tet10. [Past runs](benchmark-history.md) |
| An implicit field, grading or periodic caps | [Seven-solid coverage and replay](https://github.com/kmarchais/meshers/blob/main/demos/solver/implicit-coverage/README.md) | Independent geometry admission for gyroid, Schwarz P, Neovius, graded and thin sheets, BCC lattice and seeded spinodoid. The spinodoid witness isolates cap-distance errors. CI runs selected passing recipes for six solids. |
| Elasticity convergence | [Displacement and energy convergence](solver-verification.md) · [Manufactured-energy target](benchmark-gaps.md#where-to-work-next) | Cube, curved cylinder and torus under refinement, with independent assembly and Fedoo. A separate gyroid and BCC test scores the worst energy error across 18 displacement fields. |
| Periodic pairs, export or effective stiffness | [Gyroid homogenization](periodic-homogenization.md) · [Candidate-wheel gate](performance.md#candidate-wheel-fea-gate) | Periodic strain and stress loads, patch and energy checks. CI reopens the Schwarz P export, audits its geometry and compares stiffness with a pinned reference. |

## Comparative studies

| Study | What the comparison answers |
| --- | --- |
| [Historical Microgen benchmark](https://github.com/kmarchais/meshers/blob/602dc900/benchmarks/microgen/README.md) | Which native, VTK and MMG surface or volume routes reach the declared geometry and FEA targets? The pinned checkout supplies the independent validator used in the coverage study. |
| [Frozen wheel replay](https://github.com/kmarchais/meshers/blob/main/demos/solver/release-wheel/README.md) · [Geometry audit](https://github.com/kmarchais/meshers/blob/main/demos/solver/release-wheel/geometry-admission.md) | Do the installed wheel and exported meshes reproduce boundary distances, volume and load responses? [CPU and platform results](https://github.com/kmarchais/meshers/blob/main/demos/solver/release-wheel/configurations-20261118/README.md) |
| [Schwarz P refinement series](https://github.com/kmarchais/meshers/blob/main/demos/solver/release-wheel/schwarz-series.md) · [Stage diagnostics](https://github.com/kmarchais/meshers/blob/main/demos/solver/release-wheel/schwarz-diagnostics-20261118/README.md) | Where does a route first miss geometry or stiffness tolerance? Includes the VTK seed, individual MMG passes, repairs and a curved Tet10 control. |
| [Schwarz P matched-cost trials](https://github.com/kmarchais/meshers/blob/main/demos/solver/release-wheel/schwarz-matched-cost-20261118/README.md) | Once both routes meet 1% stiffness discrepancy, compare complete workflow time and mesh size. A separate shifted specimen tests solver reordering. |

## Candidate benchmarks

- **Periodic geometry stress suite ([PR #21](https://github.com/kmarchais/meshers/blob/b721a33f/demos/solver/periodic-bench/README.md)).** Seventeen fields, 33 settings, uniform and adaptive requests. Measures both wall distances, volume, dihedral angles, cell shape and cap-node matching. Count a uniform fallback separately when testing adaptation. [Outstanding targets](benchmark-gaps.md#periodicity-candidate)
- **Performance gate calibration ([PR #22](https://github.com/kmarchais/meshers/blob/fe3bbf31/benchmarks/gate-calibration-20261121/README.md)).** Same-build repeats measure noise. Injected time and memory costs check whether the gate catches a known regression.

## Where coverage could improve

1. **Feature loss.** Pair every thin-feature success with a missing-feature control. Probe both material and void, including an enclosed cavity whose loss leaves the solid component count unchanged.
2. **Spinodoid generalization.** Keep test phases separate from tuning phases. Report geometry admission, adaptive fallback and stiffness error at a fixed target for each realization.
3. **Curved elements.** Measure distance on the exported Tet10 walls and caps before comparing their solver cost with admitted Tet4 meshes. The release gate currently reports boundary error on the linear parent.

[Open benchmark targets and pass conditions](benchmark-gaps.md)
