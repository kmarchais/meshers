---
icon: lucide/chart-no-axes-combined
---

# Decide whether the mesh solves the physical problem

A refinement study needs a quantity with a known answer. Here the displacement is prescribed analytically. The corresponding body force loads linear elastic Tet4 elements on three successively finer meshes.

<div class="plot-pair">
  <figure><img src="assets/playground/displacement.svg" alt="Displacement error on the cube and mapped gyroids falls approximately as h squared."><figcaption>Displacement L2 error. The finest cube refinement has observed order 2.01.</figcaption></figure>
  <figure><img src="assets/playground/energy.svg" alt="Energy-norm errors decrease approximately in proportion to the maximum edge length."><figcaption>Energy error uses displacement gradients. The same cube refinement has order 1.00.</figcaption></figure>
</div>

**Read the slopes.** Follow the cube through a twofold reduction in maximum edge length h. At order two, displacement error falls by about a factor of four. At order one, the energy error halves. The horizontal axis is the measured physical edge length, including on the mapped cylinder and torus. [Specimens, loads and independent assemblies](solver-verification.md)

## Check geometry before timing the solve

| Measurement | A failure it can reveal |
| --- | --- |
| Sample points on mesh faces and measure distance to the intended wall. Then sample the wall and measure distance back to the mesh. | Mesh faces can all fit well while an entire wall is missing. Both directions are sampled. |
| Compare volume by component and probe solids and voids. | Two volume errors can cancel in a global sum. |
| Check signed orientation, shape, boundary labels and periodic pairs. | A solver can converge on a wrongly connected domain. |

For a periodic cell without an exact displacement solution, compare effective stiffness under independent loading controls and refinement. Hill–Mandel energy consistency compares macroscopic work with volume-averaged microscopic work. A small refinement change alone does not establish the continuum error. [Periodic stiffness study](periodic-homogenization.md)

Compare generation, assembly and solve time only after the candidate meets the chosen geometry and physical error budgets. Keep rejected cases visible. [Benchmark index](benchmark-index.md) · [Open targets](benchmark-gaps.md)
