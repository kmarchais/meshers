---
icon: lucide/chart-no-axes-combined
---

# Decide whether the mesh solves the physical problem

1. Prescribe the analytic displacement on the boundary.
2. Apply the matching body force to linear elastic Tet4 elements.
3. Refine the same specimen and compare both errors below.

<div class="plot-pair">
  <figure><img src="assets/playground/displacement.svg" alt="Displacement error on the cube and mapped gyroids falls approximately as h squared."><figcaption>Displacement L2 error. The finest cube refinement has observed order 2.01.</figcaption></figure>
  <figure><img src="assets/playground/energy.svg" alt="Energy-norm errors decrease approximately in proportion to the maximum edge length."><figcaption>Energy error uses displacement gradients. The same cube refinement has order 1.00.</figcaption></figure>
</div>

<details><summary>If you halve the physical edge length, which error falls faster?</summary>
<ul><li>Displacement error scales approximately as h². Expect about one quarter of the error.</li><li>Energy-norm error scales approximately as h. Expect about one half.</li></ul>
</details>

[Compute the observed order](https://github.com/kmarchais/meshers/blob/fc54032b/demos/solver/convergence.py#L318-L325) · [Specimens and independent assemblies](solver-verification.md)

## Check geometry before timing the solve

| Measurement | A failure it can reveal |
| --- | --- |
| Sample points on mesh faces and measure distance to the intended wall. Then sample the wall and measure distance back to the mesh. | Mesh faces can all fit well while an entire wall is missing. Both directions are sampled. |
| Compare volume by component and probe solids and voids. | Two volume errors can cancel in a global sum. |
| Check signed orientation, shape, boundary labels and periodic pairs. | A solver can converge on a wrongly connected domain. |

- **Periodic stiffness:** compare independent loading controls and successive refinements.
- **Hill–Mandel check:** compare macroscopic work with volume-averaged microscopic work.
- **Accuracy claim:** use a reference. Refinement stability alone leaves continuum error unknown. [Periodic stiffness study](periodic-homogenization.md)

**Compare cost at the target error:** retain rejected cases and record generation, assembly and solve times separately. [Benchmark index](benchmark-index.md) · [Open targets](benchmark-gaps.md)
