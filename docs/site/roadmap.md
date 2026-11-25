---
icon: lucide/route
---

# Meshers roadmap

Last reviewed **22 November 2026**. [Available and experimental features](index.md#project-scope) · [Benchmark index](benchmark-index.md)

## Current work on implicit solids

| Open target | What would count as progress |
| --- | --- |
| Preserve thin walls, gaps and cavities | Recover the feature or reject the unresolved mesh. Check material **and void** at known points, then measure volume and wall distance in both directions. [Lost-cavity replay](benchmark-gaps.md#reproduce-the-lost-cavity) |
| Resolve harder spinodoids | Pass the eight-cell geometry target on the requested route. Then repeat on held-out phases and report nodes, tetrahedra and time. [Measured limits](benchmark-gaps.md#where-to-work-next) |
| Improve Tet4 elasticity convergence | Reach the 2% worst energy-error target on the gyroid and BCC lattice after geometry admission. [Refinement results](benchmark-index.md#geometry-and-solver-studies) |

## Next FEA capabilities

These three priorities are **planned**. Each starts with a small problem whose result can be checked independently.

| Priority | First milestone | Benchmark |
| --- | --- | --- |
| **1. Multi-material interfaces** | Mesh two bonded layers with shared interface faces and preserve material IDs through export. | Check each layer’s volume. Load equal-thickness layers in series with zero Poisson ratio and Young’s moduli 1 and 4. The effective axial modulus should be **1.6**, their harmonic mean. |
| **2. Solver-driven refinement and solution transfer** | Solve, estimate local error, mark cells, adapt and transfer the solution. Preserve loads and boundary labels across the cycle. | Transfer a linear field exactly up to roundoff. On a problem with a localized stress gradient, compare adaptive and uniform refinement at the same energy error. Count free DOFs and time the full cycle. |
| **3. Directional sizing and hybrid meshes** | Resolve a thin layer across its thickness with wedges or swept hexes. Connect quad faces to a tetrahedral region through pyramids. | Check positive orientation and matching interface faces. Run an affine patch test across the transition, then compare node count and time at matched FEA error against Tet4. |

## Longer-term geometry and element support

- **CAD/BRep input.** Start with a plate with a cylindrical bore. Preserve its named face and applied load; check bore radius and solid volume after export.
- **Production triangle and quad surfaces.** Compare direct extraction on the same wall and crease. Measure face count and time at matched boundary distance after topology checks.
- **Hex volume meshing.** Extend beyond swept layers to more general solids. Compare curved-boundary fit and solver error with Tet4 and Tet10 on the same geometry.

[Proposed benchmark designs](benchmark-index.md#benchmarks-to-add)
