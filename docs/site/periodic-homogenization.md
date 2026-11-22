---
icon: lucide/repeat-2
---

# Effective stiffness of a periodic gyroid

<figure class="guide-figure chart"><img src="assets/playground/stiffness.svg" alt="Measured axial Young’s modulus decreases from 0.420 to 0.383 as the mesh is refined from 12 to 32 cells per axis."></figure>

At 32 cells per axis the gyroid occupies **35.67% of the cell** and reaches **15.34% of the solid’s axial modulus**. The 0.86% tensor change measures refinement stability. Continuum error remains unknown.

## Unit cell and loading

- **Geometry** — period-one gyroid, band `(-0.55, 0.55)`, matching opposite faces and traction-free internal walls. One node anchors translations.
- **Material** — isotropic Tet4 elasticity, E = 2.5 and ν = 0.25. Stresses are averaged over the full cell, including void.
- **Loads** — six prescribed engineering strains and six stress-controlled solves per mesh. Voigt order is xx, yy, zz, xy, xz, yz.

## Checks and replay

All 60 solves converged. The solid cube recovered the input material constants. Strain- and stress-controlled results agreed within 4 × 10⁻¹² relative.

Checks cover paired displacements, nodal force balance, Hill–Mandel energy consistency, symmetry and positive stiffness under translational periodicity.

```sh
uv run --with fedoo==1.0.0b2 demos/solver/fedoo_homogenization.py
```

[Full tensors, mesh sizes and solve records](https://github.com/kmarchais/meshers/blob/main/docs/fedoo-homogenization-results.json) · [Study implementation](https://github.com/kmarchais/meshers/blob/main/demos/solver/fedoo_homogenization.py)
