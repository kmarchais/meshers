# Does the displacement converge?

<div class="plot-pair">
  <figure><img src="assets/playground/displacement.svg" alt="Log–log displacement error by specimen; cube order 2.01, cylinder 1.97 and torus 1.92."></figure>
  <figure><img src="assets/playground/energy.svg" alt="Log–log elastic energy error; final orders range from 1.00 to 1.05."></figure>
</div>

| Specimen | Cells per axis | Final displacement order | Final energy order |
| --- | --- | ---: | ---: |
| Cube | 6, 9, 12 | 2.01 | 1.00 |
| Cylindrical gyroid | 12, 18, 24 | 1.97 | 1.05 |
| Toroidal gyroid | 12, 18, 24 | 1.92 | 1.03 |

Prescribe `u = (sin(πx) sin(πy) sin(πz), 0, 0)` on every boundary node. Use linear elastic Tet4 elements with λ = μ = 1. Integrate the Navier body force and errors with four-point quadrature.

Fedoo 0.8.4 and an independent P1 assembly agreed within 5 × 10⁻¹³ relative on all nine meshes. Both used diagonally preconditioned SciPy CG and passed the affine patch tests.

Rates use maximum physical edge length on each discrete domain. For periodic loading, see the [unit-cell study](periodic-homogenization.md).

## Replay

```sh
uv run --with scipy demos/solver/convergence.py
uv run --with fedoo==0.8.4 demos/solver/fedoo_convergence.py
```

[Independent assembly results](https://github.com/kmarchais/meshers/blob/main/docs/solver-convergence-results.json) · [Fedoo results and environment](https://github.com/kmarchais/meshers/blob/main/docs/fedoo-convergence-results.json)
