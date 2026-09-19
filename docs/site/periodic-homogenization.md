# Periodic homogenization with fedoo 1.0 beta

This study exercises **fedoo 1.0.0b2** with Meshers 0.1.0 tetrahedra. It tests
translational periodicity on all three pairs of faces of a unit cube, first for a fully
solid cell and then for a gyroid sheet. The gyroid uses the built-in period-one
field with `band=(-0.55, 0.55)`. Internal solid/void surfaces are traction-free.

The material is isotropic linear elastic, with `E=2.5` and `nu=0.25`
(`lambda=mu=1`). Six unit engineering-strain cases determine the full effective
6 by 6 stiffness. The Voigt order is `xx, yy, zz, xy, xz, yz`; shear strains
are engineering angles, while shear stresses are physical tensor components.
All averages use the **full cell volume, including void**, not the solid volume.

## What is checked

Fedoo assembles `tet4` elasticity, enforces `PeriodicBC("small_strain")` and
solves each prescribed-strain case. One node is pinned to remove translations;
the remaining fluctuations are free. Element strain and stress are integrated
independently from the returned nodal displacements. One integration point is
exact for the constant strain/stress in each homogeneous linear tetrahedron.

The script checks:

- Exact opposite-face node correspondence and matching surface triangulations.
- `u(x+) - u(x-) = E_macro (x+ - x-)` for every Meshers periodic pair.
- Nodal-force balance within every periodic equivalence class, including edges
  and corners, and equilibrium of every interior node.
- Hill–Mandel consistency, including cross terms between the six load cases.
- Agreement between volume-averaged stress and fedoo's macroscopic reactions.
- Stiffness symmetry, positive eigenvalues and the uniform-strain upper bound.
- Recovery of the analytic solid-material stiffness.
- Agreement with fedoo's automatic **stress-controlled**
  `get_homogenized_stiffness` calculation on each same mesh.

See the pinned [fedoo periodic-constraint implementation](https://github.com/3MAH/fedoo/blob/v1.0.0b2/fedoo/constraint/periodic_bc.py)
and [homogenization implementation](https://github.com/3MAH/fedoo/blob/v1.0.0b2/fedoo/homogen/tangent_stiffness.py).
Both loading paths use fedoo's constraints and assembly. The independent checks
use Meshers pair metadata and physical tetrahedral gradients.

## Measured results

The solid-cell test recovers `E=2.5` and `nu=0.25` to numerical precision.
The gyroid refinement results are:

| Cells per axis | Tetrahedra | Solid fraction | Effective E_x | Effective nu_xy | Full-tensor change from previous level |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 12 | 10,832 | 0.359680 | 0.419886 | 0.268179 | — |
| 18 | 34,186 | 0.357501 | 0.396716 | 0.275709 | 3.91% |
| 24 | 72,757 | 0.357098 | 0.388502 | 0.278807 | 1.38% |
| 32 | 171,679 | 0.356732 | 0.383458 | 0.280675 | 0.86% |

The change is the Frobenius norm of the tensor difference divided by the newer
tensor's norm. At 32 cells, `E_x/E_solid` is approximately 0.1534. The three
axial Young's moduli agree closely, but the full tensor is retained: equal
axial moduli alone do not establish isotropy. This gyroid is not auxetic under
these axial loading conditions.

CG uses a diagonal preconditioner with `rtol=1e-12`, `atol=1e-14`, and explicit
status/residual checks. The tighter tolerance keeps the aggregate nodal-force
check accurate on the finest mesh; its acceptance threshold is unchanged.
All 60 linear systems converge: six strain-controlled and six automatic
stress-controlled cases for the solid and each of the four gyroid refinements.
Hill–Mandel errors and differences between the two homogenization paths are
below 4e-12 relative. The effective tensors are symmetric to numerical accuracy
and have positive eigenvalues.

## Scope

This validates a connected, linearly elastic gyroid cell with **translational**
periodicity. It does not test rotational constraints on torus sections,
cylindrical homogenization, nonlinear response, buckling, or equivalence with
the Microgen VTK/MMG pipeline. Refinement measurements show discretization
stability; they do not certify the exact continuum tensor. The gyroid geometry
and its measured solid fraction change slightly with resolution.

## Reproduce with uv

From a source checkout configured with Meshers and its build toolchain:

```sh
uv run --with fedoo==1.0.0b2 python demos/solver/fedoo_homogenization.py
```

To test a downloaded candidate wheel in isolation, replace the wheel path:

```sh
uv venv .venv-fedoo --python 3.12
uv pip install --python .venv-fedoo/bin/python fedoo==1.0.0b2 /path/to/meshers.whl
uv run --no-project --python .venv-fedoo/bin/python python demos/solver/fedoo_homogenization.py
```

The script writes `docs/fedoo-homogenization-results.json`, including the full
stiffness tensors, every linear solve, verification errors, dependency versions,
and source/extension hashes. CI runs the study and uploads this report alongside
the manufactured-solution reports. No fedoo dependency is added to Meshers itself.

The recorded local run uses the tested macOS ARM candidate wheel from commit
`5f499ad5fa53cf189314f2954212a07a4856ecbb`, on an Apple M4 (4 performance and
6 efficiency cores), 24 GiB RAM, macOS 26.5.1 arm64. Timings are observations,
not controlled performance comparisons.
