# Geometry coverage of the seven historical implicit solids

All seven solids from issues [#3](https://github.com/kmarchais/meshers/issues/3) and [#4](https://github.com/kmarchais/meshers/issues/4) have a Tet4 mesh that passes the original smoke limits on the tested Linux configuration. Six have a passing mesh from eight background intervals per axis. The spinodoid required a finer starting grid in this study.

| Specimen | Cells per axis | Optimization passes | Adaptive | Tetrahedra |
| --- | ---: | ---: | :---: | ---: |
| Uniform Gyroid | 8 | 4 | No | 3,493 |
| Schwarz P sheet | 8 | 4 | No | 2,568 |
| Microgen Neovius sheet | 8 | 4 | Yes | 12,522 |
| Graded Gyroid | 8 | 4 | Yes | 6,056 |
| Thin Gyroid | 8 | 4 | Yes | 9,529 |
| BCC capsule lattice | 8 | 4 | Yes | 15,303 |
| Seeded spinodoid | 20 | 4 | No | 59,363 |

Uniform Gyroid also passed with zero and twelve passes. For the graded Gyroid, thin Gyroid and BCC lattice, the four-pass uniform mesher returned elements at eight cells, but the independent target-to-mesh distances exceeded 0.01. Setting `adapt=True` in `meshers.generate` admitted all three on the same starting grid. Refinement can introduce a different failure. The spinodoid passed at 20 cells with minimum quality 0.1133, then failed at 24 cells with minimum quality 0.0510.

For the eight-cell spinodoid, tightening `adapt_tolerance` to 0.005 and 0.002 still produced native geometry rejections. The native sampled error estimate rose to 0.200 and 2.30, respectively. Both attempts stopped before independent distance validation. The eight-cell spinodoid remains open.

The table gives one observed setting per solid. Each attempt used one meshing thread, a unit box, the original material band and a budget of 250,000 tetrahedra. The pinned validator checks orientation, connectivity, overlap, sampled outward normals, wall distance in both directions, cap overshoot, minimum MMG quality of 0.1 and relative volume error plus reference uncertainty of at most 10%. The distance limit is 0.01 unit lengths. The reverse-distance samples can miss features between grid lines. The archive includes rejected settings. Timings are single observations on a Ryzen 7 7730U running Linux. Peak memory was unmeasured. [Numerical records](2026-11-18.json)

### Elasticity

Two admitted meshes also ran the archived manufactured elasticity problem. Its worst relative energy errors were **14.08%** on the uniform Gyroid and **10.94%** on the adaptive BCC lattice. Both exceeded the 2% target. The test uses linear tetrahedra, isotropic elasticity with Young’s modulus 1 and Poisson’s ratio 0.3, exact body forces and boundary tractions, and six rigid-motion constraints per connected component. It scores the worst linear combination of 18 quadratic displacement fields and checks six affine patch loads. Geometry rejection skips the solve. [Original elasticity formulation](https://github.com/kmarchais/meshers/blob/602dc900cd028d205b36e29dcdb5b4a9cf4f581a/benchmarks/microgen/microgen_benchmark/elasticity.py)

### Replay an attempt

Run from the current repository root. The historical checkout supplies the fields, reference sampler, validator and solver. The command runs the current workspace version of Meshers and verifies the historical Python sources by digest.

```sh
git fetch origin refs/pull/2/head
git worktree add --detach .local-work/coverage-reference 602dc900cd028d205b36e29dcdb5b4a9cf4f581a
uv sync --locked --all-packages
uv run --with scipy python demos/solver/implicit_coverage.py --reference .local-work/coverage-reference/benchmarks/microgen --case implicit-lattice --cells 8 --passes 4 --adapt --require-admitted --output lattice.json
```

For the spinodoid, add `--with microgen==2.0.0b1` to the `uv run` options and use `--case spinodoid --cells 20 --passes 4`. Keep adaptive meshing off for that row. The seed, spectral resolution, mode threshold and physical cell stay fixed across mesh resolutions. On Linux, Microgen’s Gmsh dependency requires the system GLU library.

Use a fresh output filename for each attempt. Add `--fea` to run the elasticity check after geometry admission. With `--require-admitted --fea`, the process also requires the 2% energy target. Without `--require-admitted`, a completed rejection is a diagnostic result. An incomplete validation or solver failure returns a failing exit code.

### Eight-cell sheet and lattice controls

Each of the five sheets had a passing eight-cell recipe at zero, four and twelve optimization passes on Linux. Uniform Gyroid passed on the original uniform grid at all three pass counts. The other sheets used adaptive meshing for some or all settings. With zero passes, the thin Gyroid required particular care. An internal sizing target of 0.008 produced 11,606 tetrahedra with minimum quality 0.100019 against the 0.1 floor. The six other internal targets tested between 0.004 and 0.009 all failed with default snapping. Reducing snapping also failed. These results apply to this specimen and build. [All follow-up attempts](followup-2026-11-18.json)

The BCC lattice also passed at all three pass counts. At zero passes the default adaptive sizing missed the reverse-wall limit (0.01135). Three tighter sizing targets passed. The CI control uses 0.006, which gave a sampled reverse distance of 0.00398, minimum quality 0.3281 and 17,230 tetrahedra. The four- and twelve-pass controls use default adaptive sizing.

### Where the spinodoid misses the cap

Every adaptive eight-cell rejection below reached its largest first-order estimate at a wall edge midpoint on a box cap. We recovered the final face samples from the frozen wheel's field calls and matched each reported maximum. On the current build, a geometry rejection includes the coordinates directly.

| Grid and internal target | Worst sample | Native estimate | Independent distance interval |
| --- | --- | ---: | ---: |
| 8, adaptive 0.01 | x-min cap | 0.12703 | 0.01949–0.02381 |
| 8, adaptive 0.005 | y-min cap | 0.19987 | 0.00912–0.01065 |
| 8, adaptive 0.002 | x-min cap | 2.30049 | 0.02097–0.02661 |
| 20, uniform control | Interior wall edge | 0.00598 | 0.00548–0.00585 |

The adaptive intervals measure distance to the zero curve **within the cap**. The control measures distance to the wall in three dimensions. At internal targets 0.01 and 0.002, even the lower bound exceeds the 0.01 limit. At 0.005, the interval straddles the limit. Reducing the lattice snap fraction from 0.2 to 0.1 or zero also produced a native rejection. The next geometry experiment is to split and refit the offending cap edge, then repeat the independent audit on the whole mesh.

The [archived spectrum](spinodoid-spectrum.csv) contains the 52 cosine modes of the frozen unit-cell specimen, with phases in radians and the original binary64 precision. Replay checks the signed field against Microgen’s public evaluator on a 17×17×17 periodic grid before using those coefficients.

The lower bound uses the local gradient and Hessian of the pinned Fourier field and a global third-derivative bound from its mode amplitudes. The upper bound comes from a bracketed zero along an in-box ray. Both describe one sample. To reproduce the first row using the historical checkout above:

```sh
uv run --with microgen==2.0.0b1 python demos/solver/spinodoid_witness.py --reference .local-work/coverage-reference/benchmarks/microgen --point -0.5 -0.1380501438064803 -0.14563324638067843 --lock-caps
```

### Refinement at the 2% manufactured-energy target

| Specimen and route | Cells per axis | Tetrahedra | Worst energy error |
| --- | ---: | ---: | ---: |
| Uniform Gyroid | 8 | 3,493 | 14.08% |
| Uniform Gyroid | 16 | 16,218 | 8.30% |
| Uniform Gyroid | 24 | 43,692 | 5.78% |
| Adaptive BCC | 8 | 15,303 | 10.94% |
| Adaptive BCC | 12 | 27,886 | 10.49% |

The BCC adaptive 20-cell and uniform 24-cell attempts failed geometry admission before the solve. Their sampled target-to-mesh distances were 0.01230 and 0.01019, respectively. The 20-cell mesh also failed the normals check. None of the admitted Tet4 meshes met 2%. These single-configuration refinement runs establish accuracy, with corresponding assembly and solve times in the records. Full refinement sweeps on the other configurations are pending. [Numerical records](followup-2026-11-18.json)
