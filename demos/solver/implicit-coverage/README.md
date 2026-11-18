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
