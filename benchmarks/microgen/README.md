# Microgen meshing benchmark

Compare Microgen's VTK tetrahedra, two VTK-to-MMG workflows, and three Meshers
optimization budgets. One command produces a standalone HTML document, raw
measurements, mesh arrays, solved displacement fields and full failure logs.

From the Meshers repository root, with [uv](https://docs.astral.sh/uv/) installed:

```sh
uv run --project benchmarks/microgen --locked microgen-benchmark
```

Open the printed `report.html` path. The default run covers eight cases, three
resolution levels, one warmup and three measured repetitions per configuration.
It runs serially and can take hours. Start with the shorter complete case matrix:

```sh
uv run --project benchmarks/microgen --locked microgen-benchmark --suite smoke
```

The smoke suite uses eight grid intervals, no warmups, one measurement and a
60-second limit per worker (300 seconds in the standard suite). Use `--timeout`
to change that limit. It
checks that the workflows run and the report contains every case. Its timings
are not sufficient for performance claims, and coarse meshes may fail.

The project pins Python 3.12 and locks its Python dependencies. It installs
public Meshers and MMG wheels; it does not require a Rust build, an external
MMG executable, a research worktree or a modified Microgen checkout. Use a
machine without other heavy work when collecting timings. On Kevin's machine,
submit the command through the serial `pueue gpu` lane:

```sh
qrun -g gpu -l "Microgen speed and FEA benchmark" --progress report uv run --project benchmarks/microgen --locked microgen-benchmark
```

The runner streams child output and reports progress every two seconds. It uses
the Workbench `pueue_progress` helper when installed, and OSC terminal progress
otherwise. A timeout kills the owned worker and its descendants. Other jobs and
the queue's group state are untouched.

Linux also needs Gmsh's system runtime libraries, typically
`sudo apt-get install libglu1-mesa libgl1` on Debian or Ubuntu. The CI workflow
installs these before running the benchmark.

## Cases and workflows

| Case | Purpose | Raw field offset |
| --- | --- | --- |
| Solid cube | Affine patch and refinement control | Full cube |
| Gyroid | Uniform sheet | 0.5 |
| Schwarz P | Sheet with narrow junctions | 0.6 |
| Diamond / Schwarz D | Different junction structure | 0.8 |
| Microgen Neovius variant | Released function with unequal cosine coefficients | 1.0 |
| Graded Gyroid | Actual spatially varying thickness | `0.5 * (1 + 0.4*x/L)` |
| Thin Gyroid | Thin-wall resolution stress test | 0.25 |
| Gyroid at scale 0.001 | Units and scaling control | 0.5 |

Every specimen occupies the cube `[-L/2,L/2]^3`. Sheets satisfy
`abs(raw_field) <= offset(x)/2`. Offsets are field values, not physical wall
thicknesses. The HTML shows the formula, offset, scale and actual boundary
previews for each workflow. Drag a preview to rotate it. Larger surfaces use
topology-preserving simplification for display; numerical checks and saved mesh
arrays use the complete mesh.

`--cells` means intervals per box axis. Microgen uses `cells+1` sample points.
Equal intervals are only a sampling schedule. Actual element counts decide
budget comparisons.

* `vtk` explicitly uses `Tpms.grid_sheet.triangulate()`. It does not call a
  version-dependent default backend.
* `mmg` converts those VTK tetrahedra to MMG3D through mmgpy, runs one volume
  adaptation pass with `hsiz=L/cells` and `hgrad=1.3`, then converts back. Box-cap
  triangles stay fixed; MMG can adapt the curved boundary. This portable binding workflow is distinct from
  Microgen's legacy two-pass command-line remesher. Its surface approximation
  starts from VTK's polyhedral boundary. `--mmg-boundary fixed` instead protects
  every surface triangle for a comparison on the same VTK polyhedral domain.
* `mmg-microgen` calls Microgen's public `remesh_keeping_boundaries_for_fem`
  function using mmgpy's installed `mmg3d_O3` command alias. It runs
  Microgen's two passes, file exchange and built-in retries with `hsiz=L/cells`,
  `hausd=geometry_tolerance*L` and `hgrad=1.3`. These costs enter the generation
  time. `periodic=False` matches this benchmark's finite specimens and varying
  thickness; periodic homogenization needs a separate test. `--mmg-boundary`
  applies only to the direct binding workflow.
  On Windows, a scoped temporary-file wrapper sets Python 3.12's
  `delete_on_close=False`, allowing Microgen to reopen a file before its context
  closes. Context cleanup remains enabled. This changes file handling only;
  MMG passes and parameters stay unchanged.
* `meshers-fast`, `meshers-balanced` and `meshers-quality` call the public
  `meshers.generate` API with 0, 4 and 12 optimization passes. All use the same
  geometry tolerance and minimum quality. A higher budget does not guarantee a
  faster run, a successful mesh or lower FEA error. Failures remain failures.

See Microgen's [meshing examples](https://microgen.readthedocs.io/en/latest/examples/mesh.html)
for the legacy workflow, and MMG's
[required-boundary implementation](https://www.mmgtools.org/documentation/mmg3d/documentation_mmg3d_mesh_analysis.html)
for the effect of `nosurf`.

The default Meshers minimum quality is 0.1 and geometry tolerance is 0.01 times
the cell size. Configure them explicitly to test another policy. The report
records the policy; no automatic relaxation or backend fallback occurs.

## What FEA quality means here

All workflows use the same small-strain isotropic linear-elastic solver with
linear tetrahedra, E=1 and Poisson ratio 0.3. Coordinates and displacements are
nondimensionalized by L. E and L therefore do not change the relative score.

The test space contains `x², y², z², xy, yz, zx` in each of the three displacement
components, for 18 quadratic fields. Each field gives exact body forces and
natural tractions on every produced boundary face. Six independently selected
scalar constraints per material component remove only rigid motion. The runner
assembles one stiffness matrix, factorizes it once and solves all right-hand
sides. Three-point triangle quadrature integrates the linear tractions exactly;
four-point tetrahedral quadrature integrates the quadratic energy exactly.

Let G be the exact strain-energy Gram matrix and H the displacement-error
energy Gram matrix. The primary score is `sqrt(lambda_max(H,G))`. This gives the
worst relative energy error over **every linear combination** of the fixed 18
fields. The report also retains each field's error and the worst combination.
Six affine fields must pass the patch test below `1e-8`; free-equilibrium
residuals must be below `1e-8`.

This is an exact-reference approximation test on each produced polyhedral
domain. It measures discretization error without trusting a competitor's mesh
as ground truth. It does not certify a shared curved-domain engineering problem,
periodic homogenization, stress concentrations, nonlinear response or every
possible load. Add a converged same-domain engineering fixture before making
those claims. Keep that additional test family fixed during optimization.

## Validity and fair comparisons

Before FEA, the runner checks positive volumes, duplicate cells, face conformity,
closed manifold boundary edges and vertex links, disconnected material contacts,
and tetrahedral interior overlap. The overlap test uses a spatial broad phase
and separating axes at a dimensionless tolerance of `1e-10`. Exceeding its pair
cap is an incomplete check, never a pass. Very small determinants below `1e-15`
are rejected. This floating-point implementation is not an exact-arithmetic
mesh certificate.

Implicit-boundary vertices and face centroids give a sampled first-order
geometry distance, `abs(abs(field)-1)/norm(gradient)`. Cap faces are excluded.
This is a diagnostic, not a Hausdorff bound. Only runs within the declared
sampled tolerance enter FEA comparison charts. Invalid meshes and geometry
warnings remain in the tables and logs. The MMG shape score q is a separate
diagnostic, with 1 for a regular tetrahedron. A control's q below 0.1 alone does
not disqualify its FEA measurement.

An independent scrambled Sobol sample of 262,144 points estimates the intended
material volume. The report records its difference from the first 131,072
samples. This difference is a convergence diagnostic, not a confidence bound.
Meshes more than 10% from that volume estimate do not enter comparisons. Use
`--volume-tolerance` to declare another threshold. Volume agreement and sampled
boundary agreement still cannot certify every missing feature. The target
preview uses a separate 32-interval Microgen surface for display only.

The report provides:

* Generation time versus worst FEA error, and elements versus worst FEA error.
* Measured Meshers/control pairs within 5% of the same element count. Ratios
  compare Meshers to the control; below 1 is better. No match means no claim.
* The lowest observed element count reaching 5%, 10% and 20% error targets.
  These are observed configurations, without interpolation or extrapolation.
* Full timings and status for every repeat, including unsuccessful generations.

A configuration enters comparisons only if every measured repeat passes.
The largest observed element count and worst FEA error are used, so count
variation is handled conservatively. Every repetition of a matched-count pair
must stay within the declared 5% tolerance. Warmups are retained but excluded.
Generation medians include available
generation timings even when a later independent check fails; those runs never
enter the comparison charts.

## Timing and reproducibility

Each repetition runs in a fresh worker process. Computational thread counts are
fixed, and backend order rotates between repetitions. Generation time includes
Microgen case construction, backend extraction or remeshing, and copies into
owned arrays. Imports, process startup, independent validation, FEA, compression
and file writing are excluded. Assembly, solve and error-evaluation times are
saved separately. These timings measure the declared Python API workflow, not
only its native kernel or a whole Microgen application.

`manifest.json` records machine details, configuration, package versions,
Meshers/MMG native-binary SHA-256 hashes, benchmark source hashes and the git
commit. `results.json` records a hash of that manifest and each completed mesh
and field archive. A research build with the same version number can still
differ; compare the binary hashes before combining results. No run is resumed
implicitly and an existing nonempty output directory is never overwritten.
The runner stops if its benchmark source files change during the run.

Examples:

```sh
# A narrower repeated comparison with more budget points.
uv run --project benchmarks/microgen --locked microgen-benchmark --cases gyroid,schwarz-p --cells 12,16,20,24,28 --repeats 5 --warmups 1 --output benchmarks/microgen/output/repeated

# Inspect just the fast and quality profiles at a declared relaxed policy.
uv run --project benchmarks/microgen --locked microgen-benchmark --suite smoke --workflows vtk,mmg,meshers-fast,meshers-quality --minimum-quality 0.05

# Regenerate a report after a template change, without running meshes or FEA.
uv run --project benchmarks/microgen --locked microgen-benchmark --render benchmarks/microgen/output/repeated/results.json

# Numerical and failure-path tests.
uv run --project benchmarks/microgen --locked --group dev python -m pytest benchmarks/microgen/tests
```

To benchmark a local Meshers implementation, replace the wheel in the benchmark
environment with an editable native build, then use `uv run --no-sync --project
benchmarks/microgen microgen-benchmark ...`. Build in the heavy queue when
required. Keep the locked-wheel result as the baseline and retain both manifests.
Do not label a local build as the release wheel.

The `Microgen benchmark` GitHub Actions workflow runs numerical tests and the
smoke matrix on Windows, Linux and macOS for benchmark pull requests. It also
accepts a manual standard run. Download its `microgen-<OS>` artifact and open
`report.html` after extracting the complete directory. Shared CI hardware is
useful for portability checks, but a dedicated machine is better for speed claims.
