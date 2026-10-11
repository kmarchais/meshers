# Microgen meshing benchmark

Compare Microgen's VTK tetrahedra, two VTK-to-MMG workflows, and three Meshers
optimization budgets. One command produces a standalone HTML document, raw
measurements, mesh arrays, solved displacement fields and full failure logs.

From the Meshers repository root, with [uv](https://docs.astral.sh/uv/) installed:

```sh
uv run --project benchmarks/microgen --locked microgen-benchmark
```

Open the printed `report.html` path. The default run covers ten cases at 16, 24, 36 and 54 intervals, one warmup and
three measured repetitions per configuration.
It runs serially and can take hours. Start with the shorter complete case matrix:

```sh
uv run --project benchmarks/microgen --locked microgen-benchmark --suite smoke
```

The smoke suite uses eight grid intervals, no warmups, one measurement and a
60-second limit per worker and 1,800 seconds in the standard suite. Use `--timeout`
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

## Application targets

The standard suite requires sampled boundary distances below `0.002 L` in both
directions and relative material-volume error below 1%. Volume admission includes
the empirical spread of four independently scrambled Sobol estimates and a
nested-sample change. This uncertainty estimate is not a confidence bound.
Analytic grid-edge roots and material points on the box faces provide a fixed
target-to-mesh probe. Wall vertices, edge midpoints and face centers provide the
mesh-to-field probe. Cap samples must stay inside the intended material band,
and triangle normals must agree with the outward box or field direction.
Small features between probe points can still be missed;
these checks do not certify Hausdorff distance.

FEA targets are 1%, 2% and 5% worst relative energy error. `--fea-target` defaults
to 5%. A 20% result remains diagnostic data and cannot qualify as a low-error FEA
result. Every workflow must meet the same minimum element quality for FEA target
comparisons. `qualification.json` records the fastest qualifying configuration
and the smallest qualifying configuration independently. All requested repeats
must pass. It also records observed convergence using actual tetrahedron counts.
There is no interpolation or extrapolation to claim an unmeasured target.

For visualization or printing, use `--purpose surface`. Microgen VTK calls its
public `generate_surface_mesh`; Meshers and the two volume MMG pipelines include
volume generation and boundary extraction. The released Meshers API has no direct
surface generator. Surface rows retain triangle counts, volume and boundary
errors, closed-manifold and winding checks, and nonadjacent triangle intersection
checks. The intersection test excludes pairs sharing a vertex and uses floating
point tolerances, so it is not a complete self-intersection certificate.
No elasticity solve or minimum element-shape gate is required for this use.
`--purpose geometry` measures volume generation under the same geometry gates
while skipping FEA. Both modes default to native minimum quality zero, explicitly
declared in the manifest; the FEA mode still defaults to 0.1.

```sh
microgen-benchmark --purpose surface --cases gyroid --cells 16,24,36
microgen-benchmark --purpose fea --cases gyroid --cells 24,36,54 --error-targets .01,.02,.05
```

Smoke defaults retain the original `.01 L` boundary tolerance, 10% volume
tolerance, 250,000-tetrahedron budget and 2,000,000 validation-pair cap, so issues
#3 and #4 can still be reproduced. Standard defaults use a 1,000,000-tetrahedron
budget and 20,000,000 validation-pair cap. Neither cap means that a mesh passes
when the budget runs out. Override limits explicitly and retain the manifest.

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
| Implicit BCC lattice | Eight capped struts join the center to the cube corners | Radius `0.12 L` |
| Seeded spinodoid | Fixed Microgen Gaussian random field | Iso-value 0, seed 20261010 |

Every specimen occupies the cube `[-L/2,L/2]^3`. Sheets satisfy
`abs(raw_field) <= offset(x)/2`. Offsets are field values, not physical wall
thicknesses. The HTML shows the formula, offset, scale and actual boundary
previews for each workflow. Drag a preview to rotate it. Larger surfaces use
topology-preserving simplification for display; numerical checks and saved mesh
arrays use the complete mesh.

The implicit lattice uses a benchmark-defined union of eight capsule distance
fields through Microgen's public `from_field` API. It is an implicit BCC test,
not the CAD/Gmsh `BodyCenteredCubic` workflow. The spinodoid uses the public
Microgen `Spinodoid` field, with spectral resolution fixed at 16, `k0=4π/L`,
bandwidth `3/L`, and relative mode threshold 0.05. The seed and spectral grid
stay fixed across all meshing resolutions and workflows. See the
[Microgen Spinodoid API](https://3mah.github.io/microgen-docs/microgen.shape.html#microgen.shape.spinodoid.Spinodoid). Both signed fields
use `tanh` normalization and material `field <= 0`; Meshers receives band
`[-1,0]`, whose lower boundary is inactive. Microgen samples the same field
through `from_field(...).generate_volume_mesh(resolution=cells+1)`.

Charts support workflow toggles, wheel zoom, drag pan, point selection and a
reset button. Lower and left is better. Hollow points enabled by the geometry
warning switch are diagnostic scores excluded from comparisons. Table colors
compare admissible workflows at the same intervals; those colors do not prove
a matched-element win. Minimum element quality uses the configured threshold.
Element counts only become better when accuracy is comparable.

Every unavailable value has a reason. Failed generation gives no accepted mesh
or FEA score. Invalid meshes retain completed timing and mesh arrays but block
FEA. Timeouts retain completed phases and name the interrupted phase. A warning
keeps its FEA score visibly marked as diagnostic. Progress checkpoint I/O is excluded from successful generation timings.
Failed-attempt stage duration is retained separately from successful meshing time.

Failure ownership names the pipeline requiring investigation: Meshers,
Microgen/VTK, Microgen/MMG, the direct MMG binding, or the benchmark adapter and
FEA solver. Ownership is not a confirmed root cause. Native Meshers rejections
are capability failures to report with exact fields, parameters and logs;
timeouts alone are not proof of a library defect. Keep failures in published
benchmark results. `failures.json` exports every failed, timed-out and geometry-warning
measurement with its exact reproduction command and provenance. Link any confirmed issue in the row's `issue_url` field.

The retained eight-interval Meshers native failures are tracked in
[issue #3](https://github.com/kmarchais/meshers/issues/3) for the original sheet
corpus and [issue #4](https://github.com/kmarchais/meshers/issues/4) for the two
implicit solids. Matching smoke rows link these issues automatically. Rejection
at a coarse resolution does not establish failure at finer resolutions.

`--cells` means intervals per box axis. Microgen uses `cells+1` sample points.
Equal intervals are only a sampling schedule. Actual element counts decide
budget comparisons.

* `vtk` explicitly uses `Tpms.grid_sheet.triangulate()` for sheets and the public
  implicit volume generator followed by `triangulate()` for the other fields. It does not call a
  version-dependent default backend.
* `vtk-tetra` builds a Cartesian background with VTK, triangulates it first,
  samples the shared field and clips the tetrahedra at the band bounds. This is
  an additional control for extraction order, not Microgen's current volume
  workflow. `mmg-tetra` applies the same one-pass MMG binding policy to that
  seed. Original failed VTK and MMG results remain visible under their own names.
* `vtk-tetra-clean` merges nodes in normalized coordinate bins of width `1e-12`
  and discards tetrahedra with absolute determinant at most `1e-15 L³`.
  `mmg-tetra-clean` remeshes that cleaned seed. Both record merged-node counts,
  discarded-cell counts and discarded volume. They still undergo every independent
  topology, overlap and fidelity check. These four experimental controls require
  an explicit `--workflows` selection and do not change the default six workflows.
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

The FEA mode requests minimum Meshers quality 0.1. Standard geometry tolerance
is `0.002 L`; smoke retains `0.01 L`. Configure these values explicitly to test
another policy. The manifest records the policy.

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

Implicit-boundary vertices, edge midpoints and face centers give a sampled
first-order geometry distance. Reverse probes measure the distance from analytic
wall roots and material box-face samples to the generated triangle boundary.
Both directions must meet the declared tolerance. These samples are diagnostics,
not a Hausdorff bound. Invalid meshes and geometry warnings remain in the tables
and logs. MMG shape q equals 1 for a regular tetrahedron. Low-q control meshes
retain their FEA score, but cannot qualify for the FEA application target.

Four scrambled Sobol samples of 262,144 points each estimate material volume.
The report retains the individual estimates, their spread and the nested half-
sample difference. Admission includes this empirical uncertainty, with a 1%
standard volume tolerance and 10% smoke tolerance. Volume and boundary agreement
cannot certify every missing feature. The target preview uses a separate
32-interval Microgen surface for display only.

The report provides:

* Generation time versus worst FEA error, and elements versus worst FEA error.
* Measured Meshers/control pairs within 5% of the same element count. Ratios
  compare Meshers to the control; below 1 is better. No match means no claim.
* The lowest observed element count reaching 1%, 2% and 5% error targets.
  These are observed configurations, without interpolation or extrapolation.
* Full timings and status for every repeat, including unsuccessful generations.
* Continuous green-to-red table colormaps, with a separate scale for each metric
  at the same grid resolution. Diagnostic observations keep their colors and
  explicit labels; colors alone do not establish qualification. Higher q is
  better, while lower time and FEA error are better. Missing values stay gray.
  Matched-count ratios use a fixed scale centered on 1, and error-target tables
  compare element counts and timings only among qualifying observations.

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

The standard suite uses PyAMG smoothed aggregation with six rigid-body candidate
modes and conjugate gradients at relative tolerance `1e-11`. Nonconvergence fails
the worker. The smoke suite uses sparse direct LU. Select either with
`--linear-solver direct` or `--linear-solver amg`; both solve the same 24 right-hand
sides and run the same patch, equilibrium and energy checks. The solver is common
to every meshing workflow and its time is excluded from meshing time.

`manifest.json` records machine details, configuration, package versions,
Meshers/MMG native-binary SHA-256 hashes, benchmark source hashes and the git
commit. `results.json` records a hash of that manifest and each completed mesh
and field archive. A research build with the same version number can still
differ; compare the binary hashes before combining results. No run is resumed
implicitly and an existing nonempty output directory is never overwritten.
The runner stops if its benchmark source files change during the run. Each output
also retains the exact package source, `pyproject.toml` and `uv.lock` under `source/`,
with SHA-256 hashes. Run that snapshot with `uv run --project <output>/source --locked
microgen-benchmark ...` to reproduce an experiment made before the next commit.

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
