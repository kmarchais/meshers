# Real Microgen workflow on G3

This experiment starts from Microgen 2.0.0b1's own volumetric VTK output.
It does not remesh a Meshers seed. Native meshing code is unchanged.

## Fixed specimens and defaults

All specimens occupy `[-0.5, 0.5]^3`. The isotropic material has E = 2.5 and
nu = 0.25, so lambda = mu = 1. Six load cases include engineering shear.

- Gyroid sheet has raw band `[-0.55, 0.55]`. Microgen's offset is 1.1.
  `Tpms` keeps its default resolution of 20 grid points per axis.
- Schwarz P solid is the negative sublevel set of
  `cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*z)`. This is not the earlier
  Schwarz P sheet fixture. `Tpms` uses offset zero, `lower skeletal`, and
  its default resolution of 20 points.
- Spinodoid keeps the documented spectral realization fixed at resolution 16,
  seed 20261010, k0 = 4*pi, bandwidth 3, offset zero, and retained-mode
  threshold 0.05. Its native `[0, 1]^3` coordinates are translated by -0.5.
  These spectral parameters define the specimen, not Microgen defaults.
- The smooth lattice has eight capsule struts, radius 0.12, from the origin
  to the cube corners. Microgen's sequential cubic `smooth_union` uses
  k = 0.08 and the existing `scorecard.CORNERS` order. Its inherited F-rep
  volume generator keeps the default resolution of 50 points. This route
  is not the CAD strut-lattice STEP/gmsh route.

The VTK control calls `generate_volume_mesh` and `triangulate`. The MMG
control passes this fresh output to `remesh_keeping_boundaries_for_fem`
without sizing overrides. Its defaults are `periodic=True`,
`mesh_version=2`, `dimension=3`, `tol=1e-8`, and all sizing parameters
unset. Microgen calls MMG first with `nofem=True`, then with `ls=True`,
`nr=True`, and the unset size parameters. Its own retry limit is five.
The runner records actual calls and preserves installed source hashes.
The native CLI comes from mmgpy 0.17.1. No binding replacement is used.

Meshers uses a periodic 24-cell linear mesh and curved quadratic meshes
at 12, 16 and 20 cells. One optimization pass, minimum corner quality
0.1 and one thread are fixed. Its internal linear-facet geometry request
is 1.0, matching the historical quadratic convergence experiment.
That is not the admission tolerance. Independent admission stays at
0.002 boundary distance, 1% relative volume error and quality 0.1.
Unqualified solves are retained as diagnostics, never as accepted results.

## Bounded execution

The runner builds the extension with Rust 1.95.0, two compile jobs and the
locked Cargo dependencies. Each later worker uses one numerical thread,
an aggregate process-tree RSS ceiling of 8 GiB checked every 0.25 seconds,
a per-process address-space ceiling of 8 GiB,
a 250,000-tetrahedron cap and a 100,000-node solve cap. A VTK seed over the
element cap stops before MMG. Meshing and solve workers have a 1,200-second
timeout. Sampling workers have 600 seconds. No Windows binaries are used.

Use the shared `gpu` queue from the repository root:

```sh
qrun -l 'Real Microgen defaults on G3' --progress report \
  --interrupt none --interrupt-reason 'Immutable run, inspect before retrying' \
  --check-start --wait .venv/bin/python experiments/adaptive/microgen_workflow.py \
  --run experiments/adaptive/runs/g3-microgen-defaults-next
```

Choose a new run directory. Existing runs are not overwritten. Queue only
one experiment, read its results, then decide what to run next.
The Queue card uses `Job` with a bar per worker. Within a worker its bar
shows elapsed timeout budget, not solver convergence. The message states
that cap. Completed workers advance the experiment bar.

The runtime needs Python 3.12, the run's exact `requirements.txt`, and
user-local GLU/OpenGL runtime libraries. On this G3, Ubuntu packages
`libglu1-mesa` 9.0.2-1.1build2 and `libopengl0` 1.7.0-3 are extracted under
`.local-work/g3-runtime/root`. The runner sets `LD_LIBRARY_PATH` for them.
These libraries satisfy gmsh's import dependency without a system install.

## Measurements and limits

Generation time includes specimen construction, VTK extraction and MMG
where applicable. Imports are excluded. Audit time and six-load solve
time are recorded separately. Total workflow time includes those audits.
Mesh serialization is also included in total time. These are single runs
on G3, not repeatable performance rankings or comparisons with OMEN.

Boundary distances use the analytic field, independent projected Sobol
wall samples, planar cap samples, and point-to-triangle VTK queries.
The forward direction samples 13 parametric points per original wall
face and locks coordinates on box planes. Quadratic faces are evaluated
with their six-node interpolation, not replaced by their corner plane.
Reverse distance uses four then eight subdivisions per quadratic edge.
The largest change in reverse distance is recorded. Finite sampled maxima
and tessellation checks are not bounds on the continuous Hausdorff distance.

Volume comes from signed linear determinants or 15-point Tet10 Jacobian
quadrature. Four independent Sobol scrambles of 262,144 points estimate
material volume, except the Schwarz solid's exact volume of 0.5. Their
spread is recorded, not presented as a confidence interval. Outward and
inward forward deviations and signed volume bias distinguish sampled
material-domain changes.

Quality and dihedral angles refer to corner tetrahedra. Quadratic
validity also requires positive sampled Jacobians at 15 Gauss points,
all ten reference nodes, and the centroid. This is not a positivity proof.
Periodic checks include all nodes, mid-edge nodes and corner cap triangles.
No duplicate-node cleanup, zero-volume deletion, boundary repair or
nonmatching-periodic substitution is allowed in a default control.
Tet4 orientation is converted losslessly by swapping two local node
indices where needed. This changes neither the material domain nor any
boundary triangle. Negative input orientation is recorded separately
from tiny absolute volume. A signed sum before this conversion is only
an ordering diagnostic, not a physical volume.

Full tensor disagreement is the largest absolute generalized eigenvalue
of `C_candidate - C_reference` against `C_reference`. This includes mixed
loads and cross terms. Two finer curved meshes provide an assumed-order
p=4 extrapolation. The 12/16 versus 16/20 extrapolation difference and
historical Gyroid and Schwarz-solid tensor checks test its sensitivity.
References with failed geometry remain unqualified. Reference disagreement
is not certified continuum error. No local stress accuracy is measured.

The correction run adds two saved-mesh controls for each surviving native
12-cell curved mesh. The first retains its corner geometry and uses Tet4.
The second retains exactly that piecewise planar domain and uses Tet10
with straight midpoint nodes. Their tensor difference measures a change
in approximation order on the same domain. The difference between straight
and curved Tet10 measures a geometry-mapping change at fixed order and
connectivity. The two tensor differences add to the full change. Their
worst-energy scalar scores need not add. These deltas are not separate
certified continuum errors. No live generation time is assigned to replay
controls.

The runner tests triangle interpolation, tessellation, quadrature weights,
mixed-load sensitivity, and a real native Tet10 affine solid-cube solve
before the specimens. Run records contain source and native hashes,
dependency versions, worker logs, meshes, tensors, solver residuals and
the summary. Historical field archives were not transferred or revalidated.

## Runner failures

Job 2 stopped before specimen measurements because the installed `Job`
helper does not accept an explicit `eta_s` value through `Part.update`.
The parent exited while its separate build process group survived.
That build was identified and terminated before submitting job 3.
The runner now leaves ETA calculation to `Job` and terminates its child
if progress reporting raises. The first directory remains unchanged.
Job 3 uses the corrected runner and a fresh directory.

Job 3 built the Linux extension but its affine self-check requested two
cells, below Meshers' minimum of four. The test input was corrected.
Job 4 passed the affine check and completed the first specimen pilot.
Its progress audit reported no findings. Peak sampled worker RSS was
2.37 GiB. All four native spinodoid configurations failed the unchanged
quality floor. Microgen's spinodoid output failed the independent cap
triangulation match.

Two measurement checks in job 4 needed correction. VTK's mixed local
node orientations were misclassified as invalid cells. A lossless
conversion now precedes volume integration and FEA. The smooth lattice
preflight also required its scalar field to be periodic everywhere.
Sequential cubic unions differ in void, though 131,072 independent face
probes per axis found no material mismatch and no field difference within
0.03 of the boundary. The revised preflight checks material membership
and near-wall values. Final node and cap-triangle admission is unchanged.
No specimen formula was altered.

Job 5 runs fresh VTK controls, the previously skipped lattice and the
same-skeleton controls. It inherits the measured native and MMG results
from job 4 by SHA-256, checks native binary and dependency identities, and
does not repeat the fine Gyroid or Schwarz solves. Both earlier runs stay
unchanged. The `inherited.json` record points to each original artifact.
