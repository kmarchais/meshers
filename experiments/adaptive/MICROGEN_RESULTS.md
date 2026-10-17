# G3 Microgen workflow results, 17 October 2026

The real Microgen VTK and CLI MMG workflows ran on all four fixed
specimens. This is not the earlier MMG-on-a-Meshers-seed comparison.
There is no demonstrated four-fixture win for Meshers.

The useful negative results are:

- Meshers rejects all tested periodic spinodoids at quality 0.1.
- Meshers rejects the smooth lattice's raw-field periodicity, although
  Microgen produces matching discrete caps and completes diagnostic FEA.
- Curved Tet10 is much more accurate than linear Tet4 on the two admitted
  fixtures, but its strict single-thread CG solve costs more at the tested
  node counts. No matched-accuracy workflow speed comparison was run.
- Default Microgen MMG coarsens the default VTK meshes and worsens boundary
  distance on these fixtures. Every default control fails admission.
  This does not establish superiority over a tuned Microgen workflow.
- Spinodoid and smooth-lattice stiffness errors remain unproven. Neither
  has a justified, admitted fine reference in this experiment.

Native meshing code is unchanged. The branch started at `2e35c7e441b5ccaaa42ea324502ca1a2afabef1b`.
The support checkouts remain at `2b1ea80b36d4ba03842a7d9087638abf3868fe3e`
and `602dc900cd028d205b36e29dcdb5b4a9cf4f581a`.
The exact specimens, defaults, resource limits and metric definitions are
in `MICROGEN.md`. Versions and source/native hashes are in the run manifest.

## Counts, geometry and full-tensor disagreement

H is the independently sampled two-sided maximum in a unit box. Tet10
reverse queries use eight subdivisions per curved edge. The subsequent
verification also checks sixteen. These are finite samples, not continuous
Hausdorff bounds. Volume error is signed. Quality is minimum MMG corner
quality, with separate sampled Jacobian checks for Tet10.

The fixed gates are H <= 0.002, absolute relative volume error <= 1%,
quality >= 0.1, positive sampled Jacobians, and matching periodic nodes
and cap triangles. Only the Meshers Gyroid and Schwarz-solid rows below
pass. The reported Microgen stiffness disagreements are diagnostics.

| Specimen and route | Nodes | Tets | Sampled H | Volume error | q_min | Stiffness disagreement |
|---|---:|---:|---:|---:|---:|---:|
| Gyroid, VTK default | 7,808 | 30,642 | 0.00561 | +1.5190% | 0.0000270 | 14.7909% |
| Gyroid, MMG default | 3,527 | 12,282 | 0.00867 | +2.1768% | 0.000483 | 16.9218% |
| Gyroid, Meshers linear 24 | 17,647 | 72,757 | 0.00158 | +0.1891% | 0.1503 | 3.7152% |
| Gyroid, Meshers curved 12 | 19,931 | 10,832 | 0.000573 | -0.0280% | 0.3714 | 0.2951% |
| Schwarz solid, VTK default | 6,288 | 26,636 | 0.00219 | -0.0029% | 0.006997 | 3.6002% |
| Schwarz solid, MMG default | 3,151 | 10,829 | 0.00846 | -0.0551% | 0.05484 | 5.9304% |
| Schwarz solid, Meshers linear 24 | 19,004 | 90,792 | 0.00189 | +0.0086% | 0.3997 | 1.3458% |
| Schwarz solid, Meshers curved 12 | 19,628 | 11,976 | 0.000242 | -0.0002% | 0.4657 | 0.0738% |
| Spinodoid, VTK default | 5,268 | 20,702 | 0.0301* | -0.0053% | 0.00000000118 | No justified reference |
| Spinodoid, MMG default | 5,509 | 22,907 | 0.0354* | +0.0160% | 0.001835 | Solve withheld |
| Smooth lattice, VTK default | 48,353 | 233,532 | 0.00265* | -0.4017% | 0.000000232 | No justified reference |
| Smooth lattice, MMG default | 6,709 | 27,717 | 0.0106 | -2.8544% | 0.009052 | No justified reference |

An asterisk marks unresolved forward projections. There are 99 for
spinodoid VTK, 23 for spinodoid MMG and 288 for lattice VTK. Their H values
are diagnostics, not trustworthy complete maxima. The geometry gate fails.
Four Sobol volume scrambles have full spreads of 0.000191 for Gyroid,
0.000427 for spinodoid and 0.000156 for the lattice. The Schwarz-solid
reference volume is exactly 0.5. These spreads are not confidence intervals.
In particular, the tiny reported fine Gyroid volume bias is not resolved
against the sampling uncertainty.

Meshers spinodoid rejects linear 24 and curved 12/16/20 before returning
meshes. Reported minimum qualities are 0.03438, 0.03344, 0.03983 and
0.01583. The unchanged floor remains 0.1. Microgen's MMG result passes its
own node-periodicity check, but the independent cap-triangle check fails
on axis 1, so no homogenized tensor is admitted or solved for that mesh.

All four Meshers lattice configurations reject with
`field/grading is incompatible with periodic axis 0`. Sequential cubic
smooth unions are not associative. Opposite-face scalar differences reach
about 0.0000764 in void. A separate 131,072-probe face check found zero
material-membership differences and zero field differences within 0.03
of the material boundary. The Microgen meshes pass the actual periodic
node and cap-triangle checks. No field rewrite or gate relaxation was used
to make Meshers accept this specimen.

## Same-host wall time

Seconds are single live runs on G3, with one numerical thread. Generation
includes specimen setup and actual VTK/MMG or native generation.
Meshing-only time is `generation_s - specimen_setup_s` in each raw record.
Total time includes geometry audits, serialization, conversion, assembly
and six strict CG solves. It excludes interpreter imports. Fine reference
solves are not charged to a candidate's workflow. These are not repeated
performance estimates and must not be compared with OMEN seconds.

| Specimen and route | Generation, s | Six-load FEA, s | Total workflow, s |
|---|---:|---:|---:|
| Gyroid, VTK default | 0.128 | 12.387 | 13.478 |
| Gyroid, MMG default | 4.600 | 1.977 | 7.302 |
| Gyroid, Meshers linear 24 | 1.245 | 14.836 | 17.480 |
| Gyroid, Meshers curved 12 | 0.282 | 45.558 | 47.080 |
| Schwarz solid, VTK default | 0.068 | 4.365 | 5.236 |
| Schwarz solid, MMG default | 3.050 | 1.096 | 4.819 |
| Schwarz solid, Meshers linear 24 | 0.649 | 13.525 | 15.210 |
| Schwarz solid, Meshers curved 12 | 0.139 | 36.688 | 38.073 |
| Spinodoid, VTK default | 0.080 | 11.466 | 14.846 |
| Spinodoid, MMG default | 4.035 | Withheld | Incomplete |
| Smooth lattice, VTK default | 0.449 | 190.207 | 196.387 |
| Smooth lattice, MMG default | 12.623 | 5.086 | 18.669 |

The short default MMG workflows use fewer nodes and fail the accuracy
gates. Their times are not matched-accuracy wins. Conversely, the curved
Meshers results do not prove a faster end-to-end workflow than Microgen.
Only Jacobi-preconditioned CG was tested, not a direct or multigrid solver.

## References and geometry/order controls

The score is the maximum absolute generalized eigenvalue of
`C_candidate - C_reference` against `C_reference`, across all strain
combinations. It is not a maximum of six diagonal tests.

The reference extrapolates curved 16/20 tensors with assumed order four.
For Gyroid, the fine-to-extrapolated disagreement is 0.04273%, the
12/16 versus 16/20 extrapolation difference is 0.01665%, and the archived
linear-64/80 extrapolation differs by 0.00716%. For Schwarz solid, these
first two values are 0.01300% and 0.01296%. Its archived quadratic tensor
agrees to numerical precision, which checks reproduction but is not an
independent-method reference. Both fine inputs pass the geometry gates.
No value is a certified continuum error bound.

Saved 12-cell skeletons provide separate approximation-order and geometry
controls. They have identical connectivity and material for the linear
versus straight-quadratic comparison. Curved versus straight quadratic
then changes only the geometric mapping at fixed polynomial order.

| Specimen | Linear corners | Straight Tet10 | Curved Tet10 |
|---|---:|---:|---:|
| Gyroid, stiffness disagreement | 13.2473% | 1.6762% | 0.2951% |
| Schwarz solid, stiffness disagreement | 5.4935% | 0.2645% | 0.0738% |

Relative to the fine reference, the worst magnitudes of the tensor-change
eigenvalues are 11.5777% for Gyroid's order change and 1.3987% for its
geometry change. For Schwarz solid they are 5.3065% and 0.2256%.
The signed tensor changes add exactly. Their worst scalar magnitudes do
not generally add, and neither is a certified discretization-only error.
The straight controls fail the common geometry gate and remain diagnostic.
Their replay timings are not presented as live generation times.

## Evidence and checks

Jobs 4 and 5 finished sequentially and both passed their progress audits.
The first full pilot used at most 2.37 GiB sampled worker RSS. No worker
hit its 8 GiB address-space or aggregate RSS ceiling, 1,200-second timeout,
250,000-tet cap, or 100,000-node solve cap. No unrelated jobs were changed.

There are 28 final configurations, including controls and fine references:
eight admitted solves, eleven unqualified diagnostic solves, eight native
rejections, and one MMG topology rejection. Nineteen complete six-load
tensors are saved. Public native Tet10 ordering also passes a real solid
cube affine-patch solve. All numerical solver residuals are below 2e-9.

Job 7 rechecked all twenty saved meshes and nineteen tensors. It verified
frozen/inherited source hashes, positive linear orientation, exact counts,
periodic pairs for admitted meshes, mixed-load eigenvalues and solver
records. Independent five-point Tet10 volume integration agrees with the
original fifteen-point result within 7.8e-15. Sixteen-edge boundary
tessellations retain admission for all six fine/coarse curved meshes.
Every original artifact hash is unchanged. Historical saved displacement
archives were not transferred and were not revalidated.

The first verifier attempt, job 6, stopped after two controls because its
support import path was not initialized. It changed no experiment files.
The corrected read-only verifier passed as job 7.

Persistent local artifacts are at:

- `/home/kmarchais/dev/meshers/experiments/adaptive/runs/g3-microgen-defaults-v3`
- `/home/kmarchais/dev/meshers/experiments/adaptive/runs/g3-microgen-defaults-v4`
- `/home/kmarchais/dev/meshers/experiments/adaptive/g3-verification.json`
- `/home/kmarchais/dev/meshers/.audit/g3-microgen.tsv`

Text records and hashes are committed. Mesh archives and installed-source
snapshots remain local. No historical artifact or admission flag was edited.

To repeat the read-only verification on G3, keep the local mesh archives,
source snapshots and support checkouts, then run from the repository root:

```sh
qrun -l 'Verify saved G3 geometry tensors and hashes' --progress report \
  --interrupt none --interrupt-reason 'Read-only verification with one final report' \
  --check-start --wait env \
  LD_LIBRARY_PATH=/home/kmarchais/dev/meshers/.local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu \
  OPENBLAS_NUM_THREADS=1 timeout 300 \
  .venv/bin/python experiments/adaptive/verify_microgen_workflow.py \
  --run experiments/adaptive/runs/g3-microgen-defaults-v4 \
  --output /tmp/g3-microgen-verification.json
```

Set the library path before Python starts. A clone with only the committed
text evidence cannot replay the mesh checks without those local archives.

## What remains open

This is a reviewable diagnostic result, not completion of a four-fixture
matched-accuracy study. Spinodoid and smooth-lattice references are missing.
Projection failures also leave some geometry distances unresolved. A tuned
Microgen accuracy/cost comparison, repeated workflow timings and local
stress recovery have not been measured. The experiment does not justify
batching changes to the mesher or claiming continuum certification.

Kevin's later requested order remains unchanged. Direct fedoo Tet10 now
has G3 execution evidence, but public Microgen Tet10 conversion and its
periodic curved-node integration still need a separate integration result.

## Retained cap triangles limit quality

A subsequent read-only audit, job 10, explains part of the default Microgen
quality failure. The VTK seeds contain thin box-cap triangles that the MMG
workflow retains. Their shape alone limits the quality of any adjacent
tetrahedron, regardless of its fourth vertex or the requested MMG size.

For triangle area A and sum of its three squared edge lengths S, the quality
definition used in this study satisfies `q_tet <= 4*sqrt(3)*A/S`. To derive
the bound, place the fourth vertex at height h above the face and at planar
distance r from its centroid. The tetrahedron volume is A*h/3, and the sum
of all six squared edge lengths is `4*S/3 + 3*h*h + 3*r*r`. Quality is
maximized at r = 0 and `h*h = 2*S/9`, which gives the stated bound. Additional
constraints on the fourth vertex can only lower that maximum.

| Specimen | Retained cap triangles below the 0.1 quality ceiling | Smallest ceiling |
|---|---:|---:|
| Gyroid | 60 | 0.009126 |
| Schwarz solid | 0 | 0.138143 |
| Spinodoid | 46 | 0.009166 |
| Smooth lattice | 96 | 0.017177 |

The audit matched every seed cap triangle in the saved MMG mesh at coordinate
tolerance 1e-10. Spinodoid MMG also has one extra cap triangle, consistent with
its earlier topology rejection. The analytic bound was checked numerically
on regular and thin triangles, including a fourth vertex that attains it.
All eight input mesh hashes are unchanged. The job finished before the
progress-start audit could observe it; its five result records are available.

Thus tuning only MMG sizes while retaining these particular seed caps cannot
pass the existing quality gate for Gyroid, spinodoid or lattice. Schwarz has
no obstruction from this bound, though its actual mesh still fails quality.
The bound does not rule out a different VTK extraction grid or a separately
identified boundary-repair workflow. It says nothing about stiffness error.

The executable audit is `cap_quality_bound.py`; its output and input hashes
are in `g3-cap-quality-bound.json`. Run it through the shared queue with
`--source experiments/adaptive/runs/g3-microgen-defaults-v4` and a fresh
`--output` JSON path.

## Tuned Schwarz baseline pilot

Jobs 11 through 14 tested one candidate at a time using the actual Microgen
VTK generation and two-pass MMG function. Each completed result was read before
selecting the next candidate. The specimen remains Schwarz P solid at offset
zero. These are explicitly tuned runs, with `hausd=0.0005`, `hmax=0.06` and
`hgrad=1.3`, and extraction resolutions 28, 32 and 20 points per axis.
Microgen applies the size arguments to its second pass. The first pass still
receives only `nofem=True`, as recorded in every call list.

At resolution 28, all five built-in first-pass attempts failed. Job 11 exposed
the failure, but the Microgen wrapper hid the captured MMG text. Job 12 repeated
that candidate with exception-output capture. MMG reports a worst element of
zero quality and stops in `MMG5_minQualCheck`. No cells were deleted or repaired.

| Resolution | Output nodes | Output tets | Sampled H | Signed volume error | Minimum quality | Outcome |
|---|---:|---:|---:|---:|---:|---|
| 28 | None | None | Unavailable | Unavailable | Unavailable | MMG first-pass failure |
| 32 | 28,472 | 131,516 | 0.005800 | -0.01289% | 0.02092 | Geometry and quality rejection |
| 20 | 22,928 | 108,979 | 0.005347 | -0.07307% | 0.03994 | Geometry and quality rejection |

Both completed meshes pass periodic node and cap matching and have no unresolved
projection samples. Resolution 32 has 48 cap triangles whose adjacent-tet quality
ceiling is below 0.1. Resolution 20 has none, but its actual element quality still
fails. Their complete geometry-audit workflows take 12.98 and 8.80 seconds on G3.
FEA is withheld because admission fails, so these are not total meshing-plus-FEA
times and no new stiffness errors are reported.

This bounded pilot does not find an admissible tuned competitor. It does not
rule out other extraction grids, controls applied to the first MMG pass, or
separately identified periodic boundary repair. The actual first-pass geometry
has not been archived or audited, so its contribution to the final boundary
error is still a hypothesis.

The runner is `microgen_tuned.py`. It uses one thread, an 8 GiB worker ceiling,
a 600-second worker timeout, a 250,000-tet cap and a 100,000-node solve cap.
Its manifests link the frozen reference inputs by hash. Records are under
`runs/g3-microgen-tuned-schwarz28-v1`, `runs/g3-microgen-tuned-schwarz28-v2`,
`runs/g3-microgen-tuned-schwarz32-v1` and `runs/g3-microgen-tuned-schwarz20-v1`.
All runs preserve the same independent admission gates.
