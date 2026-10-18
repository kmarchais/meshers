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

## Audited intermediate MMG meshes

Job 18 repeats the resolution-20 tuned candidate without changing either MMG
call. It archives and independently measures the VTK seed and each MMG output.
The final mesh reproduces job 14. Intermediate audits add work, so this run is
explicitly excluded from comparable workflow timings.

| Stage | Nodes | Tets | Sampled two-sided H | Minimum quality |
|---|---:|---:|---:|---:|
| VTK seed | 6,288 | 26,636 | 0.002190 | 0.00700 |
| First MMG pass, nofem | 3,118 | 10,519 | 0.007183 | 0.06860 |
| Second MMG pass, tuned | 22,928 | 108,979 | 0.005347 | 0.03994 |

This confirms first-pass boundary degradation for this candidate. The seed
already narrowly misses the 0.002 distance gate, and the second pass does not
recover the lost geometry. It does not prove that every extraction grid behaves
the same way. Records are in `runs/g3-microgen-stage-audit-v1`.

The next diagnostic uses resolution 24 and overrides the first call with
`hausd=0.0001`, `hmax=0.06` and `hgrad=1.3`, keeping `nofem=True`. It still calls
Microgen's VTK generation and both MMG commands, but this override is not exposed
by its public remeshing function. The runner labels this modified route
`microgen-vtk-mmg-first-pass-override`, never as the default pipeline.

Job 19 completed that first pass with 346,003 vertices and 1,730,349 tetrahedra.
The stage audit stopped at the unchanged 250,000-tet gate before the second pass
or FEA. The VTK seed itself has H=0.001500 but minimum quality 0.000344.
This is an over-refinement rejection, not an MMG executable failure. The frozen
`result.json` uses the overly broad status `microgen-mmg-failed`; its error text
records the actual intermediate cap. The runner now distinguishes that status
as `intermediate-element-cap`. Raw first-pass meshes stay local and ignored.
Records are in `runs/g3-microgen-firstpass24-v1`.

Job 20 relaxes only the first-pass `hausd` override to 0.0005. The first pass
now has 129,222 tets and the second 168,621, both inside the cap. Final two-sided
distance is 0.001440 and signed volume error is -0.01052%, so the geometry loss
is resolved for this specimen and extraction grid. Minimum quality remains
0.02402. There are 72 retained cap triangles whose adjacent-tet quality ceiling
is below 0.1, making further size tuning insufficient on this seed. Periodic
checks pass. FEA remains withheld. Records are in
`runs/g3-microgen-firstpass24-v2`. These instrumented timings include intermediate
audits and are not used for a workflow speed comparison.

## Nearby VTK seed grids

Jobs 22 and 23 screen Schwarz resolutions 21, 22, 23, 25, 26 and 27 before
spending time in MMG. The first screen failed while serializing nonfinite quality
from zero-volume seed cells. The corrected screen records those cells as a
rejection instead of deleting them or admitting nonfinite metrics.

| Resolution | Zero-volume tets | Caps with quality ceiling below 0.1 |
|---|---:|---:|
| 21 | 538 | 48 |
| 22 | 245 | 33 |
| 23 | 0 | 48 |
| 25 | 459 | 81 |
| 26 | 0 | 96 |
| 27 | 0 | 48 |

None offers a feasible set of unchanged retained caps under the fixed quality
gate. This screen does not rule out other grids or explicit boundary repair.
It does rule out obtaining an admitted competitor merely by tuning MMG on
these six seeds. Records and the runnable driver are in
`runs/g3-microgen-seed-screen-v1` and `runs/g3-microgen-seed-screen-v2`.

## First-pass override on the default Schwarz seed

Job 61 applies the explicit first-pass hausd=0.0005 override to the original
20-point VTK seed, keeping the actual two Microgen MMG calls. It produces 35,084
nodes and 172,800 tets in 15.25 seconds including audit. Distance improves from
0.005347 to 0.002098, but still fails 0.002. Minimum quality is 0.02980 and
41 elements fail 0.1. Volume error is -0.02849%, periodic checks pass, and no
retained cap triangle has an adjacent-tet quality ceiling below 0.1. No FEA runs.
Records are in `runs/g3-microgen-tuned-schwarz20-firstpass-v1`.

Job 62 localizes the worst probes to box-cap curves and the ten worst elements
to tets adjoining a box cap. Its optional detail metrics reproduce every prior
aggregate geometry metric exactly. Jobs 63 and 64 match ten rejected curve
probes to edges retained from the default VTK seed. Job 64 recomputes the probes
from saved mesh connectivity and barycentric weights rather than trusting the
diagnostic coordinates. Endpoint matching errors are below 4.2e-16.

For a worst probe at (0.2479023, 0, 0.5), the Schwarz field stays positive
throughout its coordinate box of radius 0.002. The analytic lower estimate is
0.0005350, from the minimum cosine in each coordinate interval. The Euclidean
ball is contained in that coordinate box. Retaining this chord therefore leaves
an obstruction to the distance gate, regardless of interior MMG sizing. This is
a geometry argument, not FE or local stress certification. Records are in
`runs/g3-microgen-firstpass20-diagnostic-v1`, `runs/g3-microgen-fixed-curve-v1`
and `runs/g3-microgen-fixed-curve-v2`.

Before adding a boundary-repair algorithm, the next bounded experiment screens
the actual VTK seeds at resolutions 28, 30, 34, 36, 38 and 40 for degeneracy,
retained-cap quality ceilings and independent geometry. The specimen and fixed
admission gates remain unchanged. These are tuned extraction grids, not defaults.

Job 65 finds cap-feasible seeds at 30 and 38 points. Neither has a zero-volume
tetrahedron or degenerate cap. Their seed distances are 0.0009615 and 0.0005868,
and cap quality ceilings are at least 0.1448 and 0.1507. Raw seed tet quality
still fails 0.1, so neither is an admitted comparator before remeshing. The
other four grids have cap obstructions, and 28, 34 and 40 also contain
zero-volume tets. Records are in `runs/g3-microgen-seed-screen-large-v1`.

Select the smaller 30-point seed for the next actual two-pass Microgen/MMG
experiment, with the explicit first-pass 0.0005 override. Intermediate audits
enforce the 250,000-tet cap. Their extra cost makes this a diagnostic run, not
a comparable workflow timing.

Job 66's actual two-pass workflow retains the 30-point seed distance, H=0.0009615.
Its first pass has 139,788 tets and q_min=0.10623, passing geometry and quality.
The second pass has 163,887 tets and reduces q_min to 0.08494. The final mesh
therefore fails admission and receives no FEA. Intermediate diagnostics are in
`runs/g3-microgen-firstpass30-v1`.

The next diagnostic retains the measured first-pass output and omits the second
native call. This is explicitly a one-pass Microgen VTK/MMG override, not the
default or a public Microgen API option. The API's output cleanup and periodic
checks remain, followed by the fixed independent admission tests. If admitted,
its full six-load stiffness and live workflow time can expose a stronger
competitor than the failing two-pass controls.

Job 67's explicitly modified one-pass route passes all independent geometry,
quality and periodic gates with 29,713 nodes and 139,788 tets. H=0.0009615 and
q_min=0.10623. It completes six converged loads in 34.98 seconds and a full live
workflow in 44.35 seconds. Worst mixed-load disagreement is 1.13403% against the
justified Schwarz reference, so it still misses the 1% stiffness target. This
is an admitted geometry competitor, not a matched-accuracy winner. Records are
in `runs/g3-microgen-singlepass30-v1`.

Select the other screened cap-feasible seed, 38 points, for the next one-pass
candidate with unchanged MMG settings. Keep the standard two-pass failures and
the explicit omission of its second native call visible in the records.

Job 68's 38-point one-pass alternative rejects at minimum quality 0.09785,
despite H=0.0005868 and matching periodic caps. Its 31,172 nodes and 141,591
tets receive no FEA. Records are in `runs/g3-microgen-singlepass38-v1`.
The next candidate returns to the geometry-admitted 30-point route and reduces
only first-pass hmax from 0.06 to 0.055, seeking the missing stiffness accuracy
without changing the seed, quality floor or other MMG settings.

Job 69's hmax=0.055 candidate rejects quality at 0.08943. It has 31,949 nodes,
153,277 tets and unchanged H=0.0009615. No FEA runs. Records are in
`runs/g3-microgen-singlepass30-h055-v1`.

The next diagnostic replaces the original second level-set call with MMG's
interior optimization, `optim=True, nosurf=True, noinsert=True`. The installed
CLI documents those switches as mesh optimization, no surface modifications,
and no point insertion or deletion. It is an explicit modified two-pass recipe,
not a public Microgen default. The purpose is to repair measured quality defects
without coarsening the refined first-pass mesh or changing its surface. All
output admission gates and the full stiffness target remain fixed.

Job 70 fails all five API retries before optimization. MMG reports that its
optimization option cannot be used with an input metric. The first pass leaves
an automatically discovered sibling `.sol` file. This is an unsupported call
combination, not evidence about optimized mesh quality. Records, including raw
native errors, are in `runs/g3-microgen-interior-optim30-h055-v1`.
The corrected explicit optimization variant hashes and removes only that fresh
temporary metric before calling MMG. No frozen input or saved output is removed.

Job 71's corrected optimization runs but worsens minimum quality from 0.08943
to 0.08610 at the same 31,949 nodes. It remains rejected and receives no FEA.
Records are in `runs/g3-microgen-interior-optim30-h055-v2`. Discard this modified
optimization recipe rather than weaken quality admission. The next candidate
returns to the original one-pass variant at hmax=0.05, with intermediate audits
to enforce the unchanged 250,000-tet cap before any solve.

Job 72's hmax=0.05 one-pass variant passes geometry, quality and periodic gates
with 32,070 nodes and 154,012 tets. H=0.0009615, volume error=-0.01016%,
q_min=0.10701, and all six loads converge. Worst mixed-load disagreement is
0.92806%, below 1% against the justified Schwarz reference. This is the first
tested Microgen-derived variant to meet both admission and the stiffness target.
It is explicitly a modified one-pass recipe, not the default two-pass pipeline.
Records are in `runs/g3-microgen-singlepass30-h050-v1`.

Its 48.45-second total includes intermediate audits, so it is not a comparable
workflow timing. Before repeated timing, test a coarser eight-cell curved Meshers
candidate. The existing 12-cell candidate is much more accurate than the target,
and its FE solve dominates total cost. Independent geometry and stiffness
admission must still pass. Paired live timing must use the same six-load
formulation and solver, not mix the stress-controlled Microgen screen with the
strain-controlled Meshers admission runner.

## Schwarz matched-target timing protocol

Job 73's eight-cell curved Meshers candidate passes independent admission with
6,864 nodes, 3,852 Tet10 elements, H=0.0009637, volume error=-0.0008914% and
q_min=0.47840. All sampled Jacobians are positive and periodic nodes and caps
match. Its complete six-load tensor differs from the justified reference by
0.25870% in the worst mixed load. Its single pilot takes 7.84 seconds total.
Job 74 repeats reverse distance at eight and sixteen subdivisions, retaining
the same admitted maximum. The subdivision difference drops to 0.0001922.
Records are in `runs/g3-schwarz-curved8-pilot-v1` and
`runs/g3-schwarz-curved8-distance-v1`.

The paired experiment declares one warm-up per arm and five live repetitions,
alternating arm order between pairs. It regenerates both specimens and meshes,
audits geometry and periodicity, serializes meshes and displacement fields,
and solves six prescribed engineering strains with the same strict CG helper.
There is no saved-mesh timing replay. Fresh worker processes use one thread,
the hard 8 GiB ceiling, 250,000-tet and 100,000-node caps and 600-second timeout.
Imports and fine-reference computation are excluded consistently from reported
workflow time. Every repetition must meet all admission gates and 1% stiffness
disagreement before any timing summary is accepted. Failed repetitions stop the
experiment rather than disappear from the sample.

The Microgen arm is the explicit one-pass 30-point VTK recipe, hausd=0.0005,
hmax=0.05 and hgrad=1.3. The Meshers arm is eight-cell nonadaptive curved Tet10,
with one optimization pass and native corner tolerance 1.0. Independent final
distance admission remains 0.002. `schwarz_matched.py` records pinned sources,
the native extension hash, counts, per-load residuals, full tensors and paired
time ratios. A result applies only to Schwarz against this modified recipe,
not to the default two-pass pipeline or all four specimens.

Job 75 completes both warm-ups and all five interleaved pairs. Every run meets
the fixed geometry, quality, periodicity, positive sampled Jacobian and 1%
stiffness gates. Meshers has 6,864 nodes and 3,852 Tet10 elements; the explicit
one-pass Microgen variant has 32,070 nodes and 154,012 Tet4 elements. Their
worst mixed-load disagreements are 0.25870% and 0.92806%, respectively.

Median total workflow times are 7.7916 seconds for Meshers and 40.5036 seconds
for Microgen. Paired speedups range from 5.1155 to 5.2908. Median generation
times are 0.06008 and 7.72084 seconds; median six-load solve times are 6.70432
and 31.04057 seconds. Full workflow times range from 7.6838 to 7.8803 seconds
for Meshers and 40.2588 to 41.2876 seconds for Microgen. Both arms use the same
prescribed-strain solver. These are same-G3 repeated measurements, not OMEN
comparisons. Records are in `runs/g3-schwarz-matched-v1`.

This supports a Schwarz matched-target workflow advantage against the explicitly
modified one-pass comparator. It does not establish a four-fixture victory,
repair the default two-pass Microgen pipeline, or establish lattice stiffness
accuracy. Representative measured displacement fields still require the next
independent stress, energy and reaction audit.

Job 76 independently integrates measured fields from both arms. Stress, energy
and reaction tensors agree within 1.32e-11 relative. Periodic displacement jumps,
including curved mid-edge pairs, are at most 5.56e-17. The representative fields
therefore pass the independent FE consistency check. Records are in
`runs/g3-schwarz-matched-field-audit-v1`. This establishes the stated one-case
advantage, not local stress accuracy or a continuum certification.

## Gyroid retained-cap screening

Jobs 103 and 104 screen actual VTK Gyroid sheet seeds at resolutions 24, 28,
30, 32, 36, 38 and 31, 33, 35, 37, 39, 41. The specimen remains offset 1.1,
the original sheet, and no Meshers seed is involved. All twelve seeds contain
fewer than 250,000 tets and no zero-volume tetrahedra. Eleven have retained
cap triangles whose analytic tetrahedral quality ceiling is below 0.1. MMG
interior sizing cannot admit those fixed cap triangulations.

Resolution 31 has no such cap obstruction, with minimum ceiling 0.15055,
22,199 nodes and 95,776 seed tets. Its seed distance is 0.002496, so it is
not admitted as generated. Select one actual two-pass Microgen/MMG run,
explicitly overriding first-pass hausd=0.0005, hmax=0.05, hgrad=1.3. The
second pass uses the existing tuned hausd=0.0005, hmax=0.06, hgrad=1.3.
Stage audits distinguish a retained-curve obstruction from remeshing losses.
This is an instrumented tuned diagnostic, not default Microgen or a comparable
workflow timing. Solve only if all unchanged admission gates pass. Records are
in `runs/g3-gyroid-seed-screen-v1`, `runs/g3-gyroid-seed-screen-odd-v1` and
`runs/g3-gyroid-tuned31-v1`.

Job 105 stops after the first native MMG pass produces 288,870 tetrahedra,
above the fixed 250,000 intermediate cap. No second pass or FEA runs. Peak
RSS is 0.73 GiB, so this is an element-cap rejection, not a memory failure.
Select hmax=0.06 for one successor while leaving all other settings and gates
unchanged. Do not enlarge the cap or report the aborted route as an admitted
competitor.

Job 106 also exceeds the intermediate cap, at 289,317 tets after the first
pass. Increasing hmax from 0.05 to 0.06 does not reduce this count. The next
single-parameter test changes only first-pass hausd from 0.0005 to 0.001.
The independent 0.002 distance gate, second-pass settings and resource caps
remain unchanged. This tests sizing sensitivity, not a relaxed admission gate.
