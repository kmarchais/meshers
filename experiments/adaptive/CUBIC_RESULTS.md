# Reference-only cubic enrichment on G3

## Cubic geometry sensitivity, 20 October 2026

Job 309 separates a change in boundary mapping from another equilibrium
enrichment. The existing cubic displacement space embeds a quadratic physical
domain. The probe retains its 222,331 nodes, 44,976 elements, connectivity,
corner coordinates, and all six saved displacement coefficient arrays. It
projects only 30,972 higher-order wall nodes onto the fixed implicit boundary,
locking box-plane coordinates. It introduces no native Meshers Tet20 API.

Three focused controls pass: cubic basis identities, independent cubic surface
interpolation against the registered reference element, and rejection of an
implicit linear-geometry fallback for nonstandard elements. The experimental
geometry auditor now accepts an explicitly registered cubic element, cubic
boundary interpolation, and positive-rule volume quadrature. Default Tet4 and
Tet10 paths remain unchanged. Baseline distance and volume reproduce the saved
quadratic-domain report within 1e-10.

Both old and projected domains pass the unchanged sampled admission gates:

| Measurement | Embedded quadratic domain | Projected cubic domain |
|---|---:|---:|
| Two-sided sampled distance | 0.00054354 | 0.00013101 |
| Volume | 0.24823394598 | 0.24823372775 |
| Minimum corner quality | 0.34102 | 0.34102 |
| Minimum sampled Jacobian | 3.6811e-6 | 2.2608e-6 |
| Failed forward projections | 0 | 0 |

Periodic cubic nodes and full ten-node cap faces match. Maximum coordinate
mismatch is below 1e-15. The largest nodal movement is 0.00039995; median
movement is 3.24e-7. This is not an equilibrium or accuracy result.

The saved fields are then integrated without changing their nodal coefficients.
Positive 64-point and 125-point rules pass affine elasticity controls. The
125-point mixed-load energy change between domains is 0.0105922%, below the
predeclared 0.05% geometry-signal screen. All six signed changes are negative.
The largest 64-to-125-point change is 9.965e-7 relative. Baseline energy
reproduces the earlier saved-equilibrium quadrature difference, 6.869e-6.
Projected-field stress and energy disagree by 1.11e-5 relative; those fields
are not newly equilibrated.

Geometry takes 18.64 seconds at 0.475 GiB, and reintegration takes 105.57
seconds at 0.531 GiB. Do not run a costly projected-domain equilibrium solve
from this sub-threshold signal alone. This does not bound re-equilibrated
geometry error or continuum error. The original cubic result and all 24
empirical models remain unchanged. Before another reference solve, inspect
one finer native parent and its exact cubic node budget for an independent
higher-order refinement level.

Records are in `runs/g3-lattice-cubic-geometry-sensitivity-v1`; the maintained
runner is `cubic_geometry_probe.py`.

## Finer reference resource preflight, 20 October 2026

Job 310 generates one nonadaptive curved25 parent with the same symmetric
lattice representation and positive void clamp used by the original cubic24
parent. It passes admission at 75,113 nodes and 48,696 elements: sampled
distance 0.00038198, minimum corner quality 0.32248, positive sampled
Jacobians, no unresolved projections, and matching periodic caps.

Exact topological counting predicts 241,555 cubic nodes, inside the existing
250,000-node stage. The count is corner vertices plus twice the unique edge
count plus the unique triangular face count. It reproduces the old cubic24
count of 222,331 and passes single-tet and shared-face checks. The assembly
runner now checks and records this budget before constructing the enrichment.

Job 311 builds the actual cubic25 mesh and matrix without solving. Its
241,555 nodes match the prediction. Matrix dimensions are 724,665 squared,
with 99,084,987 nonzeros and 1,191,918,508 CSR bytes. The worker takes 58.35
seconds and peaks at 5.818 GiB, below the same 8 GiB ceiling. This is resource
feasibility, not convergence or a workflow performance claim.

One fresh six-load parent25 control is selected next. Its fields are needed
to verify exact quadratic-to-cubic embedding and the later enrichment energy
identity. Read and independently audit that parent before selecting a cubic25
equilibrium solve. No reference model is removed or replaced.

Records are in `runs/g3-lattice-cubic-parent25-preflight-v1`,
`runs/g3-lattice-cubic25-assembly-memory-v1`, and
`runs/g3-lattice-cubic-parent25-fields-v1`.

Job 312 completes all six fresh parent loads in a 128.64-second workflow.
Its arrays reproduce the preflight exactly and reduced residuals are below
9.35e-13. Independent audit 313 reproduces reactions within 7.68e-12 relative
with zero periodic jumps. The source is ready for an enrichment check.

The first fresh cubic25 solve, job 314, fails before completing any load.
Exact geometry embedding agrees within 1.67e-16 and embedded parent energy
within 9.08e-12 relative, but the full `A - A.T` symmetry check requests a
1.47 GiB temporary and hits the unchanged 8 GiB address-space ceiling. Peak
RSS is 6.068 GiB. Assembly fitting in memory did not prove solve feasibility.
The failed record remains in `runs/g3-lattice-cubic25-reference-v1`.

The research AMG helper now computes the same complete relative Frobenius
symmetry norm in blocks of at most 32,768 rows. The threshold stays 1e-12,
every row is checked, and complete matrix hashes still guard reuse. Ten
focused tests compare against the full norm, detect off-diagonal asymmetry,
and prohibit whole-matrix subtraction on the forced-block path. Job 315
passes twelve integrated cube and Schwarz loads with the same equations and
solver tolerances. Cube disagreement is 3.45e-12; Schwarz reproduces its saved
tensor within 4.42e-11 worst mixed-load disagreement. Neither is a large-case
memory proof.

Only after reading those controls, one fresh cubic25 retry is selected in
`runs/g3-lattice-cubic25-reference-blocked-symmetry-v2`. No prefix is replayed.
The 8 GiB ceiling, 250,000-node stage, 240-second first-load deadline, and
1,800-second worker timeout stay fixed. Do not infer convergence or reference
accuracy until all six loads, own-rule checks, and positive-rule sensitivity
are complete. Control records are in `runs/g3-blocked-symmetry-amg-controls-v1`.

## Completed cubic25 refinement and unresolved models

Job 316 completes all six fresh cubic25 loads, without prefix replay. The
first recorded solver call takes 185.43 seconds, and the full first-load
continuation gate passes within 240 seconds.
Final solve time is 1,020.87 seconds; worker time is 1,071.18 seconds and
peak RSS is 6.753 GiB. The original 8 GiB ceiling is unchanged. This confirms
large-case feasibility of blocked symmetry verification, not a workflow-speed
comparison with the cubic24 run that replayed four verified fields.

Own-rule stress and energy reproduce reactions within 1.03e-11 relative.
Periodic displacement jumps are below 2.23e-16. The 125-point energy shift
is 2.6354e-5 relative, below the unchanged 1e-4 limit. The minimum eigenvalue
of the parent-minus-cubic tensor is +0.000140792, consistent with same-domain
variational enrichment. Status remains `solved-needs-reference`.

Fresh-process audit 317 passes in 62.11 seconds at 0.555 GiB. Its own-rule
reaction disagreement is 9.12e-12, and the 64-to-125-point energy change is
6.8313e-6 relative. All sampled Jacobians and affine controls pass. Every
earlier empirical model is retained.

| Compared tensor | Worst mixed-load difference relative to cubic25 |
|---|---:|
| Cubic24 | 0.248177% |
| Finest adaptive Meshers Tet10 | 0.651535% |
| Eight-round curve-graded Meshers Tet10 | 0.483742% |
| Strongest Microgen-derived straight Tet10 | 0.917416% |

Cubic24 is softer than cubic25 in all six generalized modes. Two neighboring
levels do not establish an asymptotic reference. Candidate ranking also depends
on reference choice: the curve pilot is closer to cubic25 than the finest
adaptive mesh, despite being worse across the original model envelope.
Do not select cubic25 alone to manufacture a target-qualified timing comparison.

The separate model-normalized score uses the original benchmark denominator,
not the reverse model-versus-cubic diagnostic. Cubic25 itself spans 0.114924%
to 1.920819% across all 24 models and fails their 1% envelope. The earlier
candidate model envelopes remain exactly unchanged. No model is deleted or
reclassified as a certified continuum tensor.

Before a third higher-order solve, job 318 preflights one neighboring parent23.
It admits 68,095 nodes and 43,844 Tet10 elements, with distance 0.00059463,
minimum quality 0.15825, positive sampled Jacobians, no unresolved projections,
and matching periodic caps. Exact topology counting gives 218,508 cubic nodes,
inside the existing stage. A fresh parent23 six-load control is selected next.
This permits a third refinement point; it does not guarantee a reliable
three-point extrapolation or settle sensitivity to nonnested meshes.

Records are in `runs/g3-lattice-cubic25-reference-blocked-symmetry-v2`,
`runs/g3-lattice-cubic25-assessment-v1`,
`runs/g3-lattice-cubic25-reference-envelope-v1`, and
`runs/g3-lattice-cubic-parent23-preflight-v1`.

The fresh parent23 control, job 319, completes all six loads in a 119.43-second
workflow and exactly reproduces the preflight arrays. Independent audit 320
has reaction disagreement 6.57e-12 and zero periodic jumps. Assembly-only
job 321 confirms 218,508 cubic nodes, 43,844 elements, and 89,362,584 matrix
nonzeros. Its worker takes 52.91 seconds and peaks at 5.014 GiB. Only after
these gates pass, one fresh cubic23 solve is selected as job 322 in
`runs/g3-lattice-cubic23-reference-v1`. No prefix is replayed and no reference
model is removed. Read its actual completion and fresh-process audit before
interpreting a three-level sequence.

This experiment addresses lattice reference uncertainty. It does not add a
Meshers Tet20 product API, change the Microgen pipeline or establish a timing
advantage. All existing uniform and adaptive reference models remain retained.

## Analytic feasibility

Job 178 registers a degree-three tetrahedral Lagrange element through fedoo's
existing element interface. It checks all cubic monomials and their derivatives,
nodal interpolation, partition of unity and degree-four quadrature moments.
The cubic basis reproduces a nonaffine quadratic geometry mapping within
5.56e-16, and its derivatives within 1.34e-15.

A periodic solid patch contains 768 tets. The existing Tet10 control has
1,241 nodes and the enriched mesh 3,925. Both solve six prescribed strains.
Worst mixed-load stiffness errors against the analytic tensor are respectively
1.65e-12 and 4.75e-12. Independent positive 125-point reintegration agrees with
reactions and energy within 1.46e-11 absolute. Cubic periodic jumps are below
2.23e-16. This validates the experimental element and the shared helper's
unchanged default Tet10 route, not a lattice reference.

Records are in `runs/g3-cubic-reference-patch-v1`. Nine existing comparison
admission regression checks also pass after the helper extension. Unsupported
node counts no longer silently select Tet4 or Tet10 in these helpers; an
experimental element must be explicitly supplied.

## Measured lattice assembly stage

The admitted `g3-lattice-curved24-symmetric-v1` mesh has 9,907 corner nodes,
59,118 edges and 94,188 faces. Cubic enrichment therefore needs 222,331 nodes
on the same 44,976 tets, inside the existing 250,000-node stage.

Job 179 constructs this mesh, verifies periodicity and assembles its matrix
without solving. Peak RSS is 5.893 GiB under the unchanged 8 GiB ceiling.
The 666,993-square CSR matrix has 91,401,219 stored entries and occupies
1,099,482,604 bytes. This resource stage takes 54.81 seconds. It is not an
end-to-end meshing benchmark. Records are in
`runs/g3-lattice-cubic-assembly-memory-v1`.

## Initial reference solve stopped

Job 180 uses that prepared mesh. Before solving, it verifies the source
quadratic mapping at every cubic node and embeds the saved parent fields.
Mapping disagreement is 1.67e-16; integrated parent stress and energy reproduce
the archived tensor within 2.17e-11 relative. Periodic checks pass.

The run has an 8 GiB ceiling and a hard 1,800-second worker timeout. Its
240-second post-first-load continuation gate becomes unreachable while the
first load is still running. At 328 seconds after the preceding status write,
no load has completed. Stop the run rather than wait for a result that cannot
pass that gate. Final worker wall time is 381.19 seconds and peak RSS is
5.742 GiB. No lattice reference tensor or displacement archive is obtained.

Pueue cancellation initially leaves this known controller and worker alive.
Their exact argv and process group are checked, the worker group is terminated
and both known PIDs are confirmed gone before the next heavy job starts. No
other jobs are touched. Frozen partial results remain unchanged; the separate
stop record is `runs/g3-lattice-cubic-reference-jacobi-stop-v1.json`.

## Cached AMG validation and retry

PyAMG 5.3.0 is already installed and pinned in the G3 runtime requirements.
Job 181 tests smoothed-aggregation V-cycle preconditioning with three reduced
translational candidates. The original CG tolerances remain `rtol=1e-12` and
`atol=1e-14`, with independently recomputed residuals below 2e-9. Reuse is
allowed only when the complete reduced CSR matrix hashes identically.

Six-load cubic cube and saved curved Schwarz checks both pass. Cube analytic
stiffness disagreement is 4.42e-12; Schwarz reproduces its archived tensor
within 4.74e-11. Independent own-rule stress and energy disagreement remains
below 3.33e-12, and periodic jumps below 2.23e-16. Each case builds one
hierarchy and reuses it for the remaining five loads. Records are in
`runs/g3-amg-reference-validation-v1`.

The small Schwarz solve takes 6.80 seconds, versus the archived diagonal-CG
control's 6.74 seconds. This is not evidence of a speed advantage. The purpose
is to test whether the large reference becomes feasible. Benchmark defaults
are unchanged, and any future workflow comparison must use the same solver
policy on both meshers and charge hierarchy setup.

Job 182 tests the same admitted cubic lattice mesh with this validated solver.
Its interruptible 240-second first-load deadline fires inside a multigrid
coarse correction. Peak RSS is 6.976 GiB. No load completes, so there is still
no new lattice tensor or accuracy claim. The worker exits normally through
the deadline exception, without an orphaned process. Records are in
`runs/g3-lattice-cubic-reference-amg-v1`.

## Isolated system and node-block diagnosis

PyAMG's installed documentation states that its default block Gauss-Seidel
smoother becomes point-wise Gauss-Seidel on CSR input. The prior solver feeds
CSR in component-major order. Test actual node-block coupling before changing
tolerances or increasing the reference budget.

Job 183 captures the actual first-load reduced matrix and RHS through fedoo's
solver callback. It raises a dedicated capture exception rather than returning
a fake solution. The 661,029-square matrix has 91,117,159 entries. Its free
unknowns form complete three-component blocks on 220,343 nodes, verified from
fedoo's free-DOF metadata. Capture peaks at 5.850 GiB and takes 65.76 seconds.
The matrix and vectors remain local ignored archives, with hashes in
`runs/g3-lattice-cubic-reduced-system-v1/result.json`.

Job 184 permutes this unchanged system to node-major order, uses BSR blocks of
size three and retains the three translational candidates. It reaches 185 CG
iterations but does not converge before the same 240-second algorithm cap.
The hierarchy has levels of 661,029, 2,160 and 30 unknowns; setup takes 9.90
seconds. Peak RSS is only 3.048 GiB outside the FE assembly process. Reordering
alone does not establish reference feasibility. Records are in
`runs/g3-lattice-cubic-block-amg-v1`.

Job 185 adds the three rotational candidate fields used by PyAMG's own
`gallery.linear_elasticity_p1` example, retaining the original matrix, RHS,
CG tolerances and deadline. These are algebraic near-nullspace candidates,
not exact periodic rigid modes or modified physical constraints. The original
CSR residual passes after undoing the permutation. The first system converges
in 117 iterations and 165.05 seconds of algorithm wall time. Its original-system
relative residual is 7.72e-13. Hierarchy setup takes 16.19 seconds, CG 146.85
seconds and peak RSS is 2.909 GiB. This is still one load, not a full mixed-load
reference or a live workflow speed result. Records are in
`runs/g3-lattice-cubic-block-amg-rigid-v1`.

## Six-load block integration and full reference retry

Job 186 verifies the integrated block solver on the same cubic cube and curved
Schwarz controls. It passes all twelve solves and reuses the hierarchy for five
loads in each case, conditional on an identical reduced matrix hash. Cube
analytic disagreement is 3.45e-12 and Schwarz reproduces its archived tensor
within 4.42e-11. Independent own-rule stress and energy error remains below
3.16e-12, with periodic jumps below 2.23e-16. No small-case speed advantage is
established. Records are in `runs/g3-elastic-block-amg-validation-v1`.

That validation's manifest correctly records `elastic_blocks=true`, as do its
solver records, but inherits a stale scalar-method description. The frozen
description is not silently rewritten. Its correction is recorded separately
in `runs/g3-elastic-block-amg-validation-v1-method-clarification.json`, and the
maintained controller now writes the selected method explicitly.

Job 187 retries all six lattice loads through the validated integrated solver.
The actual mesh, coefficients, boundary conditions, residual tolerances,
240-second first-load deadline, 1,800-second worker cap and 8 GiB ceiling remain
unchanged. Completed fields are saved before stress/energy and 125-point audits.
The parent embedding check passes again. The first load also passes, with the
same 117 iterations and 7.72e-13 residual as the isolated system. The second
load fails before solving, while allocating a repeated `A - A.T` symmetry
check under the 8 GiB address-space cap. The rejected allocation is 695 MiB
for sparse indices; observed peak RSS is 7.105 GiB. There is no full tensor or
reference admission. Records are in
`runs/g3-lattice-cubic-reference-elastic-block-v1`.

## Hash-guarded symmetry reuse

The matrix is already hashed completely before hierarchy reuse. Check symmetry
once for the first matrix; on subsequent loads, require an identical shape,
CSR structure and data hash before reusing that verified property. A changed
matrix still fails rather than silently reusing the certificate. Independent
residual and stress/energy checks remain per load. This avoids a redundant
large sparse subtraction without changing any physical or numerical gate.

Job 188 passes all twelve integrated control solves with five verified symmetry
reuses in each case. The tensors and iterations are unchanged. Three focused
reference tests plus the nine comparison-admission regressions also pass.
These cover polynomial identities, exact quadratic-map reproduction, explicit
element selection and preserving all six solved field checkpoints in order.
Records are in `runs/g3-elastic-block-symmetry-cache-validation-v1`.

Job 189 retries the full reference with this allocation fix. Each completed
displacement field is now archived immediately, so a later allocation or audit
failure does not destroy earlier solved fields. The complete reference still
requires all six loads, own-rule stress/energy, periodic jumps, Galerkin
monotonicity and the retained 125-point sensitivity check. It solves and saves
two displacement fields, then fails during load-two reaction extraction. Fedoo
evaluates `M @ A @ X` left-to-right and attempts another large sparse matrix
product; a 349 MiB index allocation fails under the same address-space cap.
Peak RSS is 7.216 GiB. Only one force column reaches the partial result file,
so there is no complete tensor. Both field archives are intact and finite.
Records are in
`runs/g3-lattice-cubic-reference-elastic-cache-v1`.

## Vector-first reactions and verified prefix replay

The reference path now computes the equivalent `M @ (A @ X - D)` using vector
operations. It reproduces fedoo's stored MPC transpose and free/blocked identity
contributions, without forming `M @ A`. The installed dependency and benchmark
default reaction path are not patched. An explicit experimental reaction reader
selects this path for the cubic reference.

Job 190 compares native and vector-first reactions on all six cube and Schwarz
loads. Maximum relative differences are 5.62e-15 and 7.34e-14. Its two saved
prefix fields per case are accepted only after checking the current reduced
equations and reconstructing the same boundary-constrained displacements. The
remaining four loads are solved normally. Tensors, independent stress/energy
and periodic checks pass. These are twelve load checks with four verified
replays, not twelve fresh solves or a live-workflow timing. Records are in
`runs/g3-vector-reaction-prefix-validation-v1`.

Fourteen focused and comparison tests pass, including rejection of a wrong RHS
and rejection of a changed matrix even when the saved vector solves it. The
controller also now hashes the source field archive whenever a saved prefix is
used. Frozen earlier validation manifests are not rewritten.

Job 191 resumes the lattice reference from the two retained fields rather than
pay for their solves again. It verifies the identical prepared-mesh hash, each
current reduced residual and the reconstructed boundary values before accepting
them. Both retained fields pass, including a 9.56e-13 residual for load two.
The first new solve, load three, passes in 116 iterations at 9.37e-13 residual.
The 240-second deadline applies to this first newly solved load, after verified
prefix preparation. The full worker remains capped at 1,800 seconds and 8 GiB.
No field replay time is presented as a fresh meshing or workflow speed result.

It also solves and saves load four. Before load five, another reduced-matrix
allocation fails under the unchanged address-space cap. Four force columns and
four intact fields remain, not a complete reference. Records are in
`runs/g3-lattice-cubic-reference-vector-resume-v1`.

## Completed-problem cache lifetime

The replay wrapper updates its delegated problem only when the solver callback
starts, which is too late to release the previous problem before the next
reduction. Forwarding that binding immediately passes a weak-reference test
and the four-prefix control cases in job 192. It is a real retention fix, but
job 193 shows it is insufficient: replay reaches four fields and the fifth
reduction still fails. Do not claim that this binding fix alone solves the
large-case memory problem.

Inspection of pinned fedoo source also finds a problem/boundary-condition cycle,
including `problem.bc._problem = problem`. Completed problems can therefore
retain their private reduced-matrix cache until cyclic garbage collection runs.
Each cached reduction is about 1.1 GB. The reference helper now explicitly drops
that unused cache after copying the field and extracting reactions. It still
constructs a new problem for the next load, and never clears a matrix while
that problem is solving. Benchmark helper defaults remain unchanged.

A regression using actual fedoo solves and weak references verifies release of
all six completed CSR matrix objects. Sixteen focused and admission tests pass.
Job 194 also passes the two six-load control cases with four verified prefix
fields each. Numeric tensors remain unchanged. Records are in
`runs/g3-reduced-cache-release-validation-v1`.

Job 195 resumes from the four retained fields with explicit per-load cache
release. All four are revalidated against current equations and boundary values
before fresh solves. The 8 GiB ceiling, 1,800-second worker cap, first-new-load
240-second deadline and final independent audits remain unchanged. There is no
complete reference or four-fixture accuracy claim while it is pending. Records
are in `runs/g3-lattice-cubic-reference-cache-release-v1`.

## Complete same-domain cubic tensor

Job 195 completes all six loads. It revalidates four retained fields and solves
the final two in 119 iterations each. All current reduced residuals are below
9.57e-13. Explicit completed-cache release brings peak RSS to 6.341 GiB under
the unchanged 8 GiB ceiling. Worker wall time is 369.81 seconds, including
verified replay. This is not a fresh workflow timing.

Independent own-rule stress and energy reproduce reactions within 7.10e-12
relative. Periodic displacement jumps are below 2.23e-16. Positive 125-point
reintegration changes worst mixed-load energy by 6.87e-6 relative, below the
unchanged 1e-4 quadrature-sensitivity limit. The smallest eigenvalue of the
parent-minus-cubic tensor is +0.000140585, consistent with variational enrichment
on the identical quadratic physical domain. All six fields remain archived.

The tensor passes FE consistency, geometry and quadrature gates. It still has
status `solved-needs-reference`: one higher-order same-domain result does not
establish convergence to the fixed implicit specimen. All 24 earlier empirical
reference models remain retained. The next bounded stage audits these saved
fields in a fresh process and compares every retained model without re-solving.

Job 196 passes that fresh-process check at 0.514 GiB peak RSS in 53.94 seconds.
It uses explicit Tet20 registration rather than a Tet4 fallback. Affine analytic
patches pass for every rule, all sampled Jacobians are positive, and 64-to-125
point fixed-field energy changes by only 9.54e-7 relative. Records are in
`runs/g3-lattice-cubic-assessment-v1`.

Against cubic24, uniform Tet10 meshes at 24, 32 and 34 cells differ by 0.87761%,
0.56789% and 0.49251%. The finest adaptive tensor differs by 0.52086%, in the
opposite energy direction. The actual tuned Microgen Tet4 tensor differs by
3.49233%. These pairwise values do not resolve specimen accuracy. The cubic
tensor normalizes retained-model disagreements of 0.15415% to 1.75050% across
all 24 models. No tested Tet10 candidate passes every retained model at 1%.
Neither reference selection
nor a four-fixture timing victory follows from this result.

Before another reference solve, test a stronger actual Microgen control.
Its previously admitted lattice parent needs 326,716 straight Tet10 nodes.
Second-pass hausd 0.0009 reduces the lift to 267,550 but fails quality at
0.097405. One geometry-only stage changes that second-pass value to 0.001,
retaining VTK51, first-pass 0.0015 and all admission gates. Check exact lift
size and admission before any solve. This is explicitly tuned Microgen,
not its recorded defaults and not MMG on a Meshers seed.

## Worst-mode discrete error localization

Job 202 uses the worst retained-model strain of the finest adaptive Meshers
candidate, dominated by deviatoric normal strain. On the unchanged curved24
physical domain it embeds the parent Tet10 field exactly into the cubic space
and integrates the energy of its difference from the cubic solution. The
own-rule Galerkin identity reproduces the corresponding parent-minus-cubic
quadratic form within 9.50e-9 relative. This is a discrete enrichment check,
not a local stress or continuum certificate.

The highest-energy 20% of coarse tets contain 98.394% of the correction energy.
The largest-volume 20% contain only 48.515%. Elements incident on cap vertices
are 6.990% of the mesh but contain 77.461% of the correction energy. Wall-vertex
incident elements contain 95.416%. Positive 125-point integration changes total
correction energy by 0.05995% relative to the own rule; the 64-to-125 change is
0.00800%. These percentages concern the small correction, not the full tensor's
separate quadrature admission check. Records are in
`runs/g3-lattice-cubic-error-localization-v1`.

This supports one cap-focused sizing pilot, not an inference about the error
location on every refined mesh. The current finest adaptive cap edges still
reach 0.06241 despite median length 0.01080. Start from the admitted 24-cell,
adapt-tolerance 0.004 recipe and add a periodic size callback with target 0.015
on box planes. A smoothstep transitions to maximum length 0.1 over distance
0.05. Keep native budgets, geometry, quality and periodicity gates unchanged.
Check geometry and counts before selecting any solve. No native API change or
oracle field is added to the workflow.

Job 203 passes geometry at 159,147 nodes and 99,604 tets. Job 204 completes
six fresh 15-point loads with the validated explicit AMG preset. Fresh audit
205 reproduces reactions within 7.08e-12 and Hill-Mandel energy within
2.73e-12; periodic jumps are zero. Positive 125-point reintegration changes
full mixed-load energy by 4.25e-7 relative. The cap-sizing callback is positive
and has bitwise-equal opposite-plane traces in its recorded preflight.

Its worst retained-model error is 1.411382%, versus 1.427185% for the
123,570-node untargeted baseline and 1.269882% for the finest Meshers result.
Adding 35,577 nodes produces little accuracy gain and does not meet 1%.
Do not promote this as a workflow improvement or keep tuning cap lengths.
Records are in `runs/g3-lattice-cap-size015-band005-comparison-v1`.

The host suspends for 56 minutes during job 204, confirmed by the system
journal. Its monotonic timer excludes sleep, and qrun records a matching long
reporting gap. All timing claims for that run are excluded by the separate
`runs/g3-lattice-cap-size015-band005-timing-exclusion-v1.json`. Frozen numerical
results remain unchanged. The standard Job reporter is not replaced or exempted
to hide the host interruption.
