# Reference-only cubic enrichment on G3

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
