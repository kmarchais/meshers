# Difficult periodic specimens on G3

## Selecting the next change

The saved finer sizing pilots put relocation at 0.69 to 1.23 seconds within
26.22 to 37.77 second total workflows. Eliminating that phase entirely would
save about 3%. The current Gyroid volume bias is +0.20% to +0.26%, and Schwarz
volume error is below 0.004%. All pass the fixed 1% volume gate. These costs and
errors remain optimization opportunities, but neither explains the missing
spinodoid and lattice coverage.

## Adaptive spinodoid

Job 26 tests the existing adaptive route on the fixed spectral spinodoid from
the real Microgen comparison. It starts from 16 cells, uses four adaptation
rounds, geometry tolerance 0.002, minimum quality 0.1, minimum length 0.01 and
maximum length 0.1. It fails after 38.09 seconds while unfolding the torus into
the box, reporting `a tetrahedron crosses a periodic plane`. No mesh is returned
or admitted. Records are in `runs/g3-spinodoid-adaptive16-v1`.

The admission runner keeps the 8 GiB worker ceiling, 600-second meshing timeout
and 250,000-tet output cap. It measures geometry and periodic matching before
admission. An admitted result still needs a justified stiffness reference and
the full six-load comparison.

Job 27 rebuilds with a detailed unfolding error and reproduces the failure.
The offending vertex has wall mask 1 and canonical y=0.4999788240285261. Its
tetrahedron uses the image y=-0.5000211759714739 beside three vertices just
inside the lower face. This is a finite crossing, not roundoff at the face.
The rebuild passes 20 core tests. Records are in
`runs/g3-spinodoid-adaptive16-v2`.

The projection helper currently wraps a projected root into the cell before
checking bounds. A periodic sine-wave regression tests whether this permits
a wall vertex to cross a constrained box plane. It also checks that projection
to a nearby root in the same cell remains available.

Job 28 confirms the regression fails under the original helper while 20 other
core tests pass. The fix checks projected coordinates before any periodic wrap,
and applies the same rule to interior relocation targets. Job 29 passes all 21
core tests and reruns the original spinodoid. Unfolding now completes. The next
failure is sampled geometry error 0.007612 above tolerance 0.002, after 44.87
seconds. Admission still fails. Records are in
`runs/g3-periodic-projection-red-v1` and `runs/g3-spinodoid-adaptive16-v3`.

Job 30 passes all 35 focused Python checks and tests a 24-cell seed with the
same adaptation settings. Boundary error falls to 0.002570, still above the
0.002 gate, after 36.84 seconds. Unfolding succeeds again. Records are in
`runs/g3-spinodoid-adaptive24-v1`. The next candidate keeps this seed and lowers
only the minimum adaptation length from 0.01 to 0.002 to test whether the
length floor stops the remaining curve refinement.

Job 31 lowers that floor to 0.002 and reproduces the same error to all reported
digits, so the floor is not the cause. The attempt takes 37.28 seconds.
Records are in `runs/g3-spinodoid-adaptive24-v2`.

Job 32 expands the same 52 Microgen Fourier modes into a traceable scalar
expression and uses the existing embedded compiler. It checks values against
the original Microgen evaluator at 4,096 fixed random points, with maximum
difference 2.23e-16 and identical material classification. Compiled gradients
agree with central differences at 32 points within 1.86e-10 absolute. The
989-node expression takes 0.037 seconds to prepare and verify. Its full attempt
takes 6.46 seconds and reproduces the 0.002570 geometry error. This is a faster
evaluation of the same field, with the same admission rejection. Coefficients
and verification records are in `runs/g3-spinodoid-compiled24-v1`.

## First admitted spinodoid

Job 33 keeps that verified field, minimum length 0.002 and the fixed admission
gates, and starts adaptation from 32 cells. It returns 16,721 nodes and 72,468
tets in 8.75 seconds of generation. Independent two-sided distance is 0.001864,
signed volume error is +0.04458%, minimum quality is 0.2543 and minimum sampled
Jacobian is positive. Periodic node coordinates and cap triangles match exactly.
The full generation and geometry audit takes 14.08 seconds. These are not
meshing-plus-FEA timings. Records are in `runs/g3-spinodoid-compiled32-v1`.

This establishes a working adaptive route for the fixed specimen. A justified
fine stiffness reference and a matched-accuracy Microgen comparison remain open.

Job 34 regenerates that admitted linear mesh and completes six prescribed-strain
loads. The solve takes 13.78 seconds and the full generation, verification,
geometry audit and FEA takes 27.71 seconds. The full anisotropic 6-by-6 tensor
and solver residuals are saved in `runs/g3-spinodoid-linear32-solve-v1`.
There is not yet a justified reference for reporting its stiffness error.

Job 35 tests a curved Tet10 candidate from 24 cells. Its native corner-surface
tolerance is 0.01 while adaptation still targets 0.002. Final admission uses
the original independent curved-distance gate 0.002. The candidate has 105,646
nodes and 66,176 elements, positive sampled Jacobians and minimum corner quality
0.1895. Independent curved distance is 0.002591, so it is rejected before FEA.
Records are in `runs/g3-spinodoid-curved24-v1`.

## Smooth lattice input representation

Job 36 uses the original smooth-union lattice field with its positive void
values capped at 0.02. For any positive constant c, min(phi,c) has exactly the
same negative material domain and zero set as phi. The clamp leaves the field
unchanged near its wall. A 4,096-point executable check also confirms unchanged
material classification and exact equality of nonpositive values. Independent
geometry auditing still uses the original unclamped field and frozen samples.

This representation passes Meshers' existing periodic preflight. From 24 cells,
adaptation then rejects two tetrahedra below quality 0.1, with minimum 0.08209.
No mesh is admitted or solved. Records are in `runs/g3-lattice-void-clamp24-v1`.
The raw field remains unsupported by the strict raw-value preflight.

Job 37 keeps that exact material-domain representation and increases only the
seed to 32 cells. It passes the fixed admission gates and completes all six
loads. The mesh has 48,761 nodes and 222,145 tetrahedra, H=0.001773, signed
volume error -0.43152%, minimum quality 0.10073 and positive sampled Jacobians.
Periodic node coordinates and cap triangulations match exactly. Generation
takes 11.11 seconds, FEA 67.91 seconds and the full verified workflow 84.14
seconds. The quality margin is small. Records and the full tensor are in
`runs/g3-lattice-void-clamp32-v1`.

All four study specimens now have admitted Meshers routes on G3. The lattice
route requires the explicit void clamp and spinodoid uses the verified Fourier
expression. Fine stiffness references for these two specimens remain missing,
so their stiffness errors and matched-accuracy speed rankings remain unknown.

Job 38 tests a smaller curved lattice candidate using the original lattice
meshing path, a 16-cell seed and the same void clamp. It has 43,685 nodes and
24,480 Tet10 elements, minimum corner quality 0.1794 and positive sampled
Jacobians. Independent curved distance is 0.008420, so FEA is withheld.
Records are in `runs/g3-lattice-curved16-v1`. A small volume error alone does
not qualify this candidate as a stiffness reference.

## Absolute box-plane audit tolerance

The geometry audit previously used NumPy's default relative tolerance together
with `atol=1e-12` when classifying box planes. At coordinate 0.5, that also
accepted points about 5e-6 away. `on_box_plane` now sets `rtol=0` explicitly for
both cap classification and locked projection coordinates. Its executable
regression checks a point 1e-6 away is free and a point 1e-13 away is locked.

Job 39 re-audits the saved linear spinodoid and lattice and the two rejected
curved candidates. All four distances and admission outcomes are unchanged.
Job 40 re-audits all four fine Gyroid/Schwarz reference inputs, which also retain
their distances and pass admission. Source mesh and report hashes are unchanged.
Records are in `runs/g3-plane-reaudit-v1` and
`runs/g3-reference-plane-reaudit-v1`. Historical reports remain frozen.

## Locating the curved failures

Job 41 records the ten worst forward-distance probes without changing either
saved mesh. Spinodoid's worst probe is a corner node with field residual
0.003215 and distance 0.002591. Two other corners on that face lie on the wall.
Curving mid-edge nodes cannot repair the off-wall corner left by constrained
linear refinement. The lattice's worst probe is a mid-edge node on an interior
junction face, with distance 0.008420. All three mid-edge nodes on that face are
straight, and its three corners lie on the exact wall. The quadratic validity
fallback therefore leaves the unresolved linear chord. Records are in
`runs/g3-curved-boundary-diagnostic-v1`.

The next spinodoid candidate elevates the admitted 32-cell linear skeleton.
Its 114,659 Tet10 nodes require a staged solve cap above the initial 100,000.
The runner keeps that initial cap as its default and explicitly permits a
150,000-node reference stage. This is a resource stage, not a change to geometry,
volume, quality, periodicity or sampled Jacobian admission. The 250,000-tet cap,
8 GiB worker ceiling and 600-second worker timeout remain in force. The earlier
83,244-node Schwarz Tet10 workflow peaked at 2.37 GiB, which supports this
bounded next stage without increasing the memory ceiling.

Job 42 completes all six loads on that admitted curved spinodoid. Generation
takes 8.90 seconds, FEA 474.13 seconds and the full workflow 490.86 seconds.
Peak worker RSS is 3.81 GiB. The largest recomputed solver residual is
9.95e-13. The admitted linear mesh differs from this tensor by 5.3472% in
worst mixed-load strain energy and 2.6483% in tensor Frobenius norm. This is a
pairwise comparison, not yet an error against a justified fine reference.
Records are in `runs/g3-spinodoid-curved32-v1` and
`runs/g3-spinodoid-linear-curved-pair-v1.json`. The comparison is reproducible
with `stiffness_pair.py`, which checks admission, six converged loads, finite
6-by-6 tensors and positive definiteness before comparing them.

Job 43 tightens only adaptation tolerance from 0.002 to 0.0015 and audits the
candidate before selecting a solve. It passes with 147,894 nodes and 94,041
Tet10 elements, H=0.0001295 and q_min=0.2402. Records are in
`runs/g3-spinodoid-curved32-tight-v1`.

The next reference solve retains the 150,000-node and 8 GiB ceilings and uses an
explicit 1,200-second worker stage. The earlier 114,659-node solve needed 474
seconds for FEA, so the initial 600-second stage may be too short for this finer
candidate. The runner's default remains 600 seconds. Completed load columns
and residuals are checkpointed in `result.json` with status `solving`. Only a
completed six-load result receives the full stiffness tensor and solved status.

## Spinodoid reference refinement

Job 44 completes the tighter admitted candidate with 147,894 nodes and 94,041
Tet10 elements. Generation takes 9.96 seconds, FEA 678.37 seconds and the full
workflow 697.79 seconds. Peak RSS is 4.85 GiB. All six recomputed residuals are
below 1e-12. The tighter and preceding fine tensors differ by 0.020306% in the
worst mixed load. Records are in `runs/g3-spinodoid-curved32-tight-solve-v1` and
`runs/g3-spinodoid-fine-pair-v1.json`.

Job 45 independently integrates saved displacement fields from both fine
spinodoid meshes, the admitted linear spinodoid and the admitted linear lattice.
Stress, energy and reaction tensors agree within 2.03e-11 relative. All periodic
displacement jumps, including mid-edge pairs, are zero in these checks. Input
hashes are unchanged. This verifies FE integration, not local stress accuracy.
Records are in `runs/g3-new-fixture-field-audit-v1`.

`adaptive_reference.py` applies the declared length proxies and assumed orders
to the two audited fine tensors. The finer tensor is the primary empirical
reference. Its difference from the six extrapolation models ranges from
0.02611% to 0.10714%. The coarser curved candidate's difference from those
models ranges from 0.04643% to 0.12747%, below the 1% target throughout. The
admitted linear candidate ranges from 5.39348% to 5.47426% and fails that target.
These are empirical model sensitivities, not certified continuum error bounds.
The 490.86-second curved workflow remains expensive; passing accuracy alone
does not establish a performance advantage over Microgen. Records are in
`runs/g3-spinodoid-reference-sensitivity-v1.json`.

## Partial curving at lattice junctions

Job 46 adds a concave single-tet regression. Full projection causes a sampled
validity failure, and the old fallback removes curvature from all three wall
edges. The regression fails as expected, while the other 21 core tests pass.
The replacement halves failing corrections with their periodic copies, up to
eight reductions, and then restores a straight edge if still necessary. The
sampled Jacobian acceptance ratio remains greater than 0.001. Diagnostics now
count partially curved edges separately from fully straightened edges.

Job 47 passes all 22 core and 35 focused Python tests. On the unchanged 16-cell
lattice, it retains partial curvature on 55 edges and halves the independent
maximum distance from 0.008420 to 0.004210. Node and element counts remain
43,685 and 24,480, and sampled Jacobians remain positive. The mesh still fails
the 0.002 geometry gate. Job 48's 20-cell mesh also fails, at 0.008163. Records
are in `runs/g3-curved-fallback-red-v1`, `runs/g3-lattice-curved16-partial-v1`
and `runs/g3-lattice-curved20-partial-v1`.

Job 49 locates the 20-cell maximum at a partially curved junction midpoint.
Job 50 then checks the compiled gradient at equal smooth-union branches. The
documented first-branch minimum derivative and zero derivative for abs(0) do
not cancel correctly in this smooth composite. The maximum component error
against centered differences is 0.81650. Replacing the minimum term algebraically
with `(a+b-abs(a-b))/2` reduces that error to 1.82e-10. Field values agree with
the original Microgen specimen within 1.67e-16 at 4,096 fixed random points.
The generic compiler's documented primitive conventions remain unchanged.
Records are in `runs/g3-lattice-partial20-diagnostic-v1` and
`runs/g3-lattice-gradient-ties-v1.json`.

The symmetric field expression is explicitly verified and recorded as input
preparation. Job 51 reduces the 20-cell distance to 0.002241 at unchanged counts
but still fails admission. Job 52's 16-cell control remains at 0.004210. Job 53's
22-cell seed has fewer final nodes than the 20-cell seed and rejects at 0.004416.
These results do not establish monotone refinement from seed counts. Records
are in `runs/g3-lattice-curved20-symmetric-v1`,
`runs/g3-lattice-curved16-symmetric-v1` and
`runs/g3-lattice-curved22-symmetric-v1`.

Job 54 applies existing adaptation to the 20-cell symmetric expression. It
fails before returning a mesh, with 525 elements below the unchanged 0.1 quality
gate and minimum quality 0.04051. No FE solve runs. Records are in
`runs/g3-lattice-curved-adapt20-symmetric-v1`.

Job 55 adds a regression with different tetrahedron thicknesses on opposite
periodic faces. Partial-curving reductions preserve identical paired midpoint
coordinates. All 23 core tests pass. The nonadaptive 24-cell lattice passes
independent admission with 69,025 nodes and 44,976 Tet10 elements. Its sampled
two-sided distance is 0.0005435, relative volume error is -0.03105%, and minimum
corner quality is 0.3410. All sampled Jacobians are positive, periodic pairs are
exact, and caps match. This mesh needs no partial-curving reductions. The
gradient representation fix, rather than the fallback, removes this candidate's
junction obstruction. All six loads converge below 1e-12. Generation takes
1.92 seconds, FEA 234.73 seconds, and the full workflow 239.60 seconds.
Reference sensitivity remains pending. Records are in
`runs/g3-lattice-curved24-symmetric-v1`.

The next candidate uses 32 seed cells with the same symmetric expression,
void clamp, no adaptation, and Tet10 policy. It must pass independent geometry
and periodic checks before solving, within the explicit 150,000-node and
1,200-second reference stage. Seed density alone does not establish convergence.
