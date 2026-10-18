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
CSR residual must pass after undoing the permutation. This is still one load,
not a full mixed-load reference. Records are in
`runs/g3-lattice-cubic-block-amg-rigid-v1`.
