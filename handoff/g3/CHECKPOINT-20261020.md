# G3 continuation checkpoint, 20 October 2026

The four-fixture matched-accuracy goal is active and unachieved. The source
and documentation changes through `b186d39` are pushed. This checkpoint adds
the controlled quartic edge witness and its resolved weak signal.
No delegates, unrelated jobs, timers, or remote hosts are involved.

## Current execution state

Jobs 322 through 328 all completed. No Meshers successor is
queued at this checkpoint. Recheck the authoritative shared queue before
running anything. Do not restart cubic23 or treat its old observing terminal
session 71321 as a current worker.

The completed cubic23 job uses `cubic_reference.py`, source
`runs/g3-lattice-cubic-parent23-fields-v1`, prepared mesh
`runs/g3-lattice-cubic23-assembly-memory-v1`, and the validated elastic-block
AMG preset. No saved load prefix is replayed. All paths beginning with `runs`
are relative to `experiments/adaptive`.

The source parent is independently admitted, has six fresh solved loads, and
passes field audit 320. Actual cubic assembly 321 has 218,508 nodes and
43,844 elements and peaks at 5.014 GiB. The worker retains 8 GiB, the existing
250,000-node and 250,000-tet stages, a 240-second first-load deadline, and a
1,800-second total worker timeout. One numerical thread is used.

Its sampled geometry comes from the admitted quadratic parent through an
exact cubic embedding, not a projected cubic boundary. The next selected
question is now one coarse-corner coupling check on the same quartic curve-edge
witness. The edge-only controls and lattice measurements are complete; the
joint coarse-corner check has not been implemented or run.

## What changed the next action

The eight-round curve pilot improves size realization and single-run cost but
still fails the fixed model envelope at 1.439306%. Do not tune its strength
or round count again without new evidence.

A separate projected cubic boundary probe passes geometry and improves
sampled distance fourfold. With saved displacement coefficients fixed, energy
changes by only 0.0105922%, below the predeclared 0.05% signal screen. It is
not an equilibrium solution or a bound on total geometry error. No projected
domain solve was selected. Records are in
`runs/g3-lattice-cubic-geometry-sensitivity-v1`.

The actual cubic25 assembly fits, but first solve 314 fails before completing
a load because full symmetry subtraction requests a 1.47 GiB temporary beyond
the 8 GiB address-space ceiling. The AMG helper now checks the complete
relative Frobenius symmetry norm in row blocks. Its 1e-12 threshold and
hash-guarded reuse are unchanged. Ten focused tests and twelve integrated
control solves pass. The failed run stays preserved.

Fresh retry 316 completes all six cubic25 loads without replay. Worker time
is 1,071.18 seconds, peak RSS is 6.753 GiB, and own-rule error is 1.03e-11
relative. Periodic jumps are below 2.23e-16. Positive 125-point energy shift
is 2.6354e-5, below the unchanged 1e-4 sensitivity gate. Fresh-process audit
317 passes, including 64-to-125-point change 6.8313e-6. All fields remain local.

## Reference uncertainty is still real

Cubic24 and cubic25 differ by 0.248177%. Cubic25 is stiffer in every
generalized mode. Relative to cubic25, the finest native adaptive tensor
differs by 0.651535%, the eight-round curve tensor by 0.483742%, and the
strongest actual Microgen-derived Tet10 tensor by 0.917416%.

Do not use that single favorable reference to claim a timing victory. All
24 earlier models remain retained. Cubic25 itself has a model-normalized
envelope of 0.114924% to 1.920819%, so it does not qualify as a target-passing
replacement. Neither two nonnested levels nor a higher polynomial degree
proves continuum convergence. The third level is intended to test that
uncertainty, not guarantee that an extrapolation will pass.

Read `CUBIC_RESULTS.md`, `PLAN.md`, `RESULTS.md`, `PROJECT_TRACKING.md`, and
the append-only `.audit/g3-microgen.tsv`. Detailed new records are in:

- `runs/g3-lattice-cubic25-reference-blocked-symmetry-v2`.
- `runs/g3-lattice-cubic25-assessment-v1`.
- `runs/g3-lattice-cubic25-reference-envelope-v1`.
- `runs/g3-lattice-cubic-parent23-preflight-v1`.
- `runs/g3-lattice-cubic-parent23-fields-v1`.
- `runs/g3-lattice-cubic-parent23-field-audit-v1`.
- `runs/g3-lattice-cubic23-assembly-memory-v1`.

## Completed third-level result

Cubic23 completes all six fresh loads in a 1,018.27-second worker, peaking
at 5.976 GiB. Its own-rule error is 1.29e-11 and periodic jumps are below
2.23e-16. Positive 125-point energy shift is 7.6754e-6, and the minimum
parent-minus-cubic eigenvalue is +0.000137181. Fresh-process audit 323 passes,
including 64-to-125-point change 1.1212e-6. These three eight-model reports
remain retained:

- `runs/g3-lattice-curve-grading-reference-comparison-v1.json`.
- `runs/g3-lattice-reference-sensitivity-v2.json`.
- `runs/g3-lattice-microgen-reference-comparison-v1.json`.

Before reading the full cubic23 tensor, `cubic_refinement.py` defines a
common leading-matrix model, order range 0.25 to 6, at most 10% relative
matrix-fit residual, and signal at least ten times observed quadrature shifts.
Job 324 passes seven synthetic controls, verifies all field hashes and parent
provenance, and checks exact equality of the retained 24 model tensors.

Both length proxies reject the model. Correction-direction cosine is
-0.959707 and the best proportionality coefficient is negative. Matrix-fit
residual is 28.1%. Cubic24 is softer than both neighboring levels. Do not fit
only the closer odd pair, drop the middle level, promote an observed order,
or select a favorable reference for a timing claim. This rejects the stated
local model, not every possible convergence mechanism.

Records are in `runs/g3-lattice-cubic23-reference-v1`,
`runs/g3-lattice-cubic23-assessment-v1`, and
`runs/g3-lattice-cubic-three-level-screen-v1`.

## Completed edge-only diagnostic

Read the measured protocol and limits in `CUBIC_RESULTS.md`. The witness uses
the saved cubic24 fields in `runs/g3-lattice-cubic-reference-cache-release-v1`
and the prepared mesh in `runs/g3-lattice-cubic-assembly-memory-v1`. Freeze
their geometry, topology, and source fields.

The scalar mode on each periodic wall/box-curve edge orbit is
`16*lambda_i*lambda_j*(lambda_i-lambda_j)**2`. Give each orbit one shared
three-component displacement coefficient, supported on the complete incident
element star. Full overlap terms are included. Three unit tests and analytic
periodic controls pass, including nonzero quartic face traces, explicit
off-orbit coupling, and a nonzero pin-gauge test. No regularization is added.

The actual lattice has 102 edge orbits, 306 mode DOFs, and 332 support tets.
Job 327 passes energy and trace identities but fails the 64/125-point
precision gate at 6.0323% change. Preserve that rejected run. One higher-rule
125/216-point check, job 328, passes at 0.017934% change. Its worst mixed-load
trial-energy reduction is only 6.4791e-7 relative, or 0.648 ppm, versus the
500 ppm useful-signal screen. The direct energy identity agrees within
2.04e-13; continuity, periodic quartic traces, and the actual pin offset are
zero. Worker time is 26.57 seconds and peak RSS is 0.3404 GiB.

This resolved small witness does not justify a native sizing feature. It is
not an equilibrium tensor or a bound on all FE error: the original coarse
coefficients were fixed. Inputs, geometry, and all 24 models remain unchanged.
Records are in `runs/g3-quartic-edge-controls-v1`,
`runs/g3-quartic-edge-controls-v2`,
`runs/g3-lattice-quartic-edge-witness-v1`, and
`runs/g3-lattice-quartic-edge-witness-higher-v2`.

## Next single coarse-coupling check

Before closing this polynomial route, allow only the existing cubic corner
basis functions at periodic curve vertices to relax with the same quartic
edge modes. Include complete vertex stars and all cross terms. The earlier
Tet10-to-cubic audit covered 1,408 curve-vertex-incident tets versus 332 for
the edge subset; this motivates a test, not a claim about current cubic error.
Keep geometry and source fields fixed. No new mode degree or grid is selected.

Report corner-only, edge-only, and joint trial gains. Use the incremental
quartic gain after coarse relaxation for the signal, not the joint total.
Verify it by the Schur complement, and check nonnegative gain increments
under space enlargement. Reproduce the original edge-only Gram before the
comparison. Any nonzero source-pin mode value needs an explicit uniform
translation gauge, which changes no strains or periodic jumps.

Retain the validated 125/216-point pair, 1% precision gate, 0.05% signal
screen, 8 GiB ceiling, at most 512 scalar modes and 5,000 support tets.
Run controlled coupled-space checks before a lattice measurement. If the
increment remains small, stop this polynomial route rather than escalating
degree or strength. Trial energy is still not globally equilibrated or a new
homogenized reference. No such joint lattice check has run yet.

## Remaining objective and runtime

Schwarz still has the only narrow repeated matched-target workflow win.
Gyroid and spinodoid still lack admitted actual MMG comparators. Lattice
still lacks a justified target-passing comparison across the retained models.
Local stress and continuum certification remain outside the demonstrated scope.

Use the pinned `.venv`, Rust 1.95.0 if rebuilding, and
`LD_LIBRARY_PATH=/home/kmarchais/dev/meshers/.local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu`.
The native binary SHA256 remains
`908e326211792142f83f854ad2696cd879d927bd81deacf38eb6a0cb63bc4dc1`.
Large NPZ archives and native binaries are ignored local artifacts, not backed
up by the committed text. Heavy work remains in the shared `gpu` queue through
`qrun --check-start` and standard `Job` progress. Before GitHub operations,
switch to and verify `kmarchais`. Do not touch unrelated queue entries.
