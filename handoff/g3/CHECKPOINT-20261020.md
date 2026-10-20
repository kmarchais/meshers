# G3 continuation checkpoint, 20 October 2026

The four-fixture matched-accuracy goal is active and unachieved. The source
and documentation changes through `4587cd6` are pushed. This checkpoint adds
the completed finer-reference evidence and records one running successor.
No delegates, unrelated jobs, timers, or remote hosts are involved.

## Current job to inspect first

Pueue job 322 runs the fresh cubic23 reference. Its run directory is
`experiments/adaptive/runs/g3-lattice-cubic23-reference-v1`. The observing
terminal session in this thread is 71321. Check pueue and that existing handle
before doing anything else. A missing observation or stale checkpoint is not
proof that the worker stopped. Do not launch a duplicate or queue a successor
before reading the actual result.

The job uses `cubic_reference.py`, source
`runs/g3-lattice-cubic-parent23-fields-v1`, prepared mesh
`runs/g3-lattice-cubic23-assembly-memory-v1`, and the validated elastic-block
AMG preset. No saved load prefix is replayed. All paths beginning with `runs`
are relative to `experiments/adaptive`.

The source parent is independently admitted, has six fresh solved loads, and
passes field audit 320. Actual cubic assembly 321 has 218,508 nodes and
43,844 elements and peaks at 5.014 GiB. The worker retains 8 GiB, the existing
250,000-node and 250,000-tet stages, a 240-second first-load deadline, and a
1,800-second total worker timeout. One numerical thread is used.

Do not edit hashed reference, AMG, integration, or workflow sources while
this job uses them. Documentation and independent analysis can wait until
its outcome is read. Its sampled geometry comes from the admitted quadratic
parent through an exact cubic embedding, not a projected cubic boundary.

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

## After job 322

If all six loads and the existing own-rule and quadrature gates pass, run one
fresh-process `cubic_assessment.py` audit in a new directory. Use the cubic23
prepared mesh and retain exactly these three eight-model reports:

- `runs/g3-lattice-curve-grading-reference-comparison-v1.json`.
- `runs/g3-lattice-reference-sensitivity-v2.json`.
- `runs/g3-lattice-microgen-reference-comparison-v1.json`.

Compare cubic24, cubic25, the finest adaptive native candidate, the eight-round
curve candidate, and the strongest Microgen-derived Tet10 control. Confirm all
input hashes and all 24 models. Keep candidate-to-model normalization separate
from the reverse model-to-cubic diagnostic.

Then assess the full 23/24/25 sequence without assuming monotonicity or fitting
only a convenient pair. Geometry and displacement spaces are nonnested.
If the sequence is inconsistent, preserve that negative conclusion rather
than declare a reference or tune a native mesh to the model envelope. A failed
job or first-load deadline needs inspection, not an automatic retry.

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
