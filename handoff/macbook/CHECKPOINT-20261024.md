# Native MacBook continuation, 24 October 2026

This MacBook session owns primary Meshers continuation on
`handoff/g3-microgen-20261017`, starting at `0b24b564`.
The four-case objective remains active and unachieved. G3 job 476 is complete
and has not been restarted. No delegates or cross-model reviews were used.

## Native environment

The M4 has 10 CPU cores and 24 GiB unified RAM. The native environment is
`/Users/kmarchais/dev/meshers/.local-work/macbook-runtime/venv`, using Python
3.12.9. Every requested G3 numerical package version installed from arm64
wheels. No platform version replacement was needed. Each new run records
its installed requirements. Native facts, binary hash, and the resource
probe are in [environment-20261024.json](environment-20261024.json).

Meshers was rebuilt with Rust 1.95.0, the repository lockfile, two build
workers, and `experimental-surfaces`. Mac job 0 passed in 56.26 seconds at
1.307 GiB peak process-tree RSS. Job 1 passed 54 focused tests. Job 2 passed
twelve fresh AMG loads on a cubic affine cube and the saved curved Schwarz
control. Tensor disagreements were below 3.45e-12 and 4.42e-11 respectively.
Periodic, reaction, own-rule, and higher-rule checks passed. These are
environment controls, not workflow comparisons.

The first Mac numerical runs did not explicitly control Apple Accelerate,
which is NumPy's BLAS backend here. A fresh process reports native BLAS mode
0, allowing automatic threading. Their numerical results and audits remain
valid, but numerical timings from jobs 1–8 are excluded from verified
one-thread comparison baselines. Frozen values and original manifests remain
unchanged. The correction is recorded in
[threading-clarification-20261024.json](threading-clarification-20261024.json).

The shared workflow now sets `VECLIB_MAXIMUM_THREADS=1` and calls the native
`BLASSetThreading(1)` API on its calling thread. `BLASGetThreading()` verifies
mode 1 before numerical work. Jobs 9–10 pass 55 focused controls and twelve
fresh cube/Schwarz loads with the corrected policy, without changing tensors.
New manifests and worker records include the actual policy. The Apple API
is thread local; these workers perform numerical work on their main thread.

Pueue 4.0.4 runs through its Homebrew user service. Its shared `gpu` group
has one worker. The restored workbench checkout and external environment
support `qrun --check-start`. A user launch agent serves the local Queue tab
at `http://127.0.0.1:7811/#queue`. Standard Job progress and results were
verified through actual queued stages. No old timers or unrelated jobs were
changed.

macOS rejects the former 8 GiB `RLIMIT_AS` call with
`ValueError: current limit exceeds maximum limit`. On macOS,
`microgen_workflow.bounded` now uses its existing process-tree RSS watchdog
without the Linux pre-exec address-space limit. Worker records declare this
policy and the 0.25-second polling interval. Small success, timeout, and
memory-stop controls pass. Summed RSS is conservative when processes share
pages; polling does not prevent every transient allocation. The cap remains
8 GiB. No 12 GiB protocol or mesh-cap change has been needed.

The transfer was rehashed independently: all 3,509 files, totaling
1,650,689,800 bytes, match. The new receipt is
`/Users/kmarchais/dev/meshers-handoff-20261024/macbook-primary-verification.json`.
Both support repository commits match the handoff. Frozen Linux manifests
and their hashes remain unchanged.

## Second nested level

Mac job 3 completes the selected geometry and field preflight from the
finished G3 level1 mesh. Its fourteen controls pass. The output is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-nested23-level2-preflight-v1`.

The preflight has 218,135 nodes and 141,144 Tet10 elements. It bisects
64,197 parent elements in 12,861 complete periodic edge orbits. Neither
fixed 250,000 cap is reached. The curved physical domain and the six old
displacement fields are retained, with maximum interpolation errors
2.78e-16 and 3.33e-16. Relative volume change is -2.52e-14.

Sampled distance is 0.0006492808, absolute relative volume error is
0.0310173%, minimum corner quality is 0.10010709, and minimum sampled
curved Jacobian is 4.9374e-7. Periodic nodes, cap triangles, and oriented
facets pass. Worker time is 87.72 seconds at 1.047 GiB peak RSS.

The old level1 solve lacks fixture files. The new `--fixture-source` option
reads their hash-checked copies from the original parent23 directory.
Only the new output receives copies. A regression rejects changed samples.
Future nested solves also copy their prepared fixed samples and record
their hashes, so the omission does not recur.

After reading the preflight, job 4 measures the actual first reduced system
without solving or saving large matrix archives. It has 646,029 free DOFs,
51,030,214 stored entries, and 614,946,688 CSR bytes. Complete three-component
node blocks, finite entries, and positive diagonals pass. Worker time is
9.61 seconds at 4.464 GiB peak RSS. Assembly fitting does not guarantee
hierarchy or six-load solve memory, so the same watchdog remains active.

## Selected solve and interpretation

After reading assembly, job 5 completes one fresh level2 six-load solve with
validated elastic-block AMG, 15-point curved Tet10 integration, and explicit
completed-cache release. The solve has a 1,200-second timeout. Independent
field and positive-quadrature audits and Galerkin assessment each have a
300-second timeout. Its immutable output directory is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-nested23-level2-v1`.

All six residuals pass. Solve worker time is 191.45 seconds at 4.508 GiB
peak RSS. Independent fields and positive quadrature pass. Every generalized
parent-to-level2 correction is negative, with maximum magnitude 0.328272%,
versus 0.499513% for the preceding level. The full Galerkin balance agrees
within 3.78355e-8. Worst retained-model disagreement falls to 2.135117%,
still above 1%. No reference is promoted. Nonuniform bisection and two
corrections do not establish an asymptotic order or a continuum bound.
This saved-mesh diagnostic supplies no live workflow timing or cross-host
speed ratio.

After reading that result, job 6 tests one nested preflight on the best saved
native lattice26 candidate, whose original worst model disagreement is
1.2686985%. Its worst-model generalized error is positive, so measured
same-domain FE softening is a relevant hypothesis. The edge-selection rule
is unchanged and uses lengths, periodicity, quality, and budgets. It does
not use oracle fields to tune a size function. The new preflight is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-adapt26-nested-preflight-v1`.
Job 6 passes at exactly 250,000 nodes and 160,291 tets. It accepts 7,971 edge
orbits, bisects 31,991 parents, and stops further orbits at the unchanged
node cap. Map and field errors are 3.33e-16. Volume change is -1.11e-15.
Distance 0.0003236401, quality 0.14932885, positive sampled Jacobians,
periodic caps, and oriented facets pass. Worker time is 11.27 seconds at
1.140 GiB peak RSS. The source's three independent field-audit input hashes
were rechecked.

After reading that result, job 7 measures actual reduced assembly at
726,645 DOFs, 57,890,122 stored entries, and 697,588,048 CSR bytes. It passes
in 10.85 seconds at 5.270 GiB peak RSS. Job 8 then completes six fresh loads and
the same independent audits, under unchanged timeouts and the 8 GiB watchdog.
Its directory is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-adapt26-nested-v1`.
All fields, positive-rule checks, and Galerkin assessment pass. Peak solve
RSS is 5.038 GiB. The largest FE correction is only -0.0122671%, with full
Galerkin balance error 6.07082e-10. Worst disagreement across all 24 models
is 1.2565145%, still failing 1%. Stop this longest-edge recipe. Spending the
remaining node budget there does not establish a useful candidate improvement.

After reading the negative result and passing corrected thread controls,
job 11 preflights one change to edge ordering. It prioritizes actual
wall/cap intersection edge orbits, then uses the existing length order.
Earlier same-domain correction localization places 75.02% of the correction
trace near curve vertices, whereas the longest-edge test resolves little.
There is no sizing-strength parameter or change to physical geometry.
Every periodic, quality, map-restriction, and field-nesting guard remains.
The default length-only path is unchanged. A short-edge control checks the
requested bisection and exact parent map; nonexistent preferred edges reject.

The new preflight directory is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-adapt26-curve-nested-preflight-v1`.
Job 11 passes sixteen controls. It accepts 544 of 549 preferred edge orbits
and fits 249,999 nodes and 158,967 tets. Map and old-field errors are below
3.34e-16. Relative volume change is -3.77e-15. Distance is unchanged at
0.0003236401, corner quality is 0.10271976, sampled curved Jacobians are
positive, and periodic caps and oriented facets pass. Worker time is 13.37
seconds at 1.167 GiB RSS, with the corrected numerical thread policy.

After reading that preflight, job 12 measures actual reduced assembly at
721,914 DOFs, 57,483,976 stored entries, and 692,695,372 CSR bytes. It passes
in 9.97 seconds at 6.230 GiB peak RSS. Job 13 then completes one fresh six-load
solve, independent field and positive-rule audits, and full Galerkin and
24-model assessment. The 8 GiB cap and timeouts remain unchanged. The directory is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-adapt26-curve-nested-v1`.
All six loads, independent fields, positive quadrature, and full Galerkin
assessment pass. Solve worker time is 216.33 seconds at 5.472 GiB RSS. The
largest additional FE correction is only -0.0132094%, and full energy balance
agrees within 4.82986e-10. Worst disagreement across all 24 models remains
1.2556246%, failing 1%. Prioritizing the curve edges gives little additional
gain over the 1.2565145% length-only result. Stop these one-stage bisection
recipes. No public native feature, stiffness reference, or workflow win is
promoted from either negative result.

## Verification and next work

[verify_checkpoint.py](verify_checkpoint.py) independently rechecks all
3,509 transferred files and 1,650,689,800 bytes, 150 Mac input hashes, and
72 full generalized-eigenvalue comparisons for the three completed nested
studies. All pass. Twenty-six changed source references resolve to exact
frozen driver or Git content, without relaxing any mesh or field hash.
The receipt is [verification-20261024.json](verification-20261024.json).

Mac jobs 0–13 are complete. No successor is queued at this checkpoint.
The active goal remains the full four-case objective. Return to the missing
admitted actual Microgen comparator coverage, beginning with the preserved
Gyroid and fixed-spinodoid failures in `MICROGEN_RESULTS.md` and the current
coherent seed/cap evidence. Select one bounded preflight from that evidence
before any MMG or solve. Do not repeat rejected interior moves, local cavities,
seed grids, sizing strengths, or polynomial degrees without a new reason.
Lattice reference uncertainty remains open; these negative refinements do
not justify selecting a favorable tensor or removing any retained model.

Schwarz retains its historical repeated wins against explicitly modified
actual Microgen routes. Gyroid and spinodoid still lack admitted actual
comparators in tested routes. Lattice reference convergence remains unresolved.
Fresh Mac comparisons still need same-host, admitted, matched-target workflows.
The fixed specimens, distance 0.002, volume 1%, quality 0.1, sampled positive
Jacobians, periodic gates, coherent facets, and node/tet caps are unchanged.

The append-only decision trail is
[macbook-continuation.tsv](../../.audit/macbook-continuation.tsv).
