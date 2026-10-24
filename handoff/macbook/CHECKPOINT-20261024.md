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

After reading assembly, job 5 starts one fresh level2 six-load solve with
validated elastic-block AMG, 15-point curved Tet10 integration, and explicit
completed-cache release. The solve has a 1,200-second timeout. Independent
field and positive-quadrature audits and Galerkin assessment each have a
300-second timeout. Its immutable output directory is
`/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-lattice-nested23-level2-v1`.

Read its authoritative completion, all six residual records, independent
audits, parent correction, and full Galerkin balance before selecting a
successor. Retain all 24 model tensors. Compare the two nested corrections,
but do not infer an asymptotic order from nonuniform bisection alone or
promote a continuum reference from small successive changes. This saved-mesh
diagnostic supplies no live workflow timing or cross-host speed ratio.

Schwarz retains its historical repeated wins against explicitly modified
actual Microgen routes. Gyroid and spinodoid still lack admitted actual
comparators in tested routes. Lattice reference convergence remains unresolved.
Fresh Mac comparisons still need same-host, admitted, matched-target workflows.
The fixed specimens, distance 0.002, volume 1%, quality 0.1, sampled positive
Jacobians, periodic gates, coherent facets, and node/tet caps are unchanged.

The append-only decision trail is
[macbook-continuation.tsv](../../.audit/macbook-continuation.tsv).
