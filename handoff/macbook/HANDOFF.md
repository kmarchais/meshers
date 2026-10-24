# Meshers continuation on MacBook Pro, 24 October 2026

Kevin explicitly requested moving primary work from G3 to this MacBook Pro
and starting a T3 Code thread with `gpt-6.1-sol`, reasoning effort `medium`.
The new thread owns continuation. G3's last experiment, job 476, is terminal
success. Source thread `01a14a4b-66a2-7331-aa95-c9993fe588b2` belongs to the
G3 environment. The earlier OMEN thread is historical.

## Objective and authority

Improve Meshers and demonstrate that it beats actual Microgen+MMG at matched
accuracy on the fixed Gyroid sheet, Schwarz P solid, spinodoid, and smooth
lattice. All four remain in scope and the objective is unachieved. Continue
autonomously, one bounded experiment at a time, reading results before a
successor. No delegates, cross-model reviews, or unrelated jobs/timers touched.
User-space installs, useful commits, and pushes are authorized. Ask only for
a real privilege/credential blocker. Before GitHub operations, switch gh to
`kmarchais` and verify its active account.

Read destination global instructions, computers inventory, queue/progress
and decision-log skills, then this file, `PROJECT_TRACKING.md`, and
`experiments/adaptive/PLAN.md`, `RESULTS.md`, and `CUBIC_RESULTS.md`.
Older active/queued descriptions in those documents are historical.

## Checkouts and transferred evidence

- Main `/Users/kmarchais/dev/meshers`, branch `handoff/g3-microgen-20261017`.
  Its previous release branch is preserved. This handoff is committed with
  the latest experimental code and selected text evidence.
- Sibling `/Users/kmarchais/dev/periodic-microgen`, exact commit
  `2b1ea80b36d4ba03842a7d9087638abf3868fe3e`.
- Sibling `/Users/kmarchais/dev/microgen-benchmark`, exact commit
  `602dc900cd028d205b36e29dcdb5b4a9cf4f581a`.

Support sources are required by experimental imports. Git bundles and the
SHA256 transfer manifest live outside Git at
`/Users/kmarchais/dev/meshers-handoff-20261024`. The manifest records which
ignored meshes, complete six-load fields, fixture samples, and historical
inputs actually reached this machine. Linux binaries/venvs, captured sparse
systems, and redundant load checkpoints are excluded. Some original large
OMEN artifacts never reached G3 and remain absent.

Keep frozen manifests unchanged. Resolve their old `/home/kmarchais/dev`
prefix to `/Users/kmarchais/dev` when verifying historical files; retain
hashes. New runs record real Mac paths and versions. Do not fabricate a Linux
home tree or edit old manifests to make checks pass.

## Latest completed results

All named run directories are under `experiments/adaptive/runs`.

Job 473, `g3-lattice-adapt26-field-audit-v1`, independently checks job 379's
`g3-lattice-adapt26-tol0028-amg-solve-v1` fields. Reaction error is
1.0611e-11, Hill-Mandel error 4.1723e-12, and periodic jumps zero. Positive
125-point integration shifts energy by 4.0667e-8; 64-to-125 changes by
9.9686e-9. Audit passes in 117.41 seconds at 0.7454 GiB on G3. The unchanged
24-model score remains 1.2686985%, above 1%. Numerical inconsistency and
quadrature do not explain the miss. The audit now uses `worker_part` to keep
standard progress values. No native code changed in this phase.

Job 475, `g3-lattice-nested23-preflight-v1`, starts from solved and admitted
`g3-lattice-cubic-parent23-fields-v1`, 68,095 nodes and 43,844 Tet10 elements.
Uniform eight-child subdivision would exceed 250,000 tets. The new
`nested_refinement.py` bisects complete periodic edge orbits in descending
length, using disjoint incident stars and at most one split per parent.
Each child restricts its exact parent quadratic geometry map and prolongs
the six old displacement fields. Both child corner qualities must stay >=0.1.

Ten queued controls pass. The result has 123,514 nodes and 76,947 tets,
with 33,103 parents bisected and 9,450 accepted edge orbits. Geometry/field
interpolation errors are 3.33e-16, volume change 1.40e-14, sampled distance
0.0006492808, and quality 0.10385108. Curved Jacobians, periodic caps and
oriented facets pass. Worker is 14.82 seconds at 0.5769 GiB. Changed sampling
explains the distance-score difference despite exact domain preservation.

Job 476, `g3-lattice-nested23-level1-v1`, completes six fresh loads on this
saved mesh with validated elastic-block AMG and 15-point curved integration.
Solve is 245.83 seconds, worker 247.63 seconds, peak RSS 3,304,448,000 bytes.
Independent fields and positive quadrature pass. All generalized stiffness
corrections from the parent are negative, largest magnitude 0.499513%.
The full prolonged-parent/fine/correction energy balance agrees within
6.6756e-8 relative. Worst disagreement across all 24 models is 2.4480846%,
failing 1%. This coarser-domain diagnostic is not the latest best native mesh
or a fresh meshing workflow timing. One correction does not prove convergence.
Its `assessment.json` is the decisive result. No new reference is promoted.

## Selected next experiment

Establish the native Mac environment and shared queue, verify the transfer,
and run focused controls. Do not repeat completed G3 solves for context.
Then prepare ONE second-level nested geometry/resource preflight from the
solved level1 mesh and fields. Inspect actual counts and all admission checks
before deciding on a solve. Retain the fixed physical domain and old FE
space so this measures discretization separately from background-grid changes.

Two runner details need attention before that preflight:

1. `nested_solve.py` does not copy fixture sample/specimen files into its
   solve directory. Obtain hash-identical inputs from the original parent or
   prepared directory. Restore missing copies and record hashes, or support
   source resolution prospectively. Do not overwrite frozen outputs.
2. Node/tet caps remain 250,000. The next level may approach them; count
   before assembly. Extra RAM does not silently change admission caps or
   reclassify historical cap failures.

If corrections remain large or the next level cannot fit, record that result
before choosing another reference strategy. Keep all 24 historical models
visible. Their pairwise lower bound is 0.9019313%. A trial matrix passes at
0.9084715%, proving only algebraic feasibility, not physical realizability.

## Four-case status and fixed rules

- Schwarz has repeated roughly 5.2x speedup at 1% and 7.4x at supplementary
  0.1%, against explicitly modified actual Microgen routes. Strong control
  uses manual straight Tet10. G3 seconds are historical, not Mac timings.
- Gyroid Meshers curved9 passes with 11,857 nodes and 0.72921–0.95990%
  model disagreement. Tested actual Microgen routes remain unqualified.
  Boundary-only witnesses block fixed-connectivity interior-only repair;
  tested local flips, complete-star topology, and one Steiner point fail.
- Spinodoid has an admitted Meshers route, 73,951 nodes and
  0.10333–0.18439% disagreement. Tested Microgen caps obstruct admission.
- Lattice remains unresolved. Strongest Microgen-derived Tet10 is
  1.004354%, also failing 1%. Three nonnested cubic levels do not establish
  an asymptotic reference. Recovery/energy and curve-sizing pilots failed
  to close the gap. Independent volume sampling does not justify an offset.

Retain exact specimens/material from tracker section 4, distance <=0.002,
absolute volume error <=1%, corner quality >=0.1, positive sampled curved
Jacobians, periodic nodes/cap triangles, coherent oriented facets, and full
mixed-load stiffness targets/models. Unsupported routes stay explicit.
Stiffness is not local-stress or continuum certification. Do not modify
the fixed Fourier realization or reorder nonassociative smooth unions.

Product priorities remain actual comparison, Tet10 integration, stress-recovery
sizing, justified volume/parallel relocation, branch-aware creases and box
curves, then coarser extraction. Read measured failures before reviving recipes.

## Mac environment and resource setup

Verified hardware is Apple M4, 10 CPU cores, 24 GiB unified memory, arm64,
and about 223 GiB free disk. No CUDA. These meshing/FE jobs are CPU-compatible.
Rebuild native code and dependencies; do not reuse Linux wheels or venvs.
Use repository locks, Rust 1.95.0, Python 3.12, and build Meshers with
`experimental-surfaces`. Do not export G3's LD_LIBRARY_PATH.

G3 versions are NumPy 2.5.4, SciPy 1.18.1, fedoo 1.0.0b2, Microgen 2.0.0b1,
mmgpy 0.17.1, PyAMG 5.3.0, PyVista 0.49.1, VTK 9.7.1, maturin 1.9.6,
and psutil 7.2.2. Verify arm64 availability; explicitly pin any required
platform replacements. Keep solver/affine/periodic controls before comparisons.

Login-shell PATH resolves uv, Node, gh, qrun and Codex; pueue was not found.
Install/configure pueue and its single-worker `gpu` group before heavy work.
Read destination workbench progress instructions. Every heavy build/test/
solve uses qrun, standard Job fraction/message/ETA/values, `--check-start`,
bounded stages and a truthful interruption declaration. No unrelated jobs touched.

Start with one numerical thread and modest build parallelism. Existing
`microgen_workflow.bounded` hardcodes 8 GiB and Linux-tested RLIMIT_AS setup.
Check macOS behavior before using it; its psutil RSS watchdog is separate.
A measured Mac-specific ceiling around 12 GiB can be introduced prospectively,
with recorded manifests and small timeout/memory-failure controls. Do not
inherit OMEN's ceiling or consume all unified memory. Establish fresh Mac
baselines before claiming Mac speedups. Commit useful code/docs/text evidence;
keep large NPZ and build outputs ignored.
