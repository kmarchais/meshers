# Microgen comparison on G3

Kevin authorized this transfer on 17 October 2026. Continue on G3 with GPT 6.1 Sol,
medium reasoning. This is the replacement working session. Do not delegate reviews
or create duplicate timers. Install missing dependencies within the user account;
ask Kevin only for a real missing credential or privilege.

## Source and evidence

The main branch is `handoff/g3-microgen-20261017`, based on adaptive commit
`80cf40ce`. Native meshing code last changed in `8d97a785`; the timing job froze
`9dad555f`. Read `experiments/adaptive/PLAN.md` and `RESULTS.md`, then inspect
the actual scripts and recorded data. Some older Limits paragraphs in RESULTS.md
predate later experiments. Do not treat them as current conclusions.

Two support branches preserve the earlier work without merging divergent code:

- `handoff/g3-periodic-20261017`, commit `2b1ea80b`, includes the pending fix that
  refuses to estimate a live recipe time from a geometry replay.
- `handoff/g3-benchmark-20261017`, commit `602dc900`, contains the benchmark package.

Clone the main branch at `/home/kmarchais/dev/meshers`. Check out the support
branches as sibling worktrees `/home/kmarchais/dev/periodic-microgen` and
`/home/kmarchais/dev/microgen-benchmark`. Existing adaptive imports expect this
layout. `data/` contains previously ignored tensor records and the two archived
Gyroid result records used by `fea.reference()`. These result records are useful
for historical tensor comparisons, but do not include mesh/displacement archives.
Those large archives remain on OMEN. Retrieve only specific required files later.
Never claim saved-field revalidation from these small result records alone.

Copy the three JSONL files from this directory's data folder back into
`experiments/adaptive` only if absent. Copy gyroid-reference64.json to the support
worktree at `experiments/periodic-fea/runs/reference-and-seed-v1/02-periodic-gyroid-meshers-pass1-64/result.json`
and gyroid-reference80.json to
`experiments/periodic-fea/runs/reference-80-v1/00-periodic-gyroid-meshers-pass1-80/result.json`.
Rebuild the extension on Linux. Do not transfer Windows wheels or virtualenvs.

## Next work, in Kevin's requested order

1. Run the real Microgen workflow end to end. Inspect Microgen's installed source
   and use its own VTK and MMG calls with their defaults. Cases are a Gyroid sheet,
   Schwarz P solid, spinodoid, and a lattice with smooth unions. Compare against
   Meshers linear and quadratic meshes on identical specimens and material laws.
   Measure counts, independently sampled two-sided maximum boundary distance,
   volume error, element quality, generation time, total workflow time, and full
   mixed-load homogenized stiffness error against justified finer references.
   Pin versions and record exact default arguments. Report unsupported paths and
   failures; do not substitute a different pipeline and call it Microgen's default.
   Separate one-sided material-domain changes from discretization error.
   Prepare one queued experiment, with bounded pilot stages and early resource
   checks. Inspect its actual results before selecting another experiment.
2. Get Tet10 into Microgen, through from_pyvista support or a direct fedoo path.
   Check node ordering, periodic mid-edge pairing and curved Jacobian validity.
3. For linear meshes, use a coarse FEA stress-recovery error estimate to size the
   ligament necks. Curvature sizing alone has not improved Gyroid per-node error.
4. Correct adaptive volume bias and separately parallelize vertex relocation.
5. Capture min/max branch intersection creases; investigate spinodoid box curves.
6. Reduce initial extraction density to improve single-thread fast sphere speed.

Step 1 decides which later changes matter most. Do not queue steps 2 through 6 as
speculative batches or make all changes before measuring Microgen.

## What is established and what is not

Adaptive component timing was job 1676, warmup plus five interleaved repetitions,
at one and 24 threads. Fast surfaces used fewer triangles at matched sampled
distance and beat Flying Edges for Gyroid, Schwarz P and ellipsoid. Sphere was
slower on one thread. Tet adaptation did not scale with threads. The volume MMG
comparison remeshed a Meshers seed; it was not the full Microgen workflow.

Curved Tet10 stiffness results are promising, but use extrapolated references.
The large per-node improvement is not local-stress certification. Quadratic
validity currently samples Jacobians instead of proving positivity everywhere.
Do not silently admit curved meshes based only on linear corner quality.

Sharp crease lattice and coarse spinodoid still fail geometry checks. Adaptive
volume bias and box-edge wedges remain. Microgen's existing conversion rejects
Tet10. All these are real limitations, not reasons to relax acceptance gates.

The earlier periodic study has separate matched-node stiffness AND local stress
comparisons against multiple references. Read its ACCURACY.md and ADOPTION.md.
Preserve COMMON geometry 0.002, relative volume 0.01 and quality 0.1 wherever that
study's admission rules apply. Do not rewrite the Gyroid boundary seed's failed
strict 0.0002 reverse-distance flag. Saved audits must be rechecked with
verify_run(write_report=False), and their report hashes must stay unchanged.
Never infer local accuracy from stiffness alone or claim continuum certification.

## G3 execution

G3 has 14 GiB RAM, an i7-10750H and a 6 GiB mobile GPU. At transfer it had about
11 GiB available RAM and an idle shared gpu queue. Recheck before running.
Use a worker memory ceiling around 8 GiB, conservative thread counts, explicit
timeouts and staged mesh limits. Do not reuse the OMEN 16 GiB worker ceiling.
Keep timing isolated from builds and other jobs. Compare workflows on this host;
do not compare G3 seconds directly with old OMEN seconds.

Read the local computers inventory, AGENTS instructions, gpu-queue and
show-me-your-work skills. Every heavy build, solver, meshing batch and busy-lane
test goes through qrun in the gpu group. New scripts use the standard Job progress
API with fraction, message, ETA and values; submit with qrun --check-start and fix
its audit. Queue one job, inspect its result, then choose the next. No delegates.
Log decisions truthfully and commit useful code/docs on this branch. Ignored
experiment files require git add -f. Preserve existing artifacts and frozen runs.

Before GitHub operations switch gh to kmarchais and verify auth status. Kevin
authorized pushing these handoff branches. No new PR was created for the move.
The old periodic PRs are #2 and #5; link them only when actually working on them.

Global user prose rules include unslop and complete absolute paths for shared
artifacts. Read the destination's global instructions before continuing.
