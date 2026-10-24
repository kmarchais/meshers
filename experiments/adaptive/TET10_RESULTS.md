# Public fedoo integration and Tet10 controls

## A modified actual lattice control passes, 24 October 2026

One saved actual Microgen-derived lattice control now passes all 24 retained
models at 0.988497945% worst disagreement. The new native curve-star control
still fails at 1.230455074%. These are fresh Mac solves of saved meshes, not
live regeneration or a matched-accuracy timing comparison. No lattice speed
win, new reference, or continuum/local-stress certificate is claimed.

Jobs 58–59 localize the complete ranked-space increment between the two
previously audited coupled solutions. The own-rule Galerkin identity error is
7.64e-14, and positive 125-to-216-point local-trace change is 8.08e-8.
Curve-vertex incident cells carry 66.8117% of the correction in 2.4994% of
cells. Literal-curve-edge cells carry 23.3435% in 0.4557% of cells. Groups
overlap. The top 1% of cells carry 73.2238%. This is resolved discrete
correction, not an estimate of remaining error.

The first diagnostic finishes all three FE rules, then fails its wall
projection assertion. Its source, arrays and log remain frozen. The second
run reuses those hash-identical completed arrays and reports the projection
failure explicitly. Eleven Newton projections on six faces leave the box,
by up to 6.76e-5. Surface-rule stability passes, but complete trimmed-wall
geometry interpretation does not. Neither the original geometry admission
nor its samples are silently changed. No co-location conclusion follows.

The optional `--prefer-box-curve-stars` preflight mode prioritizes all edges
incident on actual wall/cap curve vertices, including radial edges, before
the existing descending-length fallback. Complete periodic orbits, disjoint
stars, exact parent maps, q0.1, positive sampled Jacobians and 250,000-node/tet
limits remain unchanged. It is mutually exclusive with literal-edge priority;
the default selection is unchanged. A radial-edge selection control and all
18 queued controls pass. No size parameter sweep or oracle field is used.

| Saved-mesh control | Nodes | Tets | Worst error across all 24 models | 1% target |
|---|---:|---:|---:|---|
| Original live native lattice26 | 199,972 | 128,300 | 1.2686985% | Fail |
| Earlier literal-curve-edge nested native | 249,999 | 158,967 | 1.2556246% | Fail |
| New curve-vertex-star nested native | 250,000 | 159,226 | 1.2304551% | Fail |
| Original tuned actual Microgen straight Tet10 | 238,411 | 158,252 | 1.0043541% | Fail |
| New star-refined actual Microgen-derived control | 249,998 | 165,772 | 0.9884979% | Pass |

Jobs 60–61 retain the original native lattice26 domain and parent field to
3.34e-16. They bisect 78.1472% of its curve-vertex incident cells. Distance
stays 0.0003236401; quality is 0.1023041 with positive sampled Jacobians.
All six stiffness changes are negative, largest 0.03868354%. Galerkin balance
is 3.71e-10 and positive-rule energy change is 9.64e-9. Six fresh loads take
212.13 seconds at 5.758 GiB. Stop this allocation recipe without promotion.

The actual control's initial worst signed error is +1.0043541%; its most
negative error across all models is -0.8937705%. That motivates one bounded
refinement test, not rounding the original failure into admission. Jobs
62–63 retain its affine domain exactly. Distance stays 0.0016215853, volume
change is zero, quality is 0.1015540, and sampled Jacobians are positive.
They bisect 73.6880% of curve-vertex incident cells. Six fresh loads take
232.18 seconds at 6.198 GiB. Galerkin balance is 3.86e-14, the largest
stiffness decrease is 0.09635235%, and positive-rule energy change is
1.64e-14. Full field and periodic audits pass. The worst final error is
negative against sensitivity model 7. More FE refinement is not automatically
helpful to this empirical target.

The passing control retains actual tuned VTK51/two-pass MMG provenance,
manual straight Tet10 elevation/direct fedoo, and an explicitly added
periodic star-refinement step. It is not the default Microgen pipeline or
MMG on a Meshers seed. All 24 original models and original failures remain.
The full parent23 reference investigation remains unresolved at 1.7413203%.

Next prepare fresh Mac actual VTK51/two-pass MMG generation with this declared
added refinement and charge every stage. Read geometry and counts before a
new solve; the frozen G3 seed does not prove fresh Mac MMG reproduction.
The native control still needs improvement before a matched-target pair.
Schwarz retains the five repeated Mac workflow wins. Gyroid and spinodoid
still lack admitted actual comparators in tested routes. All four cases and
the unchanged objective remain active.

Records are `runs/macbook-lattice-hp-ranked-localization-v1`,
`runs/macbook-lattice-hp-ranked-localization-v2`,
`runs/macbook-lattice-adapt26-curve-star-preflight-v1`,
`runs/macbook-lattice-adapt26-curve-star-v1`,
`runs/macbook-lattice-actual-microgen-curve-star-preflight-v1` and
`runs/macbook-lattice-actual-microgen-curve-star-v1`.
Jobs through 63 are terminal. Earlier selected-next descriptions are historical.

## Verified repeated Mac Schwarz comparison, 24 October 2026

Job 33 completes one excluded warmup and five fresh interleaved pairs, with
alternating arm order. Every pair passes all eight declared full mixed-load
Schwarz models at the supplementary 0.1% target. Native errors span
0.031469% to 0.093778%; Microgen errors span 0.053480% to 0.062162%.
All geometry, quality, sampled Jacobian, coherent-facet and periodic gates pass.

| Charged cost | Modified actual Microgen median | Meshers median | Five paired ratios |
|---|---:|---:|---:|
| Science-audited workflow | 131.85894 s | 19.14923 s | 6.7632 to 6.9203 |
| Generation, repair/lift and solve | 119.40758 s | 17.48343 s | 6.6949 to 6.8605 |

Workflow ranges are 130.16436 to 132.30789 seconds for Microgen and 19.11868
to 19.54234 seconds for Meshers. Meshers is faster in every measured pair
under both cost definitions. No pair is discarded, and warmup is excluded.

The comparator is actual VTK30/MMG with explicit one-pass overrides, a
charged fixed-boundary SciPy interior repair, and manual straight Tet10/direct
fedoo. It is not Microgen's default or public Tet10 conversion. Meshers uses
fresh nonadaptive curved16. Both use the same block-AMG preset and one
numerical thread. Stage limits stay at 250,000 nodes/tets and 8 GiB; the
whole study takes 980.56 seconds within its one-hour deadline. No previous
mesh, solved prefix or preconditioner is replayed.

All six pairs pass 42 controls. Five verifier unit controls also pass.
Job 34 stops before geometry replay because the verifier omitted support-path
initialization; its exact failed driver and log remain preserved. Job 35
corrects initialization and independently verifies 262 hash records, all
96 full-matrix model comparisons, every timing statistic, and all twelve
saved mesh/field sets. Geometry and periodic field jumps are recomputed
from actual arrays. Minimum quality is 0.100001001 for Microgen and 0.4499953
for Meshers. Maximum reaction/energy error is 7.60e-12; periodic jumps are
at most 5.56e-17. The timing study is not rerun after the verifier failure.

This establishes the repeated Mac Schwarz result for the declared modified
route. It does not resolve the other three fixtures or certify continuum
or local-stress accuracy. Keep all 24 lattice models and all fixed specimens.
The next reference investigation uses the already solved nested and cubic
fields on their original parent23 domain; first prove the physical mappings
and field provenance agree before computing a joint cross-energy diagnostic.
No new solve, polynomial degree, geometry offset or sizing sweep is selected.
Earlier selected-next descriptions below are historical.


## Declared live Schwarz Tet10 pilot, 20 October 2026

After closing the lattice polynomial-mode diagnostic, return to the missing
live stronger-comparator measurement. The maintained runner is
`schwarz_tet10_live.py`. Queue one pair, then read its results before deciding
whether repeated timing is justified. This supplements the original 1% study
with the already declared 0.1% Schwarz target. It does not replace any fixture
or resolve the lattice reference models.

Both arms use fresh geometry and the same explicit elastic-block AMG preset,
CG tolerances, one numerical thread, and 8 GiB worker ceiling. Microgen uses
the admitted VTK30 one-pass override at hausd 0.0005, hmax 0.1, hgrad 1.3,
followed immediately by manual straight Tet10 elevation and direct fedoo.
This is neither the default two-pass route nor public Microgen Tet10 support.
Meshers uses the verified nonadaptive curved16 recipe and fifteen-point rule.
The planar Microgen domain uses its exact four-point rule. Neither arm replays
a saved mesh, solved load prefix, or preconditioner from another run.

Each arm must pass unchanged geometry and periodic gates, six-load residuals,
fresh-process stress/reaction/energy/jump audits, and the same eight Schwarz
reference models using both length proxies and orders 1 through 4. Failed
accuracy produces no matched-target timing ratio. Four-point sensitivity
remains at most 1e-4. Controls precede the live pair.

Record both science-audited total and generation/lift/solve subtotal. The
Microgen total charges fresh parent generation and its audit, then elevation,
final audit and solves. The subtotal uses recorded generation/lift/solve
timers, so duplicate parent/final validation cost remains visible rather
than hidden in a speed claim. Checkpoint I/O inside those phases remains
charged on both arms. Separate post-workflow validation workers are
not timing samples. One pair has no warmup or repeatability claim.

## Completed live Schwarz Tet10 pilot

Job 331 passes 33 controls and runs both complete fresh arms with the same
explicit AMG preset. The actual Microgen parent has 29,204 Tet4 nodes and
136,587 elements. Its manual straight lift has 207,200 nodes. Native curved16
has 42,928 nodes and 27,492 elements. All original geometry and periodic gates
pass. Both six-load tensors pass all eight declared Schwarz models at 0.1%.

| Measurement | Microgen one-pass plus straight Tet10 | Meshers curved Tet10 |
|---|---:|---:|
| Reference-model disagreement | 0.053985–0.062716% | 0.031469–0.093778% |
| Sampled two-sided distance | 0.000961515 | 0.000161891 |
| Minimum corner quality | 0.103871 | 0.449995 |
| Fresh meshing, before elevation | 9.7197 s | 0.2791 s |
| Manual elevation and conversion | 1.1792 s | Not used |
| Six-load solve | 377.3525 s | 52.6394 s |
| Science-audited workflow | 409.8735 s | 55.1039 s |
| Generation/lift/solve subtotal | 388.2513 s | 52.9185 s |
| Solve worker peak RSS | 3.8152 GiB | 1.4316 GiB |

The single-pair total ratio is 7.4382. The phase subtotal ratio is 7.3368,
so duplicate Microgen parent/final audits do not explain most of the difference.
The frozen pilot's label says the subtotal excludes serialization. More
precisely, checkpoint I/O inside solve and generation timers remains included.
The maintained runner now states that scope. Frozen timers are unchanged.
The Microgen tensor is more accurate against these models. Both satisfy the
same target, not identical measured errors. These are single-run measurements,
not a repeated timing result, default Microgen comparison, or four-fixture win.

Fresh-process field audits reproduce reactions within 4.17e-12 for Microgen
and 7.60e-12 for Meshers. Hill–Mandel errors are below 6.54e-12 and periodic
jumps below 5.56e-17. Microgen's exact affine four-point rule agrees with
fifteen-point reintegration within 1.24e-11. All input hashes remain unchanged.

After reading that pair, job 332 independently reintegrates the fresh curved
Meshers fields with positive 64-point and 125-point rules. The final mixed-load
energy shift is 2.89e-9 relative, and the 64-to-125 change is 5.81e-10. Both
are below 1e-4. This fixed-field audit is not a higher-rule equilibrium solve
or a continuum certificate. It takes 27.08 seconds at 0.3367 GiB and is not
added to either workflow timing sample.

Records are in `runs/g3-schwarz-tet10-live-amg-pilot-v1` and
`runs/g3-schwarz-tet10-live-amg-pilot-quadrature-v1`. This pilot justifies
selecting interleaved repetitions with the same recipes and solver. No
repeated stronger-Tet10 speed claim is made before those measurements.

## Declared interleaved repetitions

`schwarz_tet10_repeated.py` selects one fresh warmup per arm followed by five
fresh measured pairs. The first measured pair runs Meshers first, then the
order alternates. Microgen's parent generation remains immediately before its
manual elevation and solve. All recipes, eight-model 0.1% target, numerical
threads, and AMG tolerances match the qualified pilot. Each pair repeats its
independent field and reference audits. No previous mesh, solved load prefix,
or previous-run preconditioner is reused.

Warmups remain visible but do not enter timing statistics. Report median,
minimum, maximum, and each paired ratio for both science-audited workflow
and generation/lift/solve subtotal. A slower pair must remain visible even
if the median favors Meshers. Both arms must qualify in every pair before
a repeated matched-target claim. This remains a supplementary Schwarz study,
not a substitute for the missing four-fixture comparison.

The entire experiment has a 3,600-second awake-time budget, with existing
per-stage timeouts shortened if less budget remains. The 8 GiB, one-thread,
250,000-node and 250,000-tet limits remain unchanged. Preserve partial results
if any stage or the whole budget fails. Compare wall and monotonic clocks
for each pair. A discrepancy above five seconds excludes this timing
experiment rather than silently treating host suspension as faster computing.
No unrelated job, timer, or host power setting changes.

## Completed interleaved Tet10 comparison

Job 333 completes all five measured pairs after one fresh warmup per arm.
Each pair passes 38 controls, unchanged geometry and periodic admission,
twelve fresh prescribed-strain loads, independent field audits, and all eight
Schwarz models at the supplementary 0.1% target. No saved mesh, solved load,
or previous-run AMG hierarchy is replayed. Both arm orders are represented.

| Recorded cost | Microgen median, range | Meshers median, range | Paired Microgen/Meshers ratio |
|---|---:|---:|---:|
| Science-audited workflow | 409.2799 s, 408.8994–412.9563 s | 55.1914 s, 54.7677–55.4881 s | 7.3858–7.5401 |
| Generation/lift/solve subtotal | 387.6108 s, 387.2127–391.2574 s | 53.0263 s, 52.6370–53.3121 s | 7.2806–7.4331 |

Meshers is faster in every measured pair under both cost definitions.
Microgen uses 207,200 straight-Tet10 nodes and 136,587 elements; Meshers uses
42,928 curved-Tet10 nodes and 27,492 elements. Their reference-model ranges
remain 0.053985–0.062716% and 0.031469–0.093778%, respectively. Microgen has
the smaller disagreement. This is a common-target comparison, not identical
measured errors or a default Microgen workflow comparison.

Maximum worker RSS is 3.8143 GiB for the Microgen solve and 1.4571 GiB for
Meshers. The complete experiment takes 2,923.13 seconds, below the 3,600-second
budget. The maximum wall/monotonic discrepancy is below 8.1e-7 seconds.
All source hashes remain unchanged during the run. Warmup data stay visible
but are excluded from every measured statistic.

The separate `verify_schwarz_tet10_repeated.py` recomputes all 96 full
generalized-eigenvalue comparisons, source-derived workflow totals, phase
subtotals, medians, ranges, ratios and warmup exclusion. It verifies 205 hash
records and reads all geometry and field gates. Two verifier controls include
an off-diagonal mixed-load error invisible to diagonal-only comparisons.
Maximum reaction error is 7.60e-12, Hill–Mandel error 6.54e-12, and periodic
jump 5.56e-17. This verifies frozen records, not a new equilibrium solve or
continuum/local-stress accuracy.

The initial progress audit passes, but the later audit records `no_chips`
at a stage transition after 16 minutes. A direct capture reproduces the cause:
the standard `Part.__exit__()` calls `done()` without values, clearing chips
after a populated worker report. Numerical results and timings are unaffected.
Preserve that warning. Retain values across stage transitions before submitting
another experiment; no shared workbench helper was modified here.
The captured reports and helper hash are in
`runs/g3-schwarz-tet10-progress-warning-v1.json`.

The subsequent reporting fix uses `job_parts.worker_part` and public standard
Job methods. Both success and failure regressions fail before the fix and
pass afterward. Stage completion retains elapsed time and the 8 GiB memory
ceiling instead of an empty values object. Selected controllers use it;
older frozen drivers and the shared helper remain unchanged. Queued checks
334 and 335 pass. The latter records six progress reports with live values.
This addresses worker-stage transitions, not a claim about every historical
controller or prolonged reporting under every failure mode.

The repeated-study verifier also accepts the exact archived driver when
maintained source later changes. It requires the same recorded SHA256 and
never substitutes a source archive for changed mesh or field data. Four
verifier controls pass. Fresh verification v2 still checks all 205 records
and 96 comparisons, with seven source records resolved to measured archives.
Original reports and verification v1 remain frozen.

Records are in `runs/g3-schwarz-tet10-live-amg-repeated-v1` and
`runs/g3-schwarz-tet10-live-amg-repeated-verification-v1.json`. This establishes
a repeated advantage for this stronger, explicitly modified Schwarz comparator.
It does not settle Gyroid or spinodoid comparison coverage, the lattice model
envelope, or the full four-fixture objective. Stop tightening Schwarz from
these results alone and return to those unresolved comparisons.

## Original public integration

`Mesh.to_fedoo()` now copies Tet4 or curved Tet10 coordinates and connectivity
into an independent fedoo mesh. The direct route passes the integration checks
on fedoo 1.0.0b2. Microgen 2.0.0b1 still rejects quadratic tetrahedra in both
`SingleMesh.from_pyvista` and `BoxMesh.from_pyvista`.

Job 9 ran the conversion tests, two solid-cell patches and the saved admitted
12-cell curved Gyroid and Schwarz P solid meshes. It used one numerical thread,
an 8 GiB worker ceiling and a 600-second timeout per stage. Peak sampled worker
RSS was 0.851 GiB. The progress audit passed. All frozen step-1 inputs retained
their hashes.

## Results

The integration uses six prescribed engineering strains and independently
integrates their stresses and mixed-load energies. The specimen comparison
uses the saved stress-controlled tensors from the Microgen workflow study.
The cube reference is the analytic isotropic material tensor, E = 2.5 and
nu = 0.25.

| Case | Nodes | Tetrahedra | Maximum relative tensor disagreement | Periodic mid-edge pairs, x/y/z |
|---|---:|---:|---:|---|
| Solid Tet4 | 189 | 768 | 1.33e-15 | 0 / 0 / 0 |
| Solid Tet10 | 1,241 | 768 | 4.44e-16 | 56 / 56 / 56 |
| Gyroid Tet10 | 19,931 | 10,832 | 4.05e-12 | 236 / 248 / 236 |
| Schwarz solid Tet10 | 19,628 | 11,976 | 1.20e-11 | 420 / 420 / 420 |

All periodic displacement jumps, including mid-edge jumps, are below 6e-17.
Periodic equivalence-class force imbalance is below 1.1e-10 relative.
Integrated stress agrees with macroscopic reactions within 1.1e-11 relative.
Mixed-load Hill–Mandel disagreement is below 4.9e-12. All 24 reduced solves
have relative residuals below 1e-12.

VTK and fedoo shape functions agree within 3.4e-16 at corresponding reference
points. Their node order agrees, including all six edge nodes. Their reference
coordinate conventions differ: a fedoo point `(xi, eta, zeta)` maps to the VTK
point `(zeta, 1-xi-eta-zeta, xi)`. Job 8 first compared the same numerical
coordinates in both conventions and failed the solid-patch audit. The test
mapping was corrected, and job 9 used a fresh run directory. Public conversion
code did not change in response to this test failure.

The experiment verifies unique edge-node ownership, matches full six-node
periodic cap faces and checks positive Jacobians at quadrature points, reference
nodes and the centroid. The smallest sampled Jacobians are 2.28e-5 for Gyroid
and 1.85e-5 for Schwarz. Sampling is not a proof of positivity everywhere.

## Reproduce

Use the existing G3 environment and local saved mesh archives from step 1.
From the repository root, select a new output directory:

```sh
qrun -l 'Tet10 public fedoo integration' --progress report \
  --interrupt none --interrupt-reason 'Immutable staged experiment' \
  --check-start --wait env \
  LD_LIBRARY_PATH=/home/kmarchais/dev/meshers/.local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu \
  OPENBLAS_NUM_THREADS=1 timeout 2400 \
  .venv/bin/python experiments/adaptive/tet10_integration.py \
  --run experiments/adaptive/runs/g3-tet10-integration-next
```

The runner records source and input hashes, package versions, worker limits,
full tensors and all six displacement fields. Text records are committed under
`experiments/adaptive/runs/g3-tet10-integration-v2`. Field archives stay local.
The three focused public conversion tests passed. Ruff and formatting checks
passed for the public method and tests.

## Limits and next step

Conversion does not create boundary conditions, transfer boundary tags or
periodic transforms, or validate an arbitrary input mesh. Rotational periodicity
is outside this result. The native mesher and the step-1 admission gates are
unchanged. These tests establish integration and tensor reproduction; local
stress accuracy and continuum error remain unproven.

## Stronger genuine Microgen control, 18 October

Job 152 lifts the admitted actual Microgen VTK/MMG Schwarz mesh to straight
Tet10. Every added node is the exact midpoint of one unique original edge.
The original planar boundary and domain are unchanged. The parent uses the
explicit one-pass override from the matched Schwarz study, not Microgen's
default two-pass route. Its source and MMG-call hashes are preserved.

The lift has 230,337 nodes and 154,012 elements, within the explicitly staged
250,000-node solve limit and unchanged 8 GiB worker ceiling. Independent
geometry and periodic checks pass. Arrays returned by `Mesh.to_fedoo()` feed
all six prescribed-strain solves; no native Meshers generation is called.
This is manual elevation plus direct fedoo, not support added to Microgen's
public Tet10 converters.

Four-point integration is checked against fifteen-point fixed-field energy,
which differs by less than 9.57e-11 relative. Independent reaction disagreement
is 3.27e-11, Hill-Mandel disagreement is 2.81e-11 and periodic mid-edge jumps
are below 5.56e-17. All six reduced solves meet the strict residual check.

Worst mixed-load stiffness disagreement is 0.062583% against the primary
Schwarz reference, and 0.054006% to 0.062702% across all declared fine-window
orders and length proxies. The original linear Microgen parent spans
0.92778% to 0.99022%. Tet10 materially strengthens the competitor.

The solve takes 1,128.92 seconds; the complete worker uses 1,152.73 seconds
and peaks at 6.05 GiB. These are saved-mesh diagnostic costs, not a live
Microgen extraction/MMG/conversion/solve workflow. Do not use them to claim
a new matched-workflow speed ratio.

The stronger control also exposes a Meshers limit. Curved12 reaches 0.13559%
across the same reference models, so it does not meet a stricter 0.1% target
that this Microgen control passes. Select one current-native curved14 pilot
against that same target. The original four-fixture 1% objective remains open.
Text records are in `runs/g3-schwarz-microgen-straight-tet10-v1` and
`runs/g3-schwarz-microgen-straight-tet10-reference-v1.json`.

## Supplementary 0.1% target

The stricter target supplements the original 1% objective; it does not replace
that objective or settle any other fixture. Curved14 and curved15 both pass
geometry and independent field checks, but their model envelopes reach
0.101239% and 0.103035%. Both strictly miss 0.1%, despite passing 1%.
Do not round either into a stricter-target pass. Seven comparison-admission
and target regression cases pass, including the genuine Microgen positive
control and the curved14 negative control.

The existing curved16 reference data indicates a passing 0.093778% envelope.
Job 157 refreshes that known control with the current verified native build
and the prescribed-strain solver. It has 42,928 nodes and 27,492 elements,
H=0.000162, and an 88.01-second live workflow. Its tensor reproduces the
archived sixteen-cell input within 2.67e-11 relative. The declared envelope is
0.031469% to 0.093778%. Job 158 independently checks its fields with reaction
disagreement 1.18e-11 and periodic jumps below 5.56e-17. This also checks the
numerical reproduction of one input to the reference models, not independent
continuum convergence.

No new strict-target workflow speed ratio is claimed against job 152's
saved-mesh solve. Before a live comparison, job 159 tests whether that
Microgen Tet10 parent overresolves its interior. Change only first-pass hmax
from 0.05 to 0.1, preserving the VTK30 seed, hausd=0.0005, hgrad=1.3, explicit
one-pass route and independent admission gates. Audit the actual returned mesh
before selecting any Tet10 lift and solve.

Job 159 passes admission with 29,204 Tet4 nodes, 136,587 elements,
H=0.0009615 and q_min=0.10387. Its exact lift has 207,200 Tet10 nodes,
only 10.0% fewer than the original 230,337-node control. Job 160 runs that
coarser straight-Tet10 control with the same six-load solver and independent
audits. No accuracy or timing conclusion is assigned before reading it.

Job 160 completes the coarser control. Its declared model range is 0.053985%
to 0.062716%, so it still passes the supplementary 0.1% target. The solve
takes 1,026.88 seconds; the worker takes 1,049.97 seconds and peaks at 5.40 GiB.
Reaction disagreement is 1.47e-11, Hill-Mandel disagreement is 1.09e-11 and
periodic jumps remain below 5.56e-17. This modest coarsening reduces nodes
without losing the declared target. It remains a saved-mesh control, not a
paired live workflow measurement. Return to the unresolved lattice recovery
validation rather than spending further runs tightening Schwarz alone.

The next requested priority is a coarse FEA stress-recovery estimate for linear
mesh sizing. The new displacement archives provide integration evidence, but
they are not qualified local-stress references for that study.

## Admitted lattice parent, 18 October 2026

After the cubic lattice reference passes, job 197 tests one more actual
Microgen VTK51/MMG seed. The first call remains `nofem=True`, hausd 0.0015,
hmax 0.06 and hgrad 1.3. Only second-pass hausd changes to 0.001, retaining
the level-set call, required boundaries and the same size controls. This is
an explicitly tuned two-pass recipe, not Microgen defaults or a Meshers seed.

The parent has 33,903 nodes and 158,252 tets. Independent sampled two-sided
distance is 0.00162159, volume error -0.48770% and minimum quality 0.109983.
Sampled Jacobians and periodic nodes and caps pass. Its exact straight Tet10
lift has 238,411 nodes, within the existing 250,000-node reference stage.
Records are in `runs/g3-lattice-tuned51-tet10-seed-hausd001-v1`.

Job 198 uses that unchanged physical mesh for a six-load direct-fedoo control,
with the existing 8 GiB ceiling, 1,800-second timeout, four-point rule exact for
affine Tet10 elasticity and independent higher-rule and field audits. No solve
result, 1% stiffness claim or live workflow timing is inferred from admission.

Job 198 reaches the unchanged 1,800-second timeout after five completed loads.
All five reduced residuals are below 9.94e-13, but there is no complete tensor
or saved field archive. Peak RSS is 6.083 GiB, so this is a time failure, not
an observed memory-ceiling failure. Its partial row remains frozen and excluded
from stiffness ranking. Records are in
`runs/g3-lattice-microgen-straight-tet10-hausd001-v1`.

The experimental controller now offers the already validated elastic-block AMG
preset explicitly, without changing default CG behavior or any benchmark policy.
A new four-point affine Tet10 control compares all six fields and tensors with
the original solver and the analytic material, then independently reintegrates
with fifteen points. That control must pass before a large retry. Every newly
completed large field is archived immediately. AMG uses identical residual
tolerances, vector-first equivalent reactions and explicit completed-cache
release. The large worker still has the same 8 GiB and 1,800-second limits.
Job 199 stops at the new small test because its fixture requests three cells,
below the documented minimum of four. No large solve runs. Correct only that
fixture and retain its failed output. Job 200 repeats the validation-gated
retry with four cells; it is not a live speed comparison.

Job 200 passes all eight small checks and completes six fresh large solves.
The solve takes 530.80 seconds; total saved-mesh study time is 555.78 seconds,
and worker peak RSS is 4.378 GiB. Every field is retained. Original CG and AMG
agree on the five available reaction columns within 2.01e-11 relative. Own-rule
reaction disagreement is 5.89e-12, Hill-Mandel disagreement is 1.94e-12 and
periodic mid-edge jumps are zero. Fifteen-point reintegration changes mixed-load
energy by 2.72e-11 relative. Job 201 reproduces those checks in a fresh process.

The full tensor differs by 0.77217% from cubic24 and ranges from 0.41254% to
1.004354% across all 24 retained models. It narrowly fails the 1% target.
Do not round that into admission. This is stronger accuracy evidence than the
tested Meshers candidates: the finest adaptive Tet10 mesh has worst retained
model disagreement 1.269882%. The new Microgen control has 238,411 nodes,
versus 193,242 for that Meshers mesh. Neither is a target-qualified timing
competitor, and their saved solves use different solver policies.

Records are in `runs/g3-lattice-microgen-straight-tet10-hausd001-amg-v2`,
`runs/g3-lattice-microgen-tet10-field-audit-v1` and
`runs/g3-lattice-tet10-retained-comparison-v1`. The comparison verifies every
candidate's archived field-audit hashes and retains all models. Its worst-load
witnesses are predominantly deviatoric normal strains. Every candidate's worst
model is the retained uniform-seed first-order extrapolation. This does not
justify removing it or distorting geometry to compensate for its bias.

Before changing Meshers, localize same-domain Tet10-to-cubic enrichment energy
in that worst strain. Verify the Galerkin energy identity and quadrature
sensitivity first. This diagnoses discrete error on an unchanged physical
domain, not local stress uncertainty or continuum accuracy.
