# Periodic stress-recovery sizing on G3

The stopping condition is a faster Meshers workflow at matched accuracy against
the actual Microgen+MMG route. The comparison still needs an admitted Microgen
baseline. Failed competitors are recorded in `MICROGEN_RESULTS.md` and do not
provide a matched-accuracy timing comparison.

## Recovery diagnostic

Job 15 solves the saved 12-cell Tet4 and straight Tet10 controls on identical
piecewise-planar domains. Both controls fail the original geometry admission,
so this is an estimator diagnostic. No local-stress continuum reference is
assumed. Their corner connectivity and all straight mid-edge positions are
checked before comparing fields.

Six prescribed engineering strains produce the full displacement response.
Element Tet4 stresses are recovered at nodes by volume-weighted averaging over
periodic equivalence classes. Four-point quadrature measures the difference
between the recovered field and the original element stress. The same rule
measures the Tet4-to-Tet10 enrichment difference on that fixed domain.

The local six-by-six error matrices use the elastic compliance inner product.
Their traces after normalization by the coarse stiffness tensor give a scalar
ranking that accounts for all six loads. The largest generalized eigenvalue of
their global sum supplies a worst mixed-load energy-norm diagnostic.

| Case | Density rank correlation | Enrichment energy in top 20% by recovery | By element volume | Best possible selection |
|---|---:|---:|---:|---:|
| Gyroid | 0.8853 | 39.88% | 35.28% | 41.07% |
| Schwarz solid | 0.9838 | 52.32% | 23.77% | 53.13% |

Ranking by stress magnitude captures 37.61% for Gyroid and 40.24% for Schwarz.
Recovery therefore adds useful localization information, especially for Schwarz.
This does not establish that remeshing with it will reduce total workflow cost.

The integrated enrichment-error matrix agrees with the difference of the two
stiffness tensors within 9.36e-11 relative for Gyroid and 4.88e-10 for Schwarz.
This checks the fixed-domain variational energy identity across mixed loads.
The six-load solves also reproduce the previously saved control tensors.
Recovery preserves constant stress in an executable check. All original control
hashes are unchanged.

Source and records are `stress_recovery.py` and
`runs/g3-stress-recovery-v1`. Displacement and indicator archives remain local.

## Spatial size input

The candidate API adds `size_field` to `meshers.generate(..., adapt=True)`.
It is a separate vectorized callback returning positive physical edge lengths.
The adapter takes the smaller of that target and its existing curvature and
thickness limits, then applies the existing bounds and gradation. Opposite
periodic traces are sampled and must match. Invalid callback values and Python
exceptions are propagated to the caller.

The sizing pilot derives an element coefficient from recovered error density
divided by the square of the coarse characteristic element length. The target
length is proportional to that coefficient to power -1/5, the three-dimensional
linear-element allocation exponent. Adjacent elements contribute the smallest
target to each periodic nodal class. Targets are rescaled to a chosen median,
clipped to [0.03, 0.1], and interpolated using four periodic nearest nodes with
inverse-squared-distance weights. This is a heuristic size field.

The first pilot compares uniform and recovery-driven median target 0.065 on
Schwarz solid. Both start adaptation from 16 cells with four rounds, unchanged
geometry tolerance 0.002, minimum quality 0.1 and quality guard 0.3. A live
12-cell coarse Tet4 generation and solve are charged to the recovery workflow.
Final admission uses the same independent geometry and periodic checks as step 1.
Jobs 16 and 17 completed in `runs/g3-sizing-schwarz-v1` and
`runs/g3-sizing-schwarz-v2`. The second uses median target 0.04. All four
meshes pass the unchanged geometry, quality and periodic admission checks.

| Target | Sizing | Nodes | Tets | Worst mixed-load disagreement | Total workflow seconds |
|---|---|---:|---:|---:|---:|
| 0.065 | Uniform | 7,104 | 29,139 | 2.6516% | 4.09 |
| 0.065 | Recovery | 5,513 | 23,463 | 2.7578% | 4.92 |
| 0.04 | Uniform | 33,657 | 179,861 | 1.0401% | 36.51 |
| 0.04 | Recovery | 22,664 | 112,134 | 1.0207% | 25.87 |

At the coarse target, recovery does not pay for its extra solve. At the finer
target, recovery uses 32.7% fewer nodes and takes 29.2% less total time with
slightly lower disagreement. These are single-run observations, not repeated
matched-accuracy timings. Neither finer candidate meets the 1% stiffness target.
They do not establish an advantage over Microgen, which still lacks an admitted
comparison mesh. Recovery does not certify local stress accuracy.

The queued Linux rebuild passed 20 core unit tests and 35 focused Python tests.
Tests cover unchanged default adaptation with a constant target, spatial
refinement, periodic matching, invalid lengths and callback error propagation.

## The 1% target

Job 21 tightens the nominal target to 0.038 without changing the estimator,
admission gates or solver. Both candidates now pass the 1% stiffness target.

| Sizing | Nodes | Tets | Worst mixed-load disagreement | Total workflow seconds |
|---|---:|---:|---:|---:|
| Uniform | 35,393 | 189,199 | 0.9492% | 37.77 |
| Recovery | 25,854 | 130,338 | 0.9385% | 31.76 |

Recovery uses 27.0% fewer nodes and takes 15.9% less total time at this common
target. The live coarse solve and size construction remain charged. This is
still one run per configuration and only a comparison against uniform Meshers.
Records are in `runs/g3-sizing-schwarz-v3`. The Microgen comparison is unresolved.

## Gyroid pilot

Job 24 applies the same sizing method at nominal 0.065 to Gyroid. Uniform sizing
uses 9,258 nodes and 34,410 tets, with 6.0894% stiffness disagreement in 7.55 seconds.
Recovery uses 9,238 nodes and 34,343 tets, with 6.0992% disagreement in 9.10 seconds.
Both pass geometry, volume, quality and periodic checks. Recovery provides no
benefit at this target. Records are in `runs/g3-sizing-gyroid-v1`.

Job 25 tightens the target to 0.04. Uniform sizing uses 21,742 nodes and 106,037
tets, with 3.9870% disagreement in 27.58 seconds. Recovery uses 18,944 nodes and
89,694 tets, with 4.2085% disagreement in 26.22 seconds. Both pass admission.
The 12.9% node reduction comes with worse accuracy and only a 4.9% time saving.
This does not establish a matched-accuracy improvement on Gyroid and remains far
from the 1% stiffness target. Records are in `runs/g3-sizing-gyroid-v2`.

The current recovery heuristic helps Schwarz at a common 1% target but has not
shown a useful Gyroid workflow gain. Keep it optional. These results do not
justify replacing curvature sizing or claiming general adaptive superiority.

## Lattice ranking check, 18 October

The saved current-native linear32 lattice field has a recovery indicator, but
no observed same-domain enrichment error yet. Job 136's density-ranked top
20% of elements contain 15.70% of recovery mass and 0.758% of material volume.
That result alone does not validate or invalidate recovery sizing.

An immutable re-ranking of the same indicator gives a different question.
The top 20% by element mass contain 78.12% of recovery mass and 67.32% of
volume. The largest-volume 20% contain 63.01% of recovery mass. Density
ranking and mass ranking are not interchangeable under an element budget.
These are indicator-concentration measurements, not observed error captures.
The script and input hashes are in `runs/g3-lattice-recovery-rank-check-v1`.

Select one fixed-domain enrichment validation before a live sizing pilot.
Use the admitted uniform24 curved mesh only as a source of corner topology.
Derive Tet4 and straight Tet10 controls with the same corners and exact edge
midpoints, keeping their polyhedral domains identical. Do not reuse the
curved mesh's tensor as a straight-control reference. The straight controls
may fail the original geometry gate and remain diagnostic-only.

Retain the earlier six-load compliance-weighted recovery and enrichment
identity. Require independent field consistency, periodic jumps and the
energy identity within their existing tolerances. Before observing this
control, choose the successor rule: recovery's mass-ranked top 20% must
capture more actual enrichment energy than both volume and stress-magnitude
ranking, and density Spearman correlation must exceed 0.5. This rule only
selects whether a live pilot is worthwhile; it does not alter benchmark
admission or certify local stress accuracy. A live pilot, if selected, must
charge its coarse meshing, six-load solve and recovery cost.

Job 161 completes the fixed-domain lattice validation on the uniform24
corner skeleton, with 44,976 elements and 69,025 straight Tet10 nodes.
Recovery's top 20% capture 72.503% of observed enrichment energy, versus
34.845% for volume and 65.201% for stress magnitude. The enrichment oracle
captures 72.815%. Density Spearman correlation is 0.992656. The summed
enrichment matrix matches the Tet4-minus-Tet10 stiffness change within
1.08e-9 relative. Independent stress/energy checks are below 1.82e-11, and
saved-field periodic displacement jump errors are zero for both controls.

These controls still fail the benchmark distance gate at H=0.005344. Their
enrichment and ranking checks validate the diagnostic selector, not physical
continuum error or an admitted workflow. The saved-field jump audit is in
`runs/g3-lattice-stress-enrichment-field-audit-v1`.

## First live lattice sizing pilot

Job 162 uses a 24-cell initial mesh, a live 24-cell coarse solve for recovery,
median target 0.065, native minimum length 0.002, four adaptation rounds and
the unchanged 0.002 distance and 0.1 quality gates. Both arms use the explicit
equivalent symmetric-union representation and positive void clamp 0.02.
No raw-field periodicity fix is claimed. Coarse generation, all six coarse
loads and recovery remain charged to the recovery workflow.

Both candidates pass admission. Uniform sizing has 41,063 nodes and 184,203
Tet4 elements, takes 63.76 seconds total and differs by 3.93949% from the
finer admitted adaptive tensor. Recovery has 41,054 nodes and 184,167 elements,
takes 73.20 seconds and differs by 3.94864%. The recovery coarse preparation
cost is 9.03 seconds. It does not pay for itself here. Job 163's independent
field audits pass, with zero jumps and reaction disagreement below 2.27e-11.

Across all retained reference families, the respective upper model differences
are 5.24158% and 5.24447%. Neither supports the original 1% stiffness target.
The near-equal counts also show that oracle ranking accuracy alone does not
establish a useful size field or a matched-workflow improvement.

An immutable target reconstruction finds that the old hardcoded floor 0.03
clips 4.37% of coarse nodes. The smallest unclipped target is 0.01720. The
final mesh's volume-per-element characteristic length is 0.02003. These
measurements suggest a possible scale mismatch, but do not identify every
active native sizing constraint. Records are in `runs/g3-lattice-sizing-target-check-v1`.

Select one bound-consistency ablation. The pilot's recovery field now uses
the already configured native minimum length, 0.002, rather than an additional
hardcoded 0.03 floor. Keep nominal length, coarse and final cells, gradation,
geometry, quality, thread count and all mesh/resource caps unchanged. Repeat
the uniform arm as a control. The generic helper's default remains 0.03 for
other callers, and frozen earlier runs retain their original protocol. Do not
claim that the floor caused the negative result until reading the ablation.

## Lattice sizing follow-up, 2026-10-18

The minimum-bound ablation does not help. Uniform takes 63.39 seconds and
differs by 3.93949% from the finest adaptive tensor. Recovery takes 72.35
seconds and differs by 3.95706%. Independent saved-field checks pass. Neither
meets 1% across the retained reference models. Removing the extra floor does
not explain a useful sizing improvement.

The final activation trial selects nominal length 0.0224765 from mean element
volume, without looking at stiffness. A count-weighted edge median was rejected
before execution because tiny cut edges dominate that statistic and its
regular-cell budget estimate exceeded the cap. These plans remain archived.

Job 166 exposes a protocol defect. Native adaptive generation does not enforce
the supplied `max_tetrahedra` budget, and the pilot omitted an independent
post-generation check. Uniform returns 424,754 tets and recovery 567,821,
both above the fixed 250,000 cap. Uniform was incorrectly solved and labeled
qualified. Recovery stops at the separate solve-node cap. Both are excluded
from benchmark admission regardless of their geometry metrics. Frozen records
are not rewritten. A new explicit check rejects over-cap coarse and final
meshes before auditing or solving, with regression cases using these counts.
This check does not fix native in-adaptation budget enforcement.

Job 169 subsequently fixes the native defect. Adaptation tracks live tets and
checks the caller's budget after each count-increasing operation, including
pending-wall repair, error splits and topology swaps. It fails at the first
operation crossing the budget rather than finishing an oversized mesh. The
exact uniform replay stops at 250,002 tets before audit or solve. The old
binary fails the new budget regression. The rebuilt code passes 24 core and
45 Python checks, using the public `MeshingError` exception contract. The
initial validation run's only failure was its mistaken `ValueError`
expectation. The runner and tensor-comparison helpers also independently
reject oversized final meshes. Earlier archives and timing claims remain
tied to their recorded native versions.

Stop this sizing-scale tuning path. Strong same-domain recovery localization
did not translate into a demonstrated lattice workflow advantage. The
reference pool also remains uncertain. Its pairwise necessary common tolerance
lower bound is 0.901931%, which does not prove a contradiction at 1% or prove
that any candidate satisfies all models. No model is discarded.
