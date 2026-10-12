# Periodic meshing investigation

The completed first pilot is pueue task 1204, retained in `runs/coarse-v2`.
It uses the packaged PR #5 build, native MMG CLI through Microgen's supported
public two-pass remesher, one thread and no timing warmup. Its five regression
checks passed. Both solid-cell controls recovered the analytic stiffness tensor
to relative error below 7e-16 and passed the periodic and equilibrium audits.

## First periodic Gyroid comparison

The material is a raw Gyroid band `[-0.55, 0.55]` in a unit cube. Geometry
admission requires both sampled boundary distances at most 0.002, volume error
plus empirical uncertainty at most 1%, and minimum tetrahedron quality 0.1.

| Workflow | Grid | Tetrahedra | Generation sample | Min q | Geometry / shape / periodic |
| --- | ---: | ---: | ---: | ---: | --- |
| Meshers, four passes | 16 | 26,008 | 1.103 s | 0.1604 | Pass / pass / pass |
| Meshers, one pass | 24 | 72,757 | 0.686 s | 0.1503 | Pass / pass / pass |
| Microgen periodic VTK + MMG | 24 | 73,049 | 5.766 s | 0.0294 | Fail / fail / pass |

The grid-24 native mesh has forward distance 0.001422, reverse distance 0.001257
and volume error plus uncertainty 0.284%. All six strain-controlled and six
stress-controlled solves pass their physical audits. Its audited homogenization
takes 28.996 s, including assembly and verification.

The grid-24 Microgen/MMG mesh is valid and periodic but has distances 0.008344
and 0.007982, volume error plus uncertainty 1.014%, and 248 cells below quality
0.1. Its 34.325 s homogenization result is diagnostic. The raw generation samples
differ by 8.4 times, but there is no qualified speedup at a common FEA target:
the MMG mesh is excluded, timings are single samples, and no converged reference
tensor has been established.

Native zero and one pass reject at grid 16 with geometry errors 0.002225 and
0.002472. The public MMG grid-16 output fails independent volume validity:
adjacent tetrahedra lie on the same side of a shared face. The input/output
archives and logs remain available; these failures are not repaired silently.

The native grid-16 four-pass tensor differs from the grid-24 one-pass tensor by
3.72% in the worst macroscopic strain-energy combination. The two grid-24
pipelines disagree by 1.45%, but that is not either mesh's FEA error. Their
geometry differs and neither tensor is a converged continuum reference.

## Frozen cap quality

The Microgen/MMG output retains all 3,072 original required cap triangles,
including their coordinates, to 12 decimal places. Its minimum normalized cap
triangle quality is 0.127919; the native grid-24 mesh's minimum cap quality is
0.354760.

For a fixed face of area `A` and sum of squared edges `S`, any attached
tetrahedron has MMG quality at most `4 sqrt(3) A/S`. Its optimal apex is above
the face centroid at height squared `2S/9`. This ceiling equals the normalized
triangle quality. The derivation and actual saved-face comparison are retained
in `cap-analysis.json` and `analyze_caps.py`.

These cap triangles do not force failure of the 0.1 gate. They do prevent the
tested fixed-boundary pipeline from matching the native mesh's minimum quality
0.1503, regardless of interior apex placement. A stronger MMG control needs a
better seed boundary or coupled periodic boundary changes. The present
measurement does not imply an intrinsic limit of MMG with different input.

## Refinement and first-pass distance

Task 1210 completed the zero/four-pass trial at grid 24 and one pass at grids
32 and 48. It also compared the public periodic MMG baseline with the
explicit first-pass-distance control at grid 24. Both MMG configurations use
`hausd=0.001` for the second pass, so the extra first-pass setting is the only
parameter difference between those two controls. Geometry admission stayed
0.002 and all other common gates stay unchanged.

Its six regression checks passed. Five complete meshes have all six saved
displacement fields, with independent stiffness reconstruction differences
below 5e-16, Hill–Mandel errors below 6e-12, force-balance errors below 1.1e-10
and periodic displacement jumps below 1e-15.

| Native configuration | Tetrahedra | Generation sample | Min q | Admission |
| --- | ---: | ---: | ---: | --- |
| Zero passes, grid 24 | No output | No timing claim | 0.0513 at rejection | 58 cells below 0.1 |
| Four passes, grid 24 | 72,757 | 2.130 s | 0.1622 | Geometry, shape and periodic pass |
| One pass, grid 32 | 171,679 | 1.475 s | 0.2737 | Geometry, shape and periodic pass |
| One pass, grid 48 | 543,496 | 4.431 s | 0.1520 | Geometry, shape and periodic pass |

The grid-32 tensor differs from grid 48 by 1.158% in the worst macroscopic
strain-energy combination, or 0.620% in the Frobenius norm. Grid 48 is still
an observed refinement reference, not a certified continuum solution. A further
level is required before interpreting a 2% reference target as a useful accuracy
measurement.

The first-pass-distance control reduces forward/reverse errors from
0.009233/0.008763 to 0.002568/0.002296. Both still fail the unchanged 0.002
boundary gate and the 0.1 shape gate. The modified output has 244 poor cells
and minimum quality 0.0312. Its generation sample costs 11.811 s, versus
7.852 s for its matched baseline. This is a geometry improvement to investigate,
not an admitted FEA or speed result.

## Fine reference and standard MMG adaptation

Task 1211 completed successfully despite the cancelled T3 background wait.
Its analytic solid-cell control, standard MMG adaptation of a native grid-32
seed, and direct grid-64 reference pass all common gates. Six regression checks
pass. Independent reconstruction of all three sets of saved displacement fields
agrees with the recorded tensors below 7e-16. Hill–Mandel differences remain
below 6e-12 and periodic force-balance errors below 1.35e-10.

The standard adaptation uses `hausd=0.0002`, native MMG CLI and Microgen's
required periodic caps. It has no `nofem` or level-set pass. Its total generation
time includes creating the native seed. This control shows that MMG can produce
an admitted periodic mesh from a better seed; the earlier failed Microgen
outputs do not establish an intrinsic failure of MMG.

| Workflow | Grid | Tetrahedra | Generation sample | Audited FEA | Min q | Worst energy difference to grid 64 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Meshers, one pass | 24 | 72,757 | 0.686 s | 28.996 s | 0.1503 | 3.166% |
| Meshers, four passes | 24 | 72,757 | 2.130 s | 28.296 s | 0.1622 | 3.176% |
| Meshers, one pass | 32 | 171,679 | 1.475 s | 75.522 s | 0.2737 | 1.581% |
| Meshers, one pass | 48 | 543,496 | 4.431 s | 401.505 s | 0.1520 | 0.418% |
| Meshers seed + standard MMG | 32 | 481,253 | 35.095 s | 424.589 s | 0.2039 | 0.777% |
| Meshers, one pass reference | 64 | 1,251,854 | 9.425 s | 1,033.062 s | 0.2927 | Reference |

Grid 48 and grid 64 disagree by 0.418% in the worst strain-energy combination.
This is observed reference stability, not a continuum error bound. Grid 32 is
below a 2% difference to this reference; grid 48 and the MMG seed control are
below 1%. All timings remain single samples without warmup. Audited FEA includes
six strain-controlled solves and six stress-controlled cross-checks.

The MMG control has forward/reverse distances 0.000769/0.000491 and volume error
plus empirical uncertainty 0.093%. Its native seed costs 1.466 s and standard
adaptation 33.629 s after diagnostic capture overhead. The grid-64 reference
has distances 0.000297/0.000224 and volume error plus uncertainty 0.094%.

More smoothing at grid 24 improves minimum shape quality modestly but does not
improve the measured tensor. Reducing the elements needed for accurate periodic
stiffness is a more promising target than adding smoothing passes alone.
Grid 48 already has lower stiffness disagreement and lower measured total cost
than this tightly constrained MMG control. MMG sizing still needs tuning before
making a general algorithm comparison.

Reproduce the cross-run comparison with `compare_reference.py`. The retained
`runs/reference-and-seed-v1/cross-run-comparison.json` records artifact hashes,
full relative energy eigenvalues, admission status and timing scope.

## Historical cost comparison against grid 64

The accuracy interpretation in this section is superseded by the
[completed grid-80 comparison](#completed-grid-80-reference-check). Generation
measurements remain valid, but grid 40 no longer passes the strict 1% target
against the finest measured tensor.

Task 1212 completed both controls and all independent saved-field checks.
All common geometry, shape, periodic and physical gates pass. Reconstruction
differences are below 8e-16, energy consistency errors below 5e-12 and periodic
force-balance errors below 9e-11. Displacement jump errors are zero.

| Workflow | Tetrahedra | Generation sample | Audited FEA | Min q | Worst energy difference to grid 64 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Native grid-32 seed + MMG, hausd 0.0005 | 198,635 | 11.880 s | 107.910 s | 0.1499 | 1.339% |
| Direct Meshers, one pass, grid 40 | 326,818 | 2.919 s | 216.870 s | 0.1144 | 0.817% |
| Native grid-32 seed + MMG, hausd 0.0002 | 481,253 | 35.095 s | 424.589 s | 0.2039 | 0.777% |

At roughly 0.8% disagreement to the same reference, direct Meshers uses 32.1%
fewer tetrahedra than the tight MMG seed control. Its measured generation is
12.0 times faster and audited FEA 2.0 times faster. These are single-sample
ratios for this Gyroid and these configurations. MMG has the higher minimum
quality in this pair. Meshers' worst cell still passes the unchanged 0.1 gate.

Loosening the MMG setting reduces generation and FEA cost sharply, but misses
the 1% reference target. This confirms that tuning affects the comparison.
The grid-40 native mesh has distances 0.000550/0.000468 and volume error plus
empirical uncertainty 0.137%. The full comparison is retained in
`runs/accuracy-cost-v1/cross-run-comparison.json`.

Task 1213 completed one warmup and three generation measurements of the roughly
0.8% pair against grid 64. This historical accuracy pairing is superseded by
the grid-80 comparison below. Every canonical mesh exactly matches its saved,
physically audited mesh. Native measured samples are 3.183, 3.239 and 3.242 s; MMG seed samples are
34.554, 34.755 and 34.495 s. Median generation is 3.239 s versus 34.554 s,
a 10.67 times ratio for these configurations. The MMG total includes its seed.
FEA timing was not repeated. The retained `generation-repeats.json` describes
that scope and records the matched mesh hashes.

## Completed grid-80 reference check

Task 1214 completed successfully in 2,813 s. Six regression checks pass. The
reference has 2,372,255 tetrahedra, minimum quality 0.1608 and no cells below
0.1. Its generation sample is 17.239 s and audited FEA takes 2,563.862 s.
Peak worker memory is 6,637 MiB. The declared four-million-element budget and
90-minute timeout allowed this reference; every admission tolerance stayed
unchanged.

Forward/reverse distances are 0.000159/0.000093 and volume error plus empirical
uncertainty is 0.085%. All periodic faces and physics audits pass. Independent
reconstruction of the six saved displacement fields agrees below 3e-16; the
energy consistency error is below 8e-12, periodic force-balance error below
1.7e-10 and displacement jump error zero.

| Candidate | Worst energy difference to grid 64 | To grid 80 | Below 1% of grid-80 tensor? |
| --- | ---: | ---: | --- |
| Direct Meshers, grid 32 | 1.581% | 1.773% | No |
| Direct Meshers, grid 40 | 0.817% | 1.00772% | No |
| Direct Meshers, grid 48 | 0.418% | 0.608% | Yes |
| Direct Meshers, grid 64 | Reference | 0.189% | Yes |
| Native grid-32 seed + MMG, hausd 0.0002 | 0.777% | 0.966% | Yes |
| Native grid-32 seed + MMG, hausd 0.0005 | 1.339% | 1.531% | No |

Grid 64 differs from grid 80 by 0.189%, compared with 0.418% for grids 48 and
64. This supports further observed convergence within the native method, but
does not remove the possibility of shared discretization errors. The retained
`runs/reference-80-v1/cross-run-comparison.json` reproduces all comparisons and
records the reference and candidate hashes.

The grid-40 configuration misses the strict 1% reference target. It must not
be rounded into a pass. The warmed 10.67 times generation ratio remains valid
for the two tested meshes, but no longer establishes a speedup for two meshes
both below 1% of the finest measured tensor. Grid 48 passes that target, with
13% more tetrahedra than the tight MMG control and a lower single generation
sample. Its single audited FEA sample is about 5% faster, and generation plus
audited FEA about 12% faster. These are unrepeated total-cost samples, not a
large established solver or overall speed advantage.

The tight MMG control's margin below 1% is only 0.034 percentage points, while
grid 40 exceeds it by 0.00772 points. The observed 0.189% reference-level shift
is larger than either margin. This is a reason to treat those verdicts as
reference-specific; it is not a quantified bound on reference error. The gate
uses the worst relative strain-energy difference over all macroscopic strain
combinations, not the tensor Frobenius norm. Grid 48 has more margin against
the measured target, but still requires an independent accuracy reference.

## Completed neighbouring-grid checks

Task 1220 completed in 482 s with six regression checks passing. Both grids
pass every geometry, shape, periodicity and physics gate. Independent stiffness
reconstruction agrees below 7e-16, energy consistency errors are below 6e-12,
periodic force-balance errors below 1e-10 and displacement jump errors zero.

| Candidate | Tetrahedra | Generation sample | Audited FEA | Min q | Worst energy difference to grid 80 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Direct Meshers, grid 38 | 285,074 | 2.736 s | 170.223 s | 0.1498 | 1.146% |
| Direct Meshers, grid 42 | 369,590 | 2.985 s | 234.483 s | 0.2548 | 0.894% |
| Native grid-32 seed + MMG, hausd 0.0002 | 481,253 | 35.095 s | 424.589 s | 0.2039 | 0.966% |

Grid 38 misses the strict 1% worst-energy target. Grid 42 passes and has
0.1548 margin above the 0.1 shape gate, substantially more than grid 40. Its
forward/reverse distances are 0.000531/0.000392 and volume error plus empirical
uncertainty is 0.132%. The neighbouring grids confirm that minimum quality
changes non-monotonically, while both tested neighbours pass the declared gate.
They do not establish robustness on other Gyroid widths or other geometries.

Grid 42 has fewer tetrahedra and a higher minimum quality than the tight MMG
control, while both are below 1% of the same grid-80 tensor. Its observed
generation and audited FEA costs are lower, but the native timing samples have
not been repeated. A follow-up repeats one warmup and three generation samples
for this pair, requiring exact identity with the saved audited meshes.
The reference remains native; continuum accuracy and broader superiority
remain unproven. The retained `runs/native-neighbours-v1/cross-run-comparison.json`
records the exact scores, hashes and admission status.

Task 1226 is queued for the grid-42 generation repeats. The non-Claude
Nemotron 3 Ultra artifact review confirmed the neighbouring-grid numbers and
retained three evidence limits: unrepeated timings, a same-method reference,
and one tested Gyroid width. It did not rerun numerical solves.

The next algorithm control is prepared in `runs/mmg-fine-seed-coarsening-v1`
but has not been queued. It generates a native grid-48 seed, then requests MMG
edge sizing `1/32`, with `hausd=0.0005` and the existing required periodic caps.
This separates seed geometry resolution from volume sizing. The unchanged
geometry, shape, periodicity and physics gates will decide whether it reduces
FEA cost while meeting the same grid-80 worst-energy target. Seed generation
remains part of its measured cost. The generation repeat result will be read
before this control enters the queue.

## Independent review and remaining limits

The resumed Claude Sonnet 5.5 review confirmed the recorded comparisons and
artifact-verification results. It inspected retained artifacts without rerunning
numerics; no transcript was available in this worktree. The earlier log entry
about an unavailable review records the initial provider failure. The resumed
review supersedes that status.

The stiffness reference uses Meshers' own discretization. Native candidates
may share discretization errors with that reference, so closeness to its tensor
does not establish better continuum FEA accuracy than MMG. Grid 80 adds a
convergence check within the same method. An independently converged reference
is still needed to test possible shared errors. Reconstruction, equilibrium
and energy audits establish numerical consistency, not continuum accuracy.

The warmed 10.67 times generation ratio compares direct native grid 40 with
standard MMG adapting a native grid-32 seed at `hausd=0.0002`. It is not a
qualified speed ratio against the current Microgen VTK-to-MMG pipeline, whose
tested meshes fail common admission gates. This result covers one Gyroid band,
one thread and one machine. FEA timings use the common audited iterative solver
and remain single samples; they do not predict timings with another solver.

Minimum native quality is not monotonic with resolution: 0.2737 at grid 32,
0.1144 at grid 40 and 0.1520 at grid 48. Grid 40 has only 0.0144 margin above
the 0.1 gate. Neighbouring grids and other geometries must pass before treating
that configuration as a robust default. The MMG control has the higher minimum
quality in the roughly 0.8% pair. These observations support a promising direct
generation result, not a general claim that Meshers has better shape quality.

## Algorithm hypotheses to test

The current native generator clips a regular lattice, moves vertices, and
refines surface edges according to geometry error. `implicit.rs` calls
`adaptive::refine_for`; that routine scores implicit boundary faces and splits
entire periodic edge orbits. It does not allocate interior resolution using the
six elasticity fields. The current experiment supports testing accuracy per
element as a separate objective from minimum element quality.

MMG's documented adaptation combines a geometry/user size map with edge
splitting, collapse and swaps in the surface and volume. A useful next control
would prescribe a size map from the saved elasticity solutions and measure
whether it reaches the same tensor target with fewer unknowns. That result
would guide a native periodic volume-adaptation implementation. It is a
hypothesis, not a measured speed or accuracy advantage.
[MMG algorithm overview](https://mmgtools.org/tutorials/tutorials_algorithm_overview.html)

Periodic restricted Delaunay refinement is another established comparison.
CGAL generates implicit domains on a flat torus, tunes surface/volume criteria,
and offers vertex optimization and sliver removal. A benchmark would need to
export a conforming cut-cell representation and pass the same independent
periodic checks. No CGAL runtime or FEA result is claimed here.
[CGAL periodic mesh manual](https://doc.cgal.org/latest/Periodic_3_mesh_3/index.html)
