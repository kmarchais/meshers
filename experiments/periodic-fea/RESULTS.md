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

Task 1226 completed in 158 s. Each of its eight outputs, including warmups,
exactly matches its saved audited mesh. Native measured generation samples are
3.587, 3.632 and 3.488 s; the MMG control samples are 34.676, 34.949 and
34.663 s. Their medians are 3.587 and 34.676 s, a 9.67 times generation ratio
for two admitted configurations below the measured grid-80 1% target. Native
grid 42 uses 23.2% fewer tetrahedra and has the higher minimum shape quality in
this pair. FEA timing has not been repeated.

The non-Claude Nemotron 3 Ultra artifact review confirmed the neighbouring-grid
numbers before those timing repeats completed. Its remaining limits were
unrepeated timings, a same-method reference and one tested Gyroid width. The
generation repeat now addresses the first limit for generation only. The
review did not rerun numerical solves.

## Completed finer-seed coarsening control

Task 1229 completed in 265 s. A native grid-48 seed is adapted with MMG sizing
`1/32`, `hausd=0.0005` and the existing required periodic caps. This separates
seed geometry resolution from volume sizing. All geometry, shape, periodicity
and physical audits pass, as do six regression checks. Saved-field stiffness
reconstruction agrees below 5e-16, energy consistency below 4e-12, periodic
force balance below 1.1e-10 and displacement jumps are zero.

| Configuration | Tetrahedra | Generation sample | Audited FEA | Min q | Worst energy difference to grid 80 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Direct Meshers, grid 48 seed | 543,496 | 4.431 s | 401.505 s | 0.1520 | 0.608% |
| Grid-48 seed + MMG sizing 32 | 275,637 | 27.578 s | 193.630 s | 0.1552 | 1.440% |
| Direct Meshers, grid 42 | 369,590 | 2.985 s | 234.483 s | 0.2548 | 0.894% |

The coarsened output has forward/reverse distances 0.000535/0.000377 and volume
error plus empirical uncertainty 0.130%. Its seed cost is 4.557 s and standard
adaptation 23.021 s, included in total generation. The output nearly halves
the seed's tetrahedra and measured audited FEA cost, but fails the unchanged
1% reference target. Preserving a finer seed boundary is insufficient in this
control; interior discretization also matters. This is not a qualified speed
advantage at the 1% target. Full scores remain in
`runs/mmg-fine-seed-coarsening-v1/cross-run-comparison.json`.

## Completed elasticity-residual diagnostic

Task 1241 completed in 10.4 s. The diagnostic reads the saved fields from the analytic solid,
direct grids 42 and 48, and the coarsened output. It builds face adjacency in
the periodic quotient mesh, so caps become interior interfaces. It measures
normal stress jumps across interior faces and traction residuals on physical
free surfaces, with face-diameter/area weights and the existing material law.
The largest local indicator over six energy-normalized unit loads supplies a
concentration map. This map is a heuristic for refinement, not a continuum
error bound or a bound for all mixed strains, and omits geometry error.

The diagnostic passed its checks of zero residual on the periodic affine solid,
detected a
periodic non-affine perturbation, and checks nonzero traction on a free isolated
tetrahedron. It verifies saved mesh/field hashes and reconstructs the reported
stiffness tensors before reporting concentration. Its retained results are
`runs/residual-location-v1/residual-analysis.json` and the hashed cell maps.

| Mesh | Top 20% of cells: indicator share | Cells containing 80% of indicator | Boundary-adjacent cells: indicator share |
| --- | ---: | ---: | ---: |
| Direct grid 42 | 42.1% | 59.3% | 10.7% |
| Direct grid 48 | 42.3% | 59.1% | 9.5% |
| Fine-seed MMG coarsened output | 44.7% | 54.5% | 9.3% |

The indicator is broadly distributed rather than concentrated in a few cells.
These are cell-based shares; the boundary-adjacent category includes all
residual contributions assigned to such cells, not just free-face tractions.
The indicator's magnitude is not a measured FEA error percentage.

The review's proposed interior sign change is rejected. The equilibrium
residual is the sum of outward tractions, `sigma_L n + sigma_R (-n)`, giving
the existing difference `(sigma_L - sigma_R) n`. A plus sign with the same
normal would fail the verified affine-solid control, whose current indicator
sum is below 1e-20. The normalization is explicitly `a(u,u) = integral sigma:epsilon`,
the squared energy norm, twice the physical strain energy. No arithmetic change
is needed. The per-cell boundary flag is intentional and named accordingly.
The face-residual approach follows the integration-by-parts reasoning used in
[the residual estimator tutorial](https://jschoeberl.github.io/iFEM/aposteriori/residualEE.html);
that scalar tutorial does not certify this periodic-elasticity indicator.

## Uniform versus residual-guided sizing

The paired test starts from the same native grid-48 seed and keeps
`hausd=0.0005`, gradation 1.3 and all admission/physics gates unchanged. One
worker requests uniform edge length `1/36`. The other supplies a scalar nodal
map with the same base length. Its heuristic coefficient comes from the saved
unit-load indicators divided by element volume and equivalent-edge length
squared; sizes vary with coefficient to the power -1/5 and are clamped between
`1/64` and `1/24`. Volume-weighted nodal averaging combines periodic equivalence
classes, giving identical desired sizes on opposite faces.

The map is frozen with the run. Its mesh must exactly match the generated seed.
The MEDIT vertex mapping is checked bijectively before writing the scalar map,
because Microgen's boundary merge can reorder vertices. The MMG metric flag is
used instead of uniform sizing. Final stiffness is measured with all six loads,
including the usual audits; the indicator never substitutes for that target.

This is an offline algorithm control. The map required the earlier grid-48
audited FEA, which took 401.5 s; that preliminary solve is not included in the
new generation time. The experiment cannot establish a faster one-shot pipeline.
It tests whether targeted volume sizing merits a native implementation or a
cheaper geometry-based predictor. MMG combines supplied isotropic sizes with its
geometry controls as described in its
[metric documentation](https://mmgtools.org/documentation/mmg3d/documentation_mmg3d_metric.html).

Task 1284 completed both workers and reconstructed their tensors and periodic
equilibrium from all six saved fields. All geometry, quality and periodicity
gates passed; the six regression checks passed too. The results are retained in
`runs/mmg-adaptive-sizing-v1`. Comparison against the actual grid-80 tensor gives:

| Sizing from native grid-48 seed | Tetrahedra | Generation, s | Audited 12-solve FEA, s | Minimum quality | Worst relative energy difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| Uniform, h = 1/36 | 422,965 | 25.890 | 305.569 | 0.1654 | 0.89483% |
| Residual-guided, base h = 1/36 | 286,066 | 27.920 | 199.329 | 0.1646 | 1.37853% |
| Direct Meshers grid 42, earlier run | 369,590 | 2.985 | 234.483 | 0.2548 | 0.89433% |

The guided map uses 32.4% fewer cells and the audited FEA takes 34.8% less time
than the uniform control, but it misses the strict 1% observed-reference target.
The different element counts prevent a conclusion about accuracy at equal
complexity. Uniform sizing reaches essentially the same tensor disagreement as
direct grid 42 with 14.4% more tetrahedra. These generation and FEA times are
single samples; only the earlier direct-grid-42 versus tight-MMG generation
comparison has warmups and repeated measurements.

The next test keeps the seed, indicator, geometry tolerance, gradation and
admission gates fixed, and decreases the guided base length to 1/42. Unclamped
three-dimensional size scaling predicts about 454,000 cells, near the uniform
control's complexity. Clamping, fixed caps and MMG geometry controls can change
that count, so the next measured result decides whether further tuning is
useful. The map still requires the earlier 401.5 s offline FEA, and comparison
with grid 80 remains a same-method convergence observation, not a continuum
accuracy certificate.

Task 1295 completed the guided base-42 control in 462.6 s. It produced 472,846
tetrahedra and 92,723 nodes, with minimum quality 0.1476. Generation took
27.881 s and the audited 12-solve FEA took 350.722 s. All admission gates and
six regression checks passed. Independent reconstruction from all saved fields
agreed with the stiffness tensor to 8.4e-16, and periodic displacement jumps
were zero. The actual grid-80 comparison gives 0.73897% worst relative energy
disagreement. Artifacts are in `runs/mmg-adaptive-sizing-42-v1`.

This is below the observed-reference target, but it uses 11.8% more tetrahedra
than uniform base-36 sizing and takes 14.8% longer in the single-sample audited
FEA. It also needs the prior offline FEA that supplied the size map. The result
does not yet establish a benefit from residual-guided allocation at equal
complexity, or a faster direct pipeline. Direct grid 42 remains a competitive
candidate at 0.89433% disagreement with 369,590 tetrahedra.

One final complexity calibration uses guided base length 1/40. Scaling from the
measured base-42 count predicts about 408,500 tetrahedra, within 4% of the
uniform control. The seed and indicator are unchanged, along with all geometry,
quality and physical audits. This tests the map near the uniform control's
complexity before deciding whether native volume adaptation deserves a prototype.

Task 1296 completed guided base-40 sizing in 420.4 s. Its 416,540 tetrahedra
are 1.5% fewer than the uniform control's 422,965. Worst relative energy
disagreement is 0.89435%, essentially equal to uniform sizing's 0.89483%.
Generation took 29.267 s and audited FEA took 313.588 s, both slightly longer
than the uniform control in these single samples. Minimum quality was 0.1080,
just above the unchanged 0.1 gate. Geometry, periodicity, saved-field physics
audits and all six regression checks passed. Artifacts and the hashed tensor
comparison are retained in `runs/mmg-adaptive-sizing-40-v1`.

The scalar residual map has not shown a useful accuracy-per-element advantage
near the target. Its offline FEA cost and weaker minimum-quality margin further
limit its appeal for direct generation. Stop calibrating this map rather than
implementing it in the native mesher on the strength of these results. This
does not rule out anisotropic metrics or a different estimator; neither has
been tested here.

## Fine MMG interior on grid-80 geometry

The remaining reference bias matters more than another scalar-sizing trial.
The next bounded run uses the native grid-80 seed and standard MMG uniform
sizing at h = 1/48, with hausd = 0.0002 and gradation 1.3. Its box-cap triangles
stay required and independently checked. The seed generation remains part of
the reported generation cost. It receives the same six strain loads, six
stress-controlled audits, saved-field reconstruction and admission gates.
The larger seed requires a declared four-million-tetrahedron budget and a
90-minute worker limit; these change resource bounds only.

This tests stiffness sensitivity to a remeshed volume topology while retaining
a much finer periodic seed boundary than earlier MMG controls. It shares the
native seed geometry and is not an independent geometry reference or a
converged continuum certificate. Read the measured count, tensor and audits
before deciding whether a finer MMG reference is warranted.

Task 1297 completed in 1,069.8 s despite the T3 completion wait being lost at
server restart. The retained run is complete and independently verified. MMG
produced 801,119 tetrahedra and 166,641 nodes. Generation took 104.916 s,
including 17.920 s for the native seed, and audited FEA took 820.707 s. Minimum
quality was 0.1314; sampled forward and reverse distances were 0.0001946 and
0.0001535. The volume discrepancy plus empirical uncertainty was 0.09275%.
All admission gates and six regression checks passed. Saved-field stiffness
reconstruction agreed to 7.6e-16, with zero periodic displacement jumps.

The full tensor differs from native grid 80 by 0.55719% in the worst relative
energy direction. All six generalized energy differences are positive: this
MMG mesh is stiffer than grid 80 in every strain combination. Against this
coarser MMG tensor, direct grid 42 disagrees by 0.34143%, compared with 0.89433%
against native grid 80. Direct grid 48 differs by only 0.07051% from the MMG
tensor. These numbers show reference sensitivity; they do not qualify a new
candidate or establish which tensor is closer to the continuum limit.
Retained results are in `runs/mmg-reference48-seed80-v1`, including comparisons
against both tensors and their artifact hashes.

The next run uses MMG h = 1/64 on the same grid-80 seed, with the same hausd,
required caps and all checks. A refinement step is needed before using a
remeshed tensor as an alternative reference. Neither the grid-80 baseline nor
the strict 1% target is replaced while this convergence check is pending.

Task 1381 completed the MMG base-64 reference in 3,105.2 s. Its 2,132,476
tetrahedra and 389,248 nodes passed all admission gates, six regression checks
and saved-field audits. Minimum quality was 0.1838. Generation took 105.720 s,
including 17.762 s for the native grid-80 seed; audited FEA took 2,668.245 s.
Both geometry-distance directions stayed below 0.000186, and volume discrepancy
plus empirical uncertainty was 0.09095%. Retained artifacts are in
`runs/mmg-reference64-seed80-v1`.

MMG base 48 to base 64 changes the worst strain energy by 0.49408%. The finer
MMG tensor differs from native grid 80 by only 0.06301%, with all six energy
eigenvalues still positive. Direct grid 42 differs from it by 0.83250%, and
direct grid 48 by 0.54651%. The tight native-seed MMG control differs by 0.90293%.
These comparisons support the earlier sub-1% observations on this Gyroid across
two fine mesh topologies. They do not provide a continuum error bound: both
fine methods share the native seed geometry, and neither refinement sequence
has a certified remaining error.

## Next algorithm control: HXT from the periodic surface

The scalar residual map did not improve accuracy per element enough to justify
a native implementation. The next control discards the seed's interior and
constructs tetrahedra from its existing boundary with Gmsh's HXT algorithm.
Gmsh 4.15.2 is already installed; its Python adapter and native library hashes
and version are frozen with the run. This experimental workflow does not change
production Meshers or add a package to the benchmark environment.

The cube control runs first, followed by native grid-42 Gyroid surface input.
Both use one thread, uniform target size 1/cells, linear tetrahedra and the
default tetrahedral optimizer. Existing surface triangles must retain a
bijective vertex mapping within 1e-12 L and exactly the same connectivity.
Independent periodic, geometry, overlap and FEA checks then run as usual. A
failed analytic cube control stops the batch before the Gyroid solve.

Timing includes full native seed generation, surface extraction, HXT generation
and boundary-preservation checks. It cannot claim the speed of a hypothetical
surface-only Meshers generator. This tests whether rebuilding the interior can
improve the practical accuracy/cost tradeoff before writing a native equivalent.
HXT is documented as an efficient Delaunay implementation in the
[Gmsh manual](https://gmsh.info/doc/texinfo/gmsh.html#Choosing-the-right-unstructured-algorithm).

Task 1387 finished the experimental batch, but did not produce a Gyroid result.
The cube passed all audits with 292 tetrahedra, unchanged boundary and analytic
stiffness error 6.8e-16. HXT failed during Gyroid volume reconstruction with
`duplicate volume 1 and 2 in reconstructed Brep` followed by `HXT 3D mesh failed`.
The batch's successful process exit reports retained case outcomes, not success
of every meshing algorithm. The Gyroid failure has no FEA or performance result.

The initial adapter described all triangles as one discrete surface, which
does not distinguish the Gyroid's two void networks in HXT's boundary labeling.
The next test registers the negative and positive implicit walls separately,
along with the six cap patches. Node coordinates and triangles stay unchanged;
only boundary entity identifiers differ. The input surface and patch labels
are archived before meshing, including on failure. Rerun the cube and Gyroid
controls to test this explanation, with the same sizing and admission gates.

Task 1401 passed the cube again and reconstructed the main Gyroid volume, but
the unchanged-boundary check rejected its output. Short single-thread diagnosis
runs were performed while the heavy queue was idle. The archived output had
153,473 tetrahedra and omitted 112 prescribed surface vertices; independent
periodic matching also failed because opposite cap node counts differed.
Adding explicit patch-interface curves did not change this outcome.

The retained native grid-42 mesh has three material components in the cut cube:
369,314 tetrahedra in the main component, plus 135 and 141 in two small corner
components. Their corner surface vertices account for the omitted 112 nodes.
The importer had declared only one volume. This is missing material in the
experimental importer, not evidence that HXT optimization removes periodic
boundary vertices. The corrected adapter declares each validated native
component as a separate volume, with its own wall and cap surface entities,
and extracts tetrahedra from every volume. It retains the exact boundary and
all common admission criteria. Numerical verification of this correction is
pending the next queued cube/Gyroid batch.

Task 1413 verified the corrected component importer. Both cases passed
saved-field reconstruction and physical audits. The Gyroid retains every
surface triangle and all three material components; its opposite cap nodes
and triangles are exactly periodic. Geometry distances and volume match the
native grid-42 seed within the unchanged admission criteria.

The HXT output has 153,704 tetrahedra and 43,237 nodes. Generation took 9.982 s,
including the full native seed and experimental boundary/component checks;
audited FEA took 100.800 s. Minimum quality is 0.07282, below 0.1, and worst
relative energy disagreement is 3.31491% against native grid 80 or 3.25139%
against fine MMG. It is rejected for the FEA target despite valid geometry,
periodicity and solver audits. Native grid 42 has 369,590 tetrahedra and
0.89433% disagreement against grid 80. A smaller mesh is not a qualified win.
Retained artifacts are in `runs/hxt-native-surface-components-v3`.

The next control keeps the native grid-42 boundary and requests HXT target
length 1/64 instead of 1/42. This isolates interior resolution. It keeps the
same optimizer, one thread, exact-boundary constraint and common audits.
The cube control still uses a grid-4 seed and target. Measure element count,
quality and tensor before deciding whether HXT optimization merits another
test. Generation timing remains a single sample including the full seed.

Task 1418 completed the finer HXT interior control in 424.1 s. The output
retains all grid-42 boundary triangles and periodic pairs, and passes geometry
and saved-field physics audits. It has 390,910 tetrahedra and 79,812 nodes,
with generation 11.479 s and audited FEA 339.333 s. Minimum quality is 0.007545;
six elements fall below 0.1. Worst relative energy disagreement is 1.34190%
against native grid 80 or 1.27931% against fine MMG. It fails both the minimum
quality and observed 1% accuracy targets. Artifacts are in
`runs/hxt-interior64-seed42-v1`.

At comparable element/node counts, direct native grid 42 has better measured
accuracy and minimum quality, and lower generation and audited solver costs
in the retained samples. Stop calibrating this uniform HXT construction rather
than adopting it as a native replacement. This result covers the tested fixed
surface, sizing and optimizer; it does not exclude other Delaunay algorithms
or optimization settings. There is no qualified HXT speed/accuracy win here.

## Next geometry: periodic Schwarz P

The measured native advantage still covers one Gyroid band. The next pilot
uses the analytic Schwarz P field `cos(2*pi*x)+cos(2*pi*y)+cos(2*pi*z)` in the
same unit cube with raw band [-0.6, 0.6], represented by offset 1.2. It compares
native one-pass grid 32 with standard MMG from that same seed at target 1/32
and hausd 0.0002. All geometry, quality, periodicity, material, solver and
saved-field audits remain unchanged. The analytic volume and boundary samples
come from this new field, not the Gyroid reference.

This first job is an admission and cost pilot. It has no fine Schwarz P tensor
reference and cannot qualify either workflow at 1%. Inspect its actual meshes
before choosing a finer reference or changing resolution. Existing default
Gyroid task selection remains unchanged; the new specimen is selected explicitly.

Task 1421 completed the Schwarz P pilot in 389.1 s. Both workflows passed all
admission gates, six regression checks and saved-field physical audits, with
zero periodic displacement jumps. The analytic sampled solid fraction is
0.343058 with empirical uncertainty 0.0001812. Artifacts are retained in
`runs/schwarz-p-pilot32-v1`.

| Schwarz P pilot, grid 32 | Tetrahedra | Nodes | Generation, s | Audited FEA, s | Minimum quality |
| --- | ---: | ---: | ---: | ---: | ---: |
| Direct one-pass Meshers | 161,328 | 35,380 | 1.107 | 65.562 | 0.3679 |
| Standard MMG from native seed | 319,114 | 69,142 | 19.166 | 223.723 | 0.2007 |

Both retain the seed caps, with maximum sampled geometry distance about
0.0017644 L and reverse distance about 0.00105 L. Relative volume discrepancy
plus empirical uncertainty is 0.1294% for native and 0.0653% for MMG. These
meshes satisfy the unchanged geometry gates, but MMG uses nearly twice as
many elements. Faster single-sample generation and solver timings alone cannot
establish a win at matched FEA accuracy.

The next job calculates native Schwarz P grid 64 with the same field and all
admission/solver settings. It uses the default two-million-element budget and
a 40-minute worker limit. Compare both pilot tensors against its actual saved
tensor before choosing a native candidate or another reference refinement.
As with Gyroid, a same-method fine tensor is a convergence observation, not
a certified continuum reference.

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
