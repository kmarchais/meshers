# Accuracy work after the first two periodic fixtures

The next objective is lower FEA error at comparable complexity and cost.
Passing a 1% effective-stiffness comparison does not mean local stresses are
accurate to 1%. Keep geometry, periodic conformity and element validity as
separate admission checks.

## Metrics

`accuracy_metrics.py` uses the six independently audited displacement fields.
It reports the full stiffness tensor's signed generalized eigenvalues against
the reference. The largest absolute eigenvalue is the worst relative strain
energy disagreement over all macroscopic strain combinations. Engineering
strains use gamma = 2 epsilon for shear. Young, shear and bulk moduli are
computed from compliance; Poisson differences are absolute because relative
errors are unstable near zero.
The full compliance tensor also gives the worst relative complementary-energy
disagreement under macroscopic stress loading.

The second metric samples local stresses at independent Sobol points. For
stress matrices S(x), with one column per macroscopic load, it constructs

```
B = mean(S_reference.T @ material_compliance @ S_reference)
A = mean((S_candidate - S_reference).T @ material_compliance
         @ (S_candidate - S_reference))
local stress score = sqrt(largest generalized eigenvalue(A, B))
```

This compares the complementary-energy norm of local stress differences over
every mixed load, including cross terms. All candidates use the same material
intersection. The report lists excluded material and both coverage fractions.
It does not set stress to zero in omitted material. Four independent scrambled
Sobol sets expose sampling sensitivity. Their spread is a diagnostic, not a
confidence interval. Analytic tests check zero error, multiplicative scaling,
mixed-load sensitivity, invariance to load basis and point location in void.

The measured constitutive law is lambda = mu = 1. Material compliance uses
physical shear stress and engineering shear strain, so the matrix product
equals sigma:epsilon. Physical strain energy includes a factor 1/2; this
factor cancels in the relative metrics. The interior traction residual remains
`(sigma_left - sigma_right) @ n`, the sum of outward tractions on opposite
normals. Element quality and solver residuals are not substitutes for FEA error.

Area-weighted wall-normal disagreement adds a geometry diagnostic relevant to
traction-free walls. It excludes artificial periodic caps and uses the
analytic field's finite-difference gradient at triangle centroids.

## First controlled result

At Gyroid grid 48, four optimization passes retain 543,496 tetrahedra and pass
the existing geometry, quality, periodic and independent saved-field audits.
Minimum quality improves from 0.15198 to 0.16326, but worst stiffness
disagreement against grid 80 changes from 0.607923% to 0.608351%. Generation
takes 14.07 s in this single run. More shape smoothing does not improve the
effective stiffness metric at this resolution.

An initial 131,072-probe comparison on a common domain covering 99.94% of the
reference material gives local stress disagreements of 10.74% for native48,
10.49% for MMG40seed48 and 8.88% for native64, against native80. This initial
run used one Sobol scramble. Its interleaved blocks condition the spatial
distribution and must not be used as sampling-uncertainty estimates. The
current script replaces them with four independent scrambles. Ranking a small
difference also needs a reference-topology sensitivity check. None of these
values certifies continuum stress accuracy.

## Experimental optimizer

Isolated branch `codex/fea-accuracy-optimizer`, initially commit `114f9145`,
changes the vertex cost from a shape-only inverse-quality penalty to

```
sum(det(tet) * sum(squared edge lengths) / h^5 + 0.05 / squared_quality^2)
```

The first term is proportional to quadratic interpolation energy. The shape
barrier, minimum-quality acceptance floor, periodic group moves, analytic wall
projection and PR #5 geometry guard remain. The normalized cost is invariant
to uniform scaling. This is a geometry-only surrogate, not an elasticity error
estimator, and the prototype performs no Delaunay connectivity optimization.

The motivation is [Chen's ODT smoothing work](https://www.math.uci.edu/~chenlong/C2004.html).
The literature establishes an interpolation objective, not a guaranteed gain
on these elasticity fixtures. Test that gain with actual FEA before adopting
the optimizer. `build_candidate.py` records committed source hashes, runs the
implicit core tests, builds an isolated wheel and records its binary hash. It
does not install over the released package or the frozen reference binary.

The completed fixed-grid prototype gives 0.603119% stiffness disagreement and
10.7469% local stress disagreement against native80. The original one-pass
mesh gives 0.607923% and 10.7489% on the same common domain. Against the fine
MMG64 reference, the prototype's local score is slightly worse. These changes
do not justify adopting the interpolation objective. It is removed from the
current prototype source, with its commit and frozen wheel retained for audit.

## Boundary refinement experiment

The original grid48 generator rejects requested tolerance 0.0002 because its
sampled error remains 0.000284176. Curvature edge splits require child quality
0.25 internally, although the profile's final minimum quality is 0.1.

Prototype `fc947b9c` makes the implicit refinement floor follow the requested
minimum quality within 0.1 to 0.25. Specialized legacy geometry keeps its
historical floor. The original shape objective is restored. The isolated build
passes all 20 implicit core tests. At grid48 and requested tolerance 0.0002,
generation succeeds in 4.84 s with about 4.7% more reduced degrees of freedom.

The completed audit has 557,583 tetrahedra, 120,248 nodes and minimum quality
0.15714. Worst stiffness disagreement to native80 is 0.556455%, an 8.5%
reduction from the original grid48 value with 2.6% more tetrahedra. The audited
FEA sample takes 431.27 s. MMG40seed48 remains closer at 0.495875%.

The independent forward boundary distance is 0.00019996; reverse distance is
0.00023631. The stricter 0.0002 experiment gate therefore fails in the reverse
direction. The artifact remains a diagnostic at that threshold. The field
comparison explicitly rechecks all meshes against the original common 0.002
gate and records the override; it does not change the saved artifact's flag.
Future runs separate `--native-geometry-tolerance` from the shared
`--geometry-tolerance` admission gate. Report both measured directions rather
than inferring them from the generator's request.

On the shared common domain, local stress disagreement drops from 10.7497% to
10.7007% against native80. MMG is still closer at 10.4843%. The alternate fine
MMG reference gives the same ranking. This is a small field gain, not an MMG
accuracy win. The next fixed-grid point requested native tolerance 0.00005,
while retaining the original common 0.002 admission threshold. Generation
failed at sampled error 0.0001682935, so it supplies no admitted FEA result.

## Interior refinement experiment

`body_refine.py` tests whether more interior resolution improves the same
grid48 boundary-refined seed. It reads the independently audited grid16 guide
in `runs/accuracy-coarse-guide16-v1`, including all six saved displacement
fields. The older `coarse-v2` artifact has no saved fields and cannot supply
this guide. Rechecking the regenerated fields reconstructs stiffness to
3.38e-16 relative difference, with periodic force imbalance 3.92e-11.

The guide combines traction residuals over all mixed macroscopic loads,
normalizes by the guide stiffness, and transfers the resulting spatial density
by periodic nearest centroid. This is a refinement heuristic, not an error
bound. The prototype bisects complete interior edge stars in disjoint cavities
and rejects children below quality 0.1. It preserves the seed boundary and
uses the original common geometry gate 0.002 and volume gate 0.01.

`runs/accuracy-body-refine48-v1` freezes the first test with at most 15,000 extra
nodes. Admission and six-load field audits must finish before comparing its
stiffness and local stresses with the seed and MMG control against both native80
and MMG64 references. `verify_preserved_boundary` additionally compares every
retained vertex and oriented boundary triangle on the saved meshes.

The prototype times adaptation of archived artifacts. Its estimated live
recipe adds archived seed generation, guide generation and guide FEA times;
this estimate is not a measured generation speed or an adoption result.

The 15,000-node run completes with 623,665 tetrahedra, 135,248 nodes and minimum
quality 0.119638. Independent checks retain all 120,248 seed coordinates and
97,818 oriented boundary triangles exactly. Volume is unchanged. All six-load
mixed energy changes are negative, and saved fields pass reconstruction,
energy, periodic displacement and force-balance audits.

On identical 131,072 probes, the results are:

| Reference | Metric | Boundary seed | Body refinement | MMG40seed48 |
| --- | --- | ---: | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.556455% | 0.501757% | 0.495875% |
| Native80 | Worst local stress disagreement | 10.700738% | 10.459140% | 10.484251% |
| MMG64 | Worst stiffness disagreement | 0.495120% | 0.440262% | 0.433647% |
| MMG64 | Worst local stress disagreement | 10.936880% | 10.695891% | 10.707194% |

Common material coverage exceeds 99.95% for both comparisons. Body refinement
improves both metrics against the seed. Its aggregate local score is slightly
lower than MMG, but their ordering varies among the four independent scrambles.
This small difference does not establish a robust local-stress advantage.
Stiffness remains slightly behind MMG against both reference topologies.
Artifacts are in `runs/accuracy-body-fields-v1`; 27 regression checks and the
optional Pardiso library check pass. The archived CG FEA sample takes 522.46 s.

## Solver calibration

The fedoo startup warning concerns its default solver selection. These
experiments explicitly select `AuditedCG`, so installing a direct solver does
not change their solve time. An isolated environment in `solver-env/pardiso`
contains PyPardiso 0.4.7 and MKL 2026.1.0; the original benchmark environment
and retained CG artifacts remain available.

`pardiso_solver.py` checks every reduced-system residual against the same
2e-9 gate. It hashes the canonical CSR matrix and reuses a factorization only
when the complete matrix matches. Its real-library check covers reuse,
changed-matrix invalidation and reuse after memory release. The first full
periodic calibration replays the admitted grid16 mesh with CG and Pardiso in
`runs/accuracy-pardiso-replay16-v1`. Both workflows retain the independent
geometry and saved-field checks, with one thread and bounded worker resources.
The first replay exposes a fedoo behavior: its automatic API enables its own
default direct backend when factor reuse is available, bypassing the supplied
callback. Its automatic solve records are empty. That run is retained as a
diagnostic. `homogenization.py` disables this built-in shortcut around the
automatic API and requires six callback records; Pardiso's callback still
reuses its own factors. The accuracy reader also rejects missing solve records.

The corrected `runs/accuracy-pardiso-replay16-v2` has twelve audited solves per
workflow and passes 29 tests. CG and Pardiso reuse the exact admitted mesh.
Total homogenization takes 5.2804 s and 2.6429 s respectively; their reduced
system solves take 4.4864 s and 1.7426 s. Peak worker RSS is 278.4 and 369.7 MiB.
Stiffness disagreement is 2.53e-12; volume-weighted local stress disagreement
over all tetrahedra is 3.79e-11. Each Pardiso formulation uses one factorization
for six right-hand sides, with the first factor released before the second.
These are single calibration samples, not repeated solver performance results.
The larger paired replay in `runs/accuracy-pardiso-replay48-v1` completes on the
same 623,665-tetrahedron body mesh. Both workflows retain all twelve solver
records and pass their individual saved-field audits. CG takes 533.10 s for
homogenization and Pardiso 258.44 s. Peak worker RSS increases from 2.39 GiB to
8.38 GiB, within the unchanged 16 GiB worker limit. Pardiso uses one factor for
each six-load formulation; maximum reduced-system residual is 7.40e-13.
The direct paired comparison confirms stiffness agreement to 9.91e-12 and
volume-weighted local stress agreement to 1.19e-10. Revalidation leaves the
saved audit unchanged. These remain single calibration samples.
`--fea-solver pardiso` explicitly selects the callback for future accuracy
experiments; CG remains the default.

## The 20,000-node result

`runs/accuracy-body-refine48-20k-v1` retains the same boundary seed and fixed
grid16 guide, with a 20,000-node budget and audited Pardiso. It completes with
645,319 tetrahedra, 140,248 nodes and minimum quality 0.119638. Volume and all
120,248 seed coordinates and 97,818 oriented boundary triangles are unchanged.
All physical and saved-field audits pass. FEA takes 282.97 s in this sample,
with peak worker RSS 8.79 GiB. Adaptation of saved artifacts takes 5.93 s;
that time still excludes live seed and guide preparation.

`runs/accuracy-body-fields20k-v1` compares all candidates on identical common
material using 1,048,576 probes per reference, split into four independent
scrambles. This finer sampling gives:

| Reference | Metric | Guided body20k | MMG40seed48 |
| --- | --- | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.486935% | 0.495875% |
| Native80 | Worst local stress disagreement | 10.352727% | 10.456603% |
| MMG64 | Worst stiffness disagreement | 0.425450% | 0.433647% |
| MMG64 | Worst local stress disagreement | 10.590793% | 10.664993% |

The body mesh's local score is lower in every independent scramble against
both references. Common material coverage exceeds 99.94%. It uses 0.74% more
tetrahedra and 3.55% more nodes than MMG, while MMG retains the higher minimum
shape quality (0.137901). This is an accuracy lead on the studied fixture at
similar complexity, with the original gates preserved. It does not certify
continuum stress error or establish a production generation speed.

## Length-only control

`runs/accuracy-body-length48-20k-v1` ranks eligible interior edges by length
alone at the same 20,000-node budget. It performs no guide FEA or residual
transfer. The admitted mesh has 642,072 tetrahedra, 140,248 nodes and minimum
quality 0.136380. It preserves every original coordinate and oriented boundary
triangle, with unchanged volume and nonincreasing energy for every mixed load.
Adaptation of the saved seed takes 5.05 s; audited Pardiso FEA takes 285.19 s
with peak worker RSS 8.21 GiB. These are single samples. The adaptation time
excludes live seed generation.

`runs/accuracy-body-marker-fields-v1` compares the length-only control with
the guided meshes and MMG on 1,048,576 identical probes per reference, split
into four independent Sobol scrambles:

| Reference | Metric | Length-only body20k | Guided body20k | MMG40seed48 |
| --- | --- | ---: | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.465206% | 0.486935% | 0.495875% |
| Native80 | Worst local stress disagreement | 10.247238% | 10.352727% | 10.456603% |
| MMG64 | Worst stiffness disagreement | 0.403727% | 0.425450% | 0.433647% |
| MMG64 | Worst local stress disagreement | 10.489613% | 10.590793% | 10.664993% |

Length-only refinement beats both controls in every local-stress scramble
against both references. Common material coverage exceeds 99.94%. It uses
0.23% more tetrahedra and 3.55% more nodes than MMG, whose minimum quality is
0.137901. Thirty regressions pass in the original environment, one optional
library test skips there and passes in the isolated Pardiso environment.
Read-only revalidation preserves the saved audit digest and checks all six
fields. All 52 source hashes and 167 input hashes match the frozen manifest.

Carry forward length-only refinement for this fixture. The coarse FEA guide
does not earn its extra preparation cost. A native implementation and live
generation measurements remain necessary before adopting a production recipe.
The local scores remain reference disagreements, not continuum error bounds.
These archived scores use VTK's original cell order for shared-face probes.
The canonical cell-order comparison below supersedes that sampling policy;
the original reports remain unchanged.

## Native recipe check

The isolated optimizer branch prototypes `interior_refine_nodes`, default zero,
in the native generator. It ranks interior edges by squared length, splits
complete disjoint edge stars and protects every boundary edge. It inserts at
most the requested number of nodes and respects the tetrahedron budget and
the quality floor of max(minimum_quality, 0.1). No optimizer runs after this
nested pass. Core tests check deterministic output, retained coordinates,
oriented boundaries, periodic metadata, volume and budget enforcement. Zero
nodes preserves the existing generation path.

`runs/accuracy-native-body-build-v1` freezes optimizer commit `83be74cb` and
binary hash `5c01b691254222832043fcbaa26daefdff333b54c2436034fea28887f66bcd93`.
Fifteen library tests, 22 implicit tests and eight public Python quality/API
checks pass. `meshers-native-body:48` in `run.py` measures live native generation
with a 20,000-node budget in `runs/accuracy-native-body48-20k-v1`. It requests
native boundary tolerance 0.0002, with the unchanged common admission tolerance
0.002. The archived boundary seed supplies an exact coordinate/boundary check
and the tensor for checking nested energy; it does not supply the generated
mesh. The live run passes independent admission and all six saved-field
checks, with unchanged volume and nonincreasing energy for every mixed load.
It produces exactly the Python control's 140,248 coordinates and 642,072
tetrahedra, allowing for cell order. The six displacement fields agree to
2.67e-13 in relative L2 norm; worst mixed-load stiffness disagreement between
the implementations is 1.79e-15. FEA takes 294.85 s with peak RSS 8.62 GiB in
this sample. The production package and default generator remain unchanged.

### Cell-order-independent stress probing

The first native comparison, `runs/accuracy-native-body-fields-v1`, exposed a
small sampling dependence. VTK chooses different adjacent cells for shared-face
or shared-edge probes when cell rows are permuted. P1 stress jumps across these
faces, so equivalent meshes could get slightly different sampled scores.
A two-tetrahedron regression reproduces the problem. `probe_cells` now sorts
cells by canonical vertex sets before building the VTK locator and preserves
original element IDs. The probe tolerance remains 1e-12. Interior membership
and all admission gates retain their definitions.

`runs/accuracy-native-body-fields-v2` repeats the million-probe comparison with
this policy. Thirty-one regressions pass, one optional library check skips in
the original environment and passes in the isolated Pardiso environment.
All 46 source hashes and 183 input hashes match. Read-only saved-field
revalidation preserves the audit digest.

| Reference | Metric | Native body20k | MMG40seed48 |
| --- | --- | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.465206% | 0.495875% |
| Native80 | Worst local stress disagreement | 10.245754% | 10.456424% |
| MMG64 | Worst stiffness disagreement | 0.403727% | 0.433647% |
| MMG64 | Worst local stress disagreement | 10.488190% | 10.664809% |

Native refinement wins in every scramble against both references. Its stress
scores now match the equivalent Python mesh to roundoff. Coverage remains
above 99.94%. Minimum shape quality is slightly below MMG, 0.136380 versus
0.137901; wall-normal RMS angle is also higher, 0.9966 versus 0.8544 degrees.
Both meshes pass the same common geometry, volume, shape and periodic gates.
The two-metric FEA lead does not mean every mesh diagnostic improves.

`runs/accuracy-native-body-generation-v2` checks live generation time with one
warmup and three measured samples per recipe in alternating order. Each recipe
retains its admitted binary, source and configuration; imports, process startup
and the independent identity check are outside the measured scope. All eight
runs must reproduce their admitted coordinates and connectivity exactly. This
comparison uses standard MMG adaptation of an admitted Meshers seed, rather
than the public Microgen two-pass remeshing path. The first timing attempt,
`runs/accuracy-native-body-generation-v1`, stopped before MMG generation because
the older archive predates `fixtures.py`. The corrected harness uses that
archive's original direct `Case` construction. Its failed output is retained.

All eight runs in v2 reproduce their admitted meshes. All 108 source hashes
and ten input hashes remain unchanged. Three measured samples give:

| Live recipe | Generation samples, s | Median, s | Maximum measured worker RSS, MiB |
| --- | --- | ---: | ---: |
| Native body20k | 5.5775, 5.4360, 5.4352 | 5.4360 | 264.8 |
| Standard MMG40seed48 | 38.5330, 37.8352, 37.9647 | 37.9647 | 575.6 |

Native generation is 6.98 times faster in this one-thread Gyroid calibration,
with the two-metric accuracy lead above at 0.23% more tetrahedra and 3.55% more
nodes. This measures the full live generation recipes, including native
interior refinement, rather than adaptation of saved artifacts. It excludes
independent validation and FEA and does not establish speed on other fixtures.

## Matched total node count

`runs/accuracy-native-body48-matched-nodes-v1` inserts 15,196 interior nodes.
Added to the seed's 120,248 nodes, that exactly matches MMG's 135,444 total
nodes. It has 621,403 tetrahedra, 2.99% fewer than MMG's 640,585, and retains
the same binary, boundary request, common gates and audited Pardiso solver.
All coordinates and oriented boundary triangles from the seed are unchanged.
Volume is unchanged and every mixed-load energy change is negative. Minimum
quality remains 0.136380. The live generation sample takes 5.61 s and FEA
303.71 s, with peak RSS 8.59 GiB. These timings are single samples; the repeated
6.98-times generation calibration above concerns the 20,000-node recipe.

`runs/accuracy-native-body-matched-fields-v1` independently revalidates the
six saved fields without rewriting their audit. Thirty-one regressions pass,
one optional library check skips in the original environment and passes in
the isolated Pardiso environment. All 47 source hashes and 199 input hashes
match. Canonical probing on one million points per reference gives:

| Reference | Metric | Native matched nodes | MMG40seed48 |
| --- | --- | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.482144% | 0.495875% |
| Native80 | Worst local stress disagreement | 10.322556% | 10.456424% |
| MMG64 | Worst stiffness disagreement | 0.420711% | 0.433647% |
| MMG64 | Worst local stress disagreement | 10.562808% | 10.664809% |

The native local score is lower in every scramble against both references.
The accuracy lead therefore survives equal total node count with fewer
tetrahedra. Reduced-system degree counts can still differ because the meshes
have different periodic boundary discretizations. Local disagreement remains
about 10%; this does not certify continuum accuracy.

`runs/accuracy-native-matched-generation-v1` repeats live generation at this
matched node count. Each recipe uses one warmup and three measured runs in
alternating order. All eight meshes match their admitted artifacts exactly;
all 112 frozen source hashes and 21 input hashes remain unchanged.

| Live recipe | Generation samples, s | Median, s | Maximum measured worker RSS, MiB |
| --- | --- | ---: | ---: |
| Native matched nodes | 5.5405, 5.3591, 5.6225 | 5.5405 | 252.2 |
| Standard MMG40seed48 | 39.5061, 39.9121, 37.1645 | 39.5061 | 566.9 |

The native recipe is 7.13 times faster in this one-thread calibration, while
retaining the two-metric lead at equal total nodes. Timing excludes process
startup, imports, independent validation and FEA. The same run verifies exact
public Microgen `BoxMesh` roundtrips and its periodic-node check for both new
Gyroid profiles and the existing Schwarz P grid52 mesh. The Schwarz P analytic
fixture matches Microgen's installed callable exactly on 2,048 fixed probes at
raw band [-0.6, 0.6]. These application checks perform no remeshing.

## Broader fixtures

`prepare_catalogue.py` freezes ten 3D TPMS callables from the installed Microgen
source, choosing sheet widths at an approximately 35% sampled solid fraction.
It also freezes a 52-mode periodic spinodoid from the local Microgen GRF source.
The spectral grid is fixed at 16, seed 20261010, k0 = 4 pi, bandwidth 3 and
amplitude threshold 0.05. Meshing refinement does not regenerate its geometry.
The spinodoid uses Microgen's negative-inside convention and raw offset zero.
Formula and spectral-source hashes distinguish these source versions.

This catalogue prepares coverage; it does not establish FEA convergence or
superiority on those shapes. Every added case still needs admitted meshes,
finer-reference sensitivity and a standard MMG control. Keep the Gyroid
accuracy investigation first before committing the shared queue to a sweep.

After the matched-node Gyroid result, `runs/accuracy-schwarz-p-baseline-fields-v1`
audits the existing Schwarz P grid52 native and standard MMG40seed52 fields.
It compares both on identical common material against the admitted native
grid80 and grid64 references. This tests resolution sensitivity, not alternate
MMG topology sensitivity. The older raw band [-0.6, 0.6] remains fixed; the
catalogue's density-adjusted width is a different specimen. Measure this
baseline before choosing a new native refinement budget or a finer MMG
reference. All four existing meshes pass read-only six-field revalidation.
Their saved audit digests and all 47 source and 157 input hashes are unchanged.

| Native reference | Metric | Native grid52 | MMG40seed52 |
| --- | --- | ---: | ---: |
| Grid80 | Worst stiffness disagreement | 0.510065% | 0.540755% |
| Grid80 | Worst local stress disagreement | 10.759541% | 10.902912% |
| Grid64 | Worst stiffness disagreement | 0.295973% | 0.327001% |
| Grid64 | Worst local stress disagreement | 11.501795% | 11.645746% |

Native grid52 has the lower local score in every scramble against both native
reference levels. It uses 129,672 nodes versus MMG's 123,172, a 5.28% excess,
and 645,684 tetrahedra versus 629,597. Coverage exceeds 99.97%. This baseline
does not yet establish a lead at matched nodal complexity or against another
reference topology.

`runs/schwarz-p-mmg-reference64-seed80-v1` generates one standard MMG64
reference from the same admitted native80 recipe and raw band. Its common
geometry 0.002, volume 0.01 and minimum quality 0.1 gates remain fixed. It uses
audited CG, with the existing finer-reference policy of 16 GiB, 5,400 seconds
per worker and a four-million-tetrahedron generation budget, also used by the
retained grid80 reference. The admitted result has 346,001 nodes, 1,956,659
tetrahedra and minimum quality 0.189098. All common gates and six-field audits
pass. CG FEA takes 2,047.21 s with peak RSS 5.87 GiB. This is one reference
sample, not a repeated solver calibration.

`runs/accuracy-schwarz-p-reference-fields-v1` compares the candidates against
native80 and this alternate MMG64 topology. All four meshes pass read-only
saved-field revalidation; audit digests and all 48 source and 225 input hashes
remain unchanged. Thirty-one regressions and the isolated Pardiso smoke check
pass, with one optional test skipped in the original environment.

| Reference | Metric | Native grid52 | MMG40seed52 |
| --- | --- | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.510065% | 0.540755% |
| Native80 | Worst local stress disagreement | 10.759541% | 10.902912% |
| MMG64 | Worst stiffness disagreement | 0.379784% | 0.406088% |
| MMG64 | Worst local stress disagreement | 11.147708% | 11.251003% |

Native grid52 has the lower local score in every scramble against both
topologies. Common material coverage exceeds 99.96%. Its 5.28% node excess
remains a limitation of this comparison.

`runs/accuracy-schwarz-p-seed50-v1` measures a grid50 seed with the native
accuracy prototype. It requests native boundary tolerance 0.001, with the
same common geometry 0.002, volume 0.01 and quality 0.1 gates, raw band [-0.6,
0.6], one optimization pass and one thread. Audited Pardiso uses the existing
16 GiB and 1,200-second worker bounds. The admitted seed has 115,004 nodes,
571,224 tetrahedra and minimum quality 0.315404. Sampled forward and reverse
boundary distances are 0.0007033 and 0.0003790. All common gates and six-field
audits pass. Its generation sample takes 3.80 s and FEA 275.24 s, with peak RSS
8.32 GiB. These are single samples, not a new speed calibration.

`runs/accuracy-schwarz-p-body50-matched-nodes-v1` inserts 8,168 length-ranked
interior nodes into this seed. The admitted mesh matches MMG's 123,172 nodes
and has 604,833 tetrahedra versus MMG's 629,597, a 3.93% reduction. Minimum
quality is 0.304329. All 115,004 original coordinates and 71,616 oriented
boundary triangles remain exact. Volume changes by -1.11e-16, and all six
eigenvalues of the mixed-load energy change are negative. All common gates,
six saved fields and independent periodic equilibrium audits pass. One live
generation sample takes 4.29 s and Pardiso FEA takes 285.35 s, with peak RSS
8.42 GiB. These single samples do not establish a speed advantage.

`runs/accuracy-schwarz-p-body-matched-fields-v1` compares this body mesh, its
grid50 seed, the original native52 mesh and MMG40seed52 on identical common
material. It uses 1,048,576 probes split across four Sobol scrambles against
each reference. Canonical cell ordering makes point location independent of
cell row order. Read-only body-field revalidation leaves the original audit
digest unchanged; exact retained-boundary checks pass. All 48 source and 375
input hashes remain unchanged. Thirty-one regressions pass, with one optional
test skipped, and the isolated real Pardiso smoke test passes.

| Reference | Metric | Native grid50 seed | Native grid50 body | MMG40seed52 |
| --- | --- | ---: | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.588839% | 0.533844% | 0.540755% |
| Native80 | Worst local stress disagreement | 10.957566% | 10.724960% | 10.902769% |
| MMG64 | Worst stiffness disagreement | 0.458454% | 0.402684% | 0.406088% |
| MMG64 | Worst local stress disagreement | 11.371331% | 11.150560% | 11.251075% |

The matched-node native body mesh has lower stiffness disagreement and lower
local stress disagreement in every scramble against both reference topologies.
Common material coverage exceeds 99.94%. Adding the grid50 candidates changes
the common intersection slightly, so local scores here differ slightly from
the earlier two-candidate comparison. Stiffness scores use the full tensors.
The stiffness margins are small, and matching total nodes does not match
reduced periodic degrees of freedom. MMG still has better wall-normal RMS,
0.7272 degrees versus 1.0415 degrees. Local stress disagreement remains about
11%; this is reference sensitivity evidence, not a continuum certificate.

`runs/accuracy-schwarz-p-matched-generation-v1` measures these admitted native
body and standard MMG recipes with one warmup and three measured live generations
each, in alternating order on one thread. All eight runs reproduce their admitted
coordinates and tetrahedra exactly. Native measured samples are 4.0341, 4.0303
and 4.0702 s; MMG samples are 31.1285, 30.5556 and 30.2592 s. Medians are 4.03 s
and 30.56 s, a 7.57-times ratio. This comparison includes complete live native
generation versus standard MMG adaptation of its freshly regenerated native
seed. Imports, process startup, independent validation and FEA are outside
timing. Three samples on one machine do not establish a generic MMG speed ratio.

Both the new body and seed pass exact public Microgen BoxMesh roundtrips and
periodic-node checks at 1e-10. The installed Microgen Schwarz P callable agrees
exactly with the fixed-width analytic recipe on 2,048 probes. All 154 source
and 20 input hashes remain unchanged. The 3.93% tetrahedron reduction and lower
sampled stress disagreement apply to these admitted meshes at matched total
node counts, with the wall-normal and reference-accuracy limitations above.

`runs/accuracy-schwarz-d-seed48-v1` starts the next case with one grid48 native
seed and audited Pardiso. It uses the actual installed Microgen Schwarz D
callable and the catalogue's fixed raw band of half-width 0.42513263140719755.
The catalogue chooses this width once from a 35% sampled density target; the
field is unchanged across meshing resolutions. Native boundary tolerance is
0.001, with unchanged common geometry 0.002, volume 0.01 and quality 0.1 gates.
The one-thread worker keeps the existing 16 GiB, 1,200-second and two-million-
tetrahedron bounds. The admitted seed has 127,943 nodes, 568,584 tetrahedra and
minimum quality 0.183251. Both sampled boundary distances pass, with forward
0.0004082 and reverse 0.0001482. Solid fraction is 0.349609, and relative volume
error including reference uncertainty is 0.0013474. Shape, geometry, manifold,
nonoverlap and periodic triangle checks pass.

All six field hashes and twelve solver residuals pass. Independent stiffness
reconstruction differs by 4.46e-16, Hill-Mandel energy by 2.71e-14, periodic
force balance by at most 1.27e-12, and displacement jumps are zero. Single-sample
generation takes 5.71 s and Pardiso FEA takes 249.59 s, with peak RSS 7.83 GiB.
These values establish admission, not reference accuracy or repeated speed.

`runs/accuracy-schwarz-d-mmg40-seed48-v1` admits one standard MMG40 control
from the identical grid48 native recipe. Its regenerated seed coordinates and
tetrahedra match the admitted native seed exactly. The native build, frozen
Microgen specimen, native tolerance, common geometry, volume and quality gates
remain the same. MMG uses required box caps, Hausdorff parameter 0.002, hgrad
1.3 and nr enabled, with audited Pardiso and the existing worker bounds.

The admitted control has 85,347 nodes, 428,865 tetrahedra and minimum quality
0.108527. Sampled forward and reverse boundary distances are 0.0012378 and
0.0009738, both below the common 0.002 gate. Volume error including uncertainty
is 0.0024850. All geometry, shape and periodic checks pass. Six field hashes
and twelve actual solver residuals pass, with independent stiffness difference
6.09e-16, Hill-Mandel difference 1.75e-14, force balance at most 2.29e-13 and
zero displacement jumps. Single-sample generation takes 32.10 s and Pardiso
FEA takes 141.76 s, with peak RSS 5.71 GiB. These samples are not repeated
speed measurements.

The native grid48 seed exceeds the control's node count by 49.91%, so adding
interior nodes to it cannot produce a matched-node comparison.
`runs/accuracy-schwarz-d-seed40-v1` admits a coarser native seed with the same
field, build, geometry request, admission gates and solver bounds. It has
73,379 nodes, 326,112 tetrahedra and minimum quality 0.129951. Sampled forward
and reverse boundary distances are 0.0006361 and 0.0004950. Volume error
including uncertainty is 0.0018131. All common geometry, quality and periodic
triangle gates pass, as do six field hashes and twelve actual solver residuals.
Independent stiffness reconstruction differs by 5.10e-16, Hill-Mandel energy
by 1.67e-14, force balance by at most 5.16e-13 and displacement jumps are zero.
Single-sample generation takes 3.01 s and Pardiso FEA takes 86.76 s, with peak
RSS 4.16 GiB. This seed establishes admission and a node budget, not accuracy
against a finer reference.

`runs/accuracy-schwarz-d-body40-matched-nodes-v1` admits all 11,968 requested
length-ranked interior insertions. It has exactly MMG's 85,347 nodes and
375,660 tetrahedra, 12.40% fewer than MMG. Minimum quality remains 0.129951.
All 73,379 seed coordinates and 68,592 oriented boundary triangles remain
exact; volume is unchanged. All six mixed-load relative energy-change
eigenvalues are negative, from -0.00281315 to -0.00049891. The fixed grid40
field, geometry request and common admission gates remain unchanged.

All six field hashes and twelve actual solver residuals pass. Independent
stiffness reconstruction differs by 5.09e-16, Hill-Mandel energy by 8.96e-15,
force balance by at most 4.13e-13, and displacement jumps are zero. One live
generation sample takes 3.42 s and Pardiso FEA takes 112.94 s, with peak RSS
4.51 GiB. The retained boundary and nonincreasing energy establish nested
refinement; they do not establish better local stress accuracy than MMG.

`runs/schwarz-d-reference80-v1` admits a finer native reference with the
same Microgen field, native build, geometry request and common admission gates.
It has 483,391 nodes, 2,394,888 tetrahedra and minimum quality 0.177350. Sampled
forward and reverse boundary distances are 0.0001446 and 0.0000810; volume error
including uncertainty is 0.0008967. All six field hashes and twelve actual
solver residuals pass. Independent stiffness reconstruction differs by
4.53e-16, Hill-Mandel energy by 4.32e-12, force balance by at most 1.60e-10,
and displacement jumps are zero. Audited CG FEA takes 2,555.98 s, with peak RSS
6.49 GiB. The one-thread worker uses the established finer-reference bounds
of 16 GiB, 5,400 seconds and four million tetrahedra.

`runs/schwarz-d-mmg-reference64-seed80-v1` admits an alternate standard MMG64
reference from this native80 recipe. Its regenerated seed coordinates and
tetrahedra match the admitted native80 artifact exactly. The specimen, native
build, geometry request, common gates and CG resource bounds remain fixed.
The alternate mesh has 350,178 nodes, 1,905,153 tetrahedra and minimum quality
0.132911. Sampled forward and reverse distances are 0.0006885 and 0.0004797;
volume error including uncertainty is 0.0013311. All six field hashes and
twelve actual solver residuals pass. Independent stiffness reconstruction
differs by 1.15e-15, Hill-Mandel energy by 2.92e-12, force balance by at most
1.70e-10, and displacement jumps are zero. CG FEA takes 2,405.42 s, with peak
RSS 5.25 GiB. These are admitted reference samples, not repeated solver timings
or continuum error certification.

`runs/accuracy-schwarz-d-body-matched-fields-v1` compares both references and
the native48, MMG40seed48, native40 seed and native40 body candidates. All 48
analysis source hashes and 585 input hashes remain unchanged, including
original sources, fields, audits and dependency sources. The comparison uses
1,048,576 probes split across four Sobol scrambles on identical common material
against each reference. Thirty-one regressions pass, with one optional test
skipped, and the isolated real Pardiso smoke check passes. Read-only body
revalidation preserves the original audit digest, and exact retained-coordinate
and oriented-boundary checks pass.

| Reference | Metric | Native grid48 | Native grid40 seed | Native grid40 body | MMG40seed48 |
| --- | --- | ---: | ---: | ---: | ---: |
| Native80 | Worst stiffness disagreement | 1.214093% | 2.061919% | 1.777094% | 2.040836% |
| Native80 | Worst local stress disagreement | 15.436120% | 15.888108% | 15.109141% | 17.478612% |
| MMG64 | Worst stiffness disagreement | 0.900925% | 1.745647% | 1.461034% | 1.722684% |
| MMG64 | Worst local stress disagreement | 16.130619% | 17.011495% | 16.239412% | 17.974925% |

The matched-node body mesh has lower stiffness and local stress disagreement
than MMG against both reference topologies, including every stress scramble.
Common reference material coverage exceeds 99.82%; geometry omissions remain
separate from the stress comparison. Minimum quality is 0.129951 versus MMG's
0.108527, and wall-normal RMS is 1.5913 degrees versus 2.2331 degrees. Matching
total nodes does not match reduced periodic degrees of freedom. Local stress
disagreement remains 15-16% and stiffness disagreement remains above 1% for
this budget, so these results do not certify low continuum error. The larger
native48 control has lower stiffness disagreement and slightly lower local
stress disagreement against MMG64, with 49.91% more total nodes.

`runs/accuracy-schwarz-d-matched-generation-v1` measures the admitted body
and MMG recipes with one warmup and three measured live generations each, in
alternating order on one thread. All eight generated meshes match their admitted
coordinates and tetrahedra exactly. Native measured samples are 3.1924, 3.2419
and 3.1672 s; MMG samples are 30.5765, 32.4619 and 32.4546 s. Medians are 3.19 s
and 32.45 s, a 10.17-times ratio for these recipes. The MMG timing includes
fresh native48 seed generation and standard adaptation. Imports, process
startup, independent validation and FEA remain outside timing.

Both the body and native40 seed pass exact public Microgen BoxMesh roundtrips
and periodic-node checks at 1e-10. The timing wrapper accepts an explicit
fixed-width field-equivalence definition for Schwarz D; older manifests retain
the Schwarz P check. The independent benchmark formula and installed Microgen
Schwarz D callable agree exactly on 2,048 probes at the fixed width. All 177
source and 20 input hashes remain unchanged. The speed ratio describes this
coarse admitted comparison; it is not a 1% stiffness-qualified adoption result.

`runs/accuracy-schwarz-d-mmg52-seed80-v1` admits standard MMG52 from the
native80 recipe, with 136,124 nodes, 670,453 tetrahedra and minimum quality
0.120484. The regenerated seed coordinates and tetrahedra match native80
exactly. Both common boundary-distance directions pass, at 0.0011168 and
0.0006674; volume error including uncertainty is 0.0017893. Six field hashes
and twelve actual solver residuals pass. Independent stiffness reconstruction
differs by 7.18e-16, Hill-Mandel energy by 1.69e-12, force balance by at most
1.62e-10 and displacement jumps are zero. Single-sample CG FEA takes 678.43 s,
with peak RSS 2.51 GiB; live generation takes 106.63 s. These samples are not
repeated speed measurements.

Its full mixed-load stiffness disagreement is 1.355278% against native80 and
1.038340% against alternate MMG64. It therefore misses the observed 1% target;
local stress has not yet been compared for this control. The fixed field,
native build and common geometry, volume and quality gates remain unchanged.

`runs/accuracy-schwarz-d-mmg62-seed80-v1` admits standard MMG62 from the
identical native80 seed recipe, changing only target sizing from 52 to 62.
It has 324,554 nodes, 1,759,577 tetrahedra and minimum quality 0.143481. The
regenerated native80 seed coordinates and tetrahedra remain exact. Both common
boundary-distance directions pass, at 0.0006991 and 0.0005141; volume error
including uncertainty is 0.0013741. Six field hashes and twelve actual solver
residuals pass. Independent stiffness reconstruction differs by 1.06e-15,
Hill-Mandel energy by 2.40e-12, force balance by at most 1.66e-10 and displacement
jumps are zero. Single-sample generation takes 111.12 s and CG FEA takes
2,104.36 s, with peak RSS 5.18 GiB.

The full mixed-load stiffness disagreement is 0.417868% against native80 and
0.103583% against alternate MMG64. This control passes the observed 1% stiffness
comparison target against both references. Local stress remains unmeasured for
this control, and this result does not certify continuum error. The rough
inverse-resolution-squared estimate used to select 62 was only a sizing
heuristic; the admitted fields determine the result.

`runs/accuracy-schwarz-d-seed64-v1` admits a finer native seed with the same
field, native build, geometry request and common admission gates. It has
267,123 nodes, 1,269,384 tetrahedra and minimum quality 0.180697. Sampled forward
and reverse boundary distances are 0.0002153 and 0.0001628; volume error
including uncertainty is 0.0010532. All common gates pass, as do six field hashes
and twelve actual solver residuals. Independent stiffness reconstruction
differs by 3.72e-16, Hill-Mandel energy by 2.09e-12, force balance by at most
1.40e-10 and displacement jumps are zero. One live generation sample takes
11.22 s and audited CG FEA takes 1,097.22 s, with peak RSS 4.76 GiB. The existing
four-million-tetrahedron, 16 GiB and 5,400-second bounds remain fixed.

The seed's full mixed-load stiffness disagreement is 0.388498% against native80
and 0.077814% against alternate MMG64. Both are below MMG62's disagreement,
with 17.70% fewer total nodes and 27.86% fewer tetrahedra. This is already below
the observed 1% stiffness target at lower complexity. The count leaves 57,431
nodes within MMG62's budget, but adding them is deferred until the existing
local-stress fields are compared.

`runs/accuracy-schwarz-d-low-error-fields-v1` compares native64 and MMG62seed80
against native80 and alternate MMG64. All 48 analysis source hashes and 435 input
hashes remain unchanged, including original fields, audits and dependency
sources. All four meshes pass read-only field revalidation, with unchanged saved
audit digests. Thirty-one regressions pass with one optional test skipped, and
the isolated real Pardiso smoke check passes. The comparison uses 1,048,576
probes split across four Sobol scrambles on the same material intersection.

| Reference | Metric | Native grid64 | MMG62seed80 |
| --- | --- | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.388498% | 0.417868% |
| Native80 | Worst local stress disagreement | 12.785651% | 10.513607% |
| MMG64 | Worst stiffness disagreement | 0.077814% | 0.103583% |
| MMG64 | Worst local stress disagreement | 13.704605% | 11.181109% |

Native64 has lower stiffness disagreement, but MMG62 has lower local stress
disagreement against both references in every scramble. Common reference
material coverage exceeds 99.89%. Native wall-normal RMS is better, 0.9627
degrees versus 1.4462 degrees, and minimum quality is higher. Those geometry
and shape diagnostics do not establish lower FEA field error. The low-error
native profile therefore does not yet meet the combined stiffness and local-
stress objective, despite its lower complexity.

`runs/accuracy-schwarz-d-body64-matched-nodes-v1` admits all 57,431 requested
length-ranked native interior insertions. It matches MMG62's 324,554 nodes,
with 1,512,918 tetrahedra, 14.02% fewer than MMG. Minimum quality remains
0.180697. All 267,123 seed coordinates and 196,416 oriented boundary triangles
remain exact. Volume changes by 5.56e-17, and all six mixed-load energy-change
eigenvalues are negative, from -0.00126415 to -0.00021471. The fixed field,
build, geometry request, common gates and CG bounds remain unchanged.

Six field hashes and twelve actual solver residuals pass. Independent stiffness
reconstruction differs by 3.37e-16, Hill-Mandel energy by 2.91e-12, force balance
by at most 1.68e-10 and displacement jumps are zero. One live generation sample
takes 12.95 s and CG FEA takes 1,764.04 s, with peak RSS 4.82 GiB. The full
stiffness disagreement is 0.262676% against native80 and 0.054269% against
alternate MMG64, both lower than the seed and MMG62. Local stress remains to
be measured; lower mixed energies do not establish a field-accuracy win.

`runs/accuracy-schwarz-d-body64-fields-v1` compares native64, MMG62seed80 and
this refined native body on identical common material against native80 and
alternate MMG64. All 48 analysis source and 510 input hashes remain unchanged,
including the previous comparison's original fields and reports. Four Sobol
scrambles use 1,048,576 probes per reference. Read-only body-field revalidation
preserves the saved audit digest, and exact retained-boundary, volume and
mixed-energy checks pass. Thirty-one regressions pass with one optional test
skipped, and the real Pardiso smoke check passes.

| Reference | Metric | Native grid64 seed | Native grid64 body | MMG62seed80 |
| --- | --- | ---: | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.388498% | 0.262676% | 0.417868% |
| Native80 | Worst local stress disagreement | 12.785678% | 12.329905% | 10.513631% |
| MMG64 | Worst stiffness disagreement | 0.077814% | 0.054269% | 0.103583% |
| MMG64 | Worst local stress disagreement | 13.704631% | 13.275278% | 11.181124% |

Interior refinement improves native local stress disagreement in every scramble,
but MMG remains lower against both references in every scramble. Common reference
material coverage exceeds 99.89%. The three-candidate intersection changes a
few point-location classifications, so the seed and MMG stress values differ
slightly from the prior two-candidate report. The fixed tolerance and all
archived reports remain unchanged. Lower stiffness disagreement and nested
energy improvements do not satisfy the combined local-stress objective.

`runs/accuracy-schwarz-d-guided-body64-v1` admits a residual-guided Python
prototype with all 57,431 requested nodes. It matches the 324,554-node ceiling,
with 1,527,204 tetrahedra and minimum quality 0.139258. All 267,123 seed
coordinates and 196,416 oriented boundary triangles remain exact; volume
changes by 5.56e-17. All six mixed-load energy-change eigenvalues are negative,
from -0.00107963 to -0.00022297. The marker uses the admitted native40 six-load
fields and transfers the worst mixed-load traction residual by periodic nearest
centroid. The field, build, common gates and CG resource bounds remain fixed.

All 74 source hashes and 24 original input hashes remain unchanged.
`run_guided_body.py` revalidates the original guide read-only before the standard
regression, admission and six-load experiment, preserving its audit digest.
Thirty-one regressions pass with one optional test skipped. The wrapper verifies
the actual oriented boundary and rejects a failed row even if the underlying
runner returns zero. Its real completed integration is retained and the helper
is committed as `3f743af3`.

Six guided-body field hashes and twelve actual solver residuals pass. Independent
stiffness reconstruction differs by 3.33e-16, Hill-Mandel energy by 2.06e-12,
force balance by at most 1.46e-10 and displacement jumps are zero. Full stiffness
disagreement is 0.282939% against native80 and 0.051593% against alternate MMG64.
CG FEA takes 1,840.00 s, with peak RSS 4.27 GiB. Python adaptation of archived
meshes takes 18.36 s. Adding archived seed and guide generation and FEA timings
gives an estimated recipe time of 119.35 s, which is not a measured live recipe
and cannot establish generation speed.

`runs/accuracy-schwarz-d-guided-fields-v1` compares this guide-refined mesh,
the native64 seed, native length refinement and standard MMG62seed80 on
identical common material against native80 and alternate MMG64. All 49 source
and 611 input hashes remain unchanged, including original guide and body
reports. Four Sobol scrambles use 1,048,576 probes per reference. Read-only
guided-body revalidation preserves its audit digest; exact retained-boundary,
volume and mixed-energy checks pass. Thirty-one regressions pass with one
optional test skipped, and the real Pardiso smoke check passes.

| Reference | Metric | Native length refinement | Native guided refinement | MMG62seed80 |
| --- | --- | ---: | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.262676% | 0.282939% | 0.417868% |
| Native80 | Worst local stress disagreement | 12.329883% | 12.404998% | 10.513604% |
| MMG64 | Worst stiffness disagreement | 0.054269% | 0.051593% | 0.103583% |
| MMG64 | Worst local stress disagreement | 13.275251% | 13.344148% | 11.181098% |

Guided placement has slightly larger local stress disagreement than length
refinement in every scramble against both references. MMG remains lower than
both in every scramble. Common reference material coverage exceeds 99.89%.
This coarse residual marker does not justify a native port or another budget
tuning run on these results. The combined accuracy objective remains open.

`runs/accuracy-schwarz-d-error-localization-v1` completes a read-only diagnostic
of these same four candidates, references and probes. All 49 source and 612
input hashes remain unchanged, as do all eight per-cell error artifacts.
Thirty-one regressions pass with one optional test skipped. `sample_fields`
optionally returns original located cell IDs; its default result remains
unchanged. The diagnostic reproduces every retained global error Gram matrix
and worst score. Per-probe energy integrates to the native seed's worst mixed-
load quadratic form, and wall bins and original-cell sums conserve total error.
The completed helpers are committed as `32c72978`.

Wall distance uses periodic reference triangle centroids as a proxy, applied
identically to all candidates. Within proxy distance 0.01 of native80 walls,
length refinement reduces anchor error density to 0.02314 versus MMG's 0.02725.
The native excess is in the interior. At distances 0.02-0.04, density is 0.01211
versus MMG's 0.00384; beyond 0.04 it is 0.01011 versus 0.000354. Median native
interior tetrahedron volume stays 3.1789e-7, nearly twice the retained fine-grid
volume 1.6276e-7. Both refinement markers mostly reduce near-wall cell volume
and leave this interior volume unchanged.

At 1e-12 coordinate quantization, MMG62 shares native80 tetrahedra at 24.61% of
the common probes; those probes contribute only 0.159% of its anchor error.
Against alternate MMG64, 22.19% of probes lie in shared tetrahedra, contributing
0.114% of anchor error. Exact-coordinate sharing is 0.034% against native80
and 22.19% against MMG64, reflecting the different coordinate representations.
None of the three native64 candidates shares reference tetrahedra by either
check. Sharing is a topology diagnostic, not a changed point-location or
geometry gate. These measurements show correlated reference regions and finer
MMG interior resolution; they do not prove that its continuum stress ranking
is biased or superior.

`runs/schwarz-d-reference92-v1` admits one finer native reference on a fresh
grid that was not used to seed the MMG candidates. It has 717,331 nodes,
3,619,944 tetrahedra and minimum quality 0.184818. The fixed Microgen field,
native build, geometry request and common admission gates remain unchanged.
Audited CG completes within the existing 16 GiB, 5,400-second and four-million-
tetrahedron bounds. FEA takes 4,643.45 s and peak RSS is 9.76 GiB. Generation
takes 27.93 s as one calibration sample, which does not establish speed.

All 61 frozen source hashes, six saved field hashes and twelve actual solver
residuals pass. Independent stiffness reconstruction differs by 6.93e-16,
Hill-Mandel energy by 3.43e-12, force balance by at most 1.85e-10 and
displacement jumps are zero. These checks admit the reference for comparison;
they do not certify continuum stress accuracy.

`runs/accuracy-schwarz-d-fresh-reference-fields-v1` freezes 50 source files and
687 inputs for the next comparison. It retains native80 and alternate MMG64,
adds fresh native92, and compares the same native64 seed, native length
refinement, guided prototype and MMG62seed80. Four Sobol scrambles use 1,048,576
probes per reference on the same common material. The job also runs regressions,
the real Pardiso smoke check and read-only fresh-reference and guided-body
audits. Exact retained boundary, volume and mixed-energy checks remain required.
The comparison completes with all source and input hashes unchanged. Thirty-one
regressions pass with one optional test skipped, and the real Pardiso smoke
check passes. Read-only fresh-reference and guided-body audits preserve both
original report digests. Both refined candidates retain all 267,123 seed
coordinates and 196,416 oriented boundary triangles. Volume and mixed-energy
checks pass. The two older comparison reports reproduce their prior scores.

| Native92 comparison | Native64 seed | Length refinement | Guided refinement | MMG62seed80 |
| --- | ---: | ---: | ---: | ---: |
| Worst stiffness disagreement | 0.571010% | 0.444732% | 0.464015% | 0.598833% |
| Worst local stress disagreement | 12.423745% | 11.930038% | 11.998160% | 12.317755% |

Both refined native candidates now have lower local stress disagreement than
MMG in every scramble against fresh native92. MMG remains lower against both
older references in every scramble. Reference material coverage ranges from
99.89% to 99.95%. All candidates remain below the observed 1% stiffness target
against all three references, but the local ranking depends on the reference.
This sensitivity supports retaining multiple reference topologies. It does not
establish a combined accuracy win or a continuum stress certificate.

`runs/accuracy-schwarz-d-volume-body64-v1` freezes one bounded diagnostic that
ranks interior edge stars by mean incident tetrahedron volume, at the same
57,431-node insertion ceiling. The localization
showed that both existing markers leave the larger interior cells untouched.
This geometry-only Python prototype uses the same archived native64 seed,
complete disjoint cavities, boundary protection and quality floor. Equal-volume
stars prefer longer edges; a fixed-seed shuffle breaks exact ties to avoid
concentrating insertions by vertex numbering. All 66 source files and 223
original inputs remain unchanged after completion. Read-only seed revalidation
preserves its audit digest. Thirty-three regressions pass with one optional
test skipped, including the two-block case where the larger cells differ from
the longest eligible edge.

All 57,431 insertions are admitted, giving 324,554 nodes and 1,518,981 tetrahedra
at minimum quality 0.180697. All 267,123 seed coordinates and 196,416 oriented
boundary triangles remain exact. Volume changes by 5.56e-17. All six mixed-load
energy-change eigenvalues are negative, from -0.00085492 to -0.00019366. Six
saved field hashes and twelve actual solver residuals pass. Independent
stiffness reconstruction differs by 4.96e-16, Hill-Mandel energy by 3.01e-12,
force balance by at most 1.92e-10 and displacement jumps are zero.

CG FEA takes 1,746.66 s with peak RSS 4.56 GiB, within the existing resource
bounds. Python adaptation of archives takes 14.21 s; this supplies no generation
speed claim. Admission and energy reduction do not establish local accuracy.

`runs/accuracy-schwarz-d-volume-fields-v1` freezes 50 source files and 823 inputs
for one completed field comparison. It adds the volume prototype to the four
retained candidates and compares all five against native80, alternate MMG64
and fresh native92. Four Sobol scrambles use 1,048,576 probes per reference on
identical common material. All source and input hashes remain unchanged.
Thirty-three regressions pass with one optional test skipped, and the real
Pardiso smoke check passes. Read-only volume-body and fresh-reference audits
preserve both original digests. All three refined candidates retain the exact
seed coordinates and oriented boundary triangles. Volume and mixed-energy
checks pass.

| Reference | Metric | Length refinement | Volume refinement | MMG62seed80 |
| --- | --- | ---: | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.262676% | 0.303476% | 0.417868% |
| Native80 | Worst local stress disagreement | 12.329883% | 12.482168% | 10.513604% |
| MMG64 | Worst stiffness disagreement | 0.054269% | 0.035441% | 0.103583% |
| MMG64 | Worst local stress disagreement | 13.275251% | 13.418484% | 11.181098% |
| Native92 | Worst stiffness disagreement | 0.444732% | 0.485790% | 0.598833% |
| Native92 | Worst local stress disagreement | 11.930038% | 12.090420% | 12.317755% |

Volume ranking has larger local stress disagreement than length refinement in
every scramble against all three references. It remains lower than MMG against
native92 but higher against native80 and MMG64. Common reference material
coverage exceeds 99.89%, and the four older candidate scores reproduce exactly.
This marker does not justify a native port or another budget tuning run.

`runs/accuracy-schwarz-d-reference-pair-fields-v1` freezes 50 source files and
879 inputs for a read-only comparison of the three retained reference fields
with one another on the same common material and four Sobol scrambles. The
candidate ranking changes with the reference, so reference-to-reference
stiffness and local stress disagreement must be measured before another mesh
recipe is tuned. All source and input hashes remain unchanged after completion.
Read-only revalidation preserves all three original audit digests. Thirty-three
regressions pass with one optional test skipped, and the real Pardiso smoke
check passes. All three self-comparisons give exactly zero stiffness and stress
disagreement, including every scramble. Common reference material coverage
exceeds 99.89%.

| Reference | Compared field | Worst stiffness disagreement | Worst local stress disagreement |
| --- | --- | ---: | ---: |
| Native80 | MMG64seed80 | 0.315476% | 9.425041% |
| Native80 | Native92 | 0.181483% | 10.897170% |
| MMG64seed80 | Native80 | 0.314484% | 9.412874% |
| MMG64seed80 | Native92 | 0.493963% | 11.912303% |
| Native92 | Native80 | 0.181813% | 10.906078% |
| Native92 | MMG64seed80 | 0.496416% | 11.938862% |

The score normalizes by the named reference, so reversing a pair changes its
value slightly. Reference local stress variation is larger than the measured
candidate margins. Small stiffness differences do not establish low absolute
stress error, and none of these finite meshes supplies a continuum certificate.
All retained reference rankings must accompany later candidate results.

`runs/accuracy-schwarz-d-seed68-v1` freezes 61 source files for one admitted native
grid68 seed at the same fixed field, native build, geometry request, common
gates and CG bounds. It increases base-grid
resolution throughout the domain rather than selecting interior edge stars.
The actual mesh has 296,787 nodes and 1,465,896 tetrahedra at minimum quality
0.134168. It leaves 27,767 nodes below the MMG62seed80 324,554-node budget.
Forward distance is 0.00018866 and reverse distance is 0.00011807. Volume error
including reference uncertainty is 0.00103837; all common and periodic gates
pass. Its boundary differs from the native64 seed, so the original-coordinate
and oriented-boundary invariants apply to later nested refinement of this new
seed, not to comparison between the two base grids.

All 61 source hashes and six field hashes remain unchanged. Thirty-three
regressions pass with one optional test skipped, and all twelve actual solver
residuals pass. Independent stiffness reconstruction differs by 3.65e-16,
Hill-Mandel energy by 2.41e-12, force balance by at most 2.13e-9 and displacement
jumps are zero. CG FEA takes 1,309.29 s with peak RSS 4.23 GiB. Generation takes
12.40 s as one calibration sample, which does not establish speed.

`runs/accuracy-schwarz-d-seed68-fields-v1` freezes 50 source files and 954 inputs
for one completed comparison. It adds native68 to the five retained candidates
and compares all six against native80, alternate MMG64 and native92. Four Sobol
scrambles use 1,048,576 probes per reference on identical common material.
All source and input hashes remain unchanged. Thirty-three regressions pass
with one optional test skipped, and the real Pardiso smoke check passes.
Read-only native68 and retained volume-body audits preserve both original
digests. The existing exact body-boundary, volume and mixed-energy checks pass.

| Reference | Metric | Native64 seed | Native68 seed | MMG62seed80 |
| --- | --- | ---: | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.388498% | 0.285422% | 0.417868% |
| Native80 | Worst local stress disagreement | 12.786201% | 12.618752% | 10.513646% |
| MMG64 | Worst stiffness disagreement | 0.077814% | 0.054859% | 0.103583% |
| MMG64 | Worst local stress disagreement | 13.705123% | 13.511935% | 11.180800% |
| Native92 | Worst stiffness disagreement | 0.571010% | 0.467659% | 0.598833% |
| Native92 | Worst local stress disagreement | 12.424362% | 11.979720% | 12.317939% |

Native68 reduces stress disagreement versus native64 in every scramble against
all three references. It still loses to MMG against both older references in
every scramble and wins against native92 in every scramble. Common reference
material coverage exceeds 99.87%. Adding the new boundary changes the common
material intersection, so retained candidates' stress scores shift slightly;
the tolerance and original reports remain unchanged. The measured reference-
to-reference variation still limits any absolute accuracy claim.

`runs/accuracy-schwarz-d-body68-matched-nodes-v1` freezes 66 source files and
225 original inputs for one admitted native length-ranked interior refinement
of the admitted grid68 seed, using at most the remaining 27,767 nodes. Length ranking
performed best among the tested native64 markers. The same native build, field,
common gates and CG bounds remain fixed. All 27,767 insertions are admitted,
giving exactly 324,554 nodes and 1,579,040 tetrahedra at minimum quality 0.134168.
All 296,787 seed coordinates and 190,512 oriented boundary triangles remain
exact. Volume changes by 1.11e-16. All six mixed-load energy-change eigenvalues
are negative, from -0.00068993 to -0.00010433.

All 66 source hashes and 225 original input hashes remain unchanged. Thirty-
three regressions pass with one optional test skipped. The wrapper revalidates
the original seed read-only, preserving its audit digest, and rejects a failed
row even if the underlying runner returns zero. Completion records the actual
live native generation scope. Six saved field hashes and twelve actual solver
residuals pass. Independent stiffness reconstruction differs by 4.99e-16,
Hill-Mandel energy by 2.10e-12, force balance by at most 1.42e-10 and displacement
jumps are zero. CG FEA takes 1,685.32 s with peak RSS 4.54 GiB. Live native
generation including refinement takes 13.17 s as one calibration sample, which
does not establish generation speed.

`runs/accuracy-schwarz-d-body68-fields-v1` freezes 50 source files and 1,089
inputs for one completed comparison. It adds the refined grid68 mesh to the six
retained candidates and compares all seven against native80, alternate MMG64
and native92. Four Sobol scrambles use 1,048,576 probes per reference on
identical common material. All source and input hashes remain unchanged.
Thirty-three regressions pass with one optional test skipped, and the real
Pardiso smoke check passes. Read-only body68 and original-seed audits preserve
both original digests. All exact retained body boundaries, volume and mixed-
energy checks pass.

| Reference | Metric | Native68 seed | Native68 refinement | MMG62seed80 |
| --- | --- | ---: | ---: | ---: |
| Native80 | Worst stiffness disagreement | 0.285422% | 0.218305% | 0.417868% |
| Native80 | Worst local stress disagreement | 12.618752% | 12.378616% | 10.513646% |
| MMG64 | Worst stiffness disagreement | 0.054859% | 0.099028% | 0.103583% |
| MMG64 | Worst local stress disagreement | 13.511935% | 13.284333% | 11.180800% |
| Native92 | Worst stiffness disagreement | 0.467659% | 0.400472% | 0.598833% |
| Native92 | Worst local stress disagreement | 11.979720% | 11.724949% | 12.317939% |

Refinement reduces native68 stress disagreement in every scramble against all
three references. At matched total nodes, it still loses to MMG against both
older references in every scramble, while winning against native92 in every
scramble. Common reference material coverage exceeds 99.87%. Lower stiffness
disagreement against all three references does not establish a combined win.
The reference-to-reference stress variation and every older result remain
retained. Another grid or insertion-budget adjustment is not justified by
these repeated older-reference misses.

The first source inspection covered the vertex optimizer and adaptive edge
bisection and missed `implicit/intersection/repair.rs`. The broader scan finds
constraint-aware collapse in intersection quality repair. It selects short
edges from tetrahedra below quality 0.1 and checks periodic orbits and feature
strata. This is distinct from general size-based coarsening of an already
admitted mesh. Its existing guards should be considered for reuse. Before
introducing a topology-changing algorithm,
`runs/accuracy-schwarz-d-latest-localization-v1` freezes 49 source files and
1,146 inputs for one read-only localization of the latest seven candidates
against all three references. It reproduces every retained global error Gram
matrix and worst score. Per-probe error integrates to the anchor mixed-load
quadratic form; wall bins and per-cell sums conserve total error. All source
and input hashes and 21 per-cell artifact hashes remain unchanged. Thirty-three
regressions pass with one optional test skipped. The diagnostic completes in
144.02 s. Progress
counts all references and candidates and reports coverage and stress values.
No meshing or solver job is queued at this checkpoint.

The host choice follows the updated computer inventory and verified SSH
procedure. OMEN's shared lane is idle with 26.99 GiB free RAM. Lenovo is
reachable with 10.29 GiB free RAM, but its Scoop pueue shim cannot launch the
queue client, so its queue is unverified. The analysis stays on OMEN with its
local frozen inputs and benchmark environment. `host-choice.json` records
these checks and is included in the pinned inputs. No new session is created.

Against native80, refined grid68 has lower near-wall anchor error density than
MMG, 0.02329 versus 0.02725 within centroid-proxy distance 0.01. At distances
0.02-0.04 it is 0.01201 versus 0.00384, and beyond 0.04 it is 0.01001 versus
0.000353. MMG shares reference cells at 24.61% of these probes; shared cells
contribute only 0.159% of its anchor error. Against alternate MMG64, sharing is
22.20% with 0.114% of anchor error. The newer base grid reduces native bulk
median volume from 3.1789e-7 to 2.6503e-7, but selective refinement still leaves
that bulk median unchanged.

Against fresh native92, no candidate shares reference tetrahedra. The dramatic
MMG bulk advantage shrinks: at distances 0.02-0.04, refined grid68 error density
is 0.01080 versus 0.01001; beyond 0.04 it is 0.00907 versus 0.00743. Its near-wall
value is 0.02082 versus MMG's 0.02763. MMG's bulk median volume is 1.6276e-7,
versus native68's 2.6503e-7. These diagnostics support investigating node
allocation between wall and bulk, while retaining reference correlation as a
limit on the observed rankings.

The next algorithm prototype will test a small, conservative surface-edge
coarsening pass as a geometry-only operation. [MMG's algorithm overview](https://mmgtools.org/tutorials/tutorials_algorithm_overview.html)
includes collapse of short surface edges before adaptation to its size map.
[Dey et al.](https://graphics.stanford.edu/courses/cs448b-00-winter/papers/dey_edge_contract.pdf)
describe link conditions using boundary cones for topology-preserving edge
contraction. The prototype will check surface and boundary-coned volume links,
protect periodic cap vertices, use existing endpoints, reject inverted or
low-quality cavities, and independently validate the resulting geometry. It
will qualify as a new seed through the unchanged common gates.

`runs/accuracy-schwarz-d-surface-coarsen-geometry-v1` freezes 52 source files and
1,221 inputs for one completed geometry-only trial on the audited native68 seed.
It removes at most 2,000 vertices on same-isovalue wall edges shorter than half
the base-grid spacing, 0.5/68. A fixed-seed shuffle breaks exact length ties.
The operation keeps an existing endpoint and uses disjoint vertex cavities.
Local quality must meet the seed minimum, 0.134168, and sampled new-wall forward
distance must remain below 0.001. Periodic cap vertices are protected; their
oriented triangles must remain exact. Boundary components and Euler
characteristic must match the original. Global topology, orientation, overlap,
both boundary distance directions, volume, quality and periodic matching use
the existing common checks. An extra relative volume-change cap of 0.001
constrains this first trial.

Five direct unit checks pass while OMEN's shared lane is idle. They restore an
exact volume mesh after removing a flat boundary-edge midpoint, protect cap
nodes, prevent deleting a material component, and reject disallowed or
misoriented wall geometry. Thirty-eight regressions pass with one optional test
skipped in the real trial. All 2,000 removals are admitted, giving 294,787 nodes
and 1,461,856 tetrahedra. Minimum quality stays 0.134168. Forward distance is
0.00024047 and reverse distance is 0.00017692; all common gates pass. Relative
volume change is 8.22e-6. All oriented cap triangles and remaining coordinates
are exact. Boundary components remain three and Euler characteristic -18.
All 52 source and 1,221 input hashes and the actual mesh hash remain unchanged.

The configured memory bound is 16 GiB and timeout 900 s, but the first RSS
monitor tracks only the Windows Python launcher. Its recorded 7.49 MiB omits
the worker child and cannot establish worker memory usage or enforcement.
The valid geometry result remains archived, and this resource-monitor defect
requires a corrected rerun before advancement. `geometry-result.json`
and `completion.json` describe only geometry qualification; no six-load FEA
result is implied. Remaining coordinates must equal the indexed original
coordinates. All original archives remain pinned. No full FEA solve or native
port is queued before this actual geometry result, and archive adaptation
supplies no generation speed claim.

`runs/accuracy-schwarz-d-surface-coarsen-geometry-v2` freezes 52 source files and
1,279 inputs for an identical geometry rerun. The algorithm and configuration
remain exact; point, tetrahedron, surface and retained-ID arrays must reproduce
v1. The corrected monitor sums the owned worker and recursive children every
two seconds, matching the original benchmark guard. A regression allocates
64 MiB in a real child process and requires the monitor to count it. OMEN's
shared lane is busy with another job, so this check and the trial are queued
together. The original launcher-only report remains unchanged. No full FEA
solve or native port is queued before the corrected result.

The corrected rerun completes in 144.61 s with peak owned process-tree RSS
1,428.61 MiB, below the 16 GiB cap and 900 s timeout. Thirty-nine regressions
pass with one optional test skipped, including the real allocating-child
memory check. All four saved mesh arrays reproduce v1 exactly. All 52 source
and 1,279 input hashes and the actual mesh digest pass. Geometry qualification
remains admitted with the values above. These results justify one bounded
six-load FEA run on this saved geometry before reusing its freed nodes.
The coarsened wall changes material slightly, so nested-energy monotonicity
against the original native68 mesh does not apply. No FEA accuracy or native
generation speed claim follows from this archive operation.
