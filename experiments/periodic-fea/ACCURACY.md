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

`runs/accuracy-schwarz-d-body64-fields-v1` freezes native64, MMG62seed80 and
this refined native body on identical common material against native80 and
alternate MMG64. It pins 48 analysis source and 510 input hashes, including the
previous comparison's original fields and reports. Four Sobol scrambles use
1,048,576 probes per reference. Read-only body-field revalidation must preserve
the saved audit digest, with exact retained-boundary, volume and mixed-energy
checks, regression tests and the real Pardiso smoke check. Local-stress results
remain pending; no alternative refinement is queued ahead of this comparison.
