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

On identical probes, the results are:

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
A larger saved-mesh replay is the next check.

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
