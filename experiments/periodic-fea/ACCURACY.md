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
generation succeeds in 4.84 s and the FEA audit starts with about 4.7% more
reduced degrees of freedom. Stiffness and local-stress comparisons remain
pending; this is not yet an accuracy or release recommendation.

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
