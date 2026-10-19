# Periodic box curves and lattice FE error

## Measured periodic-extension creases

Job 271 checks the fixed lattice rather than changing its definition. At 72
paired wall/box-curve nodes per axis, scalar traces match exactly. One-sided
wall normals differ by 48.19 to 109.47 degrees, with a 75.74-degree median.
Moving the inward probes from 2e-6 to 1e-6 changes angles by at most 0.000501
degrees. Projected curve nodes remain within 4.31e-14 of their original
positions and have maximum field-distance residual 2.78e-17.

These measured normal jumps belong to the periodic extension of the fixed
specimen. They are distinct from the compiler's corrected min/max tie
derivative. Do not smooth them away or replace the specimen with a union
that includes additional neighboring struts.

The saved full six-load Tet10-to-cubic correction is on one identical physical
quadratic domain. Positive 125-point integration gives:

| Tetrahedron group | Cell fraction | Correction trace captured | Minimum capture over every mixed load |
|---|---:|---:|---:|
| Incident on a curve corner vertex | 3.1306% | 75.0159% | 70.5900% |
| Incident on a curve edge | 0.7382% | 35.5464% | 25.9636% |
| Cap vertex, but no curve corner vertex | 3.8598% | 1.4813% | 1.1502% |

Trace uses the inverse parent stiffness to normalize all loads. Mixed-load
capture is a generalized eigenvalue of the selected and total error Gram
matrices. Curve-edge cells are a subset of curve-vertex cells, so the table
is not a partition. Own-rule and higher-rule results agree closely.

The predeclared curve-localization screen passes. This identifies one pure
FE-enrichment component near the curves. It does not prove a continuum stress
singularity, eliminate geometric interpolation error, or certify local stress.
It also explains why refining entire cap neighborhoods can waste work: cap
cells away from the actual curves carry little of this measured correction.
Records are in `runs/g3-lattice-periodic-crease-audit-v1`.

## One curve-local graded metric

The coarse source has 204 curve edges, median chord 0.03607 and maximum
0.05161. The previous finest adaptive mesh has 1,024, median 0.006735 and
maximum 0.01591. Select one stronger, geometric curve-local policy:

`target = min(0.1, 0.0025 + 0.8 * periodic_distance_to_curve_samples)`.

Native gradation is also 0.8, matching this metric's slope. This is an explicit
tuned route, not the historical 0.3 default. A same-gradation control omits
the curve field. Geometry remains the raw smooth union with the existing
positive void clamp 0.02 and corrected automatic gradients. No algebraic
symmetric-minimum rewrite or physical crease smoothing is used.

The live curve locator needs only a fresh nonadaptive curved24 mesh. It
extracts wall/cap intersection edges, samples their quadratic interpolants,
and projects samples onto the exact implicit wall while holding box planes.
Its 24,640 samples reduce to 12,212 periodic unique points. Maximum projected
sample spacing is 0.0003122, and maximum projection distance is 0.0001321.
The size callback is periodic and positive. It uses no coarse FE fields.

Job 274's geometry-only controls both pass all unchanged gates:

| Route | Nodes | Tet10 elements | Two-sided sampled distance | Minimum corner quality |
|---|---:|---:|---:|---:|
| Gradation 0.8, no curve field | 50,380 | 27,103 | 0.00033811 | 0.14465 |
| Curve-local metric | 118,475 | 68,415 | 0.00032671 | 0.16959 |

Both have positive sampled Jacobians, matching periodic caps, and no failed
boundary projections. The curve route's minimum Jacobian ratio is 0.26694.
Its geometry-only preparation takes 22.05 seconds, including coarse meshing,
coarse admission, feature projection, final meshing, and final admission.
That is not a full workflow timing or a stiffness-accuracy result. Records
are in `runs/g3-lattice-curve-grading-geometry-v1`.

## Fresh full-tensor pilot

Job 276 regenerates the curve route and solves all six prescribed strains
with the validated explicit elastic-block AMG preset. Fresh coarse meshing,
feature preparation, final generation, independent geometry, hierarchy setup,
and all six final solves are charged. There is no coarse FE solve or free
oracle field. The worker retains the hard 8 GiB ceiling, one numerical thread,
1,200-second timeout, 250,000-tet cap, and existing 250,000-node solve stage.

The completed workflow takes 342.19 seconds, including 319.08 seconds for
the six final loads, and peaks at 3.032 GiB. Job 278 independently reproduces
the reaction tensor within 1.60e-11 relative and Hill-Mandel energy within
3.20e-12, with zero periodic displacement jumps.

Job 279 reintegrates the saved fields with positive 64-point and 125-point
rules, including affine elasticity and monomial controls. The 125-point
energy tensor differs from the solved tensor by only 3.8591e-7 relative.
This is far smaller than the stiffness-target miss. The audit takes 74.34
seconds and peaks at 0.526 GiB. It is not a higher-rule equilibrium re-solve.

The candidate fails the unchanged stiffness target. Its worst disagreement
across all 24 retained models is 1.456109%, versus 1.269882% for the previous
best Meshers candidate. Difference from the actual finest saved tensor is
0.186599%. All six signed generalized energy differences are positive, so
the curve candidate is stiffer in every mixed mode than that finest tensor.
The fine-family models alone give 0.43254% to 1.31617% disagreement. Do not
discard the other model families or claim a timing win from this failed
accuracy result.

Records are in `runs/g3-lattice-curve-grading-solve-v1`,
`runs/g3-lattice-curve-grading-field-audit-v1`, and
`runs/g3-lattice-curve-grading-reference-comparison-v1.json`. True crease
constraints are not added to the native API by this geometric sizing callback.
The positive-rule records are in
`runs/g3-lattice-curve-grading-positive-rule-audit-v1`.

## Size realization and next decision

The saved live candidate has 2,822 curve edges. Median chord is 0.00260035
and maximum is 0.00727168. Median chord divided by the average requested
endpoint target is 1.01675, but 288 edges, or 10.2055%, remain above the native
4/3 split threshold. The worst ratio is 2.85078. The measured point-cloud
spacing is smaller than the requested target, so sparse curve sampling alone
does not explain these oversized edges.

This check does not prove that every remaining split is admissible. Four
adaptation rounds, quality guards, geometry constraints, and subsequent
operations can limit realization. Do not lower the target blindly. After
restart, first inspect native splitting and test one geometry-only control
with more rounds, preserving the metric, quality gates, and mesh cap. Compare
oversized-edge counts, global mesh sizes, and actual admission before deciding
whether another six-load solve is justified. A changed target and changed
bulk grading in the same experiment would not isolate the cause.

Records are in `runs/g3-lattice-curve-size-realization-v1`. The four-specimen
matched-accuracy objective remains unachieved. This is a negative sizing pilot
with useful geometric localization evidence, not a production feature win.
