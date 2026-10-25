# Physical periodic edge normals and the compact lattice

The native mesher now uses the correct physical endpoint normal when a
periodic edge crosses a cell boundary. The fixed lattice recipe produces
41,131 tets instead of 123,694, with better sampled distance and corner quality.
Its smaller quadratic displacement space fails the stiffness target. This is
a native geometry meshing improvement, not a matched-accuracy victory.

## Native defect and verification

Periodic scalar traces can agree at opposite caps while their physical normal
derivatives differ. `Context::gradient` and `Context::gradients` previously
canonicalized the upper cap to the lower cap. The adaptive mesher also reused
the lower-cap cached normal for an edge ending at the upper physical image.
That comparison invented curvature on an exactly straight wall.

The correction keeps scalar values and node identities canonical, evaluates
gradients on their physical side, and compares endpoint normals in the same
physical edge image. No geometry, sizing threshold, quality threshold, or
periodic admission rule changes.

Job 379 reproduces the gradient sign defect. Job 387 reproduces false turning
on a straight wall. Jobs 388 and 389 verify both edge orientations, real crease
turning, the actual sizing path, and the complete native release suite.
There are 101 passing core tests and one ignored research test. Job 391 has
182 passing Python tests and one optional skip. Native Rust 1.95.0 arm64 build
390 enables `experimental-surfaces`. Its binary SHA256 is
`eab1823f1234fd39223a89ec103df7d87db90c2418f59dceb7cfeebf59ceabd7`.

Jobs 392 through 398 evaluate all four fixed cases. Independent replay 403
checks saved coordinates, connectivity, facet orientation, periodic nodes and
caps, corner quality, and byte-identical fixed specimen inputs. Gyroid,
Schwarz, and spinodoid keep identical connectivity and coordinates within
4e-14. The lattice changes as follows, using exactly the same recipe.

| Metric | Previous native | Corrected native |
|---|---:|---:|
| Tets | 123,694 | 41,131 |
| Nodes | 189,180 | 65,382 |
| Minimum corner quality | 0.207354269 | 0.361611654 |
| Sampled two-sided distance | 0.000343162248 | 0.000335879905 |
| Absolute relative volume error | 0.031403518% | 0.031386755% |
| Worst disagreement across all 24 models | 1.283822914% | 2.427923243% |

The fresh corrected quadratic worker takes 46.3401 seconds, including geometry
generation, checks, and six AMG loads. Independent field audits are separate.
Both meshes fail the 1% target. These times do not establish a matched-accuracy
comparison with the admitted actual Microgen route.

## Complete cubic reference experiment still fails

Jobs 399 and 400 verify the complete existing cubic displacement space on the
compact quadratic geometry. The 208,436-node, 41,131-tet mesh stays under both
250,000 admission caps. Geometry, exact old-space embedding, and actual reduced
AMG setup pass. Setup peaks at 5.920 GiB under the unchanged 8 GiB watchdog.

Job 401 completes six fresh loads. Independent job 402 verifies reactions,
energy, periodicity, and positive 64-point and 125-point integration. Job 404
also checks the full Galerkin correction and compares the cubic tensor with
each of the unchanged 24 model tensors, normalized by that model. Worst
disagreement is **1.647040948%**, still failing 1%. The worker takes 354.056
seconds at 5.219 GiB. This is a saved-domain reference solve, not a live native
workflow. No cubic API or reference is promoted.

The `difference_from_cubic` diagnostic in frozen job 402 uses the cubic tensor
as denominator. Job 404 explicitly computes candidate-relative admission
scores and checks the reciprocal generalized eigenvalues. Frozen output 402
is preserved unchanged.

Of the resolved cubic correction energy, 75.40% lies in the 3.45% of cells
touching periodic curves. This localizes a measured correction. It does not
bound remaining FE error or certify continuum or local stress accuracy.

## Variational split preflight

The next quadratic allocation uses exact fixed-old-node energy gains of
complete periodic edge-star splits. All regions remain eligible. The ranking
balances all six parent loads and divides by the physical added-node upper
bound. It does not use the empirical reference tensors or modify sizing.

Three controls in job 405 compare the patch prediction against an independently
assembled Fedoo operator on curved geometry. They include periodic edge copies,
zero affine residual, unchanged old nodes, and direct full energy reduction.
Job 406 selects 8,237 disjoint edge-star orbits at 116,569 nodes and 73,294 tets.
Geometry and parent fields reproduce within 3.4e-16. The volume changes by
2.6e-15 relative, and minimum corner quality is 0.220317784. Positive quadrature,
periodicity, coherent facets, and reduced AMG setup pass at 3.425 GiB peak.

The reconstructed trial energy agrees with the summed patch gain to 3.43e-12
relative. Its largest directional energy reduction is 0.306408133%. This trial
holds the old coefficients fixed. It is not globally equilibrated and is not
an admitted stiffness result. One fresh globally coupled six-load solve is
selected from these gates. Its independent audits and all-model assessment
must be read before a successor is chosen.

All four cases and all 24 model tensors remain in scope. Schwarz retains its
measured comparison win against explicitly modified actual Microgen routes.
Gyroid and spinodoid still lack admitted actual comparators in the tested
routes. Lattice reference convergence and a native all-model accuracy pass
remain unresolved.

Source hashes and the four-case replay are recorded in
[the native evidence receipt](../../handoff/macbook/periodic-edge-normal-evidence-20261025.json).
Frozen run artifacts remain under `experiments/adaptive/runs`.
