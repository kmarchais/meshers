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
selected from these gates. Job 407 completes six fresh loads in a 97.150-second
worker at 3.148 GiB peak. Its independent field audit takes 27.237 seconds,
and its Galerkin and all-model assessment takes 34.516 seconds.

The largest negative stiffness correction is 0.462598229%. Galerkin energy
balance is within 2.493e-9, and all positive quadrature and periodic audits
pass. Worst disagreement against all 24 unchanged models improves from
2.427923243% to **1.962342587%**, still failing 1%. No native candidate, reference,
or comparative timing victory is promoted. The full result is saved in
[the first split evidence receipt](../../handoff/macbook/split-patch-level1-evidence-20261025.json).

Second-level preflight 408 passes at 202,738 nodes and 130,306 tets. Its
largest trial gain is 0.202002161%, and peak RAM is 4.430 GiB. No solve is
selected for that allocation. Of its 12,651 accepted disjoint stars, 613 carry
90% of the normalized gain when sorted by gain alone. The first level has
similar concentration, with 446 of 8,237 stars carrying 90%.

The next preflight retains a fixed 90% of the admissible normalized gain in
gain-per-added-node order. This is an allocation rule, not a change to any
admission threshold. It retains the complete first-level approximation space
and original curved domain. All regions remain eligible, and no reference
tensor enters selection. Actual trial energy and geometry/resource gates must
pass before a coupled solve. No sizing sweep or empirical reference fitting
is selected.

Preflight 409 retains 618 stars and 90.0004509% of the admissible balanced
gain. It adds 4,224 nodes, compared with 86,169 for the full second level.
The smaller mesh has 120,793 nodes and 76,027 tets. Minimum corner quality is
0.151628329. Its largest trial gain is 0.183889578%, energy balance is within
2.012e-12 relative, and all geometry, periodic, and resource gates pass.

Job 410 completes six fresh coupled loads. Independent field and positive
quadrature checks pass. The largest negative stiffness correction is
0.253313240%, and Galerkin balance is within 1.064e-9. Worst disagreement
against all 24 models improves to **1.723104796%**, still failing 1%. The solve
worker takes 102.836 seconds at 3.827 GiB peak, with separate 26.631-second
field audit and 34.195-second assessment. This is not a live workflow timing
or accuracy admission. See [the allocation receipt](../../handoff/macbook/split-patch-allocation-20261025.json)
and [the second-level evidence](../../handoff/macbook/split-patch-bulk90-level2-evidence-20261025.json).
One third-level preflight is selected from these fresh fields using the same
fixed 90% rule. Its actual gates must be read before another solve.

Third-level preflight 411 passes at 129,269 nodes and 81,751 tets. It retains
1,210 stars and 90.007396% of the admissible balanced gain. Its largest trial
gain is 0.131802511%, and minimum corner quality is 0.102818534. The domain,
parent fields, positive quadrature, periodicity, and actual AMG resource gates
pass. Preflight worker time is 109.672 seconds at 3.717 GiB peak.

Fresh coupled job 412 improves worst all-model disagreement to **1.549508124%**,
still failing 1%. The largest negative stiffness correction is 0.175968727%,
and Galerkin balance is within 1.600e-9. All independent field and positive
quadrature checks pass. Solve, audit, and assessment workers take 115.365,
28.850, and 34.983 seconds respectively. Solve peak RAM is 4.115 GiB. See
[the third-level evidence](../../handoff/macbook/split-patch-bulk90-level3-evidence-20261025.json).
A fourth-level preflight continues the same fixed rule. No further solve,
reference promotion, or native workflow timing win follows automatically.

Fourth-level preflight 413 passes at 144,106 nodes and 91,996 tets. It retains
2,101 stars and 90.000915% of the admissible balanced gain. Its largest trial
gain is 0.093758658%, above the existing 0.05% signal threshold. Minimum corner
quality is 0.100152291, above the unchanged 0.1 gate. Original domain and
parent-space checks, positive quadrature, periodicity, and actual reduced AMG
setup pass at 4.287 GiB peak.

Fresh coupled job 414 improves worst all-model disagreement to **1.415999636%**,
still failing 1%. The largest negative stiffness correction is 0.132695399%,
and Galerkin balance is within 1.601e-9. All independent field and quadrature
checks pass. Solve, audit, and assessment workers take 146.280, 31.676, and
38.656 seconds respectively. Solve peak RAM is 3.536 GiB. See
[the fourth-level evidence](../../handoff/macbook/split-patch-bulk90-level4-evidence-20261025.json).
One fifth-level preflight is selected from these fresh fields, with the same
fixed rule and all geometry, count, quality, resource, and accuracy gates.

Fifth-level preflight 415 is stopped during sparse-operator or AMG setup by
the unchanged process-tree RSS watchdog, at 8.008621 GiB. It had written the
mesh, prolongation, gains, and trial fields, and passed its earlier geometry
and trial assertions. It has no final admitted preflight result. No CG solve
starts. The failed result is preserved, not reclassified.

Resource-only job 416 loads exactly that saved mesh in a fresh process. It
passes reduced operator and AMG setup in 11.212 seconds at 5.040 GiB. The mesh
has 168,615 nodes and 109,059 tets. Its reduced operator has 499,914 rows and
39,464,528 nonzeros, with SHA256
`717f1d1ee1c9db7e6a04b9245a87b0c61a464f013202ee767ac94910e3b02d8e`.
There are three AMG levels and operator complexity 1.049191. CG is stopped
before it begins. This resource result does not admit the incomplete preflight.

The replay saves a full trial checkpoint before resource setup, then uses
`os.execv` to replace the worker process. The queue wrapper and 8 GiB watchdog
stay active, while patch construction and integration allocations are released.
Job 417 reruns the complete preflight and requires exact point, Tet10, and Tri6
array equality with failed 415 before setup. No new specimen samples or
geometry are substituted. The failure and isolated diagnosis are saved in
[the resource receipt](../../handoff/macbook/split-patch-level5-resource-20261025.json).
No further solve follows before every replay gate is read.

Complete replay 417 passes at 4.590 GiB peak under the same 8 GiB watchdog.
The original failed 415 remains failed. Point, Tet10, and Tri6 arrays match
415 exactly. The reduced-operator hash matches independent 416 exactly. All
trial, fixed-domain, parent-space, positive quadrature, facet, and periodic
checks pass. Minimum quality is 0.100032836, and the largest trial gain is
0.065557630%, above the existing 0.05% signal threshold. See
[the complete replay receipt](../../handoff/macbook/split-patch-level5-replay-20261025.json).
Fresh coupled job 418 is selected from this complete replay. Its six loads
and independent all-model assessment must pass before another decision.

All four cases and all 24 model tensors remain in scope. Schwarz retains its
measured comparison win against explicitly modified actual Microgen routes.
Gyroid and spinodoid still lack admitted actual comparators in the tested
routes. Lattice reference convergence and a native all-model accuracy pass
remain unresolved.

Source hashes and the four-case replay are recorded in
[the native evidence receipt](../../handoff/macbook/periodic-edge-normal-evidence-20261025.json).
Frozen run artifacts remain under `experiments/adaptive/runs`.
