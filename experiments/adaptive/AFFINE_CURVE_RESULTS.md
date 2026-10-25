# Affine curve mobility, 25 October 2026

The periodic prism investigation exposed a constraint classification bug.
The adapter froze every vertex touching two affine planes. Two independent
planes leave a movable line. A third independent constraint fixes a point.
The worst cells in the previous saved-seed replay touched short feature edges
whose endpoints were both frozen.

Job 294 fails the focused line-mobility regression before the correction.
Job 295 passes afterward, including the cap endpoint and joint projection
checks. The change uses the rank of active affine normals and held box axes.
Redundant planes and box locks do not add rank. Scalar and nonlinear
constraint behavior is unchanged. Link, geometry and quality guards remain.
The red source was archived after the failed run and before the native edit;
the supplemental source receipt does not replace the original run manifest.

## Exact saved-seed replay

Job 296 reuses the exact floating-point seed, options and field from job 282.
It performs four adaptation rounds and one polish pass, with q >= 0.1,
distance <= 0.002, absolute volume error <= 1%, and 250,000 node/tet caps.

| Metric | Before, run v2 | Rank-aware mobility, run v3 |
| --- | ---: | ---: |
| Tetrahedra | 5,547 | 3,770 |
| Nodes | 1,551 | 838 |
| Minimum quality | 0.062001989 | 0.304656138 |
| Cells below q = 0.1 | 10 | 0 |
| False error splits | 0 | 0 |
| Native final gate | Reject quality | Pass |
| Independent geometry gates | Reject quality | Pass |

The recovered mesh has 42 periodic node pairs and 61 matching cap triangles
per end, positive determinants and coherent oriented facets. Sampled distance
is below 2e-16 and relative volume error is at roundoff. Probe values and
gradients match the compiled fixture at 1,028 points. The physical geometry
and the 33.4 degree periodic wall-normal jump are unchanged.

This is a private-core saved-seed recovery. It is not a public generator
workflow or a timing comparison. The public scalar prism recipe remains
rejected. The public periodic affine path still requires each plane to be
invariant along the enabled periodic axes. Supporting this prism requires a
justified complete-field trace correspondence and a verified periodic seed;
relaxing individual named-field periodic admission is not authorized by this
result. No FE solve, stiffness claim or reference promotion follows.

The four fixed study cases and all 24 lattice model tensors remain in scope.
Their earlier FE results apply to their recorded mesh and native image.

## Final validation

The complete native release suite, job 297, passes 97 tests with one explicit
research test ignored. Core and binding Clippy pass with warnings denied.
Clippy requested an equivalent iterator expression; final build 300 contains
that correction. Job 309 reproduces the recovered mesh byte for byte from
the final source. Independent replay 310 verifies its archived source hashes,
identical input seed, geometry, quality and periodic correspondence.

Final arm64 image SHA256 is
`99817e40c74341717ea58cc1bc5a54cacb5a1305e64b8b38763e5e5963492e06`.
Python suite 301 passes 206 tests with one optional skip. Fresh controls
302, 303, 305 and 306 admit Gyroid, Schwarz, spinodoid and lattice.
Independent replay 308 verifies all four exactly match the previous native
image's canonical arrays and preserve the fixed sample/specimen hashes.
These are single admission controls, not repeated timing comparisons.

Spinodoid launch 304 omitted the recipe-selection flag and stopped before
generation. Its failed directory is retained. Run 305 uses the previously
recorded `--strengthened --select meshers-tet4-historical32` configuration;
the audit explicitly identifies the replacement directory.

The public periodic affine column, job 307, remains admitted with 2,614 tets
and 604 nodes, compared with 2,663 and 655. Minimum quality rises from
0.302251 to 0.321550, with roundoff geometry/volume error and matching caps.
No all-metric or comparative speed claim follows from this control.

Receipts are in `handoff/macbook/affine-curve-native-candidate-20261025.json`,
`affine-curve-study-audit-20261025.json`,
`affine-curve-research-final-audit-20261025.json` and
`affine-curve-local-artifacts-20261025.json`.
