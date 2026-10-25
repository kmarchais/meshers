# Periodic image error-query correction, 25 October 2026

The focused regression is red before correction and green afterward. Full
native/Python suites and independent four-study replay pass. The lattice uses
fewer elements and improves its reported quality quantiles at unchanged sampled
distance. No new stiffness reference or comparative timing claim is promoted.

## Seed and label diagnosis

The physical six-plane tapered-prism seed exposes two separate stages:

| Diagnostic | Construction only | One independent nonperiodic cleanup round |
|---|---:|---:|
| Tetrahedra / nodes | 46084 /10255 | 45849 /10126 |
| Minimum quality | 0.00922081 | 0.09790673 |
| Low / high end nodes | 139 /139 | 130 /122 |
| Cap triangles per end | 218 /218, matching | 209 /201, incompatible |
| Sampled distance | 1.24e-16 | 1.24e-16 |

Neither is admitted as a periodic candidate. The declared diagnostic turns off
generator periodicity and uses a zero generation-quality floor to expose the
seed, while retaining the original external quality0.1 and periodic admission
requirements. Jobs277–278 and their preflight manifests record this distinction.
The actual original raw recipe272 remains rejected at geometry error0.0106425.

The construction-only seed has31 pairs with different raw constraint masks.
Applying its exact six-plane end-trace permutation `[1,0,3,2,4,5]` reconciles all
of them. Independent cleanup destroys cap correspondence; merely turning off
the periodic field guard cannot be a correct native fix.

## Fixed-seed private-core control

An explicit ignored research test loads that saved seed through integer float
bits, preserving every coordinate exactly. It uses the actual torus adapter
and original final core gate with four adaptation rounds,one polish round,
quality0.1,distance0.002 and250000node/tet limits. The original scalar field
and its gradient agree with the native compiled evaluator on1028 probes,
including nonsmooth ties. This bypasses no production guard: it is an isolated
private-core experiment, not a public generator workflow or fresh timing.

Job 279 produces5797tets /1662nodes with exact geometry and matching53 node pairs
and77 cap triangles per end, but ten cells miss quality0.1; qmin is0.06360266.
The real core final gate rejects them. Its adaptation error statistic reaches
0.0152681 despite the final exact planar boundary, prompting a query audit.

## Proven query bug and narrow correction

The adapter's error enforcement forms local minimal-image simplex coordinates.
A local point can lie outside the fundamental box. The old distance query
passed that point directly to the field and gradient. A field need only match
its requested box traces; its unrestricted continuation need not represent the
periodic extension, especially when the physical extension has a normal jump.

For `x + 0.3*abs(z) - 0.35`, the point whose local z is-0.51 can be exactly on
the periodically repeated wall at physical z0.49. The old query assigns it a
nonzero residual. Ordinary regression 280 fails before correction. All original
core sources used by it were copied before the native edit, with a supplemental
receipt; no earlier manifest is claimed to have already captured those sources.

Distance queries now wrap only coordinates outside enabled periodic bounds.
Points already inside and all nonperiodic coordinates retain their original
values. Field and one-sided gradient are evaluated at the correct physical
image; no normal is averaged or smoothed. Regression 281 passes, including
locked-coordinate and unchanged nonperiodic controls.

On exactly the same saved seed and options, job 282 removes all111 spurious
error splits; its maximum error is1.06e-16. Counts fall to5547tets /1551nodes.
Exact geometry and periodic caps remain, but ten cells still fail quality0.1
with qmin0.06200199. The real final gate continues to reject them. Do not claim
a quality improvement or an admitted public tapered-prism route from this.

## Current verification boundary

Core and binding all-target/all-feature Clippy pass. Native release suite 283
passes 95 tests with the explicit saved-seed research test ignored by default.
Build 284 produces native arm64 image 30506b13. Python/fixture/tools suite 286
passes 206 tests with one optional skip. Public tapered-prism control 285 still
rejects sampled geometry error0.0106425; the query fix does not close its
feature-preservation failure.

Fresh controls 287–290 admit all four study specimens. Independent replay 292
proves Gyroid,Schwarz P and spinodoid exactly unchanged after canonical
relabeling, while recording the intentionally changed adaptive lattice:

| Lattice Tet4 metric | Before | Corrected periodic queries |
|---|---:|---:|
| Tetrahedra / nodes | 218239 /48077 | 213800 /46131 |
| Error-driven splits | 1960 | 24 |
| Sampled distance | 0.001768700533 | 0.001768700533 |
| Minimum quality | 0.171901839 | 0.179051949 |
| Quality p1 /p5 | 0.377941 /0.544803 | 0.517194 /0.577340 |
| Quality median /p95 | 0.730036 /0.925150 | 0.735863 /0.926112 |
| Nominal absolute volume error | 0.4291234% | 0.4291763% |

The nominal volume error worsens slightly, although both are inside1%; no
all-metric dominance is claimed. Positive Jacobians,coherent facets,periodic
nodes/caps,counts and fixed specimen/sample hashes pass independent replay.
The positive periodic affine control 291 is exactly unchanged from its admitted
pre-fix counterpart. These are admission controls, not repeated timings.

Research audit 293 independently reads the exact float-bit seed and saved
outputs, verifies that tet_adapt.rs is the only changed input source between
the two fixed-seed runs, reproduces the 111-to-0 error-split correction, and retains
both ten-cell quality failures and the public tapered-prism rejection. All ten
bad cells touch an affine curve with at least two active physical planes. This
is a useful next diagnostic; it is not a demonstrated quality fix.

See [the four-study audit](../../handoff/macbook/periodic-distance-study-audit-20261025.json),
[the research audit](../../handoff/macbook/periodic-distance-research-audit-20261025.json)
and [the native receipt](../../handoff/macbook/periodic-distance-native-candidate-20261025.json).
The saved-seed control remains supporting correctness evidence, never a
comparative workflow timing or a replacement for public admission.

The production affine periodic eligibility rule remains unchanged. Plane
permutations and the remaining quality obstruction are not fixed by this
query correction. Preserve all 24 lattice models and the original four-case
accuracy,geometry,quality,periodicity and same-host timing requirements.
