# Mesh-generation comparison, 24 October 2026

Fresh M4 Mac measurements, one excluded warmup and five measured repetitions per setting. Meshers and the numerical libraries use one configured thread; VTK uses its sequential backend. The unit-cell side length is 1. Configuration order alternates between repetitions. No FE solve time is included. Ratios below are ratios of median generation times.

Tested packages: Meshers 0.1.0 native experimental build, Microgen 2.0.0b1, mmgpy 0.17.1, NumPy 2.5.4 and VTK 9.7.1. Both the Meshers library and MMG executable are native arm64. The manifests pin the actual source and executable hashes.

## What the comparison establishes

- **Schwarz, matched Tet10:** native fast generation is 45.31× faster, with 4.97× fewer tetrahedra, 4.82× fewer nodes and 5.94× lower sampled distance than the declared repaired actual control. All five reported quality quantiles improve. Four-pass native generation is 18.68× faster and slightly more faithful. Both arms pass every mesh gate in 5/5 runs.
- **Lattice, matched Tet10:** the geometry-only uniform25 recipe is 10.90× faster, with 3.25× fewer tetrahedra, 3.18× fewer nodes and 4.45× lower sampled distance. Four-pass native generation is 3.32× faster and has 7.44× lower sampled distance. Minimum, p1, p5 and median quality improve. The p95 quality is slightly lower, so this is not strict dominance of all reported quantiles. Both arms pass every mesh gate in 5/5 runs.
- **Gyroid:** the tested native Tet4 and Tet10 settings pass 5/5. Actual defaults fail distance, volume and quality in 5/5; stronger first-pass attempts exceed the fixed budget or fail in MMG. No admitted matched-quality speedup is claimed. Coarse native Tet10 has better minimum/low-tail quality but a lower median than the rejected default actual output.
- **Spinodoid:** native curved Tet10 passes 5/5. The tested actual first-pass controls fail distance, quality, facets/periodicity as recorded. Native adaptive Tet4 also fails distance. No admitted speedup is claimed. The actual reference-to-mesh distance alone is about 0.003506, already above 0.002; its larger forward-foot value is not interpreted as a certified Hausdorff distance.
- **Linear Tet4 remains a limitation:** Schwarz adaptive native generation is 6.65× faster with 3.78× fewer tetrahedra, but sampled distance is 3.3% larger and median quality is slightly lower. On lattice the actual Tet4 control passes while the tested native Tet4 controls miss fidelity or budget. Thus curved-Tet10 results do not establish a universal Tet4 advantage.

## Timing and admission

Generation includes fresh specimen/field preparation, native compilation where needed, VTK extraction, MMG calls, explicit repair/elevation and output conversion. Imports, archive writes, independent audits, measured observer instrumentation and all FE work are excluded. These are API generation times, not cold process startup or full validation workflows.

Every timed generator runs in a bounded worker with 8 GiB RSS, a 300-second timeout and 250,000 node/tet caps, through the shared one-worker queue. Mesh admission is sampled distance ≤ 0.002, absolute relative volume error ≤ 1%, minimum corner quality ≥ 0.1, positive sampled Jacobians, coherent facets and periodic nodes/caps. The same box-constrained observer is used for both arms.

Internal Meshers geometry_tolerance=1.0 follows the documented recipes; the independently checked output gate remains 0.002. Native adaptation tolerance is a separate parameter. This prevents an internal parameter name from being mistaken for the measured output error.

Corner quality is the common normalized tetrahedron metric. Higher is better. Curved Jacobian positivity is checked separately and is sampled, not globally certified. Surface distances are sampled; forward feasible-foot distances are conservative observations, not certified closest-point or Hausdorff bounds.

Schwarz volume is analytically 0.5. The other volumes use the frozen four-scramble estimates. Small differences below the sampling spread are unresolved, especially the spinodoid volume ranking. No stiffness reference, FE candidate or continuum certificate is promoted.

## Actual routes and native representations

All actual meshes begin with genuine Microgen extraction, never a Meshers seed. Defaults remain separate from explicit variants. Schwarz uses VTK30, first MMG nofem hausd0.0005/hmax0.1/hgrad1.3, omits the second call, and charges an interior repair. Lattice uses VTK51, first hausd0.0015/hmax0.06/hgrad1.3, second hausd0.001/hmax0.06/hgrad1.3, and charges an interior repair. Tet10 elevation is explicit and preserves the affine actual parent geometry.

Gyroid tables use actual defaults at resolution20 because stronger budget controls did not produce a qualifying output. Spinodoid retains the FFT16/52-mode realization while sampling at36 and uses a declared first-pass-only control. The ordinary two-pass spinodoid baseline also fails periodicity; tighter second-pass controls exceed the intermediate cap. Earlier stronger preclipped seed studies also failed the required gates and are not pooled into fresh timing statistics.

Native spinodoid compilation preserves the frozen Fourier coefficients. Native lattice uses the previously validated symmetric-min expression and positive-void clamp; the material domain and original smooth-union order are unchanged. Exact arguments and source/binary hashes are in each manifest and result.

## Full measured settings

Each time cell is median [minimum, maximum] in seconds. All quality columns are medians across the five repeated per-mesh measurements. Counts and geometry metrics are likewise medians; the machine-readable report retains ranges. Rejected rows remain descriptive and receive no admitted speed ratio.

### Schwarz P solid

| Setting | Gate passes | Tets / nodes | Generation s | Sampled distance | Abs. volume error % | q min | q p1 | q p5 | q median | q p95 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Actual defaults | 0/5 | 10,829 / 3,151 | 1.15352 [1.13567, 1.16181] | 0.00846258 | 0.0551498 | 0.0548389 | 0.135826 | 0.410837 | 0.679021 | 0.905416 |
| Actual tuned, declared additions | 5/5 | 136,511 / 29,192 | 3.7653 [3.71672, 3.78763] | 0.000961515 | 0.0111356 | 0.100001 | 0.460605 | 0.561414 | 0.747047 | 0.934503 |
| Actual tuned + straight Tet10 | 5/5 | 136,511 / 207,099 | 4.13033 [4.08543, 4.18417] | 0.000961515 | 0.0111356 | 0.100001 | 0.460605 | 0.561414 | 0.747047 | 0.934503 |
| Meshers Tet4, uniform24 | 5/5 | 90,792 / 19,004 | 0.223513 [0.22261, 0.224394] | 0.00189406 | 0.00857178 | 0.399714 | 0.544233 | 0.623103 | 0.929516 | 0.965765 |
| Meshers Tet4, adaptive | 5/5 | 36,132 / 8,741 | 0.566459 [0.559208, 0.574788] | 0.00099339 | 0.00192656 | 0.410699 | 0.527616 | 0.578604 | 0.743023 | 0.935526 |
| Meshers curved, fast 1 pass | 5/5 | 27,492 / 42,928 | 0.0911598 [0.0900254, 0.0927175] | 0.000161891 | 5.47504e-05 | 0.449995 | 0.557311 | 0.615507 | 0.916233 | 0.974491 |
| Meshers curved, 4 passes | 5/5 | 27,492 / 42,928 | 0.221097 [0.218018, 0.224053] | 0.000157498 | 3.99323e-05 | 0.487272 | 0.562917 | 0.620915 | 0.914944 | 0.974142 |

[Raw repetitions](/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-mesh-only-schwarz-repeated-v1/summary.json)

- Actual defaults: distance, quality.

### Smooth-union lattice

| Setting | Gate passes | Tets / nodes | Generation s | Sampled distance | Abs. volume error % | q min | q p1 | q p5 | q median | q p95 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Actual tuned, declared additions | 5/5 | 158,430 / 33,953 | 8.75784 [8.71674, 8.82532] | 0.00170124 | 0.490534 | 0.100001 | 0.489743 | 0.594626 | 0.803707 | 0.946063 |
| Actual tuned + straight Tet10 | 5/5 | 158,430 / 238,729 | 9.19878 [9.13129, 9.23043] | 0.00170124 | 0.490534 | 0.100001 | 0.489743 | 0.594626 | 0.803707 | 0.946063 |
| Meshers Tet4, uniform32 | 0/5 | 121,420 / 27,503 | 2.29609 [2.25589, 2.3092] | 0.00462883 | 0.317816 | 0.27688 | 0.488194 | 0.57494 | 0.920406 | 0.935995 |
| Meshers uniform25 curved, 1 pass | 5/5 | 48,696 / 75,113 | 0.843539 [0.8264, 0.854344] | 0.000381982 | 0.0315376 | 0.322484 | 0.492254 | 0.640447 | 0.929516 | 0.944808 |
| Meshers uniform25 curved, 4 passes | 5/5 | 48,696 / 75,113 | 2.76691 [2.72238, 2.79597] | 0.000228592 | 0.0317062 | 0.335695 | 0.503628 | 0.646851 | 0.929516 | 0.945353 |

[Raw repetitions](/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-mesh-only-lattice-repeated-v2/summary.json)

- Meshers Tet4, uniform32: distance.

### Gyroid sheet

| Setting | Gate passes | Tets / nodes | Generation s | Sampled distance | Abs. volume error % | q min | q p1 | q p5 | q median | q p95 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Actual defaults | 0/5 | 12,250 / 3,533 | 1.82458 [1.77391, 1.84626] | 0.00992548 | 2.24893 | 0.000546158 | 0.0542905 | 0.203073 | 0.715926 | 0.927329 |
| Actual defaults + straight Tet10 | 0/5 | 12,250 / 21,865 | 1.84778 [1.78455, 1.87674] | 0.00992548 | 2.24893 | 0.000546158 | 0.0542905 | 0.203073 | 0.715926 | 0.927329 |
| Meshers Tet4, uniform24 | 5/5 | 72,757 / 17,647 | 0.400699 [0.395271, 0.449157] | 0.00158399 | 0.189121 | 0.150254 | 0.513315 | 0.594989 | 0.882537 | 0.950839 |
| Meshers Tet4, adaptive | 5/5 | 103,493 / 23,717 | 1.8387 [1.79557, 1.85882] | 0.0010229 | 0.117887 | 0.249383 | 0.524888 | 0.574596 | 0.704416 | 0.909048 |
| Meshers curved, fast 1 pass | 5/5 | 6,008 / 11,857 | 0.0683243 [0.0666756, 0.0692565] | 0.00118934 | 0.0263276 | 0.14469 | 0.337904 | 0.457465 | 0.673738 | 0.924575 |
| Meshers curved, 4 passes | 5/5 | 6,008 / 11,857 | 0.201651 [0.195767, 0.204093] | 0.00118934 | 0.026778 | 0.160483 | 0.348674 | 0.467914 | 0.675667 | 0.920407 |

[Raw repetitions](/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-mesh-only-gyroid-repeated-v1/summary.json)

- Actual defaults: distance, facets, quality, volume.
- Actual defaults + straight Tet10: distance, facets, quality, volume.

### Fixed spinodoid

| Setting | Gate passes | Tets / nodes | Generation s | Sampled distance | Abs. volume error % | q min | q p1 | q p5 | q median | q p95 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Actual tuned, declared additions | 0/5 | 143,585 / 31,787 | 6.19809 [6.1424, 6.2135] | 0.0399574 | 0.0151691 | 2.32785e-07 | 0.177697 | 0.503471 | 0.723602 | 0.924259 |
| Actual tuned + straight Tet10 | 0/5 | 143,585 / 221,826 | 6.58675 [6.55105, 6.62613] | 0.0399574 | 0.0151691 | 2.32785e-07 | 0.177697 | 0.503471 | 0.723602 | 0.924259 |
| Meshers Tet4, adaptive | 0/5 | 131,598 / 30,499 | 4.20584 [4.18537, 4.2215] | 0.00483425 | 0.0378576 | 0.152884 | 0.520466 | 0.581486 | 0.76689 | 0.944459 |
| Meshers curved, fast 1 pass | 5/5 | 45,994 / 73,951 | 1.8312 [1.81822, 1.84126] | 0.00133298 | 0.0258625 | 0.268581 | 0.526352 | 0.592462 | 0.785077 | 0.939914 |
| Meshers curved, 4 passes | 5/5 | 45,994 / 73,951 | 3.57705 [3.56879, 3.57785] | 0.00133298 | 0.0258731 | 0.333453 | 0.53458 | 0.59737 | 0.780426 | 0.935125 |

[Raw repetitions](/Users/kmarchais/dev/meshers/experiments/adaptive/runs/macbook-mesh-only-spinodoid-repeated-v1/summary.json)

- Actual tuned, declared additions: distance, periodicity, quality.
- Actual tuned + straight Tet10: distance, periodicity, quality.
- Meshers Tet4, adaptive: distance.

## Verification and limits

Independent job114 recomputes mesh metrics from all 138 saved meshes, including warmups, and checks 76 input hashes. Counts, quality quantiles, geometry metrics, sampled Jacobians, facet consistency, periodicity and straight-elevation coordinates reproduce. Every node is used. Three fixed-recipe/Fourier benchmark controls pass.

The timed runs and independent replay finished before the later clamshell sleep. Pueue UTC elapsed minus summed monotonic worker durations is 1.29–1.75 seconds per final experiment, consistent with parent setup and bookkeeping. This supports clock consistency, not CPU isolation.

Pilots and all failed attempts remain frozen. Early incorrect known-recipe settings and the straight-elevation adapter error are identified as benchmark reproduction/adapter failures, not evidence of a mesher defect. They are excluded from the final timing statistics.

This compares selected reasonable settings on four fixed periodic specimens and pinned versions. It is not an exhaustive optimizer search, a general result for all geometries, or a stiffness-equivalence comparison. The four-case project objective remains unachieved.

## Evidence

- [Verified machine-readable report](/Users/kmarchais/dev/meshers/handoff/macbook/mesh-comparison-20261024.json)
- [Timing check](/Users/kmarchais/dev/meshers/handoff/macbook/mesh-comparison-timing-20261024.json)
- [Protocol](/Users/kmarchais/dev/meshers/experiments/adaptive/MESH_COMPARISON.md)
- [Generator](/Users/kmarchais/dev/meshers/experiments/adaptive/mesh_only_comparison.py)
- [Independent replay](/Users/kmarchais/dev/meshers/experiments/adaptive/mesh_comparison_report.py)
