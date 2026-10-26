# Verified Meshers versus actual Microgen+MMG, 26 October 2026

Follow-up: Schwarz D now passes after correcting material-trace matching. The original linear-thickness Tet4 case and a new calibrated linear-density case also admit with geometry adaptation. See [the follow-up results](PERIODIC_GRADING_RESULTS.md). The eight-case tables below retain the historical baseline binary and its original failures.

Meshers has admitted meshes for seven of these eight specimens with curved Tet10 and six with Tet4. The tested modified actual routes admit Schwarz P and the lattice. Default MMG on the fixed study inputs admits none. The additional public-sampling32 controls on both frozen spinodoids also fail. These are bounded results for named routes, not universal geometry support.

The strongest generation results are Schwarz P curved Tet10 and the lattice Uniform37 curved Tet10. The latter uses almost equal element counts. Meshers still has periodic-matching failures on Schwarz D, a Tet4 distance rejection on linear grading, quality-distribution tradeoffs and unresolved lattice stiffness. No native meshing code changed during this comparison.

## Coverage under the unchanged gates

The user selected linear density/thickness grading as the primary added
graded-TPMS benchmark. Periodic modulation remains a secondary control.
The existing linear case varies the level-set thickness parameter. A
prescribed linear density ramp needs a separate calibrated, frozen specimen.
This priority change leaves the measured protocol and results unchanged.

Pass means five measured successes after an excluded warmup. The linear-grade case is bounded and nonperiodic for both arms. Every other case requires full periodic nodes and cap facets. Schwarz D failures repeat in three fresh attempts per route.

| Fixed specimen | Meshers Tet4 | Meshers curved Tet10 | Actual Microgen+MMG tested routes |
|---|---|---|---|
| Gyroid sheet | Pass | Pass | Default fails quality, volume, facets and reverse distance; modified MMG fails |
| Schwarz P solid | Pass, metric tradeoffs | Pass, strong generation result | Modified route passes; default minimum quality fails |
| Schwarz D solid | Fails periodic nodes | Fails periodic nodes | Both seeds fail public periodic-input check |
| Fixed spinodoid | Pass | Pass | Default16, public sampling32 and modified route fail |
| Anisotropic spinodoid | Pass | Pass | Default16, public sampling32 and modified route fail |
| Smooth-union lattice | Pass, metric tradeoffs | Compact and Uniform37 pass | Modified route passes; default fails |
| Gyroid, periodic grading | Pass | Pass | Default fails; modified route reaches 300 s timeout |
| Gyroid, linear grading | Distance 0.002002836 exceeds 0.002 | Pass | Default fails; modified first MMG output has 277,505 tets, above cap |

## Admitted mesh-generation measurements

Times exclude imports, common fixture loading, archival writes, audits and measured MMG observer I/O. They include fresh arm-specific preparation, native compilation, extraction, MMG retries, explicit repair or straight elevation, and conversion to common arrays. They are generation times, not FE workflow timings. All measurements use the same M4 and one configured numerical thread in the shared heavy queue.

Quality q is the corner-tet mean ratio, from 0 to 1. Higher is better. Geometry distance is a finite two-sided sampled observation in unit-cell length units. Volume is absolute relative error in percent against the frozen reference. These distances are not certified Hausdorff bounds. Small volume differences below the reference scramble spread do not establish a fidelity ranking.

| Specimen | Route | Measured passes | Generation s, median [range] | Tets | Nodes | q min / p01 / median | Sampled distance | Abs. volume % |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Gyroid sheet | Meshers Tet4 | 5/5 | 0.397108 [0.395791, 0.40603] | 72,757 | 17,647 | 0.150254 / 0.513315 / 0.882537 | 0.001583989 | 0.18912144 |
| Gyroid sheet | Meshers curved Tet10 | 5/5 | 0.0663062 [0.0651256, 0.0679088] | 6,008 | 11,857 | 0.144690 / 0.337904 / 0.673738 | 0.001189342 | 0.02632761 |
| Schwarz P solid | Actual, modified Tet4 | 5/5 | 3.77608 [3.75445, 3.80937] | 136,511 | 29,192 | 0.100001 / 0.460605 / 0.747047 | 0.000961515 | 0.01113555 |
| Schwarz P solid | Meshers Tet4 | 5/5 | 0.590813 [0.58546, 0.595512] | 35,994 | 8,701 | 0.437911 / 0.528044 / 0.745663 | 0.000993390 | 0.00144543 |
| Schwarz P solid | Meshers curved Tet10 | 5/5 | 0.090221 [0.0898842, 0.0903981] | 27,492 | 42,928 | 0.449995 / 0.557311 / 0.916233 | 0.000146600 | 0.00005475 |
| Schwarz P solid | Actual, modified straight Tet10 | 5/5 | 4.11032 [4.09913, 4.14003] | 136,511 | 207,099 | 0.100001 / 0.460605 / 0.747047 | 0.000961515 | 0.01113555 |
| Fixed spinodoid | Meshers Tet4 | 5/5 | 3.28449 [3.28206, 3.30636] | 72,378 | 16,702 | 0.233170 / 0.527518 / 0.783152 | 0.001864496 | 0.04505673 |
| Fixed spinodoid | Meshers curved Tet10 | 5/5 | 1.90355 [1.89621, 1.91472] | 45,923 | 73,811 | 0.301282 / 0.528144 / 0.786726 | 0.000303984 | 0.02589182 |
| Anisotropic spinodoid | Meshers Tet4 | 5/5 | 2.54285 [2.52512, 2.6179] | 61,868 | 14,005 | 0.121849 / 0.526182 / 0.780050 | 0.001835045 | 0.05509647 |
| Anisotropic spinodoid | Meshers curved Tet10 | 5/5 | 1.40567 [1.39631, 1.40895] | 40,351 | 63,919 | 0.138481 / 0.529079 / 0.785397 | 0.001145505 | 0.00421173 |
| Smooth-union lattice | Actual, modified Tet4 | 5/5 | 8.79075 [8.75204, 8.87461] | 158,430 | 33,953 | 0.100001 / 0.489743 / 0.803707 | 0.001701239 | 0.49053370 |
| Smooth-union lattice | Meshers Tet4 | 5/5 | 2.04208 [2.02072, 2.05651] | 58,847 | 13,698 | 0.385465 / 0.528418 / 0.777657 | 0.001885109 | 0.47679962 |
| Smooth-union lattice | Meshers curved Tet10 | 5/5 | 1.32501 [1.31948, 1.33179] | 41,340 | 65,721 | 0.365931 / 0.525756 / 0.768100 | 0.000278662 | 0.03124175 |
| Smooth-union lattice | Actual, modified straight Tet10 | 10/10 | 9.19781 [9.17134, 9.27162] | 158,430 | 238,729 | 0.100001 / 0.489743 / 0.803707 | 0.001701239 | 0.49053370 |
| Smooth-union lattice | Meshers Uniform37 curved Tet10 | 5/5 | 1.88953 [1.85354, 1.90097] | 157,521 | 231,026 | 0.339372 / 0.577287 / 0.929516 | 0.000178307 | 0.03147631 |
| Gyroid, periodic grading | Meshers Tet4 | 5/5 | 0.477695 [0.471724, 0.479689] | 78,602 | 19,595 | 0.122621 / 0.472262 / 0.861128 | 0.001482919 | 0.19671321 |
| Gyroid, periodic grading | Meshers curved Tet10 | 5/5 | 0.208152 [0.205599, 0.208258] | 25,744 | 44,857 | 0.321158 / 0.471548 / 0.795067 | 0.000389127 | 0.00169024 |
| Gyroid, linear grading | Meshers curved Tet10 | 5/5 | 0.20595 [0.205841, 0.206806] | 25,624 | 44,841 | 0.159092 / 0.455714 / 0.792993 | 0.000277250 | 0.01083434 |

The lattice actual straight Tet10 row pools ten measured samples from two repeated comparisons. Each ratio below uses only its own five interleaved pairs. Warmups and pilot timings are excluded from the repeated pool. Failed routes never support a successful speed ratio. Their partial measurements and failure logs remain in the numeric ledger.

## Same-order comparisons where both routes admit

| Specimen and native route | Actual time / native time | Actual tets / native tets | Remaining tradeoff |
|---|---:|---:|---|
| Schwarz P solid, Meshers Tet4 | 6.391x | 3.7926x | Native median q 0.745663 versus 0.747047; distance 0.000993390 versus 0.000961515 is higher. |
| Schwarz P solid, Meshers curved Tet10 | 45.558x | 4.9655x | Native minimum, p01 and median q, distance and volume improve. Actual q max is slightly higher; no full quality-CDF dominance claim. |
| Smooth-union lattice, Meshers Tet4 | 4.305x | 2.6922x | Native median q and distance are worse despite much fewer tets. |
| Smooth-union lattice, Meshers curved Tet10 | 6.952x | 3.8324x | Compact native median q and upper quantiles are worse. Its stiffness disagreement fails 1%. |
| Smooth-union lattice, Meshers Uniform37 curved Tet10 | 4.857x | 1.0058x | Almost equal tets: 157,521 versus 158,430. Native q99 0.968982 versus 0.973730 is lower; stiffness still fails 1%. |

Uniform37 versus the actual straight Tet10 route has median generation times 1.889529 and 9.176655 s in its five paired measurements. It has 0.574% fewer tets and 3.227% fewer nodes. Native minimum/median q are 0.339372/0.929516 versus 0.100001/0.803707. Sampled distance is 0.000178307 versus 0.001701239. The full seven quality quantiles remain in the ledger; improved minimum and median do not imply a better entire distribution.

The same Uniform37 connectivity with straight Tet4 has sampled distance 0.003424939 and remains rejected. Its quality and timing cannot establish an admitted linear-mesh win.

## Method checks and scope

- The original four sample/specimen files retain their historical hashes. The smooth union retains its original cubic blend and order. Native positive-void clamping preserves the material domain. The added TPMS and frozen Fourier fields have independent 4,096-probe identity controls. The original TPMS and ordered lattice definitions retain their historical identity checks and fixed input hashes.
- Every actual route starts from genuine Microgen output. No actual comparator uses a Meshers seed. The installed versions are Microgen 2.0.0b1 and mmgpy 0.17.1. Modified routes explicitly name first-pass overrides, second-pass omission and any bounded SciPy repair; those costs are charged.
- "Default" means public MMG defaults on fixed study inputs. TPMS sampling is 20 and generic lattice sampling is 50. Both spinodoids require their frozen FFT16 field; sampling16 and a separate public-sampling32 control are recorded. A fresh FFT32 realization would change the specimen and was not substituted.
- Every intermediate actual seed/MMG mesh and every final mesh must fit 250,000 tets and nodes. Workers have an 8 GiB polled process-tree RSS ceiling and 300 s deadline. The native binary remains `2005c3e01892098a9d734567220a56ed3ea5d748b20a4f96754fbb55168060b1`.
- Admission requires sampled distance <=0.002, absolute relative volume <=1%, corner q>=0.1, positive sampled curved Jacobians, coherent oriented facets, and full periodic matching where required. Curved checks include a positive 125-point rule, nodes and centroid. Quadratic reverse geometry uses subdivisions 8 and 16. Straight Tet10 reduction is allowed only after exact affine midpoint and full connectivity verification.
- The old corner-face audit could miss different Tet10 midpoint IDs on a shared edge. A deliberate cracked-mesh control exposes that gap; the new audit rejects it. Analytic regular-tet quality/volume, straight Tet10, plane distance and sphere projection controls pass. The current focused control suite passes 12 tests.
- Rejected cheap gates trigger a fail-fast observer. Such rows may have only a reverse-distance observation, not a complete two-sided score. The independent linear-grade witness minimization gives 0.002002834963, still above 0.002, and does not reclassify the original rejection.
- Added grading changes the Gyroid level-set half-width to `0.55*(1+0.25*cos(2*pi*z))` or `0.55*(1+0.5*z)` on z in [-0.5,0.5]. This grades the level-set parameter, not a claimed constant physical wall thickness or prescribed density. Linear grading is nonperiodic in z, so both arms use bounded nonperiodic meshing for that case.
- Schwarz D has exact reference volume 0.5. The other added cases use four fixed Sobol scrambles of 2^20 volume probes each. Spread is empirical, not a confidence interval. The original four volume references and all mesh gates remain unchanged.

## Lattice stiffness remains unresolved

The saved six-load tensors were independently re-evaluated against all 24 unchanged retained model tensors. No new FE solve or reference was created. The statistic is worst relative mixed-load strain-energy disagreement, not certified continuum or local-stress error.

| Saved route | Tets | Nodes | Worst disagreement across all 24 | 1% target |
|---|---:|---:|---:|---|
| compact-native | 41,340 | 65,721 | 2.412306220% | Fail |
| uniform37-native | 157,521 | 231,026 | 2.085389679% | Fail |
| experimental-saved-relocation | 162,189 | 249,824 | 1.213246982% | Fail |
| modified-actual-with-star-refinement | 165,722 | 250,000 | 0.994670612% | Pass |

The admitted actual stiffness route includes explicit repair, straight elevation and star refinement to 165,722 tets and 250,000 nodes. It differs from the 158,430-tet mesh-only timing route. The mesh-only speed ratios do not apply to the complete stiffness workflow. Uniform37 is byte-identical to its earlier admitted preflight archive. The original four curved default outputs also reproduce their prior field meshes exactly. The seven Uniform37 comparisons against the compact-default archive in the ledger are expected differences, not a failed Uniform37 replay.

No matched-stiffness claim is made for the added specimens. The full four-case objective remains incomplete.

## Evidence and next work

The final audit verifies 158 attempts across 20 runs and independently replays 30 distinct meshes. It verifies measured source and artifact hashes, configs against the frozen selections, all successful worker limits, counts/order, quality, volume, geometry, quadratic conformity and required periodic caps. Every selected admitted recipe has five fresh measured samples after an excluded warmup.

- [Numeric findings and verification](../../handoff/macbook/comparison-findings-20261026.json)
- [Frozen eight-case protocol](../../handoff/macbook/comparison-protocol-20261026.json)
- [Near-equal-count lattice selection](../../handoff/macbook/comparison-matched-count-selection-20261026.json)
- [Spinodoid sampling control](../../handoff/macbook/comparison-spin-sampling-selection-20261026.json)
- [Fixture hashes and controls](../../handoff/macbook/comparison-fixture-receipt-20261026.json)
- [Decision trail](../../.audit/macbook-continuation.tsv)

Failure work should start with Schwarz D periodic-node matching. Then address the linear-grade Tet4 distance miss. Lattice work still needs better quality allocation and a candidate that passes every retained stiffness model; another nearby sizing sweep has no support from these results.
