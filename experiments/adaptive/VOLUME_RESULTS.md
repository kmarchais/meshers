# Independent lattice volume reference on G3

## Why the reference needs refinement

The frozen lattice volume is 0.2483110428, estimated from four scrambled
Sobol sets of 2^18 points each. Their absolute spread is 0.00015640. The
finest Meshers curved mesh has volume 0.2482327338, an apparent deficit of
about 0.000078 against that estimate. The reference spread is too large to
justify a geometry offset of that size.

## Nested independent audit, 19 October 2026

Job 258 extends the same four sequences to 2^24 points each, for 67,108,864
material classifications. Every original 2^18 prefix reproduces exactly. The
field is still the eight radius-0.12 capsules, joined in the frozen corner
order by Microgen's cubic smooth union with k=0.08. Another 4,096-point check
matches the explicit field against Microgen's shape evaluation.

The four final volume estimates are:

| Seed | Volume |
|---|---:|
| 0 | 0.2482439280 |
| 1 | 0.2482404113 |
| 2 | 0.2482352257 |
| 3 | 0.2482302189 |

Their mean is 0.2482374460, absolute spread is 0.0000137091, and estimated
standard error of the mean is 0.000003000. The final nested mean changes by
only 0.0000007004 from 2^23 points per seed. Both predeclared precision
screens pass: spread at most 0.00003 and mean-prefix shift at most 0.00001.
Peak RSS is 0.235 GiB, below the unchanged 8 GiB ceiling. Each replicate
has a 300-second cap and takes less than eight seconds.

These are empirical QMC stability measures, not a rigorous volume interval.
The maintained half-cube control uses a strict boundary comparison and is
checked against an unscrambled two-point net containing a boundary point.
The frozen driver's inclusive control remains unchanged; this checker fix
does not change lattice classifications or recorded volume estimates.

## Consequence for volume correction

| Saved mesh | Relative error at new mean | Range across replicate estimates |
|---|---:|---:|
| Meshers nonadaptive curved24 | -0.001410% | -0.004021% to +0.001501% |
| Meshers finest adaptive | -0.001898% | -0.004509% to +0.001013% |
| Meshers energy-sizing pilot | -0.001862% | -0.004473% to +0.001049% |
| Actual tuned Microgen-derived straight Tet10 | -0.458199% | -0.460798% to -0.455301% |

At the new mean, about 94% of the finest Meshers mesh's old apparent deficit
comes from the old reference estimate. Zero lies inside its empirical
replicate range. Do not add a Meshers volume offset based on the old mean.
This does not prove zero geometry error or explain stiffness discrepancies.
The Microgen control still passes the 1% volume gate, but its full stiffness
target remains unproven at 1.004354% worst retained-model disagreement.

No frozen result, admission threshold, or stiffness reference is overwritten.
The 24 lattice stiffness models remain retained. The failed energy-sizing
pilot remains failed. No continuum or local-stress certificate follows from
the closer volume estimate.

Records are in `runs/g3-lattice-volume-reference-audit-v1`. The reproducible
driver is `volume_reference_audit.py`.
