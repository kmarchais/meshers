# Periodic meshing investigation

The completed first pilot is pueue task 1204, retained in `runs/coarse-v2`.
It uses the packaged PR #5 build, native MMG CLI through Microgen's supported
public two-pass remesher, one thread and no timing warmup. Its five regression
checks passed. Both solid-cell controls recovered the analytic stiffness tensor
to relative error below 7e-16 and passed the periodic and equilibrium audits.

## First periodic Gyroid comparison

The material is a raw Gyroid band `[-0.55, 0.55]` in a unit cube. Geometry
admission requires both sampled boundary distances at most 0.002, volume error
plus empirical uncertainty at most 1%, and minimum tetrahedron quality 0.1.

| Workflow | Grid | Tetrahedra | Generation sample | Min q | Geometry / shape / periodic |
| --- | ---: | ---: | ---: | ---: | --- |
| Meshers, four passes | 16 | 26,008 | 1.103 s | 0.1604 | Pass / pass / pass |
| Meshers, one pass | 24 | 72,757 | 0.686 s | 0.1503 | Pass / pass / pass |
| Microgen periodic VTK + MMG | 24 | 73,049 | 5.766 s | 0.0294 | Fail / fail / pass |

The grid-24 native mesh has forward distance 0.001422, reverse distance 0.001257
and volume error plus uncertainty 0.284%. All six strain-controlled and six
stress-controlled solves pass their physical audits. Its audited homogenization
takes 28.996 s, including assembly and verification.

The grid-24 Microgen/MMG mesh is valid and periodic but has distances 0.008344
and 0.007982, volume error plus uncertainty 1.014%, and 248 cells below quality
0.1. Its 34.325 s homogenization result is diagnostic. The raw generation samples
differ by 8.4 times, but there is no qualified speedup at a common FEA target:
the MMG mesh is excluded, timings are single samples, and no converged reference
tensor has been established.

Native zero and one pass reject at grid 16 with geometry errors 0.002225 and
0.002472. The public MMG grid-16 output fails independent volume validity:
adjacent tetrahedra lie on the same side of a shared face. The input/output
archives and logs remain available; these failures are not repaired silently.

The native grid-16 four-pass tensor differs from the grid-24 one-pass tensor by
3.72% in the worst macroscopic strain-energy combination. The two grid-24
pipelines disagree by 1.45%, but that is not either mesh's FEA error. Their
geometry differs and neither tensor is a converged continuum reference.

## Frozen cap quality

The Microgen/MMG output retains all 3,072 original required cap triangles,
including their coordinates, to 12 decimal places. Its minimum normalized cap
triangle quality is 0.127919; the native grid-24 mesh's minimum cap quality is
0.354760.

For a fixed face of area `A` and sum of squared edges `S`, any attached
tetrahedron has MMG quality at most `4 sqrt(3) A/S`. Its optimal apex is above
the face centroid at height squared `2S/9`. This ceiling equals the normalized
triangle quality. The derivation and actual saved-face comparison are retained
in `cap-analysis.json` and `analyze_caps.py`.

These cap triangles do not force failure of the 0.1 gate. They do prevent the
tested fixed-boundary pipeline from matching the native mesh's minimum quality
0.1503, regardless of interior apex placement. A stronger MMG control needs a
better seed boundary or coupled periodic boundary changes. The present
measurement does not imply an intrinsic limit of MMG with different input.

## Refinement and first-pass distance

Task 1210 completed the zero/four-pass trial at grid 24 and one pass at grids
32 and 48. It also compared the public periodic MMG baseline with the
explicit first-pass-distance control at grid 24. Both MMG configurations use
`hausd=0.001` for the second pass, so the extra first-pass setting is the only
parameter difference between those two controls. Geometry admission stayed
0.002 and all other common gates stay unchanged.

Its six regression checks passed. Five complete meshes have all six saved
displacement fields, with independent stiffness reconstruction differences
below 5e-16, Hill–Mandel errors below 6e-12, force-balance errors below 1.1e-10
and periodic displacement jumps below 1e-15.

| Native configuration | Tetrahedra | Generation sample | Min q | Admission |
| --- | ---: | ---: | ---: | --- |
| Zero passes, grid 24 | No output | No timing claim | 0.0513 at rejection | 58 cells below 0.1 |
| Four passes, grid 24 | 72,757 | 2.130 s | 0.1622 | Geometry, shape and periodic pass |
| One pass, grid 32 | 171,679 | 1.475 s | 0.2737 | Geometry, shape and periodic pass |
| One pass, grid 48 | 543,496 | 4.431 s | 0.1520 | Geometry, shape and periodic pass |

The grid-32 tensor differs from grid 48 by 1.158% in the worst macroscopic
strain-energy combination, or 0.620% in the Frobenius norm. Grid 48 is still
an observed refinement reference, not a certified continuum solution. A further
level is required before interpreting a 2% reference target as a useful accuracy
measurement.

The first-pass-distance control reduces forward/reverse errors from
0.009233/0.008763 to 0.002568/0.002296. Both still fail the unchanged 0.002
boundary gate and the 0.1 shape gate. The modified output has 244 poor cells
and minimum quality 0.0312. Its generation sample costs 11.811 s, versus
7.852 s for its matched baseline. This is a geometry improvement to investigate,
not an admitted FEA or speed result.

## Next experiment

The next batch starts with an analytic solid-cell proof of single standard MMG
adaptation, then applies it to a valid native grid-32 Gyroid seed. This isolates
MMG's adaptation from the VTK seed and the public pipeline's two-pass treatment.
The native seed's required periodic caps stay constrained and every common
admission/physics check remains active. Total generation includes the seed cost.

A direct grid-64 reference follows with the same native profile, constitutive
law and solver. Its prepared worker timeout is 40 minutes. The result will
measure reference stability and compare the standard adaptation control at
the actual periodic stiffness target. No 2% or 1% periodic FEA qualification
is claimed before that comparison.
