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

## Next experiment

The next batch tests zero and four native passes at grid 24 and one pass at
grids 32 and 48. It also compares the public periodic MMG baseline with the
explicit first-pass-distance control at grid 24. Both MMG configurations use
`hausd=0.001` for the second pass, so the extra first-pass setting is the only
parameter difference between those two controls. Geometry admission stays
0.002 and all other common gates stay unchanged.

Saved displacement fields allow an independent reconstruction of the full
stiffness tensor, periodic force equilibrium and the Hill–Mandel energy matrix.
The refinement results will determine which native profile and resolution to
repeat, and whether another MMG seed control is useful. No 2% or 1% periodic
FEA qualification is claimed until a fine reference has measured convergence.
