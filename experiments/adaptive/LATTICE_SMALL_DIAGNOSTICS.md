# Small lattice diagnostics after the comparison

26 October 2026. Six bounded jobs test where saved six-load fields fail
pointwise traction equilibrium and whether stiffness-decreasing corrections
are compatible with all 24 retained models. Each job takes less than six
worker seconds and peaks below 0.5 GiB. No mesh, field, material, reference
or native implementation changes. No new FE solve runs.

The experiments support a broader error diagnostic before another allocation
change. Curve-only refinement misses almost half of the native traction
signal. A raw residual score also does not reproduce the ranking by retained
model disagreement. These results do not yet improve a mesh.

[Numeric findings](../../handoff/macbook/lattice-small-diagnostics-20261026.json)
include hashes, two samples, dense quadrature checks, all-model directional
floors and worker resource records. The fixed four-case objective remains
open; the original geometry, quality, periodic and count gates remain intact.

## Traction defects extend beyond curve patches

The sources are the best saved native relocation mesh at 162,189 tets and
249,824 nodes, and the modified actual star-refined mesh at 165,722 tets and
250,000 nodes. Both have six complete saved displacement fields. Their
existing geometry admission and retained-model roles are unchanged.

Geometry defines six face strata before stress evaluation. Each stratum
contributes at most 512 uniformly selected faces, with fixed seeds 20261026
and 20261027 in separate jobs. At interior and periodic interfaces, the
diagnostic measures the jump in traction. At physical free walls it measures
traction against zero. The sampled Gram matrices include every pair of the
six loads. Each face has a diameter weight; sampled stratum sums are expanded
by population/sample count. The trace is normalized by the route's own
stiffness tensor and divided by six.

This is a residual diagnostic. It omits volume-divergence terms and has no
equilibrated reconstruction, error constant or sampling confidence interval.
Its numeric value is not a stiffness-error percentage. Periodic traction
jumps are an FE equilibrium observation, not a failure of the admitted
periodic node/displacement constraints.

| Share of residual score | Native sample 1 | Native sample 2 | Actual sample 1 | Actual sample 2 |
|---|---:|---:|---:|---:|
| Interior interfaces incident to curve cells | 35.19% | 31.31% | 47.91% | 53.75% |
| Physical walls incident to curve cells | 17.25% | 21.16% | 23.63% | 25.79% |
| Other boundary-region interior interfaces | 27.26% | 28.50% | 20.13% | 10.93% |
| Bulk interior interfaces | 14.91% | 14.65% | 2.13% | 2.05% |
| Other physical walls | 3.09% | 3.07% | 2.53% | 3.91% |
| Periodic cap interfaces | 2.29% | 1.31% | 3.67% | 3.57% |

Combined native curve-incident share is 52.45–52.47%. Other groups account
for 47.53–47.55%, consistently including about 15% in bulk interfaces.
Individual shares vary between samples, particularly in the actual mesh.
These are stratified sample observations, not whole-face census values.

The actual residual scores are 0.08095–0.08253, versus native 0.03561–0.03642,
despite actual having better disagreement against the retained models.
The meshes represent different approximate domains and use their own tensor
normalizations. This is not evidence that actual has a larger continuum
stiffness error. It shows why minimizing this raw score cannot replace the
declared all-model stiffness comparison or justify a new allocation alone.

## Independent controls and dense integration pass

Job 684 checks affine and exact quadratic displacement fields on two affine
tets, detection of a continuous displacement kink, orientation invariance,
zero wall traction for tangential stress and detection of nonzero wall
traction. Jobs 686 and 688 also check an affine physical displacement field
on matching curved maps.

The first 9/16-point surface comparison exposes a 5.96% local Gram-vector
change in one native wall stratum. Its original result is retained rather
than treated as converged. Replaying the same faces at positive 36/64-point
rules reduces the final change to 0.00530%. Every stratum in both samples
passes the declared final 1% Gram-vector stability screen.

Independent central differences of shape values reproduce analytic stress
evaluation on saved faces within 4.60e-9 relative. The check uses values
rather than the analytic shape-derivative implementation. Matching interface
maps and periodic translations are checked at every quadrature point.

## Fixed-domain softening has limited directional room

For reference tensor M, the comparison uses eigenvalues of C−M relative to M.
If a future tensor X satisfies X≤C in energy for every strain, then along
the most negative direction of C−M its error cannot improve. Therefore
max over all models of max of zero and minus the smallest eigenvalue is a
necessary disagreement floor for that class of corrections.

This is a conditional algebraic statement for the recorded matrices. Exact
fixed-domain conforming Galerkin enrichment motivates the stiffness-decreasing
class. A changed domain or nonnested approximation is not constrained by this
test. Numerical quadrature and actual FE realizability remain separate checks.

| Saved tensor | Current worst disagreement | Necessary floor for any X≤C |
|---|---:|---:|
| Compact native | 2.412306% | 0% |
| Uniform37 native | 2.085390% | 0% |
| Best saved native relocation | 1.213247% | 0.688944% |
| Complete enriched diagnostic, positive125 | 1.091235% | 0.813475% |
| Modified actual with star refinement | 0.994671% | 0.994671% |

The enriched diagnostic's solved-rule value remains the historical
1.091221736%; the positive125 reintegration differs slightly as already
recorded. The new screen retains all three recorded integration rules.
No floor exceeds 1%, so this necessary test does not rule out further
improvement. It also does not prove that an admissible FE correction exists.
Five matrix controls and explicit witness Rayleigh quotients pass.

A separate algebraic sensitivity control asks whether multiplying an entire
tensor by a scalar no greater than one could pass all models. For generalized
energy ratios with extremes L and U across all 24 references, the allowed
scale interval is [0.99/L, 1.01/U], intersected with scales no greater than one.
Direct eigenvalue replay verifies every reported result.

| Tensor | Uniform algebraic softening can reach 1%? | Best disagreement in that scalar class |
|---|---|---:|
| Compact native | No | 1.155374% |
| Uniform37 native | No | 1.026852% |
| Best saved native relocation | Yes, scale 0.996868–0.997893 | 0.948609% |
| Complete enriched diagnostic, positive125 | Yes, scale 0.998119–0.999098 | 0.951034% |

Those scaled matrices are not FE solutions, references or material changes.
They are not promoted. Their narrow intervals show the need to reduce the
right mixed-load energy errors without worsening already softer directions.
Uniform softening cannot fix the two simpler native recipes, and excessive
softening can make a passing tensor fail.

## Next useful experiment

Develop a genuinely equilibrated stress reconstruction on a tiny controlled
patch before using an estimator to allocate mesh changes. It must preserve
normal traction continuity, zero physical-wall traction, periodic traction
relations and compatibility with rigid translations and rotations. Then check
curved maps and all six loads before trying saved lattice patches.

[Weakly symmetric stress-equilibration research](https://arxiv.org/abs/1808.02655)
supports this direction, including the need for rigid-body compatibility.
Its displacement-pressure spaces do not directly certify our curved Tet10
displacement formulation. That applicability remains
work to prove, rather than a bound to assume.

Do not repeat the closed cubic-bubble residual family, refine only curve
patches, change material stiffness, fit a geometry offset to the references,
or start another global sizing sweep from these diagnostics. No native win
or reference convergence is established. Gyroid, Schwarz P and spinodoid
remain in scope alongside lattice.
