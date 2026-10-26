# Improvement directions after the nine-case comparison

26 October 2026. These are ranked hypotheses and experiment designs. Only
the native-lookup feasibility control below has run in this investigation.
No production change, admitted new mesh, FEA result or speedup is claimed.

The largest scientific gap is lattice stiffness accuracy. The lowest-risk
implementation opportunity is the graded-field evaluator. Element economy
and bulk quality need changes to the algorithm rather than another nearby
global resolution choice.

| Direction | Measured reason | Proposed change | Confidence and first decision |
|---|---|---|---|
| Resolve lattice error sources, then allocate by energy benefit | Best saved native experiment is 1.213247% disagreement across all 24 models; actual star-refined is 0.994671%. Good Uniform37 shapes still give 2.085390%. | Separate fixed-domain discretization error from geometric/reference discrepancy. If discretization remains material, use six-load patch energy changes for both refinement and coarsening, including conformity cost. | Highest importance, substantial research risk. Validate a diagnostic before another full candidate solve. |
| Evaluate the frozen density lookup natively | Admitted density recipes take 10.92–11.39 s and use the callback path. Direct compilation rejects `np.interp`. Its share of generation time is not yet profiled. | Express the same immutable piecewise-linear table with existing comparisons and selections. Start with finite-difference gradients to isolate evaluator placement. | Feasibility verified. Next measure complete generation with unchanged recipes and independent admission. |
| Make Tet10 adaptation aware of the curved output | Density uniform Tet10 uses 25,787 tets and misses distance at 0.002020511. The admitted adaptive recipe has 141,392 tets and distance 0.000583760. | Construct and inspect curved faces, then refine their actual defects locally. Avoid forcing every quadratic recipe to satisfy the linear-face error estimator first. | Code path supports the hypothesis; element savings are unmeasured. First locate all rejected curved patches on the frozen mesh. |
| Optimize curved patches jointly at fixed counts | Native Tet4 and compact lattice meshes have quality-distribution tradeoffs. The larger saved relocation pass raises lattice median q only from 0.685789 to 0.690140. | Try a coupled patch solve over interior corner and edge nodes, using a high-order Jacobian functional and boundary/periodic constraints. | Moderate research risk. Test local feasibility and whole quality distributions before any FEA solve. |

## Lattice accuracy needs a diagnostic that distinguishes causes

The saved-mesh joint-cavity experiment improves metric excess and all seven
reported quality quantiles, yet stiffness disagreement gets slightly worse:
1.213246982% to 1.213444591%. The larger shape pass improves stiffness by
only 0.001929 percentage points. Neither edge length nor quality is a reliable
stand-in for this stiffness objective.

The previous complete fixed-domain enriched field reaches 1.091221736%.
Its next tested residual family offers only a 0.015002691% trial energy gain,
below the existing 0.05% solve screen. That signal is not an upper bound on
the remaining error. Repeating that closed enrichment family is unsupported.

The next useful diagnostic is a stress-equilibration or complementary-energy
feasibility study on a small control, then on saved lattice patches. It must
check normal-traction continuity, traction-free physical walls, periodic
tractions, rigid-body compatibility and curved mapping/quadrature. Retain a
separate geometric-error investigation against the unchanged implicit field.
Do not infer geometric bias merely because the residual screen is small.

Stress-equilibration research provides a basis for local error assessment,
but the cited formulation uses a displacement-pressure pair and specific
stress spaces. It is not a ready-made certificate for our curved, perforated
periodic Tet10 domain. Likewise, published primal/dual bounds for effective
elastic stiffness use a different discretization. Their applicability must
be established before calling any number a bound.
[Elasticity stress equilibration](https://arxiv.org/abs/1808.02655),
[effective-stiffness bounds](https://doi.org/10.1002/pamm.202300098).

If the diagnostic supports a discretization-driven change, the new allocation
experiment should estimate both the energy lost by coarsening a patch and
the energy gained by refining it. Aggregate a 6×6 correction matrix across
all six macro loads, normalized by the current tensor, and account for the
actual periodic closure and added-node cost. Do not select against a favorite
reference tensor. This is different from the previous residual-score-per-node
split followed by length-based coarsening: the net benefit of the exchange
must be measured, and patch interactions still require a fresh global solve.

## The native lookup can reuse the existing compiler

[The preflight script](graded_lookup_preflight.py) reproduces
`unsupported traced NumPy function: interp`, then encodes the same 257-knot
table as a balanced selection tree. No fitted polynomial, new compiler
opcode, external compiler or dependency is required. Endpoint clamping and
the original piecewise-linear intervals are preserved.

[The recorded control](../../handoff/macbook/graded-lookup-preflight-20261026.json)
passes 8,965 probes, including every knot, adjacent representable values,
clamped exterior points and independent random points. The maximum width
difference is 1.11e-16 and the full Gyroid-field difference is 4.44e-16.
The field has 2,071 graph nodes; preparing both native evaluators takes
0.0142 s in this control. Analytical gradients remain disabled. These are
finite numerical checks, not a bitwise mesh-equivalence or speed result.

The first implementation experiment should change only the evaluator for
the admitted density recipes. Measure stage costs and callback counts, keep
the frozen table/fixture hashes, and compare the output arrays. If arrays
change, perform full independent admission; do not transfer old FE scores.
Only after a successful pilot should five interleaved measurements establish
a generation gain. The current 11-second total does not prove that Python
overhead accounts for all or most of it.

## Curved-error adaptation could avoid expensive over-resolution

`implicit.rs` builds/adapts the linear mesh before `quadratic::elevate` runs.
`tet_adapt.rs::enforce` checks straight face centroids and edge midpoints.
The quadratic stage projects edge nodes and limits curvature to retain
positive sampled Jacobians. Therefore, the current pipeline does not choose
its refinement by observing the final curved Tri6 approximation.

The density example is a useful discriminator. Uniform Tet10 misses the
distance gate by about 1.03%; the chosen admitted recipe uses 5.48 times as
many tets. This does not prove that a small local repair exists. Start by
mapping the independently checked distance witnesses to curved boundary
patches, including any patches where curvature was reduced. Then try one
bounded local-refinement policy with complete periodic orbits where needed,
reconstructing and independently auditing the curved output.

The acceptance target is the unchanged distance 0.002 with fewer elements
than the current admitted recipe, without losing quality, volume, Jacobian,
facet, periodic or count admission. A changed curved approximation is a new
FE domain, so old displacement fields are not a nested correction on it.

High-order mesh fitting research supports treating boundary fitting and
element validity together. It does not demonstrate this proposed adaptation
policy or its performance on our fixtures.
[Implicit high-order mesh fitting](https://arxiv.org/abs/2208.05062).

## Quality needs a stronger search direction

The existing research motion moves one corner toward a centroid and carries
incident edge midpoints by half that displacement. This preserves the curved
patch trace, but it restricts the search directions. The strict sorted-quality
guard is useful: its global pass prevents distribution regression. It accepts
few operations and does not itself supply a better search direction.

A new experiment should optimize several interior corner and edge degrees
of freedom together, while fixing every physical Tri6 boundary node and
coupling periodic copies. Use a high-order shape functional with a positive
Jacobian barrier, a bounded line search and a separate full corner-quality
distribution check. Independent midpoint freedom can affect curved validity;
it does not by itself improve corner-tet q. Count gains and FEA gains must
be measured separately.

The published target-matrix approach optimizes quality within high-order
elements through node movement. It motivates the coupled objective; it does
not guarantee dominance of our recorded quality quantiles. A small patch
prototype is preferable to introducing a whole meshing library.
[High-order target-matrix optimization](https://arxiv.org/abs/1807.09807).

## Bounded sequence and stopping rules

1. Complete the lattice error-diagnostic feasibility design before selecting
   another expensive lattice solve. Keep all 24 models and the existing
   0.05% auxiliary-solve screen. Do not revive the closed residual family.
2. Run one native-density evaluator pilot, isolating this low-risk change.
   Read its geometry, quality, arrays and timings before a repeated job.
3. Run one curved-patch localization/preflight on the rejected density Tet10
   mesh before attempting a new adaptation algorithm.
4. Run one coupled-quality patch control only after specifying its objective
   and acceptance rules; do not repeat the old centroid pass with more rounds.

Every heavy stage uses the shared one-worker queue, standard progress,
explicit timeout/memory limits, and read-before-successor policy. Preserve
the 250,000 node/tet caps and every geometry/quality/periodic gate. Any
production change must retain regression coverage for the four original
cases and all nine dashboard specimens. Only fresh same-host timings and
admitted actual routes support new comparison ratios.

These directions leave the fixed field definitions, ordered smooth union,
original samples, historical failures and references unchanged. Broader
coverage and FEA evidence remain part of the objective.
