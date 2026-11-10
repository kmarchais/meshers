# Project and study reassessment, 10 November 2026

Meshers has a useful native mesher and several measured advantages. The study
does not establish that it beats actual Microgen+MMG on every metric or geometry.
The best next engineering direction is to fit and adapt the final curved mesh.
The stiffness-reference investigation needs its own acceptance plan, rather than
another round of mesh sizing selected against uncertain extrapolations.

This assessment covers the public interfaces, linear-to-quadratic pipeline,
comparison fixtures, accumulated negative experiments and new measurements. It
is not a complete release audit. The original four cases and all nine expanded
specimens remain in scope. No admission limit or reference tensor changes.

## What the project already does well

The Rust core evaluates implicit fields and generates meshes in process. The
Python layer compiles supported fields and derivatives, including the recently
added frozen interpolation tables. Periodic constraints, cap facets and
quadratic node association have explicit checks. These are concrete strengths;
there is no evidence here that a rewrite or a new meshing dependency would help.

The previous production improvements addressed specific defects: Schwarz D
periodic matching, native field evaluation and density-table compilation.
The last density implementation reduced complete meshing time by about 2.7x
under its public analytic-gradient policy. That was a Meshers evaluator comparison,
not a new Microgen comparison. The native binary remains unchanged today.

The public solver examples separately verify manufactured displacement
convergence and fedoo integration. Those checks are useful, but they do not fill
the missing effective-stiffness comparisons in this nine-specimen study.

## What the comparison actually establishes

| Area | Current evidence | Remaining gap |
|---|---|---|
| Specimen coverage | Selected native recipes admit on 9/9 fixed specimens. Admitted actual modified routes cover Schwarz P and the lattice. | Failure of tested actual routes is not a speed or accuracy victory, nor proof of universal native support. |
| Schwarz P, same-order Tet10 meshing | Prior interleaved runs show about 45.56x faster generation, 4.97x fewer tets and 6.56x smaller sampled boundary error. Most quality quantiles improve. | The maximum corner-quality quantile is slightly lower. This mesh-only ratio is not the separate FEA workflow ratio. |
| Lattice, nearly equal Tet10 counts | Prior Uniform37 comparison uses 157,521 native tets versus 158,430 actual tets. Native generation is about 4.86x faster with 9.54x smaller boundary error; most quality quantiles improve. | Native upper-tail quality is slightly worse. Uniform37 stiffness disagreement is still 2.08539%. |
| Best saved lattice mechanics | Best saved native experiment is 1.213247% across all 24 models; the separate modified actual star-refined workflow is 0.994671%. | Native still misses the retained-model criterion. The actual FE mesh differs from the actual mesh-only timing route. |
| Other TPMS, spinodoids and grading | Native meshes admit; several actual routes fail quality, topology, periodicity, count or resource checks. | Successful matched-accuracy actual comparisons and broad FEA convergence evidence are missing. |

These are existing same-host measurements, not timings rerun today. Sources are
[the comparison findings](../../handoff/macbook/comparison-findings-20261026.json),
[periodic and grading follow-up](../../handoff/macbook/periodic-grading-findings-20261026.json)
and [the existing dashboard](comparison-dashboard.html).

Selected recipes are another practical limitation. Gyroid, spinodoid, lattice
and grading use different resolutions or adaptation choices. Passing 9/9 with
those choices does not demonstrate that one automatic public policy handles
all nine. A robust product should accept accuracy and resource requirements,
inspect the result and perform a bounded repair or refinement. That policy
still needs a fresh comparison; it cannot inherit the hand-selected timings.

## The lattice target needs interpretation, not deletion

The 24 tensors include three reference families and several assumed
extrapolation orders. Different discrete domains and nonnested spaces do not
justify treating every extrapolation as a certified continuum value.

The previously verified algebraic lower bound is 0.901931% worst disagreement.
A trial SPD tensor reaches 0.908472%, so satisfying all 24 below 1% is possible.
It leaves little room. Neither observation identifies the physical continuum
tensor. The actual route's 0.994671% result passes the retained test; it is not
automatically a continuum certificate. Conversely, native failure remains a real
failure of that test and is not excused by reference uncertainty.

The fixed-domain enrichment, quadrature checks and residual diagnostics were
useful. They also show why further quality-only motion, curve-only allocation,
the closed bubble family or a global sizing sweep is a poor next experiment.
The completed enriched diagnostic is still about 1.0912%; its next tested residual
family had only a 0.0150% trial signal. Raw traction residuals rank the actual
and native routes differently from the model-disagreement test.

Keep all 24 historical scores. Define a separate reference investigation with
an explicit uncertainty target and independent refinement evidence. First do a
resource preflight. A future distinction between candidate budgets and reference
budgets may be useful, but it needs a recorded protocol decision. Today both
250,000 admission caps remain unchanged. Do not substitute a preferred tensor,
fit geometry or material to references, or call stiffness a local-stress result.

## A new measured geometry opportunity

The core adapts linear tetrahedra before `quadratic::elevate` projects edge nodes.
Linear-face errors therefore drive refinement before the final curved faces
exist. Projecting nodes onto the implicit wall also does not minimize the
distance of the entire curved face between those nodes.

We tested that second limitation directly. The experimental fitter fixes every
corner, every tetrahedron and every surface association. It moves quadratic
wall-edge nodes along constrained normal directions, locks cap coordinates and
ties periodic copies. Nodes may lie slightly off the exact implicit wall so that
the complete quadratic patch approximates it better. This changes the discrete
domain and displacement space; old FE fields and stiffness scores cannot transfer.

The objective uses the unchanged implicit function, not frozen benchmark samples,
reference volumes or stiffness tensors. Six bounded IRLS iterations approximate
an area-weighted L8 objective. Movement is limited to 2% of the shortest incident
wall edge. Cap-intersection distances use gradients restricted to allowed motion.
The existing independent observer makes the admission decision afterward.

| Fixed specimen | Tets, unchanged | Standard boundary error, before → after | Reduction | Absolute volume error %, before → after | Added fit s |
|---|---:|---:|---:|---:|---:|
| Gyroid sheet | 6,008 | 0.00118934 → 0.00117012 | 1.62% | 0.0263276 → 0.0257566 | 0.520 |
| Schwarz P solid | 27,492 | 0.00014660 → 0.00013515 | 7.81% | 0.0000548 → 0.0000250 | 0.449 |
| Schwarz D solid | 27,900 | 0.00012987 → 0.00009792 | 24.60% | 0.0002255 → 0.0002561 | 1.258 |
| Fixed spinodoid | 45,923 | 0.00030398 → 0.00029269 | 3.72% | 0.0258918 → 0.0257594 | 13.626 |
| Anisotropic spinodoid | 40,351 | 0.00114550 → 0.00067775 | 40.83% | 0.0042117 → 0.0045298 | 5.476 |
| Smooth-union lattice | 41,340 | 0.00027866 → 0.00025282 | 9.27% | 0.0312417 → 0.0315781 | 4.271 |
| Linear thickness Gyroid | 25,624 | 0.00027725 → 0.00027027 | 2.52% | 0.0108343 → 0.0109864 | 1.507 |
| Periodically graded Gyroid | 25,744 | 0.00038913 → 0.00037309 | 4.12% | 0.0016902 → 0.0017256 | 1.708 |
| Linear density Gyroid | 25,747 | 0.00209983 → 0.00170496 | 18.80% | 0.0232557 → 0.0206473 | 1.882 |

All nine final meshes pass the unchanged standard admission gates. All corner
coordinates, element/node counts and corner-quality distributions are unchanged.
The density seed changes from rejected to admitted. Its count is 81.8% below
the existing 141,437-tet adaptive analytic recipe, but those two meshes do not
have matched boundary, quality or FE accuracy. The earlier tangent-fit prototype
already admitted this coarse seed; today's new result is broader face fitting,
not the first recovery of that specimen.

Volume error increases in five cases. Many small changes are below the sampled
volume reference's uncertainty; raw values are retained without claiming physical
volume rankings. Schwarz D has an analytic half-cell volume and a real small
volume-error regression. Minimum sampled lattice Jacobian decreases by 1.94%,
remaining positive. Corner-quality invariance is not curved-quality dominance.
Added times are single saved-mesh observations including the diagnostic holdout
work, not repeated complete generation timings. This Python implementation is
too expensive to become the default unchanged.

The initial L2 Gyroid fit was a negative result: its independent maximum boundary
error worsened by 6.55%, despite lower RMS error. Correcting cap-distance scaling
alone still worsened it by 2.38%. The L8 policy then improved it by 1.62%.
These failures are preserved. A lower average fitting loss is insufficient.

Jobs 690–694 all finish below their limits. The nine-case screen takes 78.27 s
at 1.575 GiB; its audit takes 83.58 s at 1.779 GiB. The audit reproduces the
standard results and uses an independent degree13 surface grid. It confirms
improvement on all nine, with zero failed projections. Dense reductions range
from 1.66% to 40.22%. Sphere radial-distance, plane, fixed-corner and face-order
controls pass. Sampling still does not prove global positivity or a Hausdorff bound.

See [numeric findings and 82 input hashes](../../handoff/macbook/face-fit-findings-20261110.json),
[the predeclared selection](../../handoff/macbook/face-fit-selection-20261110.json),
[the experiment](curved_face_fit.py) and [the independent audit](curved_face_fit_audit.py).
The wider approach of joint implicit fitting and element-validity optimization
has a research basis, but that literature does not validate this particular
implementation: [Barrera et al., high-order mesh morphing](https://arxiv.org/abs/2208.05062).

## What to do next

1. **Test a cheap local version of curved-face fitting.** Use native field
   evaluation and select connected defective patches from the final Tri6 geometry.
   Measure complete fresh generation, including fitting and compilation, before
   claiming a speed gain. Start with the density recovery and anisotropic
   spinodoid, then replay all nine. Preserve validity and report volume tradeoffs.
   Do not port the global Python sparse solve into the core without this test.
2. **Drive Tet10 refinement from remaining curved errors.** If local fitting
   cannot meet tolerance, refine those patches with periodic closure. One
   evidence-selected coarsening control can then test fewer elements at matched
   accuracy. Keep the existing linear path for Tet4. An opt-in fitting strategy
   is appropriate while off-surface nodes and new FE domains remain experimental.
   Test one requirements-driven policy across the specimens, not nine manually
   chosen recipes presented as automatic robustness.
3. **Resolve the mechanics target separately.** Use the lattice reference
   preflight above and a small, validated error estimator. Do not run another full
   six-load candidate merely because its geometry or corner quality improved.
   For nonperiodic linear grading, define the finite specimen and loading first;
   periodic homogenization cannot be copied onto unequal opposite material traces.
4. **Complete comparisons where evidence is missing.** Investigate one actual
   workflow failure at a time. Keep defaults, explicitly modified actual routes
   and MMG-on-Meshers seeds distinct. Show same-order Pareto tradeoffs, node/DOF
   costs and reference uncertainty alongside the existing strict dominance test.

The study needs a small current decision index more than another large sweep.
Keep frozen evidence, but treat old planned next steps as historical. Today's
specific next experiment is the local native-evaluator fitting control. This
does not complete the four-case objective or establish superiority on every metric.
