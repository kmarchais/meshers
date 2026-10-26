# Schwarz D periodic matching and linear grading, 26 October 2026

The Schwarz D periodic failure is fixed in the native scalar and quadratic paths. Both original fixed Schwarz D recipes pass five measured runs after an excluded warmup. The original linear-thickness Tet4 case also admits with one geometry-adaptive policy. A separately calibrated linear-density Gyroid admits in both orders using that policy.

These are coverage improvements. There is no new admitted actual comparator for these cases, no comparative speed victory, and no new stiffness reference. The broader four-case matched-stiffness objective remains incomplete.

## Why Schwarz D failed

The minimal red case uses four cells per axis, zero snapping and zero optimization. It contains the supported material corner `(-0.5, 0.5, 0.5)`, but the opposite `(0.5, 0.5, 0.5)` has no adjacent material volume. The old matcher demanded a node there solely because it matched coordinate planes. A phantom node would be unused.

The corrected matcher keeps every available old pair, requires a bijection for every material-cap node, compares complete Tri3/Tri6 cap facets and rejects any disconnected class of existing periodic edge/corner images. It uses the same coordinate tolerance. The minimal fixture keeps all 173 original points and 408 tets; all material traces match. The matching change leaves the field, phase and mesh coordinates unchanged and adds no orphan nodes.

Linear and quadratic paths share this check. Independent Python checks use spatial searches and a separate union-find implementation. Negative controls reject unmatched cap nodes, wrong triangulations, displaced quadratic cap midpoints and disconnected existing node classes. Six affine macro modes have correct jumps through transitive classes.

The historical all-coordinate-plane rule and its rejected Schwarz D outcomes are retained in the old study and recorded alongside the corrected audit. This is an explicit matching correction, not a relaxed cap, quality or geometry tolerance.

## Admitted repeated results

All rows below pass five measured repetitions after one excluded warmup on the same M4, one configured numerical thread. Times include fresh generation and arm-specific preparation, excluding imports, common fixture preparation, archival I/O and independent audits. Density-table fields use the supported Python callback evaluator, which remains a measured cost.

| Fixed specimen | Mode | Tets | Nodes | Generation s, median [range] | q min / median | Sampled distance | Abs. volume % |
|---|---|---:|---:|---:|---:|---:|---:|
| Schwarz D solid | Tet4 | 94,332 | 20,431 | 0.321319 [0.315781, 0.327324] | 0.445715 / 0.929516 | 0.001352064 | 0.004348 |
| Schwarz D solid | Curved Tet10 | 27,900 | 44,079 | 0.120475 [0.118210, 0.123482] | 0.457266 / 0.918489 | 0.000129866 | 0.000226 |
| Linear density Gyroid | Tet4 adaptive | 132,812 | 31,318 | 10.920944 [10.718819, 10.996630] | 0.261089 / 0.753412 | 0.000992296 | 0.122415 |
| Linear density Gyroid | Curved Tet10 adaptive | 141,392 | 218,792 | 11.392458 [11.192465, 11.417527] | 0.266349 / 0.732800 | 0.000583760 | 0.020342 |
| Linear thickness-parameter Gyroid | Tet4 adaptive | 79,809 | 19,797 | 1.226776 [1.178509, 1.238990] | 0.349324 / 0.822227 | 0.000997822 | 0.129659 |

All rows retain distance <=0.002, absolute volume discrepancy <=1%, corner q>=0.1, positive sampled curved Jacobians, coherent facets, complete quadratic associations and both 250,000 caps. Schwarz D also requires full material-cap periodicity and complete existing node classes. Graded cases use bounded nonperiodic meshing for both arms because their opposite-z material traces differ.

## Linear grading is measured explicitly

The new specimen prescribes cross-sectional material area fraction `rho(z)=0.27+0.18*(z+0.5)` for z in [-0.5,0.5]. The 27% and 45% endpoints are an agent-selected baseline near the previous mean fraction; the user specified linear grading, not those endpoint values.

A fixed two-dimensional Sobol quantile calibrates the Gyroid half-width on 257 z planes. Both workflows use the exact same immutable piecewise-linear table. Four independent scrambles check 66 validation planes with 131,072 samples each. The maximum observed absolute density discrepancy is 0.000466846, or 0.046685 percentage points. This is finite validation, not a certified bound for every z.

The frozen field volume reference is 0.3599233627, with four-scramble spread 0.0002536774. The ideal integrated target 0.36 is not substituted as an exact geometric volume. Field identity against the actual Microgen callable offset passes 4,096 probes.

Initial uniform24 Tet4 and uniform16 Tet10 outputs fail distance at 0.002686906 and 0.002020511. Independent constrained minimization leaves the curved witness at 0.002018558. One declared adaptation-on trial preserves those seed resolutions, sets adaptation tolerance 0.001 and minimum length 0.002, and admits both forms. No new size sweep or fixture revision follows.

The original half-width ramp `0.55*(1+0.5*z)` remains a separate thickness-parameter fixture. Its uniform Tet4 failure at 0.002002836 is unchanged. The same one-policy adaptive trial passes at 0.000997822 with 79,809 tets, versus 78,261 in the rejected uniform mesh. Generation is slower and median q is lower than that uniform output; this is admission recovery rather than an all-metric optimization claim.

## Actual workflow controls and limits

The actual Schwarz D control uses genuine Microgen seeds, the same public MMG calls, a corrected input node check at the public input tolerance and the public bounds-snap operation. Input triangulations may be remeshed; the final audit always requires full caps at 1e-10. Both tested seeds reject ambiguous nodes before MMG. They do not become admitted competitors because the native matcher is fixed.

For the calibrated density case, actual default MMG exceeds the intermediate tet cap and the modified first-pass route fails MMG. The uniform native failures are retained. No successful timing ratio is calculated from any of these rejected routes.

Density callbacks currently cost about 10.9 s for the admitted Tet4 recipe and 11.4 s for Tet10. Compiling the frozen lookup representation is a possible later performance task. It is not implemented or claimed as a gain here.

## Regression and provenance checks

- Full native release suite across all workspace targets/features: 110 passed, one ignored research test.
- All-target/all-feature Clippy passes with warnings denied.
- Full public Python API suite: 185 passed, one skipped.
- Final independent trace and fixed-field controls: 11 passed.
- The original four curved point/tet/Tri6 arrays regenerate exactly. The anisotropic spinodoid, periodic grade and original linear-thickness curved arrays also regenerate exactly. Their prior admitted geometry and original four saved field spaces remain applicable.
- The final audit independently replays seven distinct new meshes across 41 attempts. It also verifies all 20 frozen historical runs and their 30 distinct meshes under the original gate, resolving the preserved baseline binary by hash. All 24 lattice model inputs remain unchanged.

The new native binary is `6394fdcd5e13c42955f96625c2e91800e5b4c2cb2817fb3349de3366dbcfab0a`. Baseline `2005c3e01892098a9d734567220a56ed3ea5d748b20a4f96754fbb55168060b1` remains preserved locally. No FE solve, reference or complete matched-stiffness workflow timing is added. Lattice disagreement remains 2.412306220% for compact native, 2.085389679% for Uniform37, 1.213246982% for the best saved experiment, and 0.994670612% for the separate modified actual star-refined route.

- [Audited findings](../../handoff/macbook/periodic-grading-findings-20261026.json)
- [Matching decision and preserved red fixture](../../handoff/macbook/schwarz-d-periodic-trace-decision-20261026.json)
- [Frozen density calibration table](../../handoff/macbook/linear-density-calibration-20261026.json)
- [Density specimen selection](../../handoff/macbook/linear-density-selection-20261026.json)
- [Density adaptive trial selection](../../handoff/macbook/linear-density-adaptation-selection-20261026.json)
- [Original thickness adaptive trial selection](../../handoff/macbook/linear-thickness-adaptation-selection-20261026.json)
