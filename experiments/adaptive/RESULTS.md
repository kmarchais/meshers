# Adaptive and quadratic implicit meshing: results

## Remaining hierarchical residual screen closes another solve, 25 October 2026

Four reconstruction/Schur controls pass in job 364. They check the enriched
parent gradient and scattered residual against independent quadratic/bubble
energy decomposition. Job 365 then screens six new residual directions of
the complete audited field, including the analytic bubble contribution.
The geometry, original Tet10 space, cubic degree and all physical nodes/tets
remain unchanged. This is a fixed-parent trial, not another equilibrium solve.

Largest additional trial gain is 0.017533229%, below the existing 0.05% signal
screen. Positive 64/125-point gain matrices differ by 9.703e-6 relative, well
within 1%. The trial still fails all 24-model admission at 1.119427778% worst
disagreement. Worker time is 75.50 seconds at 1.238 GiB. No further auxiliary
solve is selected from this screen, and it does not bound the full remaining
FE error or establish convergence. No reference or native candidate is promoted.

The last independently equilibrated enriched result remains 1.136813109%;
the current native workflow remains 1.283822914%. Preserve that distinction.
The next investigation should compare geometric error on the fixed native
and admitted actual domains before another allocation or solve. Keep all 24
models and all four cases. No fitted geometry offset, new sizing sweep or
change to the 1% stiffness, distance, quality, periodic or count gates follows.
See [remaining residual trial](../../handoff/macbook/lattice26-hierarchical-remaining-trial-20261025.json).
Earlier selected-next descriptions below are historical.


## Complete hierarchical fields audited; lattice target still fails, 25 October 2026

Six fresh zero-macro auxiliary responses finish in job 361, using the exact
current lattice26 sparse parent operator. The 6x6 Schur correction preserves
the curved geometry, all original Tet10 functions and the 189,180 nodes /
123,694 tets. Its six extra scalar mode amplitudes do not turn it into a native
Tet10 mesh workflow. Nodal displacements alone omit its analytic cubic bubbles.
Worker time is 123.82 seconds at 5.148 GiB under the unchanged 8 GiB watchdog.

Independent job 362 reconstructs the complete fields with positive 15, 64 and
125 point integration and freshly assembles all coarse equilibrium rows.
The operator SHA reproduces exactly. Maximum coarse residual is 2.057e-10,
auxiliary residual 1.775e-11, modal residual 6.730e-13 and reaction/energy
relative difference 3.650e-12. Galerkin correction-energy balance differs by
2.169e-11 relative. Positive 125-point fixed-field energy differs from the
solved tensor by at most 1.185e-7 across mixed loads. Complete cubic traces
and zero-macro response periodic jumps pass. Audit takes 58.66 seconds at
4.820 GiB. No additional equilibrium solve occurs during the audit.

The largest stiffness correction is 0.147582534%. Worst disagreement across
all 24 unchanged models drops from 1.283822914% for the native parent to
1.136813109% for the enriched field. This still fails 1%. No reference,
native candidate or comparative timing claim is promoted. The admitted
modified actual Microgen route remains ahead on this lattice accuracy gate.

Saved-array replay 363 checks local positive semidefiniteness, global gain
balance and rule agreement. Curve-vertex incident cells carry 54.81% of the
largest-direction correction in 6.10% of cells, down from the unrelaxed trial's
60.18%. Wall cells without curve vertices carry 25.43%; cells touching neither
wall nor cap carry 17.32%. Groups overlap. This maps the resolved correction,
not the remaining error, and does not justify another curve-only sizing sweep
or copying scores to a different domain. Before another native allocation
change, distinguish unresolved approximation error from geometric error with
a same-domain residual check of these complete equilibrated fields.

See [complete-field audit](../../handoff/macbook/lattice26-hierarchical-coupled-audit-20261025.json) and
[localization replay](../../handoff/macbook/lattice26-hierarchical-localization-20261025.json). Schwarz's repeatable modified-actual-route
wins stand. Gyroid and spinodoid still lack admitted actual comparators in
tested routes. All four specimens, all 24 lattice models, fixed geometry,
quality, periodic gates and 250,000-node/tet caps remain in scope.
Earlier selected-next descriptions below are historical.


## Current-domain hierarchical coupling preflight passes, 25 October 2026

Eight virtual-work/force-balance/Schur controls pass. Job 359 then verifies
consistent 15-point mode coupling against positive 125-point integration:
mode Gram and gain changes are 6.655e-5 and 4.023e-5 relative, within 1%.
The actual sparse parent has 545,136 free rows and 44,320,304 stored entries,
534,024,196 CSR bytes. Six coupling columns need 26,166,528 reduced bytes.
All saved-field residuals pass, and the exact parent operator SHA reproduces
the six-load run. Four-level AMG setup passes at 5.018 GiB with CG intercepted
before entry. No dense global mode rows are passed to AMG or matrix saved.

Parser failure 357 and key-comparison failure 358 remain frozen. The latter
compared the (shape,SHA) tuple with a SHA string; comparing its SHA member
passes without changing the exact identity or physics gates. No numerical
operator discrepancy was established. Worker 359 takes 43.46 seconds.
See [coupling preflight](../../handoff/macbook/lattice26-hierarchical-coupled-preflight-20261025.json).
Next permit six fresh auxiliary zero-macro response solves and a 6x6 Schur
correction, followed by independent complete-field and positive-quadrature
audits. Geometry, original Tet10 space and node/tet caps remain fixed. No
reference, native candidate or timing victory is promoted; all four cases
and all 24 model tensors remain in scope.
Earlier selected-next descriptions below are historical.


## Lattice26 same-domain residual trial resolves a correction, 25 October 2026

After three passing analytic controls, job 355 builds six global trial modes
from shared cubic edge/face bubbles on the immutable current lattice26 mesh.
All original Tet10 nodal functions and geometry maps remain available; this
trial holds the original solved coefficients fixed and optimizes only six
extra scalar mode amplitudes per load. It is not a new equilibrium solve,
full cubic nodal discretization, native mesh candidate or workflow timing.
The 411,670 support features and 7,410,060 fixed coefficient entries are
reported explicitly. No new geometric nodes/tets are created or admitted.

The largest fixed-domain energy gain is 0.117081127%, above the existing
0.05% signal screen. Positive 64/125-point gain matrices differ by 8.564e-6
relative, within the unchanged 1% gain-precision gate. The complete cubic
trace check covers 494,776 shared/periodic face records, using the six exact
P2 zeros and four complementary points; maximum jump is zero. About 60.18%
of the largest-direction correction energy lies in the 7,548 curve-vertex
incident cells, 6.10% of the mesh. This localizes this trial correction only.

Worst disagreement remains 1.167137138% over all 24 unchanged models, failing
1%. Worker time is 68.12 seconds at 1.091 GiB under the unchanged 8 GiB cap.
The failed first control 352 exposed an installed Tet4 derivative-list type;
its source is preserved, conversion is corrected, and controls 353–354 pass.
See [hierarchical trial](../../handoff/macbook/lattice26-hierarchical-trial-20261025.json).
Next check consistent 15-point coupling, actual reduced assembly and AMG
setup feasibility before any fresh auxiliary response solve. Do not feed
six dense global mode rows to AMG. No reference or native workflow is promoted;
all four cases and all 24 models remain in scope.
Earlier selected-next descriptions below are historical.


## Fresh lattice26 fields pass; all-model target still fails, 25 October 2026

The unchanged recorded curved Tet10 recipe on native 7bad4a45 admits 189,180
nodes / 123,694 tets, versus 199,972 / 128,300 historically. Minimum corner
quality rises from 0.149329 to 0.207354 and q_p01 from 0.293528 to 0.522231.
Sampled distance is slightly worse, 0.000343162 versus 0.000323640, within
0.002. Volume error is also nominally slightly worse; all original geometry,
positive Jacobian, facet and periodic gates pass. No all-metric win follows.

Fresh Mac workflow 348 completes all six AMG loads in 151.245 seconds including
generation/workflow checks, at 5.380 GiB RSS under the unchanged 8 GiB cap.
Independent audit 349 passes reactions, energy and periodic jumps; positive
125-point fixed-field reintegration shifts worst mixed-load energy by only
4.10e-8, below 1e-4. Assessment 351 preserves all 24 model tensors and reports
worst disagreement 1.283822914%, versus previous 1.268698522%. The 1% target
still fails. No reference is promoted, and no comparative timing ratio is
claimed. The existing admitted modified actual lattice route remains a case
where Meshers is behind on the retained-model target.

Assessment 350's raw-label identity requirement failed on numbering alone.
Its source/failure is preserved; 351 proves exact coordinates and all local
Tet10/Tri6 associations after node/whole-row canonicalization, without sorting
local shape nodes or changing any scientific gate. The pilot now uses the
active venv, records actual host/dependency/thread metadata, enforces both
250,000-node/tet caps and fails the queue stage on incomplete requested loads.
See [current lattice assessment](../../handoff/macbook/lattice26-current-native-assessment-20261025.json).
Next use these fully audited fields for one same-domain hierarchical residual
trial, with no sizing sweep, geometry change, new reference or native workflow
claim. Keep all four cases and all 24 models in scope; the full goal is open.
Earlier selected-next descriptions below are historical.


## Public re-entrant periodic recipe recovered, 25 October 2026

Complete region-conjunction end traces now receive a conservative native
certificate before physical seed construction and the existing periodic gate.
Red 333 fails; green 336 and independent public replay 339 pass. The unchanged
public notched recipe admits 3,896 tets / 898 nodes, qmin 0.301405 and roundoff
geometry/volume, with 54 paired nodes and 82 matching cap triangles per end.
The public arrays exactly match the earlier private recovery. Historical
failure 328 is unchanged. Full native suite 338 passes 101 tests / one explicit
research skip; Python 340 passes 211 / one optional skip. Independent audit
346 proves all four fresh fixed study meshes and the periodic column exactly
unchanged from native 02d1c618. Final arm64 image is 7bad4a45.
See [periodic region results](PERIODIC_REGION_RESULTS.md).
No FE solve, reference promotion or comparative timing victory follows.
Next replay the recorded lattice26 curved Tet10 recipe as one geometry and
resource preflight on the new native image, before considering fresh loads.
Keep its settings, 250,000-node/tet caps, fixed specimen and all 24 models.
The full objective remains open. Earlier selected-next descriptions below
are historical.


## Re-entrant periodic failure has a controlled private recovery, 25 October 2026

The new auxiliary notched prism fails the public scalar recipe 328 at error
0.0106363 against 0.002. Its fixed volume is 0.18625 and its compiled graph has
eight planes / two region masks. Manufactured 329 proves feasibility at qmin
0.111168. Physical seed 330 has 54,534 tets / 12,264 nodes and exact geometry,
160 matching end nodes and 252 cap triangles per end, but qmin 0.0062124 and
is explicitly unqualified. Private coupled adaptation 331 of the exact seed
bits passes native final and independent geometry gates at 3,896 tets / 898
nodes, qmin 0.301405, roundoff distance/volume and 54 paired nodes / 82 matching
cap triangles. Independent investigation replay 332 passes. This does not
recover the public route or promote a reference. Native 02d1c618 and all four
fixed study meshes are unchanged. See [periodic region results](PERIODIC_REGION_RESULTS.md).
Next certify complete region end traces before considering public dispatch.
All four study cases and all 24 lattice tensors remain in scope.
Earlier selected-next descriptions below are historical.


## Periodic prism orientation control passed, 25 October 2026

Job 326 admits the recovered prism rotated 0.37 radians about z and translated
(0.02,-0.03,0), without clipping the side walls or changing its volume/end
traces. It uses 3,879 tets / 865 nodes, qmin 0.304041, roundoff distance/volume,
46 node pairs and 68 matching cap triangles per end. Independent replay 327
passes at identical options/gates. Thirteen fixture controls pass. Native
02d1c618 and all four fixed study specimens are unchanged. No FE/reference
or speed claim follows. Next run one periodic re-entrant notched-prism
control to test nonconvex branch coverage; no sizing sweep is queued.
Earlier selected-next descriptions below are historical.


## Public periodic affine normal-jump recipe recovered, 25 October 2026

The unchanged public tapered-prism recipe now admits 3,811 tets / 844 nodes,
qmin 0.309003 and roundoff geometry/volume error, with 45 periodic node pairs
and 66 matching cap triangles per end. Complete normalized affine end-trace
multisets must match before physical cutting; the actual seed must then pass
the existing node/cap gate before coupled torus adaptation. Named fields keep
their individual trace requirement. Red 311 rejects; green 314 and independent
replay 317 pass. No old failure was reclassified. Full suite 316 passes 99 native
tests / one explicit research skip; Python 318 passes 208 / one optional skip.
Fresh controls 319–323 and independent replay 325 prove all four study meshes
and the periodic column exactly unchanged from native 99817e40. Audit 324
selected the earlier failed spinodoid CLI directory; its failure is retained,
and 325 explicitly selects the valid v2 predecessor. No FE solve, reference
promotion or comparative timing claim follows.
See [periodic affine trace results](PERIODIC_AFFINE_TRACE_RESULTS.md).
Next check one rotated version of the recovered prism for orientation
robustness before optimization. Original four specimens and all 24 lattice
model tensors remain fixed and in scope. The full objective remains open.
Earlier selected-next descriptions below are historical.


## Verified affine curve mobility, 25 October 2026

The adapter incorrectly froze vertices on two affine planes. Rank-aware
classification now permits their valid line motion and holds endpoints fixed
by an independent third constraint. The exact saved periodic-prism seed
recovers from ten cells below q = 0.1 to none, qmin 0.304656, 3770 tets / 838 nodes,
with geometry, positive Jacobians, periodic caps and coherent facets admitted.
This is private-core research; the public prism route is still rejected.
Full native suite 297 passes 97 tests / one explicit research skip; Python 301
passes 206 / one optional skip. Independent study audit 308 proves all four
fresh fixed study meshes exactly unchanged against native 30506b13. The
periodic affine column remains admitted, now 2614 tets / 604 nodes, qmin 0.321550
and roundoff distance. No timing, FE accuracy or reference claim follows.
See [affine curve results](AFFINE_CURVE_RESULTS.md).
The next failure is the public construction of a verified periodic seed with
permuting end traces. Keep the individual named-field periodic gate intact.
All four cases and all 24 lattice models remain in scope; the full objective
is unachieved. Earlier selected-next descriptions below are historical.


## Verified periodic image-query correction, 25 October 2026

The red/green normal-jump regression proves local out-of-cell error queries
needed the prescribed periodic image. The correction changes no field, gradient
policy or admission gate. Native suite 283 passes 95 tests / one explicit
research skip; Python suite 286 passes 206 tests / one optional skip. Fresh
controls 287–290 admit all four cases, independently replayed in 292. Gyroid,
Schwarz and spinodoid are exactly unchanged. Lattice uses 213,800 tets / 46,131
nodes instead of 218,239 / 48,077, with unchanged distance 0.001768701 and
improved reported quality quantiles. Nominal volume error is slightly worse
(0.429176% versus 0.429123%), within the same 1% gate. No timing ratio or FE
accuracy claim follows. The original periodic affine control remains exact.
See [periodic distance results](PERIODIC_DISTANCE_RESULTS.md). Saved-seed replay 293 removes 111
false error splits, but ten bad cells remain and the public tapered-prism
recipe still rejects. All ten touch affine curves; inspect their actual
constraint freedom before another native change. Original four cases and
all 24 lattice models remain in scope; the full objective is unachieved.
Earlier selected-next descriptions below are historical.


## Periodic normal-jump failure isolated to a feasible input,25 October2026

The invariant affine fix is committed/pushed as6c843b6. The next exact tapered
prism has matching periodic scalar traces and a33.4degree wall-normal jump.
Raw adaptive272 rejects error0.0106425; mapped273 rejects0.00234375; zero-round
mapped274 rejects q0.045969. None is reclassified. A separate manufactured
20tet mesh passes all gates at q0.171170,proving geometric/quality feasibility
without claiming a native workaround. Independent replay276 passes;12fixture
controls pass. No production source,native image,FE result or reference changed.
See [periodic normal-jump investigation](PERIODIC_SEAM_RESULTS.md). Next inspect one physical
six-plane seed and actual periodic end correspondence before a native change.
Keep the mapped recipe closed and preserve all four cases and24 lattice models.
Earlier selected-next descriptions below are historical.


## Explicit resume and verified periodic affine recovery, 25 October 2026

Kevin resumed from the pause checkpoint. Native test265 had finished; it was
not restarted. The periodic column now admits 2,663 tets / 655 nodes, qmin 0.302251
and distance 6.21e-17 at unchanged settings. Independent replay 271 reproduces
the rejected scalar corner witness and exact named-domain control, and proves
all four freshly generated study meshes exactly unchanged after relabeling.
Python/fixture/tools suite 264 passes 205 tests / 1 optional skip; full native
release suite 265 passes 94 tests. Final build 266 reproduces the measured image.
See [periodic affine results](PERIODIC_POLYHEDRAL_RESULTS.md). No new FE/reference/timing claim follows.
The supported periodic planes must each be invariant along requested axes;
non-affine normal jumps and branch permutations remain open. Next diagnose
one actual periodic-seam feature mechanism before any further native change.
All four cases and all 24 lattice models remain in scope. Earlier pause and
selected-next descriptions below are historical.



Eligible nonadaptive Tet4 calls now retain analytic outer-max branches.
The bored-prism scalar route that failed at matched settings now admits with
145924 tets, qmin 0.121400, distance 0.00105814 and volume error 0.05022%.
Its mesh is exactly equal after relabeling to the admitted seven-field named
control. All four original study meshes also replay exactly against their
preserved baseline. Full suite 236 passes 198 tests / 1 skip; final all-feature
binding Clippy 238, Python lint and typing pass. See
[nonlinear constraint results](NONLINEAR_CONSTRAINT_RESULTS.md).
Adaptive nonlinear postcontrol 244 still rejects at 0.003393814; no adaptive
option is bypassed. Next carry actual field identities through adaptive
splitting, collapse and projection. No FE solve, reference promotion, model
removal or new actual Microgen timing ratio. Actual Gyroid/spinodoid admission
and all 24 lattice convergence remain open. The full goal stays active.


Meshers now supports bounded nonconvex min/max combinations of convex affine
fields for nonperiodic Tet4. Shared cuts preserve one conforming union; actual
boundary facets define constraints, and smoothing accepts any complete region.
The L-prism default adaptive control admits with 7354 tets, qmin 0.269396 and
near-zero geometric error. Rotated adaptive controls pass; the rotated,
unsnapped nonadaptive recipe needs two explicit optimization rounds. Its
one-round quality failure remains rejected, and the failed reordering trial
is reverted. See [nonconvex results](NONCONVEX_RESULTS.md).
Final native 83247d4d controls admit all four original specimens and exactly
match prior meshes after canonical relabeling, with identical fixture hashes
and configurations. Full Python suite 223 passes 190 tests / 1 skip; full native
suite 202 passes 91 tests; final lint-clean native check 224 passes 62 unit/intersection tests.
All-target core and all-feature binding Clippy pass. No FE solve, reference
promotion or new actual Microgen timing ratio. The full goal stays active.
Next investigate one nonlinear wall/cap intersection failure relevant to the
lattice boundary limits, keeping all four original controls. General nonlinear
and periodic feature handling, admitted actual Gyroid/spinodoid comparators
and all 24 lattice reference convergence remain open.


## Adaptive sharp-feature continuation, 24 October 2026

Adaptive convex polyhedral Tet4 generation now preserves exact plane constraints.
The rotated control exposed a nonperiodic seed-welding panic; preserving vertex
identities and checked collapse of numerical corner copies fixes it. The final
adaptive octahedron admits with 2565 tets, qmin 0.300440 and distance 8.97e-17.
All four fresh original study controls admit and exactly match their previous
meshes after canonical relabeling. Final suite 176 pass / 1 skip, core 49 pass and
Clippy pass. See [adaptive polyhedral results](ADAPTIVE_POLYHEDRAL_RESULTS.md).
General nonconvex/nonlinear/periodic feature handling, actual Gyroid/spinodoid
comparators and lattice reference convergence remain open. No FE solve,
reference promotion or new timing ratio; the full goal stays active.


## Automatic exact polyhedral feature route, 24 October 2026

The default fused octahedron now admits without a manual plane rewrite. Native
compiler metadata recovers its eight affine constraints and reuses the existing
constrained cutter. Same-input baseline rejects at error0.003240370; candidate
has13,488tets/3,451nodes, qmin0.317878, sampled distance8.28e-17 and zero reported
relative volume error. The old package still fails the regression. Final suite
passes171tests1skip and Clippy/lint pass. All four fresh fixed-study controls
pass with exact geometric equality to saved baselines, allowing lattice node
relabeling. This does not add general adaptive/nonlinear/periodic crease support
or prove the original Microgen/MMG objective. See
[scope and evidence](POLYHEDRAL_RESULTS.md).

## Analytic topology coverage, 24 October 2026

Eight saved Tet4 cases pass fresh replay, including all four study specimens
and four new analytic geometries. A torus and separated spheres initially
fail only the 1% volume gate. One analytic volume/area-derived tolerance each
reduces deficits to 0.230699% and 0.464510%. A hollow shell passes directly.
The fused sharp octahedron fails reverse distance 0.002418269. Its exact
named-plane route passes with 13,488 tets, minimum quality 0.317878 and sampled
distance below 8.41e-17. A production regression now protects four-plane
corners. No native algorithm, study reference or comparative timing claim
changes. See [full evidence](TOPOLOGY_COVERAGE.md), including auditor/source
provenance failures and their corrections. Universal input support remains open.

## Recovered linear coverage, 24 October 2026

The resumed failure investigation corrects the native admission gap reported
below. Historical adaptive32 recipes pass three fresh generations each on
the Mac. Lattice has 218,239 Tet4, distance 0.001768701 and minimum quality
0.171902. Spinodoid has 72,528 Tet4, distance 0.001850419 and minimum quality
0.232844. Jobs 118–119 require admission, and independent job 120 reproduces
all four saved meshes' geometry, quality, facets and periodicity. No native
patch, threshold relaxation or reference promotion occurred. All four fixed
specimens have verified linear routes; arbitrary-input coverage remains open.
See [Tet4 coverage correction](TET4_COVERAGE.md).

## Four-case mesh-generation comparison, 24 October 2026

Kevin explicitly prioritized generation speed, quality and geometry fidelity
at fewer or equivalent elements. The completed comparison separates Tet4 from
matched-order Tet10, excludes FE solve time, and uses fresh same-host workers,
one warmup and five alternating-order measured repetitions per setting.

Schwarz Tet10 fast generation is 45.31 times faster than the declared repaired
actual route, with 4.97 times fewer tetrahedra, 4.82 times fewer nodes and 5.94 times
lower sampled distance. All five reported quality quantiles improve. Four-pass
native generation is 18.68 times faster. Both sides pass every mesh gate 5/5.

Lattice's geometry-only uniform 25 recipe is 10.90 times faster, with 3.25times
fewer tetrahedra, 3.18 times fewer nodes and 4.45 times lower sampled distance.
Four-pass generation is 3.32 times faster with 7.44 times lower sampled distance.
Minimum, p1, p5 and median quality improve; p95 is slightly lower. This is not
strict dominance of every listed quality quantile. Both sides pass 5/5.
The larger FE-oriented adaptive 26 recipe is retained as supporting evidence,
but it is not the appropriate mesh-only Pareto point.

The linear Tet4 comparison is weaker. Schwarz native adaptation uses fewer
cells and improves the low-quality tail, but has 3.3 percent larger sampled
surface error and slightly lower median quality. The actual lattice Tet4
control passes; the tested native Tet4 settings miss fidelity or budget.
Gyroid native Tet4/Tet10 and spinodoid native Tet10 pass 5/5, while the tested
actual routes reject. No admitted speedup is claimed from those failures.

Independent job 114 passes 138 saved-mesh replays and 76 input hash checks,
including counts, quality, geometry, Jacobians, facets, periodicity and affine
Tet10 mapping. No unused nodes are found. Three fixed-recipe/Fourier controls
pass. The native field domain and fixed specimens are unchanged. No FE solve,
reference promotion or four-case project completion follows from this benchmark.

Next use these mesh-only gaps to prioritize linear lattice/spinodoid admission
and robust quality/fidelity improvements, rather than resuming the older
stiffness-reference sequence automatically. Keep all four cases and every
original gate in scope. [Mesh comparison results](/Users/kmarchais/dev/meshers/experiments/adaptive/MESH_COMPARISON_RESULTS.md) contains all measured settings, ranges,
versions, caveats and source links.

## Complete parent25 correction localization, 24 October 2026

Saved-field job96 passes in157.422seconds at0.7823GiB, without another solve.
Both positive rules localize the complete coupled-minus-base field, including
all relaxed quadratic responses. The cell matrices sum to its global correction
energy. Galerkin identity differences against the independently audited stiffness
reduction are0.0779213% and0.0780281%, below the declared1% gate. Positive-rule
correction tensors differ by5.03e-7 relative. Every local matrix is positive
semidefinite within1e-12. On64 actual cells per rule, the original grouped
integrator reproduces the local correction matrices independently.

At216points, curve-vertex cells carry16.3654% of normalized six-load correction
energy in0.08627% of volume. Wall-vertex cells without curve vertices carry
74.6016%, cap-vertex cells without curve vertices32.9546%, and cells touching
neither wall nor cap4.86065%. Wall and cap groups overlap and are not additive.
This is the resolved correction's distribution, not remaining-error certification.
It replaces the trial-only20.6100% curve fraction for this purpose; that older
result remains correct for its different, unrelaxed trial field.

A direct independent replay reads both local archives, checks exact input
hashes, recomputes PSD, global sums, normalization, group fractions and Galerkin
identities. It passes67 input hashes. Receipt:
`handoff/macbook/verification-20261024-parent25-localization.json`. This is a
localization replay; full checkpoint95 remains the last complete replay. No
claim is made that95 already covered96. The new replay branch is available to
future full checkpoint checks.

Stop curve-only allocation. Previously rejected cap-length and energy-sizing
recipes also remain closed; this map does not revive them or justify copying
cell scores to a different physical domain. Next review the native periodic
normal/feature path against the existing normal-jump evidence and this full
correction map. Select a controlled geometry regression only if it tests a
mechanism distinct from those closed recipes. Do not infer a new native win
from reference-only auxiliary work.

The full field still misses all24 models at1.711485183%. No reference, native
candidate or workflow ratio is promoted. Schwarz's repeated modified-actual-
route wins stand; Gyroid/spinodoid actual-comparator gaps and lattice reference
convergence remain. All four specimens and original gates remain in scope.
Record: `runs/macbook-lattice-parent25-hp-localization-v1`. Jobs through96 are
terminal and no successor is queued at this checkpoint.

## Parent25 coupled field passes audits but misses the target, 24 October 2026

The user explicitly resumed work after the pause. The queue and native worker
remain healthy while another macOS user is logged in. Parent25's fresh coupled
field passes energy, stress, periodic, pin and reassembled equilibrium audits.
Worst mixed-load disagreement across all 24 retained models is 1.711485183%,
still above 1%. No reference, native candidate or timing advantage is promoted.

Assembly91 finishes in 10.093 seconds at 5.3634 GiB. Its constrained matrix has
643,332 rows and 48,227,764 nonzeros. Symmetry and constraint virtual work pass;
the saved first-load residual is 9.79e-13. Setup92 finishes in 13.364 seconds at
4.8203 GiB, with four AMG levels and operator complexity 1.055315. It stops
before CG. Each job is read before its successor is selected.

Solve93 completes all six fresh zero-macro responses in 143.600 seconds at
5.3803 GiB, under 240-second per-response and 1,600-second worker timeouts.
Every saved baseline residual is below 2.15e-10. The Schur matrix is positive
definite; response periodic jumps and pin values are zero. These are auxiliary
responses on the frozen domain, not six fresh macro workflows or a native mesh.

Audit94 finishes in 225.131 seconds at 5.1301 GiB. Its own-rule full-field energy
agrees with the prediction within 4.16e-13, and stress-energy balance within
3.96e-12. Positive125/216 energy changes are 5.82e-9. The complete fields have
periodic jump below 4.45e-16 and an exact zero pin. Fresh operator replay gives
all auxiliary residuals below 1.06e-11, complete-field residuals below 2.10e-10,
and mode stationarity within 2.10e-12. Relative to the same-domain fine quadratic
base, all six high-rule stiffness corrections are negative, the largest 0.513183%.
The own15 predicted maximum reduction is 0.513381%; the two numbers use different
quadrature and must remain distinguished.

The positive-rule audit algebraically folds the original cubic and quadratic
fields with the mode coefficients and response fields before bounded integration.
It retains the full energy and stress. Six focused controls include comparison
with the original unfolded18-field energy/stress formula on distinct child maps.
Default parent23 paths remain unchanged. Exact measured Python snapshots preserve
the source versions used by earlier runs; no frozen data or model hash is edited.

All 24 models remain in scope. The 1.711485% result is a controlled fixed-domain
correction, not a continuum certificate or proof of convergence between different
parent geometries. It does not replace native26's live 1.268699% result or saved
1.183637% control. The fresh explicitly modified actual lattice route remains at
0.994671%, with no lattice matched-target timing ratio. Schwarz retains five fresh
paired wins against its declared modified actual route. Gyroid/spinodoid actual
comparators remain inadmissible in the tested routes. All four cases and every
original geometry, quality, periodic, volume and resource gate remain unchanged.

Next localize the resolved complete correction on this same domain before choosing
another native boundary allocation or crease change. Use its Galerkin energy
identity and both positive rules; distinguish correction energy from remaining
error. Do not launch another curve-only sizing sweep or transfer cell scores
between different physical domains without a separate justification.

Independent replay95 passes3,509 transferred files,2,557 hashes and312 full
model comparisons, including both positive-rule unfolded energy/stress samples.
Receipt: `handoff/macbook/verification-20261024-parent25-coupled.json`.

Records are `runs/macbook-lattice-parent25-hp-assembly-v1`,
`runs/macbook-lattice-parent25-hp-setup-v1`,
`runs/macbook-lattice-parent25-hp-solve-v1` and
`runs/macbook-lattice-parent25-hp-audit-v1`. The prior pause and older next-action
entries below are historical.

## User-requested pause after coupling preflight90, 24 October 2026

The user requested stopping after the current job. Job90 is terminal success,
78.101 seconds at 0.50165 GiB. No successor is queued. The goal is paused,
not achieved. Resume only on explicit user instruction.

The parent25 15-point coupling preflight passes. Its mode Gram and gain differ
from the positive216 rule by 0.0265102% and 0.0288212% respectively, below the
unchanged 1% consistency gate. Virtual-work error is 6.66e-16 and all six mode
Gram eigenvalues are positive. The saved coupling has 749,994 by six entries,
35,999,712 bytes. Four focused controls pass before submission; three cross-
energy controls also pass inside the worker. Saved coupling algebra, force
balance and exact input hashes were rechecked directly after completion.

This is a force/rule preflight only. No reduced matrix assembly, AMG setup,
auxiliary solve, reference promotion or native-candidate promotion follows.
A further full checkpoint replay was not queued after the stop request.
The prior full replay89 remains the last complete evidence audit. Receipt:
`handoff/macbook/pause-20261024.json`.

On explicit resume, first read job90's saved result and source metadata, then
inspect one bounded reduced-assembly memory preflight using the existing
`hp_coupled_preflight.py --parent25 --assembly` route. Read it before choosing
AMG setup or a solve. Do not restart completed jobs or run older next-action
instructions automatically. All four cases and every gate remain in scope.
Schwarz's five modified-actual-route paired wins stand; Gyroid/spinodoid actual
comparators and lattice native accuracy/reference convergence remain unresolved.

## Parent25 saved h/p trial resolves a surface signal, 24 October 2026

The exact-domain parent25 trial passes positive quadrature checks, but worst
mixed-load disagreement across all 24 retained models is 1.859718171%, above
1%. No reference or native candidate is promoted. The six periodic cubic-minus-
quadratic modes reduce the fine quadratic field energy by at most 0.374103334%.
This is a saved-field trial with fixed coefficients outside those six modes,
not a new equilibrium solution, remaining-error bound or live workflow timing.

Preflight85 verifies transferred cubic data hashes, original parent embedding
within 1.67e-16, composed maps, complete partition, original field reproduction,
and periodic fields on both node sets. Fine quadratic counts remain 249,998
nodes and 127,549 tets; the cubic mesh has 241,555 nodes on 48,696 parent tets.
The physical domain is unchanged. Frozen G3 evidence is read through resolved
path prefixes; no source hashes or inputs are rewritten.

The grouped integrator finds 42,838 map patterns and projects 5,198 seconds for
the two positive rules. Bounded cell batches avoid those singleton map loops.
Three polynomial controls pass. On the actual fields, benchmark86 measures
30.224 seconds grouped and 2.206 seconds batched, with all 144 entries agreeing
within 7.50e-16. This 13.70x integrator ratio is a resource result, not a
meshing-workflow speedup. Four-point tensors are discarded for science claims.

Trial87 fails before integration because its provenance list names the wrong
`periodic.py` path. Its log, failed source and failure note are preserved.
Corrected trial88 finishes in 160.155 seconds at 0.8571 GiB. Positive 125- and
216-point rules change the trial tensor by 3.86e-9 and the gain tensor by
4.58e-7 relative. The high-rule quadratic/cubic diagonal blocks agree with
their saved own energies within 3.17e-7 and 2.68e-5 respectively, below the
unchanged 1e-4 cross-energy check.

The normalized six-load gain puts 20.6100% in curve-vertex cells, occupying
0.08627% of volume. Wall-vertex cells without a curve vertex carry 70.6453%;
cap-vertex cells without a curve vertex carry 36.7881%. These last two groups
overlap and must not be added. Cells touching neither wall nor cap carry
3.71554%. These fractions localize this trial's gain only. They support stopping
the coarse-bulk curve-only allocation, but do not certify the remaining error
or justify transferring a score between different physical mesh domains.

Next inspect a bounded resource preflight for re-equilibrating these six
periodic modes against the fine quadratic space, with consistent assembly
quadrature and full residual checks. Read that preflight before any solve.
Keep the 8 GiB ceiling and 250,000 physical node/tet caps. Do not repeat a
curve-only sizing sweep or treat the trial field as equilibrated. Any subsequent
native change needs its own geometry and full-tensor evidence.

The fresh explicitly modified actual lattice route remains at 0.994671%;
native26's best saved control remains at 1.183637%, and its live candidate at
1.268699%. No lattice matched-target timing advantage follows. Schwarz retains
five fresh paired wins against its declared modified actual route. Gyroid and
spinodoid still lack admitted actual comparators in tested routes. All four
fixed specimens, every admission gate and the original unachieved goal remain.

Independent replay89 passes 3,509 transferred files, 2,261 input hashes and
288 full model comparisons, including sampled grouped checks of both positive
rules. Receipt: `handoff/macbook/verification-20261024-parent25-hp.json`.

Records are `runs/macbook-lattice-parent25-hp-preflight-v1`,
`runs/macbook-lattice-parent25-hp-batched-benchmark-v1` and
`runs/macbook-lattice-parent25-hp-trial-v2`. Older next-action entries below
are historical.

## An opt-in constrained auditor admits the geometry, but parent25 still misses accuracy, 24 October 2026

The unchanged rejected parent25 mesh now passes a separately versioned
geometry audit. The subsequent cap-budget correction also passes geometry
and FE consistency, but worst disagreement across all 24 models is
2.214326639%, failing 1%. Stop the coarse-bulk, curve-only allocation recipe.
No native candidate, reference or timing advantage is promoted. Native26's
best saved control remains 1.183636899%; the fresh modified actual route
passes at 0.9946706124%.

`geometry(..., box_projection_fallback=True)` opts into
`newton-box-slsqp-v1`. The default historical auditor is unchanged. Only
failed unbounded Newton samples enter SLSQP; every attempted sample remains
in the result. A recovered target must meet the same finite level, normalized
field-residual tolerance of 1e-9, fixed box and original locked coordinates.
The forward distance is measured to that feasible point, a conservative
upper bound rather than a certified closest point. Distance 0.002, volume 1%,
quality 0.1, positive sampled Jacobians, coherent facets and periodic gates stay fixed.
Failed constrained trials remain failed. Ten controls cover a trimmed plane,
locked coordinates, an unreachable level, fixed corners and facet admission.

Job 80 re-audits the exact rejected round13 mesh, prolongation and sample
files. Their hashes are identical. The default still returns two failures;
the opt-in method recovers both. Counts, distance, volume, quality and sampled
Jacobians are unchanged. No refinement or FE occurs in this re-audit, and
job77's rejected result remains frozen. The corrected record explicitly
identifies the changed observation method before continuation is allowed.

Job 81 stops before refining because the checkpoint reader eagerly evaluates
a missing legacy-field fallback. Its exact source and failed log are retained.
The reader and independent replay are corrected; a regression covers a
record with `total_rounds` but no `completed_rounds`, legacy records, and
cumulative counts. Job 82 passes twenty queued controls and continues only
from the explicit same-mesh re-audit, using the same selection rule and caps.

| Continued round | Nodes | Tets | Curve parents bisected |
|---|---:|---:|---:|
| 14 | 177,268 | 95,236 | 7,033 |
| 15 | 196,416 | 103,617 | 8,381 |
| 16 | 218,704 | 113,308 | 9,691 |
| 17 | 244,530 | 124,797 | 11,489 |
| 18 | 249,998 | 127,549 | 2,752 |
| 19 | 249,998 | 127,549 | 0 |

Every accepted parent is curve-vertex incident; no bulk fallback is used.
The implementation finds no further accepted star in the remaining budget.
All fifty-one unbounded projection failures at the final mesh are recovered
by valid box-constrained targets. Distance stays 0.0003819821, quality is
0.1000034664 and minimum sampled Jacobian is 2.75613e-11. The original parent25
domain, approximation space and pin are retained through composed maps.
This domain is distinct from native26; their controls are not nested.

Job 83 completes six fresh loads in 162.95 seconds at 5.715 GiB peak RSS.
All reduced residuals are below 9.79e-13. Independent field, periodic and
positive-rule audits pass. The largest stiffness decrease from parent25 is
0.64684004%; all six signed changes are negative. Galerkin balance is
4.20e-8, and positive-rule energy change is 8.02e-8. These verified corrections
still do not meet the retained target. The shorter saved-mesh solve is not
a matched-accuracy or live timing win.

Next inspect the available cubic25 fields for one shared-domain bulk
approximation probe, before another sizing change. The prepared mesh is
`runs/g3-lattice-cubic25-assembly-memory-v1/mesh.npz`; all six fields are in
`runs/g3-lattice-cubic25-reference-blocked-symmetry-v2/fields.npz`, with audit
`runs/g3-lattice-cubic25-assessment-v1`. First verify exact parent25 embedding,
composed maps, periodic fields, unique-map count and resource feasibility.
Then decide on a bounded positive-rule cross-energy trial without a new
equilibrium solve. Keep every model and distinguish a trial subspace gain
from remaining-error, reference or live-candidate certification. Do not
assume that still deeper curve-only refinement fixes the unresolved bulk.

The fresh actual lattice control remains admitted with its explicit extra
repair/refinement and diagnostic cost of 314.03 seconds. Schwarz retains five repeated
Mac wins against its declared modified actual route. Gyroid and spinodoid
actual-comparator gaps remain. All four cases, fixed specimens, original
thresholds and the full unachieved objective stay in scope.

Records are `runs/macbook-lattice-parent25-box-projection-reaudit-v1`,
`runs/macbook-lattice-parent25-preferred-only-preflight-v3`,
`runs/macbook-lattice-parent25-preferred-only-preflight-v4` and
`runs/macbook-lattice-parent25-preferred-only-v1`. Jobs through 83 are terminal.
Earlier selected-next descriptions are historical.

## Preferred-only refinement improves the control but exposes a projection failure, 24 October 2026

The new native26 fixed-domain control improves worst disagreement across all
24 models to 1.183636899%, from 1.230455074% with a one-round curve-star
priority plus bulk fallback. It still fails 1%. No native candidate, reference
or workflow timing advantage is promoted. The fresh modified actual route
remains admitted at 0.9946706124%, with its separately recorded 314.03-second
diagnostic workflow cost.

The research refiner adds an optional `preferred_only` switch. It filters
complete periodic edge orbits and retains all quality, overlap, count,
geometry and periodic guards. Default and measured-score behavior stay
unchanged; preferred-only mode requires explicit edges and rejects mixed
score ranking. A disjoint-tet control proves that a remote cell remains
untouched despite spare budget, while default refinement still splits it.
All thirteen nested controls and nineteen queued controls pass.

Job 74 repeatedly refines only edges incident on actual wall/cap curve
vertices. Every accepted parent is in that neighborhood; no bulk fallback
is used. Composed maps retain the original native26 domain, old approximation
space and pin 50999. Original geometry and six displacement fields reproduce
within 4.45e-16. Distance stays 0.0003236401, volume changes by 8.88e-16,
quality is 0.10021975, and minimum sampled Jacobian is 1.52083e-9.

| Native26 round | Nodes | Tets | Curve parents bisected | Added nodes |
|---|---:|---:|---:|---:|
| 1 | 213,098 | 135,234 | 6,934 | 13,126 |
| 2 | 227,732 | 143,312 | 8,078 | 14,634 |
| 3 | 244,696 | 152,748 | 9,436 | 16,964 |
| 4 | 249,996 | 155,829 | 3,081 | 5,300 |
| 5 | 249,996 | 155,829 | 0 | 0 |

The implementation accepts no additional star in the remaining budget and
stops. Job 75 completes six fresh loads in 218.97 seconds at 5.832 GiB peak
RSS. Independent field, periodic and positive-rule audits pass. All six
stiffness changes from the original parent are negative, largest 0.08623610%.
Galerkin balance is 3.92e-11 and positive-rule energy change is 9.97e-9.
This supports the allocation change, but the target remains unmet.

The already audited nonadaptive parent25 offers comparable geometry accuracy
at 75,113 initial nodes, rather than native26's 199,972. Its distance is
0.0003819821 and initial quality is 0.322484. It is a distinct fixed domain;
no nested relation or monotone energy comparison to native26 is claimed.
One bounded preflight uses the identical preferred-only rule and 250,000
node/tet caps, without a new size parameter or seed sweep.

Job 76 reaches 108,850 nodes and 64,716 tets after eight rounds. All gates pass,
but the declared round limit stops work while stars are still being accepted.
This is not convergence. Job 77 continues that verified geometric checkpoint
for one more bounded eight-round block. Original source fields and composed
maps are retained; no equilibrium solve or solved prefix is replayed.
Rounds 9–12 pass. Round 13 stops at 161,494 nodes and 88,203 tets because two
forward projections fail. Distance, volume, quality, sampled Jacobians,
coherent facets and periodicity otherwise pass. The mesh remains rejected;
no FE is selected, and no earlier round is substituted to hide the failure.

Job 78 independently replays exactly those two forward samples on face 17087.
Unbounded Newton feet leave the box by 4.00488e-6 and 2.35737e-6. Box-constrained
SLSQP trials find valid zero-level targets at distances 2.68001e-5 and
2.41392e-5, with field residuals at roundoff. A trimmed-plane control proves
that a constrained target can exist when the unbounded normal foot leaves
the box. These feasible points provide distance upper bounds, not certified
closest points. The original rejection, samples and mesh remain unchanged.

Next validate an explicitly versioned, opt-in box-constrained projection
fallback on the same samples. Preserve the default historical auditor,
0.002 distance, 1% volume, q0.1, positive-Jacobian and periodic gates. Keep
failed samples in the calculation; a constrained candidate must satisfy the
same target level and box constraints. Re-audit the unchanged rejected mesh
before any further refinement or solve. Do not silently reclassify the old
failure. This is an observation-method investigation, not permission to
waive geometry admission or fit a mesh to empirical stiffness models.

Schwarz retains its repeated Mac wins against the declared modified actual
route. Gyroid and spinodoid comparator gaps remain open. Lattice native
accuracy and reference convergence remain unresolved. All four cases and
the unchanged objective remain active and unachieved.

Records are `runs/macbook-lattice-adapt26-preferred-only-preflight-v1`,
`runs/macbook-lattice-adapt26-preferred-only-v1`,
`runs/macbook-lattice-parent25-preferred-only-preflight-v1`,
`runs/macbook-lattice-parent25-preferred-only-preflight-v2` and
`runs/macbook-lattice-parent25-projection-failures-v1`.
Jobs through 78 are terminal. Earlier selected-next descriptions are historical.

## Fresh Mac actual lattice workflow passes, 24 October 2026

The declared actual lattice route now passes fresh Mac generation and all
24 retained stiffness models. Worst disagreement is 0.9946706124%, negative
against sensitivity model 7. This closes the frozen-seed reproduction gap
for one explicitly modified route. It does not establish a lattice Meshers
win: the native star-refined control still fails at 1.230455074%.

Job 65 regenerates the fixed smooth-union specimen through actual VTK51 and
two MMG calls. First pass uses nofem, hausd0.0015, hmax0.06 and hgrad1.3;
second pass uses hausd0.001, hmax0.06, hgrad1.3, ls and nr. The VTK seed has
51,309 nodes and 249,432 tets. First MMG has 24,936 nodes and 115,536 tets;
second has 33,953 nodes and 158,430 tets. Every intermediate count passes.
Full intermediate geometry diagnostics are saved and their cost is charged.

Fresh final quality is 0.0998809921, so the unrepaired Mac output is rejected.
Distance 0.0017012390, volume error -0.4905337%, positive Jacobians, coherent
facets and periodicity pass. Job 66 identifies one bad cell with one free
interior vertex and no fixed-face quality obstruction. No FE precedes repair.
Job 67 adds the existing guarded SLSQP interior relocation, explicitly enabled
only for this exact lattice recipe. Eight controls pass, including rejection
of a changed recipe or failed geometry. One vertex moves 1.43275e-5; complete
star quality reaches 0.100001001. Boundary coordinates and connectivity stay
exactly fixed, and signed volume changes by 3.70e-17. The initial failure is
retained rather than rounded into admission.

Job 68 counts straight Tet10 elevation before allocation: 238,729 nodes.
Shared manual elevation is factored from the existing control without changing
its midpoint construction or default behavior. A shared-edge/node-cap control
and twelve nested controls pass. The declared curve-vertex-star priority
then produces 250,000 nodes and 165,722 tets. Parent geometry and six affine
control fields reproduce exactly within recorded tolerances. These controls
are not equilibrium fields. No unnecessary parent FE solve is performed.
Distance and volume are preserved; q0.100001001, positive sampled Jacobians,
periodic nodes/caps and coherent facets pass under unchanged admission limits.

Job 69 performs six fresh 15-point loads with the same cached block-AMG preset,
new fields and a new hierarchy. Each load has a 240-second continuation gate,
the stage has a 1,200-second timeout, and the RSS watchdog stays at 8 GiB.
Solve plus final field archive takes 218.55733 seconds; peak worker RSS is
6.10088 GiB. Independent job 70 verifies reactions within 8.74e-12,
Hill-Mandel energy within 2.66e-12, and zero periodic jumps. Positive
64-to-125-point energy change is 1.79e-14. Job 71 retains every original
model and passes the unchanged 1% condition at 0.9946706124%.

| Charged fresh Mac stage | Seconds |
|---|---:|
| VTK/MMG generation and intermediate/final geometry audits | 15.290462 |
| Quality localization | 0.811610 |
| Explicit interior repair and admission | 3.422781 |
| Straight elevation, star refinement and admission | 12.487351 |
| Six-load solve and final field archive | 218.557328 |
| Independent field and positive-quadrature audit | 63.405356 |
| All-model assessment | 0.053772 |
| Total | 314.028660 |

Job 72 preserves this full cost receipt. The job 71 named-stage subtotal is
313.217050 seconds; the supplemental receipt also charges quality localization.
Timers retain diagnostics and serialization within their named boundaries.
Synthetic controls, untimed provenance checks, bookkeeping, queue and inspection
gaps are excluded. Intermediate audit overhead is deliberate, and
`timing_comparable` remains false. This is one fresh diagnostic workflow,
not production/repeated paired timing, the default Microgen pipeline, MMG on
a Meshers seed, or a continuum/local-stress certificate. No reference is promoted.

The next Meshers action is a bounded preferred-star-only refinement budget
preflight on the same native lattice26 domain. Its prior priority pass
bisects 6,934 curve-vertex parent cells, then its fallback bisects 23,992
other parents: 77.5787% of all 30,926 bisections are outside curve-vertex
cells. This is a cell-split census, not exact node-cost attribution. Inspect
node allocation and test repeated preferred-star refinement with bulk fallback
suppressed, preserving the parent space and every geometry, quality and
periodic guard. Use the unchanged 250,000-node/tet caps. No new sizing sweep,
model-fitting geometry offset, or solve before reading the preflight is selected.

Schwarz retains its five repeated Mac modified-actual-route workflow wins.
Gyroid and spinodoid still lack admitted actual comparators in tested routes.
Lattice native accuracy and reference convergence remain unresolved. All four
cases remain in scope, and the full objective stays active and unachieved.

Records are `runs/macbook-lattice-actual-vtk51-stage-audit-v1`,
`runs/macbook-lattice-actual-vtk51-quality-localization-v1`,
`runs/macbook-lattice-actual-vtk51-interior-repair-v1`,
`runs/macbook-lattice-actual-live-preflight-v1`,
`runs/macbook-lattice-actual-live-solve-v1`,
`runs/macbook-lattice-actual-live-field-audit-v1`,
`runs/macbook-lattice-actual-live-assessment-v1` and
`runs/macbook-lattice-actual-live-costs-v1`. Jobs through 72 are terminal.
Earlier selected-next descriptions are historical.

## A modified actual lattice control passes, 24 October 2026

One saved actual Microgen-derived lattice control now passes all 24 retained
models at 0.988497945% worst disagreement. The new native curve-star control
still fails at 1.230455074%. These are fresh Mac solves of saved meshes, not
live regeneration or a matched-accuracy timing comparison. No lattice speed
win, new reference, or continuum/local-stress certificate is claimed.

Jobs 58–59 localize the complete ranked-space increment between the two
previously audited coupled solutions. The own-rule Galerkin identity error is
7.64e-14, and positive 125-to-216-point local-trace change is 8.08e-8.
Curve-vertex incident cells carry 66.8117% of the correction in 2.4994% of
cells. Literal-curve-edge cells carry 23.3435% in 0.4557% of cells. Groups
overlap. The top 1% of cells carry 73.2238%. This is resolved discrete
correction, not an estimate of remaining error.

The first diagnostic finishes all three FE rules, then fails its wall
projection assertion. Its source, arrays and log remain frozen. The second
run reuses those hash-identical completed arrays and reports the projection
failure explicitly. Eleven Newton projections on six faces leave the box,
by up to 6.76e-5. Surface-rule stability passes, but complete trimmed-wall
geometry interpretation does not. Neither the original geometry admission
nor its samples are silently changed. No co-location conclusion follows.

The optional `--prefer-box-curve-stars` preflight mode prioritizes all edges
incident on actual wall/cap curve vertices, including radial edges, before
the existing descending-length fallback. Complete periodic orbits, disjoint
stars, exact parent maps, quality 0.1, positive sampled Jacobians and 250,000-node/tet
limits remain unchanged. It is mutually exclusive with literal-edge priority;
the default selection is unchanged. A radial-edge selection control and all
18 queued controls pass. No size parameter sweep or oracle field is used.

| Saved-mesh control | Nodes | Tets | Worst error across all 24 models | 1% target |
|---|---:|---:|---:|---|
| Original live native lattice26 | 199,972 | 128,300 | 1.2686985% | Fail |
| Earlier literal-curve-edge nested native | 249,999 | 158,967 | 1.2556246% | Fail |
| New curve-vertex-star nested native | 250,000 | 159,226 | 1.2304551% | Fail |
| Original tuned actual Microgen straight Tet10 | 238,411 | 158,252 | 1.0043541% | Fail |
| New star-refined actual Microgen-derived control | 249,998 | 165,772 | 0.9884979% | Pass |

Jobs 60–61 retain the original native lattice26 domain and parent field to
3.34e-16. They bisect 78.1472% of its curve-vertex incident cells. Distance
stays 0.0003236401; quality is 0.1023041 with positive sampled Jacobians.
All six stiffness changes are negative, largest 0.03868354%. Galerkin balance
is 3.71e-10 and positive-rule energy change is 9.64e-9. Six fresh loads take
212.13 seconds at 5.758 GiB. Stop this allocation recipe without promotion.

The actual control's initial worst signed error is +1.0043541%; its most
negative error across all models is -0.8937705%. That motivates one bounded
refinement test, not rounding the original failure into admission. Jobs
62–63 retain its affine domain exactly. Distance stays 0.0016215853, volume
change is zero, quality is 0.1015540, and sampled Jacobians are positive.
They bisect 73.6880% of curve-vertex incident cells. Six fresh loads take
232.18 seconds at 6.198 GiB. Galerkin balance is 3.86e-14, the largest
stiffness decrease is 0.09635235%, and positive-rule energy change is
1.64e-14. Full field and periodic audits pass. The worst final error is
negative against sensitivity model 7. More FE refinement is not automatically
helpful to this empirical target.

The passing control retains actual tuned VTK51/two-pass MMG provenance,
manual straight Tet10 elevation/direct fedoo, and an explicitly added
periodic star-refinement step. It is not the default Microgen pipeline or
MMG on a Meshers seed. All 24 original models and original failures remain.
The full parent23 reference investigation remains unresolved at 1.7413203%.

Next prepare fresh Mac actual VTK51/two-pass MMG generation with this declared
added refinement and charge every stage. Read geometry and counts before a
new solve; the frozen G3 seed does not prove fresh Mac MMG reproduction.
The native control still needs improvement before a matched-target pair.
Schwarz retains the five repeated Mac workflow wins. Gyroid and spinodoid
still lack admitted actual comparators in tested routes. All four cases and
the unchanged objective remain active.

Records are `runs/macbook-lattice-hp-ranked-localization-v1`,
`runs/macbook-lattice-hp-ranked-localization-v2`,
`runs/macbook-lattice-adapt26-curve-star-preflight-v1`,
`runs/macbook-lattice-adapt26-curve-star-v1`,
`runs/macbook-lattice-actual-microgen-curve-star-preflight-v1` and
`runs/macbook-lattice-actual-microgen-curve-star-v1`.
Jobs through 63 are terminal. Earlier selected-next descriptions are historical.

## Ranked nested coupled correction, 24 October 2026

Mac jobs 53–56 complete the measured-score nested-space investigation.
The 249,999-node, 162,599-tet mesh retains the original parent23 curved
physical domain and the previous approximation space. Root-map geometry
reproduction is within 1.67e-16. Fifteen versus positive 216-point integration
changes the mode Gram matrix by 0.020572% and the fixed-lift trial gain by
0.025114%, below the declared 1% precision gate.

The prolonged quadratic base is not equilibrated on the refined space.
Its six reduced residuals are 19.29–31.69% of their internal-force norms.
The generalized correction therefore includes six fresh lift-residual
responses as well as six fresh mode responses. A dense full augmented-system
control verifies the condensation and rejects dropping the lift residual.
Combined AMG setup fits at 5.584 GiB. All twelve responses finish in 467.17
seconds at 5.769 GiB, under the unchanged 8 GiB ceiling and explicit timeouts.
These are reference-research auxiliary responses, not twelve macro workflows.

The complete modal field audit takes 174.96 seconds at 5.580 GiB. Own-rule
energy agrees with the predicted tensor within 3.96e-14, and stress-energy
balance is 2.46e-14 relative. Positive 125-to-216-point energy change is
4.37e-9. Fresh operator replay checks all twelve auxiliary responses below
5.56e-12, all six complete equilibrium states below 1.66e-12, and mode
stationarity within 2.56e-12. Periodic jumps are below 4.45e-16 and the pin
is exactly zero. All six signed stiffness changes from the preceding coupled
space are negative; the largest decrease is 0.09312825%.

Worst disagreement across every retained model improves from 1.82618412%
to 1.74132028%. It still fails 1%. No model is deleted, reference promoted,
or native candidate promoted. This controlled fixed-domain result is not a
live Meshers regeneration, same-host competitor timing, or continuum/local
stress certificate. The best live lattice candidate remains the audited
adapt26 result at 1.2686985%. The reference problem remains open.

Read the remaining correction and geometry evidence before selecting a
native change. Do not repeat the exhausted literal-curve-edge sizing recipe
or assume this parent23 investigation resolves geometry bias. Gyroid and
spinodoid still lack admitted actual comparators in tested routes. The five
fresh Mac Schwarz pairs retain their explicitly modified actual-route win.
All four cases and their original admission gates remain in scope.

Records are `runs/macbook-lattice-hp-ranked-coupling-v1`,
`runs/macbook-lattice-hp-ranked-setup-v1`,
`runs/macbook-lattice-hp-ranked-solve-v1` and
`runs/macbook-lattice-hp-ranked-audit-v1`. Jobs through 56 are terminal.
Earlier selected-next descriptions are historical.

## Measured-score lattice allocation preflight, 24 October 2026

Job 49 localizes the complete coupled correction, using all six loads and
two declared mixed directions. Every element tensor contributes to the full
Galerkin identity; the 15-point identity error is 1.15e-13 relative to the
baseline. Positive 125-to-216-point local-vector changes are 1.066e-5 for
the six-load normalized trace, 1.364e-5 for the largest-correction direction,
and 1.019e-5 for the worst retained-model direction. All pass the 1% precision
gate. The worker takes 155.91 seconds at 0.622 GiB.

| Positive 216-point localization | Element fraction | Six-load correction fraction |
|---|---:|---:|
| Incident on a curve vertex | 1.98379% | 59.10487% |
| Containing a curve edge | 0.43006% | 17.30880% |
| Top 1% ranked by correction | about 1% | 54.1449% |
| Top 5% ranked by correction | about 5% | 87.4500% |

Geometry groups overlap; these fractions are not a partition. The trace is
normalized by the complete baseline stiffness matrix. Localization also
retains the worst correction and worst model mixed-strain witnesses. This
is resolved discrete correction, not certified remaining continuum error.
It does not justify repeating the earlier literal-curve-edge-only recipe.

The nested-refinement research helper now optionally ranks complete periodic
stars by measured element-score sum per added-node upper bound. Its default
ordering and every quality, overlap, count, geometry and periodic guard stay
unchanged. A disjoint-star budget regression proves that the measured score
can select a shorter high-score star instead of the longest edge. All eleven
nested-refinement controls pass. No product Meshers API is added.

Job 50 independently replays both positive-rule localization arrays and all
geometry-group shares, then uses the six-load normalized trace for one ranked
preflight. It bisects 21,455 parent cells through 4,723 accepted periodic
orbits. The selected 15.2008% of parents carry 88.80881% of the measured
score. Final counts are 249,999 nodes and 162,599 tets, within unchanged caps.
Quality is 0.10003495, sampled Jacobian minimum is 2.47990e-7, and distance
remains 0.0006492808. Volume changes by 6.67e-16 relative. Geometry and the
quadratic base reproduce within 4.45e-16. Periodicity, coherent facets and
the original pin node 4847 are preserved. The unchanged cubic function span
and included old quadratic space retain the old augmented space.

Job 51 captures the actual reduced operator without solving. It has 740,964
mesh unknowns and 58,642,444 nonzeros, with 706,673,188 CSR bytes. Peak RSS is
5.394 GiB under the unchanged 8 GiB watchdog. This is mesh-operator assembly,
not proof that coupling columns and AMG together fit. Job 52 independently
replays all child nodal geometry and field values, score coverage, sampled
geometry and periodic base fields from saved arrays. Maximum nodal replay
error is 1.12e-16. All new meshes and arrays remain saved, with hashes.

Next is one coupled-force and AMG resource preflight on this ranked space.
The saved prolongation is the old coupled field's quadratic base, not a new
Tet10 equilibrium solution. Account for its six-load residual through a
general coupled correction or fresh macro solves; do not reuse the earlier
formula that assumes an equilibrated base. Retain the complete modal field,
fixed cubic span, original pin, 8 GiB and 250,000-node/tet limits. Read those
checks before choosing any solve. No new accuracy result or reference is
promoted; the last audited coupled disagreement remains 1.826184% across
all 24 models. Gyroid and spinodoid actual-comparator gaps remain open, and
the repeated Mac Schwarz result stands. The full four-case goal is active.


## Audited coupled lattice correction, 24 October 2026

Jobs 41 to 43 establish a consistent quadrature and resource path before any
coupled solve. The 15-point mode Gram matrix differs from the positive
216-point result by 0.0478437%; its trial gain differs by 0.0572437%, both
below the fixed 1% precision gate. Assembled mode forces satisfy direct
virtual work within 5.03e-16 and annihilate rigid translations. Full and
reduced coupling storage is 31,411,440 and 31,009,392 bytes.

Reduced assembly has 646,029 unknowns and 51,030,214 nonzeros. The six global
scalar modes give 646,035 augmented unknowns, with the same 218,135 physical
nodes and 141,144 tets. The sparse mesh operator stays separate from those
six columns. Assembly peaks at 6.044 GiB; setup of the validated four-level
block AMG peaks at 5.359 GiB and stops before CG. The 8 GiB ceiling and all
250,000-node/tet admission caps stay unchanged.

Job 44 then solves six fresh auxiliary mode responses with a 240-second
limit per response and a 1,800-second worker timeout. It uses saved baseline
fields as explicit reference inputs; it is not six fresh macro workflows.
Every baseline load is checked against the assembled operator. The positive
six-by-six Schur system gives a 0.3254513% own-rule energy reduction from the
nested baseline. The worker takes 177.49 seconds at 5.991 GiB. No dense
augmented graph is passed to AMG, and no solved prefix is replayed.

The complete displacement representation combines quadratic response fields
with the frozen cubic-minus-nested modes. It is saved as `hp-fields.npz`;
ordinary nodal Tet10 interpolation would lose its cubic part. Job 45 checks
that representation independently, including the constant pin gauge:

| Check | Result |
|---|---:|
| Own-rule energy versus Schur prediction | 1.23e-13 maximum mixed-load difference |
| Own-rule stress/energy relative difference | 2.14e-12 |
| Positive 216-point energy shift | 3.003e-6 |
| Positive 125-to-216-point energy difference | 6.001e-9 |
| Freshly reassembled auxiliary reduced residual | at most 1.075e-11 |
| Freshly reassembled complete-field reduced residual | at most 2.038e-10 |
| Mode-equilibrium relative residual | 1.594e-12 |
| Six-load periodic jump | at most 4.45e-16 |
| Pin displacement | 0 |

The full-field audit takes 139.39 seconds at 5.720 GiB. Every sampled
Jacobian remains positive. All 24 original model tensors are retained.
Worst disagreement is 1.82618412%, so the 1% target still fails. This is an
audited equilibrium in the stated augmented space on the frozen approximate
parent23 domain. It is not a promoted reference, a best native Meshers
candidate, a live workflow timing, or continuum/local-stress certification.

Job 46 records a provenance-verifier failure: its old source lookup tried
a commit predating the h/p helper. Job 47 verifies the exact committed
source bytes. Job 48 generalizes that lookup to exact same-path Python blobs
on the current branch and records the verifier hash and repository head.
A regression checks historical lookup after HEAD changes and rejects wrong
bytes. Data hashes never use this fallback. No solve or field audit is rerun.

Next, localize the resolved complete-field correction energy, using all six
loads and the full Galerkin energy identity. Check positive-rule precision
and concentration on actual geometry groups before selecting any new
refinement or Meshers change. This is resolved correction, not a certified
estimate of remaining continuum error. Do not revive curve-only sorting,
raise caps, tune an offset or add another polynomial degree from this result.
Gyroid and spinodoid comparator coverage remain open; the verified repeated
Mac Schwarz comparison stands. The full four-fixture goal remains active.


## Shared-domain lattice field combination, 24 October 2026

Job 37 proves that the saved second nested Tet10 level and cubic23 fields
use the same original curved parent23 physical domain. Both solve-to-mesh
provenance links and eight transferred NPZ hashes pass. All ten child Tet10
nodes reproduce the composed quadratic map, and all twenty cubic nodes
reproduce the original quadratic geometry, within 1.67e-16. Reference-volume
fractions sum to one for every original parent. The two prolongation checks
cover all 76,947 and 141,144 children, reproducing fields within 3.34e-16.
Periodic caps/nodes and all six affine field jumps pass.

The nested mesh has 218,135 nodes and 141,144 tets; the cubic mesh has
218,508 nodes and 43,844 tets. Neither mesh, physical geometry, field array,
fixture, quality threshold or count limit changes. The original cubic mesh
is read from its hash-identical prepared directory, not regenerated.

Job 38 integrates all 144 cross terms between the twelve saved load fields.
It optimizes six global periodic cubic-minus-nested difference fields. Source
coefficient arrays remain immutable; mode weights vary, without independent
nodal re-equilibration. Two controls check exact mixed polynomial integrals
under nonidentity child maps and cancellation of manufactured opposing errors.

| Positive rule | Maximum relative trial-energy reduction | Direct energy-identity error |
|---|---:|---:|
| 125 points | 0.259381777% | 1.31e-16 |
| 216 points | 0.259381631% | 1.98e-16 |

The gain-matrix rule change is 6.71442e-7 relative, below the existing 1%
precision gate. The gain exceeds the existing 0.05% signal screen. Own-field
energy differences from solved tensors stay below 7.69e-6; minimum sampled
Jacobian is 4.94816e-7. The worker takes 90.98 seconds at 0.472 GiB, under
its 600-second timeout and unchanged 8 GiB watchdog. These are diagnostic
costs, not live meshing-workflow timings.

Job 39 independently solves the 18-by-18 constrained-basis system instead
of reusing the six-mode formula. Both rule results reproduce, with saddle
residual below 1.44e-16. Every original lattice model is retained. Trial
worst disagreement is 1.89983697% across all 24, still failing 1%. This is
not a new equilibrium tensor, native candidate, reference or continuum bound.

The resolved gain justifies one coupled-space feasibility preflight. First
check a consistently assembled quadrature rule against the positive-rule
screen, then measure coupling storage and reduced assembly under the existing
8 GiB and 250,000-node/tet caps. Record the six extra scalar mode unknowns
explicitly. Do not apply AMG blindly to globally dense mode rows. No full
coupled solve follows before those checks pass. This uses already solved
spaces, with no higher degree, offset or speculative sizing sweep.
Gyroid and spinodoid comparator coverage remain open; the repeated Mac
Schwarz result stands. The full four-fixture objective remains active.


## Verified repeated Mac Schwarz comparison, 24 October 2026

Job 33 completes one excluded warmup and five fresh interleaved pairs, with
alternating arm order. Every pair passes all eight declared full mixed-load
Schwarz models at the supplementary 0.1% target. Native errors span
0.031469% to 0.093778%; Microgen errors span 0.053480% to 0.062162%.
All geometry, quality, sampled Jacobian, coherent-facet and periodic gates pass.

| Charged cost | Modified actual Microgen median | Meshers median | Five paired ratios |
|---|---:|---:|---:|
| Science-audited workflow | 131.85894 s | 19.14923 s | 6.7632 to 6.9203 |
| Generation, repair/lift and solve | 119.40758 s | 17.48343 s | 6.6949 to 6.8605 |

Workflow ranges are 130.16436 to 132.30789 seconds for Microgen and 19.11868
to 19.54234 seconds for Meshers. Meshers is faster in every measured pair
under both cost definitions. No pair is discarded, and warmup is excluded.

The comparator is actual VTK30/MMG with explicit one-pass overrides, a
charged fixed-boundary SciPy interior repair, and manual straight Tet10/direct
fedoo. It is not Microgen's default or public Tet10 conversion. Meshers uses
fresh nonadaptive curved16. Both use the same block-AMG preset and one
numerical thread. Stage limits stay at 250,000 nodes/tets and 8 GiB; the
whole study takes 980.56 seconds within its one-hour deadline. No previous
mesh, solved prefix or preconditioner is replayed.

All six pairs pass 42 controls. Five verifier unit controls also pass.
Job 34 stops before geometry replay because the verifier omitted support-path
initialization; its exact failed driver and log remain preserved. Job 35
corrects initialization and independently verifies 262 hash records, all
96 full-matrix model comparisons, every timing statistic, and all twelve
saved mesh/field sets. Geometry and periodic field jumps are recomputed
from actual arrays. Minimum quality is 0.100001001 for Microgen and 0.4499953
for Meshers. Maximum reaction/energy error is 7.60e-12; periodic jumps are
at most 5.56e-17. The timing study is not rerun after the verifier failure.

This establishes the repeated Mac Schwarz result for the declared modified
route. It does not resolve the other three fixtures or certify continuum
or local-stress accuracy. Keep all 24 lattice models and all fixed specimens.
The next reference investigation uses the already solved nested and cubic
fields on their original parent23 domain; first prove the physical mappings
and field provenance agree before computing a joint cross-energy diagnostic.
No new solve, polynomial degree, geometry offset or sizing sweep is selected.
Earlier selected-next descriptions below are historical.


## First verified Mac Schwarz pair, 24 October 2026

Job 29 completes one fresh same-host pair after the separately preserved
quality rejection. The actual Microgen route is VTK30 plus one MMG call
with `nofem`, hausd 0.0005, hmax 0.1 and hgrad 1.3, explicit omission of the
second call, one fixed-boundary SciPy interior move, and manual straight
Tet10/direct fedoo. The repair is charged. It is not a Microgen default.
Meshers freshly generates curved16 with one optimization pass and no adaptation.
Both use the same explicit six-load block-AMG preset, one numerical thread,
250,000-node/tet limits and the 8 GiB process-tree RSS watchdog.

Both pass geometry, quality, positive Jacobians, coherent facets and periodic
caps/nodes. Full mixed-load error against all eight declared Schwarz16/20
models is 0.053480–0.062162% for Microgen and 0.031469–0.093778% for Meshers.
Both pass the supplementary 0.1% target. Independent reactions and energy
agree within 7.60e-12, with periodic jumps at most 5.56e-17.

Charged science-audited workflow time is 132.63564 seconds for Microgen and
19.27508 seconds for Meshers, a single-pair ratio of 6.8812. The generation,
lift/repair and solve subtotal is 119.88952 versus 17.60300 seconds, ratio
6.8107. Peak RSS is 4.910 GiB for the Microgen solve and 2.244 GiB for Meshers.
This is one pair with no excluded warmup, not a repeated timing claim.

Job 30 independently replays saved artifacts, full model errors, charged
workflow totals, the earlier geometry witnesses and unchanged repaired
boundary/connectivity. Job 31 reintegrates the fresh native six-load fields
with positive 64- and 125-point rules. Maximum energy shift from the solve
is 2.88454e-9; successive-rule shift is 5.80694e-10. These are fixed-field
checks, not new equilibrium solves or continuum certification.

The selected next experiment is one bounded repeated Mac study: use this
verified pilot and quadrature record, exclude one warmup per arm, then run
five fresh interleaved pairs with alternating arm order and the same explicit
charged repair. Port the repeated runner's historical G3 evidence lookup
without rewriting frozen records. Every pair must pass before its successor.
Gyroid and spinodoid still lack admitted actual comparators in tested routes;
lattice remains at 1.2556246% against all 24 retained models. No new reference
or four-fixture victory is promoted.


## Spinodoid geometry and fresh Mac comparison, 24 October 2026

Job 24 replays the fixed spinodoid samples against the coherent sampling36
seed. All 113 reverse samples above 0.002 are wall samples; none is a cap
sample. The worst 32 nearest faces are wall triangles, including witnesses
well inside the box. Reverse maximum is 0.0044944787. Closest-point replay
agrees with the existing VTK distance evaluator and a manufactured plane
control. Forward Newton failures remain separate diagnostics.

Job 25 tests hypothetical damped wall-node projection on the exact frozen
field. All 22,602 periodic representatives converge below 1e-9 residual
distance, with box coordinates locked and periodic copies consistent.
Reverse maximum improves to 0.0028468542, but 27 samples still fail and 146
tets invert. No projected mesh is applied or admitted. Stop this simple
projection recipe as well as the previous two-orbit cap search.

Job 26 attempts the recorded Schwarz pair afresh on the Mac. The Microgen
parent passes distance, volume, facets and periodicity, but quality 0.0956791
rejects it before FE. Historical G3 wins remain historical. Job 27 localizes
one failing tet, with one free interior vertex and no fixed-face obstruction.
Job 28 reuses the existing complete-star relocation control: one move of
0.000748602 gives global quality 0.100001001 while boundary coordinates,
connectivity and volume stay unchanged. All geometry gates pass. This is
saved-mesh feasibility, not a live workflow or a timing win.

Job 29 selects one fresh paired workflow with that explicit SciPy repair
charged to Microgen, followed by manual straight Tet10 and direct fedoo.
Keep the existing supplementary 0.1% eight-model mixed-load test, all original
admission gates, and independent six-load field audits. Read its actual
result before deciding on repetitions or another experiment. Gyroid and
spinodoid comparator coverage and lattice reference convergence remain open.


## Independent two-orbit cap search, 24 October 2026

Mac job 21 stops at its manufactured fixture: the optimized cube does not
retain the assumed quarter-grid node. Its exact test and failed log are
preserved. Job 22 selects the actual neighboring periodic orbit, passes both
analytic-derivative and known-repair controls, then tests four independent
planar unknowns on the spinodoid. All eight deterministic starts miss the
unchanged 0.1 cap gate. Direct replay finds minimum cap quality about
0.00503593, positive incident tetrahedra, and exactly matching periodic
translations across all candidates. The patch has 24 cap triangles and 51
incident tets. No move, MMG call, or FE solve follows.

This stops this two-orbit search recipe; nonlinear search failure does not
prove infeasibility. Job 23 independently replays every returned candidate
from saved coordinates, all transferred hashes, run inputs, and all 72 lattice
comparisons. No reference or comparator is promoted. Before expanding cap
motion again, localize the coherent spinodoid seed's independent reverse
distance failure (0.00449448 versus 0.002) and unresolved forward projections.
Use that evidence to decide whether the extraction geometry is worth repairing;
do not spend another cap sweep on a source that still fails geometry.
All four cases remain in scope, and fresh same-host workflows are still needed.


## Fixed-spinodoid comparator preflights, 24 October 2026

The shared experiment fixture loader now restores all 52 frozen Fourier
coefficients exactly, after a regression catches Mac FFT roundoff up to
8.88e-16. A paired actual Microgen extraction at sampling36 preserves those
coefficients and identical sampled inputs. Five-tet preclipping removes all
282 source face conflicts, but the coherent 180,277-tet seed still fails
distance, quality, and retained-cap feasibility. No MMG or FE runs.

Single flips cover only 58 of 164 blocked caps. Fixed polygon-corner bounds
find no obstruction. One complete periodic cap-interior star has a checked
fixed-neighbor obstruction, independently replayed at the original 0.1 gate
with 80-digit arithmetic. Its active support includes another free cap vertex,
so coupled motion is not excluded. A whole-cap single-vertex sweep is not
selected. One closed active patch preflight is next.

Twenty-nine focused controls pass. The failed initial manufactured control
is retained. Independent replay verifies both saved seeds, 210 input hashes,
all 3,509 transferred files, and 72 earlier lattice model comparisons.
The four-case goal remains active. Details are in `MICROGEN_RESULTS.md`.

## Curve-priority nested lattice26 remains rejected, 24 October 2026

Prioritizing actual wall/cap curve edge orbits accepts 544 of 549 preferred
orbits while retaining every guard. The mesh fits 249,999 nodes and 158,967
tets. Geometry and old fields reproduce within 3.34e-16, and all admission
checks pass. Six fresh loads and independent field, positive-rule, and
Galerkin checks pass with verified single-thread Accelerate. Peak solve RSS
is 5.472 GiB, below the unchanged 8 GiB watchdog.

The largest FE correction is only -0.0132094%, and full energy balance agrees
within 4.82986e-10. Worst disagreement across all 24 models is 1.2556246%,
still failing 1%. This is little gain over longest-edge refinement at
1.2565145%. Stop these one-stage bisection recipes. No native product feature,
reference, continuum bound, or workflow advantage is promoted.

Independent checkpoint verification passes 3,509 transferred hashes, 150
Mac input hashes, and 72 full model comparisons. The four-case objective
remains active; missing admitted actual comparators and lattice reference
uncertainty remain unresolved. See the current Mac checkpoint.

## Lattice26 bulk refinement remains rejected, 24 October 2026

The fixed-domain longest-edge refinement reaches exactly 250,000 nodes and
160,291 Tet10 elements. Geometry and old-field reproduction, periodic caps,
sampled Jacobians, oriented facets, six fresh loads, and independent audits
pass. The largest stiffness correction is only -0.0122671%, and the full
Galerkin balance agrees within 6.07082e-10. Worst retained-model disagreement
is 1.2565145%, failing 1%. Stop this recipe without a reference promotion or
workflow claim.

The initial Mac numerical timings lack explicit Apple Accelerate threading
control and are excluded from verified one-thread comparison baselines.
Numerical records remain frozen. The corrected native mode is now set and
queried as single threaded. New 55-test and twelve-load controls pass.
One preflight prioritizing actual wall/box curve edges is selected from the
earlier localization evidence and the negative bulk test. Read the current
Mac checkpoint before any solve.

## Native Mac second nested level, 24 October 2026

Native setup passes 54 focused tests and twelve fresh AMG control loads.
The second parent23 nested level preserves the curved domain and old fields
within 3.33e-16, at 218,135 nodes and 141,144 Tet10 elements. All geometry,
quality, sampled Jacobian, periodic, and oriented-facet gates pass.

Six fresh loads complete in a 191.45-second solve worker at 4.508 GiB peak
RSS. Independent fields and positive quadrature pass. The largest additional
signed stiffness correction is -0.328272%, and the Galerkin balance error is
3.78355e-8. Worst disagreement across all 24 unchanged models is 2.135117%,
failing 1%. No reference or timing advantage is promoted. Two nonuniform
nested corrections do not prove continuum convergence.

The next selected experiment is one fixed-domain nested geometry preflight
on the best saved native lattice26 candidate, under unchanged caps. Details
and native resource controls are in `handoff/macbook/CHECKPOINT-20261024.md`.

## MacBook handoff after nested correction, 2026-10-24

Job 473 completes the lattice26 field and positive-quadrature audits, with
reaction error 1.0611e-11 and 125-point energy shift 4.0667e-8. The candidate
still fails the retained 1% target at 1.2686985%.

Jobs 475–476 establish one conforming nested refinement of parent23 on its
fixed quadratic physical domain. Ten controls and independent admission pass.
The refined space has 123,514 nodes and 76,947 Tet10 cells. Six fresh loads
and saved-field/positive-rule audits pass; maximum parent-to-refined energy
correction is -0.499513% and the Galerkin balance agrees within 6.6756e-8.
Worst retained-model disagreement is 2.4480846%, so no reference is promoted.
G3's final job is terminal success. Kevin requests primary continuation on
MacBook Pro with Sol 6.1 medium. Read `handoff/macbook/HANDOFF.md` for data,
native environment setup and the next bounded preflight.

## Paused after completed lattice26 workflow, 2026-10-20

Job 379 completes all six fresh AMG loads on a geometry-admitted curved
Tet10 mesh with 199,972 nodes and 128,300 elements. Charged workflow is
736.51 seconds; bounded worker peak RSS is 5,781,733,376 bytes. Preliminary
full-matrix comparison against the unchanged 24 models gives 1.2686985%
worst disagreement, above the fixed 1% target. Independent field and
higher-quadrature audits remain pending. Geometry qualification alone is
not stiffness verification. No successor was queued after the user's pause
request, and the four-fixture objective remains unachieved. See the root
`PROJECT_TRACKING.md` and `CUBIC_RESULTS.md` for the final checkpoint.

## Lattice target feasibility and finer preflight, 2026-10-20

Assessment 377 retains all 24 model tensors and checks full-matrix target
feasibility. Their necessary pairwise error lower bound is 0.9019313%.
A harmonic trial matrix passes every model at 0.9084715%, so the 1% condition
is algebraically feasible. This trial is not a mesh or reference. No admission
gate, model, or frozen result changes.

Current-native lattice26 preflight 378 passes geometry at 199,972 nodes and
128,300 tets, distance 0.0003236401 and quality 0.149329. One fresh regenerated
six-load AMG workflow is selected under the existing 250,000-node and 8 GiB limits.
No stiffness, audited timing ratio, or workflow victory is claimed yet.

## Fixed-boundary repair does not yet qualify Gyroid, 2026-10-20

Job 366 accepts 64 interior moves and reduces the coherent actual-MMG output
from 104 bad tets to 16. Minimum quality 0.05126434 still rejects. Boundary,
connectivity, volume, periodicity, and sampled distance remain unchanged.
The second sweep accepts no moves. No FE runs.

Job 368 and independent extended-precision replay establish a fixed-neighbor
obstruction for the worst vertex at quality 0.1. Its active witness includes
three other interior vertices, so this does not rule out coupled moves.
The next diagnostic isolates fixed-boundary-only necessary constraints before
choosing another repair. Details and preserved failures are in
`MICROGEN_RESULTS.md`.

Boundary-only check 370 subsequently finds an obstruction independent of other
interior vertices, verified in extended precision. All eight 2-to-3 candidates
in 372 fail. A complete 13-vertex star solve in 374 reports infeasible after
quality filtering. Job 375 rejects one noncontained analytic Steiner point before
insertion. No source changes or FE follow. These bounded negatives do not
exclude larger remeshing or other boundary/Steiner constructions.

## Complete cap preparation still needs MMG quality, 2026-10-20

Jobs 350–354 complete the fixed Gyroid31 cap replacement on a fresh copy.
Single-neighbor growth fails, while joint one-rings supply the missing volume
triangulations. Independent auditing catches two near-flat cones in the first
complete copy. A relative determinant guard and actual-source regression
then produce 20 safe cones and four ring replacements.

The rebuilt seed preserves points, wall triangles, volume, oriented incidence,
and periodic cap matching. It has no analytically obstructing caps, minimum
ceiling 0.105301, and distance 0.00163287. Its actual minimum corner quality is
still 0.00109506, so this is not an admitted FE mesh or a workflow victory.
The next selected stage is an explicitly modified actual Microgen/MMG replay,
not a native Meshers change or a defaults claim. See `MICROGEN_RESULTS.md`.

Actual replay 359 later fits the intermediate and final tet caps at 29,876
nodes and 133,826 tets. It passes distance 0.00192726, volume, oriented facets,
and periodic cap matching, but quality 0.047711 fails 0.1. No FE runs.
Localization 361 finds 104 bad cells with no analytically obstructing boundary
faces. Free-vertex classification is next, before any relocation assumption.
The explicitly modified saved-seed diagnostic is not an admitted competitor
or a complete-workflow timing. Forty-one focused tests pass at this checkpoint.

## Gyroid cap alternatives remain incomplete, 2026-10-20

Jobs 345–349 follow the fixed-field extraction check. All 24 analytically
blocked cap faces have a valid local planar edge alternative. Complete-star
checks find 22 disjoint valid volume cavities, but two three-tet stars reject
both existing-apex cones. Independent exhaustive enumeration checks all
8,191 nonempty subsets of 13 candidate tetrahedra for each unresolved cavity.
Neither has a matching signed boundary with these six existing vertices.

No partial replacement, MMG call, or FE solve is applied. This does not rule
out larger cavities or added vertices. Next assess bounded one-ring cavity
expansion while preserving target cap polygons and periodic boundaries.
Detailed controls, failed serialization attempt, and limits are preserved in
`MICROGEN_RESULTS.md`. No comparative victory follows from these diagnostics.

## Gyroid extraction-order check, 2026-10-20

Jobs 342–344 compare fixed Gyroid31 samples and offset 1.1 with normal clipping
or an explicitly pretriangulated background. Normal extraction has 288
oriented-face conflicts. Five- and six-tet backgrounds remove those conflicts
and preserve periodic cap triangles, with identical original point/scalar
hashes. Each introduces 24 cap faces whose analytic quality ceiling is below
0.1, minimum 0.08710052. Neither is admitted and no MMG or FE solve follows.

The five-tet arm has 30,271 nodes, 101,472 tets, sampled distance 0.00163287,
and volume error +0.718972%. Six-tet distance fails at 0.00253049. The failed
first five-tet attempt is preserved, and its point-order guard receives a
regression before the successful diagnostic retry. Sixteen combined focused
tests pass. This is an extraction diagnostic, not a matched-accuracy win.
Detailed evidence and limits are in `MICROGEN_RESULTS.md`.

## Gyroid source consistency changes the repair direction, 2026-10-20

A convex single-star check finds a negative tangent bound at quality 0.1.
Independent extended-precision replay agrees. The subsequent 2-to-3 cavity
preflight stops at an invalid oriented source cavity rather than applying
a flip. A twelve-mesh lineage audit finds 684 same-orientation shared faces
in the VTK43 seed and 146 in the one-pass FEM output and its partial repair.
All have positive corner determinants. Retire that invalid repair source.
The two-pass Gyroid and all tested admitted controls have zero conflicts.

Every archived Schwarz timing arm is checked separately or exactly matched
to its audited control: 24 arm meshes and six manual lifts pass. Earlier
Schwarz results are unaffected. The common future geometry admission now
fails closed on incoherent oriented incidence, with all numerical thresholds
unchanged. Queued regression 340 reproduces legacy geometry metrics and the
old 96-comparison, 205-hash verification. Sixty-eight focused tests pass.

Next test the extraction order at fixed resolution 31 before any new MMG
setting: triangulate the same sampled background before clipping, explicitly
non-default, preserving field samples and box geometry while auditing the
new periodic caps. No such extraction trial
has run yet. See `MICROGEN_RESULTS.md` for limits and source-validity scope.

## Gyroid interior repair remains unqualified, 2026-10-20

One complete fixed-neighbor star on actual MMG output reaches quality 0.1
without changing boundary or volume. That local witness justifies one bounded
full interior-only repair on the smaller one-pass FEM precursor. Job 335
accepts 108 moves and reduces bad tets from 222 to 49, but worst quality stays
0.0079758. No move is accepted on the second sweep. Geometry and periodicity
are unchanged; the route remains unqualified and no FEA runs.

This explicitly added SciPy postprocess is not a Microgen default or a Meshers
seed. Stop its unchanged sweep recipe and diagnose one failed star before
choosing coupled movement or a topology change. See `MICROGEN_RESULTS.md`.
The progress-transition regressions and both queued checks pass using the
standard Job APIs, without changing the shared workbench helper.

## Repeated stronger Schwarz Tet10 comparison, 2026-10-20

Job 333 completes five interleaved measured pairs after warmup. All geometry,
periodic, field and eight-model supplementary 0.1% gates pass with common
AMG. Median audited workflows are 409.28 seconds for actual VTK30/one-pass
MMG plus manual straight Tet10 and 55.19 seconds for Meshers curved16.
Every paired ratio favors Meshers, 7.3858–7.5401. Generation/lift/solve-only
ratios are 7.2806–7.4331. Microgen has the smaller model disagreement; this
matches a target, not identical errors. It is not the default Microgen route.

The independent frozen-record verifier checks 96 mixed-load comparisons,
all timing statistics, 205 hashes, and geometry/field gates. Both clock and
resource limits pass. A late progress warning exposes cleared stage-end
values in the standard context manager. Fix reporting before another run;
do not erase the warning or discard scientific results because of it.
The four-fixture goal remains unachieved. See `TET10_RESULTS.md`.

## Fresh stronger Schwarz Tet10 comparison, 2026-10-20

Job 331 compares fresh actual VTK30/one-pass MMG plus manual straight Tet10
against fresh Meshers curved16. Both use the same explicit AMG preset and
pass the supplementary 0.1% target across all eight Schwarz models. Audited
workflow times are 409.87 and 55.10 seconds, a single-pair ratio of
7.4382. Generation/lift/solve subtotals are 388.25 and 52.92 seconds, ratio
7.3368. Both field audits pass; no saved mesh or solved load is replayed.

The native fresh-field positive 125-point audit 332 shifts energy by only
2.89e-9 relative. This qualifies the pilot for repeated timing, not a repeated
speed claim. The actual Microgen route remains explicitly one-pass with
manual Tet10 support, not defaults. Gyroid, spinodoid, and lattice conclusions
are unchanged. See `TET10_RESULTS.md` for recipes, errors, counts, and costs.

## Closed quartic corner-coupling diagnostic, 2026-10-20

Five algebraic tests and manufactured periodic controls pass in job 329.
Job 330 then checks 102 edge and 96 corner orbits on 1,408 frozen-domain
tets. The original edge-only gains reproduce exactly. Joint-minus-corner
gain agrees with the Schur complement within 5.37e-16 relative.

The incremental worst mixed-load reduction is only 0.676219 ppm, versus the
unchanged 500 ppm screen. Its 125/216-point change is 0.018333%, passing
precision. Trial-energy and periodic-trace checks pass. The worker takes
263.45 seconds and 0.6863 GiB. Stop this polynomial-mode route without a
native feature, a new reference, or another degree/strength trial.
These are nonequilibrated trial-energy diagnostics, not stiffness errors.
See `CUBIC_RESULTS.md` for records and limits.

## Fixed-domain quartic edge witness, 2026-10-20

A coupled residual probe on the saved cubic24 domain uses 102 periodic
curve-edge orbits and 306 mode DOFs. Derivative, affine, manufactured-energy,
pin-gauge, overlap, and quartic periodic-trace controls pass. The original
64/125-point precision check fails at 6.03% and remains rejected. A single
125/216-point follow-up passes precision at 0.017934% change.

The trial-energy reduction is only 0.648 ppm, far below the declared 500 ppm
signal screen. This is not a stiffness error or a new equilibrated reference.
No native sizing feature is justified. One possible omission remains to test:
coarse curve-corner relaxation coupled to the same edge modes, with separate
coarse-only gain and a Schur-complement check. See `CUBIC_RESULTS.md`.

## Three-level cubic consistency, 2026-10-20

Cubic23 completes six fresh loads at 5.976 GiB and passes its independent
field and positive-rule audits. The full 23/24/25 sequence does not support
the predeclared common leading-error matrix model. Correction directions
reverse with cosine -0.959707 and the best fit has 28.1% residual, above the
10% screen. Both existing length proxies fail. Seven synthetic controls pass,
and the exact 24-model set remains unchanged.

Do not fit only the closer odd levels, discard the middle tensor, or adopt
a favorable cubic reference for a timing claim. Next diagnose unresolved
curve-edge FE modes on one frozen domain, rather than another grid-count
or sizing-strength sweep. See `CUBIC_RESULTS.md` for the proposed bounded
quartic-edge residual witness and its limits.

## Completed finer cubic reference, 2026-10-20

All six fresh cubic25 loads complete within 8 GiB, peaking at 6.753 GiB.
The same geometry, quality, periodic, residual, and quadrature gates pass;
an independent fresh-process audit reproduces the fields. Its tensor differs
by 0.248177% from cubic24, so reference convergence is still unproven.
The curve pilot differs by 0.483742% from cubic25, the finest adaptive
Meshers tensor by 0.651535%, and the strongest Microgen-derived Tet10 tensor
by 0.917416%. These single-reference comparisons do not override the original
24-model failures. Cubic25 itself fails their model-normalized envelope at
1.920819%. A third neighboring parent23 passes geometry and fits the existing
cubic node stage. Read `CUBIC_RESULTS.md` before selecting further refinement.

## Cubic mapping sensitivity, 2026-10-20

Projecting higher-order cubic boundary nodes improves sampled distance from
0.00054354 to 0.00013101 while keeping topology and corners fixed. Both domains
pass admission. With the existing six displacement coefficient arrays held
fixed, positive 125-point energy changes by only 0.0105922%, below the frozen
0.05% signal screen. Quadrature sensitivity is 9.965e-7 relative. This does
not qualify the projected fields as equilibrated, establish a continuum bound,
or justify discarding any lattice model. No projected-domain solve is selected.
See `CUBIC_RESULTS.md` before choosing a higher-order refinement experiment.

## Eight-round curve follow-up, 2026-10-20

Eight adaptation rounds improve curve-target realization, quality, and cost,
but do not resolve lattice stiffness. At the unchanged geometric metric,
oversized curve edges fall from 288 to 32, minimum quality rises to 0.30065,
and the admitted mesh uses 111,605 nodes. A fresh six-load workflow takes
274.06 seconds at 2.851 GiB. Independent fields pass. Worst disagreement
across all 24 unchanged models is 1.439306%, still above 1% and worse than
the previous best Meshers result, 1.269882%. The primary fine-family envelope
itself reaches 1.301824%. Do not promote further round-count or curve-strength
tuning from this evidence. See `CREASE_RESULTS.md`.

## Resumed G3 study, 2026-10-19

The subsequent curve-local lattice pilot also fails the fixed 1% stiffness
target. It admits 118,475 nodes and 68,415 Tet10 elements, takes 342.19 seconds
for the fresh six-load workflow, and peaks at 3.032 GiB. Worst disagreement
across all 24 unchanged models is 1.456109%, worse than the previous best
1.269882%. Independent saved-field checks pass. Physical periodic-extension
normal jumps and full same-domain FE correction strongly localize near the
wall/box curves, but that evidence does not establish a successful live size
field. About 10.21% of curve edges remain above the requested split ratio
after four adaptation rounds. See `CREASE_RESULTS.md` before any successor.

The four-fixture objective remains unachieved. Tet10 recovery fails its
predeclared incremental localization screen, and required-entity grading
worsens the tested actual Gyroid MMG output. Neither negative result changes
admission or justifies a new sizing feature. See `STRESS_RESULTS.md` and
`MICROGEN_RESULTS.md`.

A focused native change has one measured geometry benefit. Quality-preserving
partial settling admits the previously rejected fixed spinodoid curved24
mesh at 105,618 nodes. Independent distance is 0.00173079, volume error is
+0.02597%, corner quality is 0.18951, sampled Jacobians are positive, and
periodic caps match. All 25 core and 37 focused Python tests pass. Its six-load
tensor passes independent stress, energy, and periodic-jump checks and differs
by 0.05430% to 0.13529% across all six unchanged spinodoid reference models.
The 640.08-second workflow is slower than the historical curved32 control,
despite 7.89% fewer nodes. This is no new workflow-speed victory. All four
admitted geometry controls pass, and positive 125-point fixed-field energy
differs by only 1.79e-8 relative from the solved tensor. Retain the native
admission fix and test one coarser geometry. See `PERIODIC_RESULTS.md`.

Further seed coarsening rejects, but a looser adaptive sizing tolerance of
0.003 on the fixed 24-cell seed passes at 73,951 nodes. All six audited loads
complete in a 364.07-second workflow with 2.61 GiB peak RSS. Disagreement is
0.10333% to 0.18439% across the unchanged reference models. Higher-rule energy
shift is 3.01e-8 relative. The archived pre-patch binary rejects the exact
recipe at distance 0.00237022, so the native admission gain is independently
confirmed. These single runs do not establish a Microgen timing victory.

The previously skipped actual Microgen spinodoid seed41 is now inspected in
a separate 300,000-tet seed-only preflight. It has 190 cap faces with quality
ceilings below 0.1, so retaining those caps cannot admit it. Final mesh and FE
caps remain unchanged. No MMG or FE successor is run on this known obstruction.
The four-case objective is still unachieved.

A non-oracle all-six-load stress-energy indicator passes a separate
fixed-domain localization screen, but its live lattice sizing pilot is
negative. After one native budget rejection and a geometry-only correction,
the fresh charged workflow uses 204,956 nodes and takes 812.17 seconds. Field
and positive-rule audits pass. Worst disagreement across all 24 unchanged
models is 1.352008%, versus 1.269882% for the previous best Meshers candidate.
Higher-rule energy changes by only 2.49e-7 relative. Stop that sizing recipe;
strong localization alone does not establish a workflow gain. See
`STRESS_RESULTS.md`.

An independent 67-million-point lattice volume audit reproduces every frozen
prefix and passes its precision screens. The mean becomes 0.2482374460, with
absolute replicate spread 0.0000137091. The finest Meshers volume error is
then -0.001898%, with zero inside the replicate range, instead of the old
nominal -0.03154%. A volume offset is not justified. The actual tuned
Microgen-derived straight Tet10 control remains about -0.4582% by volume.
No frozen errors, gates, or stiffness models change. See `VOLUME_RESULTS.md`.

The exact-tie compiler gradient defect is corrected separately. Smooth-minimum
and smooth-maximum regressions fail before the change; afterward 133 Python
tests pass with one optional wide-constraints skip, and all four geometry
controls remain admitted. The actual lattice gradient error drops to 1.81e-10.
A same-binary raw-clamped comparison reduces distance from 0.00039027 to
0.00026671 at 516 fewer nodes. Both meshes already pass geometry. Fully raw
input still fails the unchanged periodicity preflight, and the new raw-clamped
mesh has no six-load accuracy result. See `PERIODIC_RESULTS.md`.

## Current G3 outcome, 2026-10-18

One narrow matched-target result is supported: Schwarz Meshers beats the
explicitly modified actual one-pass Microgen comparator in five interleaved
pairs by 5.1155 to 5.2908 times. Both pass admission and the 1% stiffness target.
This is not the default Microgen pipeline or a four-fixture victory.

Gyroid curved9 now meets admission and the 1% target with 11,857 nodes and
a 20.56-second single-run workflow. Its fine-window order-sensitivity envelope
is 0.72921% to 0.95990%. The same-helper curved12 control takes 39.17 seconds.
Tested Microgen alternatives remain unqualified. Spinodoid has no cap-feasible
seed in the tested fixed-realization 16-through-40 range. Lattice still lacks
a cross-family reference supporting the 1% target. The actual tuned two-pass
Microgen grid51 route now passes admission and field consistency, but differs
by 3.89090% from the fine adaptive tensor. The fresh Meshers linear32 control
differs by 3.95566%. A further admitted adaptive level changes by only 0.03654%
from its predecessor but still reaches 1.26988% against retained reference
models. This does not establish 1% accuracy. See `MICROGEN_RESULTS.md` and
`PERIODIC_RESULTS.md` for the full recipes and exclusions.

A stronger genuine Microgen-derived straight-Tet10 Schwarz control reaches
0.05401% to 0.06270% under the same reference models. Fresh Meshers curved16
also meets the supplementary 0.1% target, with 42,928 nodes and an 88.01-second
live single run. No new paired workflow speed ratio is claimed from the
Microgen saved-mesh solve. Its coarser admitted variant is being tested before
selecting a live comparison. See `TET10_RESULTS.md`.

## Historical measurement protocol

The opening study used branch `claude/adaptive-implicit` and made no timing
claim. Later dated G3 sections contain same-host queued workflow measurements.
Its independent geometry metrics used `scorecard.py`:

- mesh to surface: Newton foot-point distance at 13 points per wall triangle;
- surface to mesh: exact point-to-triangle distance (VTK) from 2^17 projected samples;
- volume: analytic, or four scrambled Sobol sets of 2^23 points (uncertainty < 0.002%);
- quality: triangle shape 4√3A/Σl², minimum angle; MMG tet quality, minimum dihedral.

Cases: sphere r=0.35, ellipsoid (0.42, 0.25, 0.12), Gyroid sheet with raw band ±0.55
(the periodic-fea fixture), all in [-0.5, 0.5]^3, float64.

## Surfaces: curvature-adaptive remeshing (`remesh.rs`)

Target edge length from the chord sag of the exact surface, with gradation 0.3.
Split, collapse, flip and tangential relaxation, with exact projection onto walls,
box planes and their intersection curves. Matched Hausdorff distance H, against VTK
Flying Edges plus box clipping (the PR #2 control):

| Case | Flying Edges | Adaptive | Triangles |
|---|---|---|---|
| Sphere, H≈4e-4 | 10,616 tri, min angle 0.6° | 6,460 tri, min angle 33° | -39% |
| Ellipsoid, H≈4.5e-4 | 21,480 tri, min angle 0.2° | 7,234 tri, min angle 28° | -66% |
| Gyroid, H≈9e-4 | 44,938 tri, min angle 0.0°, ΔV +0.23% | 25,576 tri, min angle 18.7°, ΔV +0.10% | -43% |

With `balanced=True` the wall vertices move outward by half the local centroid sag,
so triangles deviate to both sides. Sizing then allows about twice the sag and the
measured deviation is still checked against the tolerance. Same independent H:

| Case | Adaptive | Balanced | Flying Edges at similar H |
|---|---:|---:|---:|
| Sphere, H≈4.2e-4 | 6,570 | 3,446 | 10,616 |
| Ellipsoid, H≈4.5e-4 | 7,194 | 4,130 | 21,480 |
| Gyroid, H≈9.8e-4 | 25,576 | 15,920 | 44,938 |

Balanced vertices are off the surface, so this mode is for printing and display;
curved quadratic FEA meshes keep vertices on the surface.

### Fast mode (`fast=True`): fewer triangles at Flying Edges speed

For printing and display, where triangle shape does not matter: collapse the 24-cell
extraction towards the size field (two passes, shape guard 0.05 instead of 0.3, no
flips or smoothing), then split only where the error exceeds the tolerance.
Supporting changes, which also apply to the quality modes:

- the error check per triangle is the exact maximum of the quadratic through its
  vertex and mid-edge samples (it was the largest sample, which read up to 2.7% low
  against the independent reverse distance on irregular triangles);
- balanced mode repairs locally: after balancing, faces still above the tolerance
  are split and balanced again (it used to tighten the tolerance everywhere, and
  in fast mode failed and fell back to unbalanced);
- lattice extraction solves all edge crossings in parallel and triangulates slabs
  in parallel with deterministic numbering; error scans, balancing and the final
  error check are parallel and incremental; manifold validation sorts instead of
  hashing.

Tolerance 1e-3, independent scorecard. Counts, H, ΔV and angles are deterministic.
Times: `timing_job.sh`, pueue job 1676 alone in the `gpu` queue on the OMEN (Intel
family 6 model 198, 24 threads), release build of 9dad555f frozen in `builds/`.
Compiled fields through Python, one warm-up, then the median of 5 runs interleaved
with all other configurations (`timing-threads1.json`, `timing-threads24.json`, with
every run). The pip VTK wheel uses the sequential SMP backend, so Flying Edges is
single-threaded in both columns; its two columns differ by run-to-run noise only.

| Case | Mesh | Triangles | H | ΔV | Min angle | 1 thread | 24 threads |
|---|---|---:|---:|---:|---:|---:|---:|
| Gyroid | Flying Edges 40 | 31,034 | 1.20e-3 | +0.34% | 0.0° | 0.107 s | 0.122 s |
| | Flying Edges 48 | 44,938 | 9.3e-4 | +0.23% | 0.0° | 0.175 s | 0.179 s |
| | fast | 26,780 | 9.9e-4 | +0.14% | 1.1° | 0.119 s | 0.086 s |
| | fast, balanced | 19,616 | 1.0e-3 | +0.08% | 1.8° | 0.125 s | 0.078 s |
| | balanced (quality) | 15,802 | 1.0e-3 | +0.05% | 18.0° | 0.226 s | 0.206 s |
| Schwarz P | Flying Edges 36 | 19,360 | 1.16e-3 | +0.26% | | 0.061 s | 0.072 s |
| | Flying Edges 40 | 24,064 | 9.5e-4 | +0.20% | 0.0° | 0.081 s | 0.091 s |
| | fast | 16,872 | 9.8e-4 | +0.19% | 0.0° | 0.076 s | 0.055 s |
| | fast, balanced | 11,880 | 1.0e-3 | +0.12% | 0.0° | 0.073 s | 0.048 s |
| | balanced (quality) | 8,946 | 1.0e-3 | +0.09% | 18.9° | 0.146 s | 0.116 s |
| Sphere | Flying Edges 32 | 4,808 | 1.01e-3 | -0.47% | | 0.016 s | 0.016 s |
| | fast, balanced | 1,954 | 9.9e-4 | -0.26% | 8.8° | 0.018 s | 0.015 s |
| Ellipsoid | Flying Edges 64 | 9,464 | 1.08e-3 | -0.41% | | 0.047 s | 0.049 s |
| | Flying Edges 80 | 14,856 | 6.9e-4 | -0.26% | | 0.077 s | 0.089 s |
| | fast, balanced | 2,950 | 1.0e-3 | -0.34% | 0.1° | 0.016 s | 0.013 s |

- Against the Flying Edges grid that meets H ≤ 1e-3 (48, 40, 32, 80), fast balanced
  has 2.0-5.0x fewer triangles and a smaller volume error except on the ellipsoid
  against grid 80. Speed: gyroid 1.4x faster on one thread and 2.3x on 24; Schwarz P
  1.1x and 1.9x; ellipsoid 4.7x and 6.8x (against grid 64, which misses 1e-3:
  2.9x and 3.8x); sphere 1.2x slower on one thread, parity on 24.
- The quality modes stay slower than Flying Edges at one thread (gyroid balanced
  0.226 s against 0.175 s) with 2.8x fewer triangles and 18° minimum angle.
- Volume block, same job (gyroid, tetrahedra): lattice-32 1.73 s / 0.58 s,
  adapt 1.81 s / 1.78 s (the tet adaptation barely uses the threads), MMG3D on the
  non-periodic seed 4.77 s / 4.65 s.
- Balancing offsets vertices by half the centroid sag, which minimizes H; the volume
  bias left is about a third (for a spherical cap the mean sag is 3/4 of the centroid
  sag). An offset of 3/4 would cancel the volume error at a larger H.
- Limits: the hard-min lattice (creases) fails the 1e-3 check in every adaptive mode,
  as before; Flying Edges 48 reaches only 1.26e-3 there, without a check. Fast mode
  keeps triangles valid (positive shape) but has slivers.

## Volumes: size-driven tet adaptation (`tet_adapt.rs`)

MMG-style split, collapse, edge removal, 2-3 swaps and relocation, with boundary
vertices projected onto the exact field. Against MMG3D (mmgpy 0.17.1) adapting the
same 24-cell seed with the same targets (hausd = tolerance, hmax 0.1, hgrad 1.3,
wall and box patches tagged):

| Case, tolerance | Adapted: tets / H / q_min / q_p01 | MMG3D: tets / H / q_min / q_p01 |
|---|---|---|
| Sphere 1e-3 | 13,689 / 7.7e-4 / 0.51 / 0.60 | 17,154 / 8.7e-4 / 0.42 / 0.52 |
| Ellipsoid 1e-3 | 16,059 / 9.6e-4 / 0.52 / 0.62 | 15,980 / 1.05e-3 / 0.36 / 0.51 |
| Gyroid 2e-3 | 42,441 / 1.5e-3 / 0.27 / 0.59 | 52,427 / 1.8e-3 / 0.14 / 0.47 |
| Gyroid 1e-3 | 74,662 / 1.0e-3 / 0.18 / 0.58 | 83,801 / 3.8e-2 / 0.15 / 0.46 |

MMG only sees the seed's facets, so its geometry cannot beat the seed (1.6e-3 for the
24-cell Gyroid); at 1e-3 it produced a 0.038 defect. Periodic cells are adapted on the
torus (copies merged, minimal-image geometry, box planes as internal constraints);
opposite faces match exactly. Periodic Gyroid at 1e-3: 74,938 tets, q_min 0.309,
against 171,679 tets, q_min 0.274 for the 32-cell lattice at the same H.

## FEA: homogenized stiffness error per DOF

Periodic Gyroid sheet, E=2.5, ν=0.25, the audited fedoo homogenization of the
periodic-fea experiment (six loads, single-thread CG). Error is the worst relative
strain-energy difference over all macroscopic strains against the quadratic
extrapolation of the archived native64/native80 linear tensors (observed order 2.1).

| Mesh | Order | Nodes | Elements | Error |
|---|---|---:|---:|---:|
| Lattice 24 | linear | 17,647 | 72,757 | 3.71% |
| Adapted 1e-3 | linear | 18,669 | 74,938 | 3.50% |
| Lattice 32 | linear | 39,717 | 171,679 | 2.12% |
| Lattice 8 | quadratic, straight | 6,914 | 3,257 | 4.24% |
| Lattice 12 | quadratic, straight | 19,931 | 10,832 | 1.67% |
| Lattice 8 | quadratic, curved | 6,914 | 3,257 | 1.31% |
| Lattice 12 | quadratic, curved | 19,931 | 10,832 | 0.29% |
| Adapted 4e-3 | quadratic, curved | 32,602 | 18,600 | 0.16% |
| Lattice 16 | quadratic, curved | 45,527 | 25,598 | 0.10% |
| Adapted 2e-3 | quadratic, curved | 70,024 | 42,768 | 0.04% |
| Lattice 20 | quadratic, curved | 79,848 | 46,881 | 0.04% |

- For linear tets in this thin sheet, adaptivity does not change error per DOF
  (3.50% vs 3.56% interpolated from the lattice): through-thickness resolution rules.
- Curved quadratic tets give 13-22x lower error at equal node count. A 1% target needs
  about 8k quadratic nodes and ~4k elements, against roughly 100k linear nodes and
  ~500k elements. Observed quadratic order is 3.7, close to the expected h^4.
- Curving the wall mid-edge nodes onto the exact field is most of the gain (3-6x over
  straight quadratic elements). Elevating a VTK or MMG mesh without the field gives
  only the straight result.
- The reference is cross-checked: extrapolating the quadratic 16- and 20-cell
  tensors with p=4 agrees with the linear 64/80 extrapolation within 0.007-0.012%,
  so quadratic errors are resolved down to ~0.04%. Quadratic lattice 20 (79,848
  nodes, 46,881 elements) is at 0.037%.

Balancing linear wall vertices (half the centroid sag outward, periodic copies
sharing one offset) cuts the volume bias about 2.5x but lowers the stiffness error
only ~5%: lattice 24 3.71% -> 3.51%, lattice 32 2.12% -> 2.01%, adapted 1e-3
3.50% -> 3.33%. Linear discretization, not geometry bias, dominates.

### Thick solid: Schwarz P network (`fea_solid.py`)

cos 2πx + cos 2πy + cos 2πz < 0, about 50% solid with thick ligaments. Reference:
p=4 extrapolation of the quadratic 16- and 20-cell tensors (they differ from it by
0.03% and 0.01%).

| Mesh | Order | Nodes | Elements | Error |
|---|---|---:|---:|---:|
| Lattice 16 | linear | 6,240 | 27,492 | 3.06% |
| Lattice 24 | linear | 19,004 | 90,792 | 1.35% |
| Lattice 32 | linear | 41,956 | 210,576 | 0.76% |
| Adapted 4e-3, hmax 0.2 | linear | 2,364 | 8,889 | 5.20% |
| Adapted 2e-3, hmax 0.2 | linear | 3,888 | 15,825 | 3.58% |
| Adapted 1e-3, hmax 0.2 | linear | 8,451 | 34,918 | 2.58% |
| Lattice 8 | quadratic | 6,864 | 3,852 | 0.26% |
| Lattice 12 | quadratic | 19,628 | 11,976 | 0.074% |
| Adapted 4e-3, hmax 0.2 | quadratic | 15,221 | 8,889 | 0.065% |
| Adapted 2e-3, hmax 0.2 | quadratic | 25,890 | 15,825 | 0.031% |

Geometric adaptivity lowers the linear error by about 20% at equal nodes on coarse
meshes and loses slightly at 8.5k nodes: FEA error follows stress resolution in the
ligament necks, which curvature sizing does not target. An error-driven size field
would be the next step for linear meshes. Curved quadratic elements again win by an
order of magnitude: 0.26% with 6.9k nodes and 3.9k elements against 0.76% with 42k
linear nodes and 211k elements.

Exact (quadrature) solid volume error of the same meshes:

| Grid | Linear | Curved quadratic |
|---:|---:|---:|
| 8 | +2.25% | +0.025% |
| 12 | +0.96% | -0.0007% |
| 16 | +0.47% | +0.0000% |

## Limits

- One FEA fixture (Gyroid sheet, one band). Thicker solids, lattices and spinodoids
  are untested; for thick solids interior grading should also reduce linear DOFs.
- Quadratic validity is checked by sampling the Jacobian at 25 points per element,
  not by Bernstein bounds; failing elements fall back to straight edges.
- No timing claims. Single-run wall clocks on a loaded machine suggest the adaptive
  surface is near Flying Edges' cost and the tet adapter near the lattice pipeline's;
  both are single-threaded.

## Issues #3 and #4 (8 cells, tolerance 0.01, minimum quality 0.1)

`issues.py`, full log in `issues.txt`. Independent two-sided distance H and volume
error of accepted meshes:

| Fixture | Lattice, 0 / 4 / 12 passes | Adapter (8-cell seed) |
|---|---|---|
| Gyroid | accepted all; H 1.0e-2, q_min 0.15 | H 9.9e-3, ΔV +1.5%, q_min 0.33 |
| Schwarz P | rejected at 0 passes | H 8.1e-3, ΔV +1.4%, q_min 0.39 |
| Neovius | rejected (2.9e-2) | H 9.8e-3, ΔV +1.8%, q_min 0.28 |
| Graded Gyroid | rejected at 0 passes (quality) | H 9.8e-3, ΔV +1.4%, q_min 0.34 |
| Thin Gyroid | rejected at 0 passes (quality) | H 9.0e-3, ΔV +1.8%, q_min 0.32 |
| Implicit BCC lattice | rejected (4.2e-2) | rejected (3.5e-2) |

Issue #3's non-monotonic Gyroid is fixed on this branch by PR #5's geometry guard.
Two adapter changes made the hard cases pass: projection falls back to bracketing
and bisection, and a wall vertex whose exact position would invert a neighbour is
inserted part of the way and moved onto the surface once refinement makes room.
That move keeps the split quality floor, so geometry is never bought with slivers.

The capsule lattice is a hard `min` of eight distances, so its surface has concave
creases. Triangles that straddle a crease keep a sag proportional to their size,
and the vertices that would resolve it cannot settle without slivers. The lattice
pipeline first passes at 24 cells. Capturing creases needs feature-aware
projection onto the curve where two `min` branches are equal; the compiled-field
JIT knows those branches. Spinodoid is untested here (needs Microgen).

### Curve chords and the spinodoid

The spinodoid (`spinodoid.py`, run in the PR #2 environment for Microgen's
`Spinodoid`) exposed a measurement gap. On a wall face, a chord that lies in a box
plane follows the wall/plane curve. The first-order check used the full field
gradient, which measures the distance to the unclipped wall; where the wall meets
the box at a shallow angle that wall passes close outside the box, while the curve
inside the plane is far. The 16-cell adapted spinodoid passed its own check (7.3e-3)
but its independent surface-to-mesh distance was 1.23e-2, at a curve on y = 0.5.

The adaptive paths now measure such chords with the in-plane gradient, size curve
edges by the curve's in-plane turning, and the surface remesher samples every edge of
every wall face (curve edges with a cap neighbour were sometimes never sampled).
The adapted spinodoid is now rejected at 16 cells (2.0e-2) instead of passing with a
hidden error. The lattice pipeline keeps its released check: with the in-plane
measure, three of its existing tests exceed their tolerance by 8% to 6x. Whether to
tighten the released check is a tolerance decision left open.

### Lattice with smooth unions

`smooth_lattice.py` rebuilds the issue #4 lattice with Microgen's cubic polynomial
smooth union (`implicit_ops.smooth_union`). Microgen's strut lattices currently go
CAD -> STEP -> gmsh; its F-rep path already offers smooth unions. A sharp concave
crease is also a stress singularity, so FEA would not converge there anyway.

| Blend k | Cells | Lattice pipeline | Adapter |
|---|---:|---|---|
| 0.08 | 12 | 5,696 tets, H 5.8e-3, ΔV -3.3%, q_min 0.32 | 6,225 tets, H 1.0e-2, ΔV -3.1%, q_min 0.48 |
| 0.08 | 16 | 24,480 tets, H 7.4e-3, ΔV -0.5%, q_min 0.20 | 8,129 tets, H 6.9e-3, ΔV -2.0%, q_min 0.34 |
| 0.03 | 8-16 | rejected | rejected (closest 1.09e-2 at 16) |

A blend much smaller than the cell size behaves like a crease. When the exact
position of a wall vertex would flatten a neighbour, the adapter now splits the
longest edge facing it in each such element and retries, which made the smooth
fixtures slightly better (Neovius q_min 0.28 -> 0.34). Equalizing sag to the
tolerance raises the mean deviation, hence the larger volume error of the adapted
16-cell mesh at equal H; balanced placement would cancel it for linear meshes.

## Real Microgen defaults on Linux G3, 2026-10-17

`MICROGEN_RESULTS.md` records the full default VTK-to-MMG workflow on the
Gyroid sheet, Schwarz P solid, fixed spinodoid and smooth-union lattice.
This replaces no earlier frozen results and does not relabel the historical
MMG-on-a-Meshers-seed measurements. Both Microgen paths now have real G3
execution evidence on all four fixtures, with unsupported FEA routes marked.

At about 20k nodes, curved Meshers has 0.2951% Gyroid and 0.0738% Schwarz-solid
full mixed-load stiffness disagreement, versus 3.7152% and 1.3458% for the
tested linear meshes. Their geometry passes the unchanged common gates,
but curved CG workflows take longer at these unmatched accuracies.
Every default Microgen control fails geometry or quality admission.

The negative findings matter. Meshers rejects periodic spinodoids at the
quality floor and rejects the smooth lattice's raw-field periodicity,
although Microgen completes diagnostic lattice FEA. Neither fixture has a
justified fine stiffness reference. Same-domain order controls separate
approximation-order changes from curved geometric-mapping changes without
claiming certified continuum error. Native meshing code remains unchanged.

## Public direct-fedoo Tet10 integration, 2026-10-17

`TET10_RESULTS.md` records the successful G3 integration through `Mesh.to_fedoo()`.
Four six-load cases pass: solid Tet4, solid Tet10, curved Gyroid and curved Schwarz
solid. Independent strain-controlled stress integration reproduces the saved
specimen tensors within 1.2e-11 relative mixed-load energy disagreement. Periodic
mid-edge displacement jumps remain below 6e-17. This validates conversion and
periodic FE integration, not local stress accuracy. Microgen's public converters
remain unsupported for Tet10.

## Spatial sizing from recovered stress, 2026-10-17

`STRESS_RESULTS.md` records an optional periodic target-length callback and two
Schwarz pilots. Recovery sizing loses at the coarse target after charging its
coarse solve. At target 0.04 it uses 22,664 rather than 33,657 nodes and takes
25.87 rather than 36.51 seconds, with 1.0207% rather than 1.0401% worst mixed-load
stiffness disagreement. All four meshes pass the unchanged admission checks.
This is a single-run Meshers comparison, not a matched-accuracy Microgen victory
or a local-stress accuracy claim.

## Matched-target comparison and remaining uncertainty, 2026-10-18

Five interleaved Schwarz pairs establish a narrow matched-target result.
Meshers curved Tet10 uses 6,864 nodes at 0.25870% mixed-load stiffness
disagreement, versus 32,070 Tet4 nodes at 0.92806% for an explicitly modified
actual Microgen VTK-to-single-MMG route. Both pass all admission gates and the
1% stiffness target. Median live total workflows are 7.7916 and 40.5036 seconds,
with Meshers faster in every pair by 5.1155 to 5.2908 times. This is not the
default two-pass Microgen pipeline and not a four-fixture victory.
See `MICROGEN_RESULTS.md` for the route, same-solver protocol and field audits.

The fixed spinodoid now has admitted curved meshes and a sensitivity-checked
empirical reference, but no admitted Microgen timing comparator. The lattice
also admits curved meshes, after explicit equivalent-field preparation and
measured native corrections. Its uniform refinement family does not have a
stable observed stiffness convergence rate. An independently audited adaptive
123,570-node tensor differs by 0.82695% from uniform32; positive higher-rule
reintegration changes fixed-field energy by at most 5.65e-7 relative. Lattice
1% stiffness accuracy remains unproven. See `PERIODIC_RESULTS.md` for rejected
changes, admission failures, native regression tests and reference limits.

The tighter adaptive lattice subsequently completes at 175,818 nodes and
5.64 GiB peak RSS. Its six-load workflow takes 989.67 seconds. Independent
field checks pass within 4.02e-11. Its adaptive-family sensitivity envelope is
0.15983% to 0.94546%, but comparison with earlier uniform-family models reaches
1.30614%. Both families remain recorded; no tested lattice candidate passes
the union of declared models at 1%. That unresolved uncertainty is not a
matched-accuracy success. Actual Gyroid VTK seed screening identifies one
cap-feasible resolution among twelve tested, and genuine tuned MMG calls are
being assessed under the unchanged element and admission caps.

## Gyroid comparison limits and cheaper Meshers pilot

Seventeen actual VTK Gyroid seeds reveal two cap-feasible resolutions, 31 and
43. Tuned genuine two-pass MMG routes on both can pass distance, volume,
periodicity and resource caps, but fail minimum quality near 0.0406. Two
explicit third-pass optimization variants do not repair that failure and are
removed from the maintained runner. Frozen drivers retain the negative tests.

Separate unranked six-load diagnostics find 2.58544% and 2.29279% stiffness
disagreement on those returned meshes. Independent stress and energy agree
with reactions within 7.26e-12. Quality repair alone therefore would not meet
the 1% target. These meshes remain unqualified and excluded from admitted
timing comparisons. See `MICROGEN_RESULTS.md` for exact recipes and limits.

A live Meshers curved10 pilot passes all gates with 13,953 nodes, 0.45586%
stiffness disagreement and 24.71 seconds total. Independent field checks pass
within 6.79e-12. A fresh same-native, same-helper curved12 control is selected
before attributing cost savings to coarsening. The earlier 47.08-second
curved12 control used a different solver formulation and native revision.
No Gyroid matched-accuracy Microgen victory is established.

## Stronger lattice control, 18 October 2026

A reference-only cubic enrichment on the admitted 24-cell Meshers geometry
now completes all six loads and independent stress, energy and positive
125-point quadrature checks under 8 GiB. It does not resolve the disagreement
among the 24 retained empirical reference models or certify continuum error.
See `CUBIC_RESULTS.md`.

An actual tuned VTK51/two-pass MMG lattice seed admits 33,903 Tet4 nodes,
158,252 tets, distance 0.00162159, volume error -0.48770% and quality 0.109983.
Manual straight Tet10 elevation and direct fedoo conversion preserve that
physical geometry and use 238,411 nodes. This is not Microgen's unsupported
public Tet10 conversion or its default meshing recipe. The original CG control
times out after five loads. An explicitly selected, validated AMG retry
completes all six and passes fresh-process field audits at 4.378 GiB peak RSS.

Worst disagreement across all retained models is 1.004354%, versus 1.269882%
for the best tested Meshers mesh. Neither passes 1%. The Microgen result is
stronger accuracy evidence and must not be rounded into admission or omitted.
Saved solves with different solver policies do not establish a timing win.
See `TET10_RESULTS.md` and `runs/g3-lattice-tet10-retained-comparison-v1`.

The remaining worst direction is mainly deviatoric normal strain. A verified
same-domain Tet10-to-cubic correction puts 77.46% of its energy in the 6.99%
of coarse elements incident on cap vertices. This motivates one cap-focused
size-field pilot rather than more global refinement. Its geometry stage passes
with 159,147 nodes; read its full solve before selecting another change.

The completed cap-focused pilot passes independent field and positive
125-point quadrature checks, but worst retained-model error is 1.411382%.
Its untargeted 123,570-node baseline has 1.427185% error, so the extra 35,577
nodes give little benefit. It does not improve the best Meshers result or beat
the stronger actual Microgen Tet10 control. A journal-confirmed host suspend
also excludes its timing from comparison. Stop blind cap-length tuning.

Kevin requests stopping G3 within the next hour. Finish and preserve this
negative bounded study, then stop research for shutdown. The four-case goal
remains unachieved. Restart instructions and exact current evidence are in
`/home/kmarchais/dev/meshers/handoff/g3/CHECKPOINT-20261018.md`.
