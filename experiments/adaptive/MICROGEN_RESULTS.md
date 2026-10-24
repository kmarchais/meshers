# G3 Microgen workflow results, 17 October 2026

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


## Fixed-spinodoid extraction and cap feasibility on Mac, 24 October 2026

Mac jobs 14–20 investigate missing actual-comparator coverage. No MMG or FE
solve runs in this sequence. All final mesh admission gates and the
250,000-node/tet caps remain unchanged.

### Exact field across FFT backends

The Mac FFT constructor differs from the frozen G3 realization in 33 of
260 coefficient entries, by at most 8.88e-16. The new exact-coefficient
regression fails before correction. `microgen_workflow.specimen` now restores
the authoritative 52-mode amplitudes and phases before any spinodoid grid
sampling. It first checks identical mode indices and threshold and rejects
differences beyond platform roundoff. The source fixture is
`runs/g3-spinodoid-compiled32-v1/fourier-coefficients.json`, unchanged.
Its hash and raw-constructor difference are recorded per new seed.

Three focused controls verify exact coefficients, rejection of a changed
mode set, and unchanged field values when sampling resolution changes.
The constructor FFT remains 16. This loads the specified fixed specimen;
it does not regenerate a different realization or alter the stored samples.

### Paired extraction at one sampling resolution

Job 14 compares original Microgen extraction with explicit five-tet
background triangulation before scalar clipping, at 36 sampling cells.
The 233,280 background tets fit the existing cap. Ten queued controls pass.
All 52 coefficients, background coordinates and scalar arrays match between
arms. Field values at 4,096 probes remain bitwise unchanged within each arm.
Both call the actual `Spinodoid.generate_volume_mesh` implementation.
The frozen-field load and non-default sampling are explicitly labeled.

Original extraction has 38,712 nodes and 182,287 tets, reproducing the old
seed36 counts and 122 obstructing cap triangles. It also has 282 oriented
shared-face conflicts. Pretriangulation produces 49,032 nodes and 180,277
tets with zero such conflicts and matching periodic cap triangles. Its
164 obstructing caps have minimum adjacent-tet quality ceiling 0.0000491890.
Actual minimum corner quality is 1.40805e-9. Both arms remain rejected.

The coherent arm's reverse sampled distance is 0.00449448, above 0.002.
It also has 148 unresolved forward projections; its raw forward maximum is
not promoted as a reliable complete distance bound. Absolute volume error
is 0.00528745%. No MMG follows on analytically blocked caps, and no stiffness
diagnostic or comparison timing is inferred.

### Necessary planar cap checks

Job 15 tests 348 individual planar alternatives incident on the 164 blocked
faces. Only 58 blocked faces have a passing single-edge alternative; 106
do not. Alternatives are not applied or claimed jointly compatible.

Job 16 scans all 2,388 cap-plane boundary corner instances. Each planar
boundary link is unambiguous, and none has the acute-corner obstruction.
For a fixed convex wedge of angle `a <= 60 degrees`, every contained
triangle has quality at most `sqrt(3)*sin(a)/(2-cos(a))`. Equal adjacent
lengths attain this bound. Combining it with the existing retained-face
tet bound gives a necessary fixed-polygon condition, not a sufficient
triangulation theorem. Five controls cover actual triangle/tet quality,
orientation, reflex corners, and 80-digit arithmetic.

Every blocked face has at least one cap-interior vertex. The per-plane
counts by one or two free vertices are 24/5 on each x face, 24/6 on each
y face, and 19/4 on each z face. No blocked face has zero free cap-interior
vertices. This motivates a planar-star test rather than another fixed-volume
interior relocation sweep.

### One periodic cap-interior star

Job 17 stops in its manufactured control because it incorrectly uses
`cells=2`; the native API requires at least four. No production diagnostic
runs. Its failed log, input hashes, and original control source remain.
Job 18 uses the corrected four-cell control, passes analytic-gradient,
concavity, periodic, feasible and obstructed-star checks, then tests the
worst real cap face, nodes 13320/37253/37254.

That face has one free cap-interior orbit, nodes 13306 and 13320. The test
allows one common planar displacement, fixes all neighbors and the outer
cap polygon, and checks all incident cap triangles plus positive volume
seed guards. These guards preserve a valid seed; they are not the final
tet-quality admission gate. No proposed move is applied.

The complete star cannot attain the cap target with fixed neighbors. Its
supporting-plane certificate uses only two cap triangles. Job 19 independently
replays their area and squared-edge derivatives and box support using
80-digit Decimal arithmetic and normalized weights. At the original 0.1
gate, the residual upper bound is -0.00063315085191554. This is a checked
fixed-neighbor floating-point witness, not interval certification.

The second active triangle includes another free cap vertex, 13894. The
witness therefore does not rule out coupled cap motion or retriangulation.
Do not apply a whole-cap single-vertex sweep from this failed local star.
The active boundary edge is only 4.03873e-6 long, whereas the edge between
the two free cap vertices is 0.0277778. A rigid shared translation preserves
that long edge and is too narrow a next test. The selected preflight gives
the two active free vertices independent planar motion, with their periodic
copies and outer curves fixed. It has four unknowns. Its constraints are
nonlinear; a failed search is inconclusive rather than an infeasibility
certificate. Read direct candidate quality, positivity, and periodic checks
before any mesh edit, MMG, or FE. Geometry rejection also remains unresolved.

Job 20 passes 29 focused extraction, fixture, corner, planar-star, witness,
progress, and resource controls. Independent checkpoint verification rechecks
both saved seed artifacts, 210 input hashes, all 3,509 transferred files,
and the previous 72 full lattice model comparisons. The failed job 17 is
explicitly retained in the receipt. Neither new spinodoid seed is an admitted
actual comparator, and the four-case objective remains unachieved.

Records are `runs/macbook-spinodoid-preclip36-five-v1`,
`runs/macbook-spinodoid-cap-flip-preflight-v1`,
`runs/macbook-spinodoid-cap-corner-bounds-v1`,
`runs/macbook-spinodoid-cap-plane-star-v1`,
`runs/macbook-spinodoid-cap-plane-star-v2`,
`runs/macbook-spinodoid-cap-plane-witness-v1`, and
`runs/macbook-spinodoid-preflight-controls-v1`.

## Valid-source interior repair remains incomplete, 20 October 2026

Classification 364 confirms all 104 bad cells have free interior vertices.
Job 366 therefore tests the existing two-sweep, 400-attempt fixed-boundary
repair on the coherent prepared two-pass output. Seven queued controls pass.
The repair accepts 64 moves in sweep zero and none in sweep one. There are
97 attempts, with 33 unsuccessful searches. Bad elements fall from 104 to 16,
but minimum quality is only 0.0512643362. The result stays `repair-incomplete`
and receives no FE solve.

Boundary coordinates and connectivity remain fixed. Signed volume changes
by -7.40e-17. Distance remains 0.0019272611, volume error +0.72004462%, and
periodic caps and oriented incidence still pass. Repair takes 8.73 seconds,
the bounded worker 14.42 seconds at 0.4054 GiB. These saved-mesh diagnostic
costs are not charged live-workflow timings.

Do not repeat the unchanged second sweep. Job 368 diagnoses the new worst
element, 62548, whose only free vertex is 19446. Its complete star has 20
tets. At quality 0.1, the checked concave-tangent upper bound is
-0.001009524655965. Independent 64-bit-mantissa replay agrees within 2e-16.
This rejects that single vertex with all neighbors fixed, not every repair.
The active witness also depends on interior vertices 805, 8895, and 10865.
It does not yet exclude moving those vertices together.

The selected next diagnostic uses only tets whose three opposite vertices
are all on the fixed boundary. It checks whether those necessary constraints
already obstruct the worst vertex before selecting a coupled relocation or
topology change. No new moves or FE runs in that diagnostic. Numerical gates,
references, and native code remain unchanged.

Job 370 performs that check on seven tets. The boundary-only tangent bound
is -0.000976627527263. Independent 64-bit-mantissa replay agrees, and every
active opposite vertex is on the fixed boundary. Moving other interior
vertices cannot change those supported constraints while the cells and
boundary remain fixed. Coupled interior-only relocation is therefore not
selected. This does not exclude changing connectivity or boundary geometry.

Job 372 tests eight adjacent two-to-three cavities on this coherent source.
Three fail oriented-boundary matching, three fail new-cell quality, and two
are degenerate. No local candidate passes and no flip is applied. This is
distinct from the earlier failed preflight on an invalid one-pass source.

Job 374 tests the complete 20-tet star, with 13 existing vertices and its
exact oriented outer boundary. The minimum fixed-face ceiling is 0.105301.
There are 477 candidates after positive-determinant and quality 0.100001
guards. The bounded binary solve reports infeasible. It does not certify
all possible larger cavities, extra vertices, or boundary changes.

Job 375 tests one additional candidate at the analytic quality-maximizing
apex of the highest-weight cap face. Its single-tet quality is 0.105301, but
its maximum minimum barycentric coordinate over the old star is -0.866917.
It is not contained, so no point is inserted and no binary solve follows.
This rejects one analytic candidate, not every Steiner construction. Ten
queued controls include apex-side, containment, and original-volume checks.
The extended solve helper treats added candidate vertices separately from
the original cavity, whose domain and volume stay fixed.

Stop these unchanged sweeps and local recipes. Gyroid still has no admitted
actual Microgen comparator. Larger-cavity or boundary-triangulation work
remains unproven, not universally impossible. Before more meshing, the next
comparison diagnostic checks whether the retained lattice model envelope
even admits a common one-percent full tensor, without changing its models.

Records are `runs/g3-gyroid-prepared-mmg-interior-repair-v1`,
`runs/g3-gyroid-prepared-repair-star-feasibility-v1`, and
`runs/g3-gyroid-prepared-repair-star-witness-audit-v1.json`, plus
`runs/g3-gyroid-prepared-repair-boundary-star-v1`,
`runs/g3-gyroid-prepared-repair-boundary-witness-audit-v1.json`,
`runs/g3-gyroid-prepared-repair-two-to-three-v1`,
`runs/g3-gyroid-prepared-repair-star-topology-v1`, and
`runs/g3-gyroid-prepared-repair-star-steiner-v1`.

## Actual MMG replay of the prepared Gyroid seed, 20 October 2026

The prepared source is `runs/g3-gyroid-cap-replacement-v2`. These are genuine
calls through `microgen.remesh.remesh_keeping_boundaries_for_fem` on an
explicitly modified Microgen-derived seed. They are not Microgen defaults,
Meshers-seeded MMG, or fresh complete-workflow timings. Preparation cost is
excluded, each intermediate tet cap remains 250,000, and no FE runs.

| Job | First-pass hausd | Second-pass hausd | Measured stop |
|---|---:|---:|---|
| 355 | 0.0005 | 0.0005, not reached | First pass returns 649,230 tets; 58.02 s diagnostic |
| 356 | 0.0015 | 0.0005, not reached | First pass returns 259,134 tets; 27.23 s diagnostic |
| 357 | 0.00175 | 0.0005 | First pass fits at 229,778 tets; second returns 446,978 and hits cap; 62.20 s diagnostic |
| 359 | 0.00175 | 0.0015 | Final 133,826 tets fit, but minimum quality 0.047711 rejects; 48.71 s diagnostic |

Every first pass retains actual `nofem=True`, with explicit hmax 0.06 and
hgrad 1.3. The second public call has `ls=True`, `nr=True`, and the recorded
sizing arguments. Captured output counts stop the worker before another pass
or solve. A rejected cap is not waived because an output is close to it.

The captured first stage of 357 has 48,565 nodes, sampled distance
0.0020306511, volume error +0.69621401%, and minimum quality 0.03323404.
It fails both distance and quality. Passing the count limit does not qualify
that intermediate mesh. The second output is not geometry-ranked above the
stage cap. All three records have `qualified=false`.

Job 359 completes that looser second-pass replay at 29,876 nodes. Sampled
distance is 0.0019272611, volume error +0.72004462%, and minimum sampled
Jacobian 2.388798e-8. Oriented incidence and periodic cap triangles pass.
No cap face is below the analytic ceiling, whose minimum remains 0.105301.
Only final corner quality rejects this candidate. Peak worker RSS is 0.814 GiB.
No FE solve runs, and none of these saved-seed controls is a timing victory.

Localization 361 finds 104 elements below 0.1. All touch both cap and wall
vertices. Seventy-five have a cap face, 26 a wall face, and none has an
obstructing fixed cap or wall face under the analytic ceiling. The sets can
overlap. This does not establish that interior relocation can repair them.
The selected next check classifies free interior vertices in each bad cell
before choosing a repair or another MMG call.

Classification 364 confirms all 104 bad cells have at least one interior
vertex, and none is all-boundary. The worst cell has two free interior
vertices. This makes interior relocation applicable as a trial, not a
guaranteed repair. The necessary face-ceiling checks still pass.

The existing fixed-boundary relocation runner now accepts an explicit source
path. Its default retains the retired one-pass source and its incidence
rejection. The additional supported source is this recorded two-pass prepared
Gyroid route. Seven focused star, whole-mesh, and progress controls pass.
The selected trial is two sweeps, at most 400 vertex attempts, complete-star
quality target 0.100001, fixed boundary/connectivity, and signed volume
preservation. It allows 600 seconds, one thread, and 8 GiB. No FE follows
before fresh full admission. This is explicitly added SciPy postprocessing,
not a Microgen default or an uncharged end-to-end comparison.

Older saved-seed rows reuse the historical `vtk_extraction_s` field for load
time, while explicitly marking `saved_seed_replay` and excluded preparation.
Future rows use `saved_seed_load_s` and `vtk_extraction_performed=false`.
Do not interpret the old field as a measured live extraction or edit frozen
rows to relabel it.

Records are `runs/g3-gyroid-prepared-two-pass-mmg-v1`,
`runs/g3-gyroid-prepared-two-pass-mmg-hausd0015-v1`, and
`runs/g3-gyroid-prepared-two-pass-mmg-hausd00175-v1`, plus
`runs/g3-gyroid-prepared-two-pass-mmg-second0015-v1` and
`runs/g3-gyroid-prepared-mmg-quality-localization-v1`, with free-vertex
classification in `runs/g3-gyroid-prepared-mmg-quality-localization-v2`.

## Complete rank-guarded Gyroid cap preparation, 20 October 2026

Jobs 350–354 resolve the previously incomplete cap replacement without moving
coordinates, wall triangles, or periodic boundary curves. This remains an
explicit modification of Microgen extraction, not its default workflow.

Job 350 tries each single face neighbor of the two unresolved three-tet stars.
Each has five neighbors. All ten seven-vertex binary models are infeasible.
The binary signed-face solver first reproduces positive exhaustive controls
and the previous six-vertex obstruction. It bounds solve time at two seconds,
branch nodes at 1,000, and threads at one. Any returned subset needs exact
integer face replay, incidence counts, positive determinants, and volume
agreement. Solver status alone cannot admit a cavity.

Job 351 considers each measured complete one-ring jointly. Both nine-vertex
rings solve, with 105 candidate tets and eight selected tets per ring. Signed
face equations replay exactly and relative volume differences are below
2.5e-16. The resulting 24-plan copy in 352 passes global oriented incidence,
periodic nodes/cap triangles, unchanged wall/point arrays, and unchanged volume.
Its cap ceiling improves to 0.105300997. However, two of the earlier accepted
cones contain roundoff-scale cells with determinants near 1e-20. The final
quality check rejects that copy at approximately 2.75e-16. It is preserved,
and no MMG call uses it.

A regression uses one actual affected six-vertex star. Preparation now has an
explicit relative determinant guard, rejecting candidate determinants at most
1e-12 times the original cavity's summed determinant. This is a numerical
rank guard, not a replacement or relaxation of final quality 0.1. Job 353
recomputes 20 safe cones and four full rings. Each ring solves from 102 retained
candidate tets, with minimum selected determinant at least 2.36e-7. Fifteen
queued controls pass, including the actual near-flat example.

Job 354 applies those 24 plans to a fresh copy and independently audits it:

| Measurement | Result |
|---|---:|
| Nodes / tets | 30,271 / 101,472 |
| Sampled two-sided distance | 0.0016328718 |
| Relative volume error | +0.71897165% |
| Relative volume change from seed | 0 |
| Minimum retained-cap quality ceiling | 0.1053009971 |
| Cap faces below ceiling 0.1 | 0 |
| Minimum replacement determinant / original cavity sum | 0.0015171466 |
| Minimum actual corner quality | 0.0010950575 |
| Minimum sampled Jacobian | 1.417393e-9 |

Oriented facets and periodic cap triangles pass. Original points and wall
triangles are identical. This is a cap-feasible seed for a bounded MMG test,
not an admitted final FE mesh. The seed still fails the actual quality gate.
No FE solve or matched-accuracy timing claim follows from preparation.

Records are `runs/g3-gyroid-cap-cavity-expansion-v1`,
`runs/g3-gyroid-cap-cavity-fullring-v1`,
`runs/g3-gyroid-cap-replacement-v1`,
`runs/g3-gyroid-cap-cavity-rankguard-v1`, and
`runs/g3-gyroid-cap-replacement-v2`. Old binary-solver source snapshots match
their recorded hashes. Worker memory remains below 0.49 GiB, with one thread,
120-second preparation stages, and unchanged 8 GiB/250,000-tet limits.

## Gyroid cap-flip feasibility, 20 October 2026

The coherent five-tet seed has 24 cap faces below the necessary analytic
quality ceiling. Jobs 345–349 test only local replacements. No MMG call or
FE solve runs, and no source or partial mesh replacement is applied.

Job 346 enumerates 48 cap-edge alternatives incident to those blocked faces.
Twenty-four pass the 0.1 ceiling on both new triangles, twelve fail the
ceiling, and twelve fail convexity or oriented patch-boundary matching. All
24 blocked faces have a passing local alternative. Each passing alternative
retains the exact oriented planar patch boundary and area. This is not yet
a compatible periodic volume mesh.

First attempt 345 passes seven controls but fails while serializing NumPy
edge IDs. Its failed log and original source/test snapshots remain frozen.
The retry converts IDs to Python integers and adds a JSON roundtrip regression.

Job 347 tests complete tetrahedral stars around the 24 selected cap edges.
Twenty have valid two-tet replacements, with positive determinants, unchanged
non-cap oriented cavity boundary, and unchanged volume. Four have three-tet
stars and are not silently truncated. No partial replacement follows.

Job 348 permits only those measured two- or three-tet stars. It tests cones
from existing non-patch vertices and finds 22 valid, mutually disjoint cavity
plans. Both existing-apex cones fail on each remaining three-tet star. Since
all 24 replacements are needed, none is applied. A manufactured three-tet
control passes before the real diagnostic.

Job 349 independently enumerates every nondegenerate positive-oriented
tetrahedron on the six existing vertices of each unresolved cavity. Each has
13 candidate tets and 8,191 nonempty subsets. None matches the target signed
boundary-face equations. Valid matches would also need volume agreement and
at most two cells per face. Four controls include reproduction of the original
cavity, a successful changed-boundary convex cavity, the archived nonconical
case, and rejection of an unclosed boundary. Six controls including progress
tests pass in the queued run.

This is a finite-precision fixed-vertex obstruction for two particular
cavities, not a general remeshing impossibility. It does not exclude larger
cavities or added vertices. The next selected check is bounded one-ring cavity
expansion around those two stars, with unchanged target cap polygons and
periodic boundaries. Do not rerun failed cones or submit the blocked seed to
MMG. Native Meshers, references, specimens, and admission gates are unchanged.

Workers remain below 0.294 GiB, one thread, 120 seconds per stage, and the
8 GiB ceiling. These are saved-seed diagnostics, not workflow timings. Records
are `runs/g3-gyroid-cap-flip-preflight-v1`,
`runs/g3-gyroid-cap-flip-preflight-v2`,
`runs/g3-gyroid-cap-cavity-preflight-v1`,
`runs/g3-gyroid-cap-cavity-preflight-v2`, and
`runs/g3-gyroid-cap-cavity-enumeration-v1`.

## Fixed-field Gyroid extraction order, 20 October 2026

Jobs 342 and 344 compare normal Microgen scalar clipping with background
tetrahedralization before its two clips, at resolution31 and offset 1.1.
Coordinates and every sampled scalar array have identical hashes across arms.
These are explicitly modified extraction variants, not default Microgen
pipelines or Meshers seeds. No MMG or FE call runs.

| Extraction | Nodes | Tets | Oriented-face conflicts | Cap faces with ceiling below 0.1 | Sampled two-sided distance | Relative volume error |
|---|---:|---:|---:|---:|---:|---:|
| Normal clips, then triangulation | 22,199 | 95,776 | 288 | 0 | 0.00249648 | +0.588923% |
| Six-tet background before clips | 34,743 | 123,344 | 0 | 24 | 0.00253049 | +0.719315% |
| Five-tet background before clips | 30,271 | 101,472 | 0 | 24 | 0.00163287 | +0.718972% |

Both variants pass periodic nodes and cap triangles, with coordinate errors
below 1.2e-16. Neither is admitted. Their retained caps have minimum analytic
quality ceiling 0.08710052, below 0.1. Minimum actual corner quality is
0.00050237 for six-tet and 0.00109506 for five-tet extraction. Positive corner
determinants alone would not distinguish the normal arm's 288 face conflicts.

The five-tet construction uses VTK's rectilinear-grid conversion. Microgen's
structured grid has a different point order. First attempt 343 stops at the
exact-coordinate guard before clipping. Its failure remains frozen. A new
regression and explicit bijective connectivity mapping preserve original point
order and samples in retry 344. VTK's float32 output coordinates are replaced
with the exact original coordinates before scalar clipping, after checking
their expected point correspondence. No field value or sample is resampled.

Job 342 passes four controls. Retry 344 passes seven, including permuted point
order and rejection of nonrectilinear backgrounds. Sixteen combined extraction,
progress, incidence, and admission tests pass after the retry. Both queued
progress audits pass. Worker peaks are at most 0.530 GiB, with 120-second stage
timeouts, one numerical thread, and unchanged 250,000-tet and 8 GiB ceilings.
These small diagnostic costs are not end-to-end workflow timings.

Do not send either fixed-cap seed to MMG. The next useful check is local
cap-triangulation feasibility on the coherent five-tet seed, preserving its
periodic cap boundaries. It must establish a valid replacement before any
volume remeshing or FE cost. No grid-resolution or sizing sweep is selected.
This experiment identifies an extraction-order effect, not the exact upstream
VTK operation responsible for all conflicts or a generally repaired pipeline.

Records are `runs/g3-gyroid-preclip31-v1`,
`runs/g3-gyroid-preclip31-five-v1`, and
`runs/g3-gyroid-preclip31-five-v2`. The first two retain exact driver and test
snapshots matching their original manifest hashes.

## Gyroid fixed-boundary interior feasibility, 20 October 2026

After the repeated Schwarz comparison, return to the missing admitted Gyroid
control. Do not repeat audit 212: it already finds no cap or wall face below
the analytic quality ceiling. Instead job 334 tests only the worst two-pass
tet's free interior vertex, node 30242, over its complete 26-tet star. All
neighbors and boundary coordinates stay fixed. The target is quality 0.1,
with a positive 1e-9 numerical guard. Fixed opposite faces give upper bound
0.14561, which is necessary but not sufficient for a feasible position.

Four controls pass. One constrained minimum-displacement solve finds a point
with minimum star quality 0.100000001 in nine iterations and 0.00842 seconds.
The move is 0.0116271. Every signed determinant stays positive, and full-cell
signed volume changes by only -1.48e-16. The source boundary and connectivity
are unchanged. Global minimum quality rises only to 0.04554, so the whole
mesh remains rejected and receives no FEA. This is one feasible local move,
not a repaired workflow or an admitted Microgen comparison.

Exact edge counting gives 325,911 nodes for a straight Tet10 lift of this
two-pass mesh, above the unchanged 250,000-node stage. The earlier one-pass
FEM43 precursor would lift to 153,010 nodes and fits that stage. Every bad
tet in both precursors has at least one interior vertex. This motivates one
bounded repair attempt on the smaller precursor, not a resource-cap increase.

## Incomplete Gyroid interior-only repair

Job 335 retains the actual VTK43/one-pass FEM precursor, its fixed scalar
specimen, all boundary coordinates and connectivity. It explicitly adds SciPy
constrained interior relocation, rather than claiming this is a Microgen
default or remeshing a Meshers seed. Each accepted move must put its complete
incident star above 0.100001 and strictly reduce the number below 0.1.
The limits are two sweeps and 400 attempted vertices. No stiffness solve is
selected before full admission.

Six controls pass, including boundary immobility and an all-boundary bad tet
that must remain rejected. The trial accepts 108 moves in its first sweep.
Bad elements fall from 222 to 49, but minimum quality remains 0.0079758.
The second sweep accepts no moves. There are 210 attempts total: 108 feasible,
92 without a feasible point found, and ten fixed-opposite-face rejections.
Attempt counts are not distinct-element counts. A failed optimization is not
an impossibility proof.

Boundary arrays remain identical, signed volume is unchanged, sampled
distance stays 0.00151229 and relative volume error stays +0.36758%.
Periodic cap triangles match, with coordinate errors below 1e-15. All signed
determinants remain positive. Repair takes 9.24 seconds; the full worker takes
12.58 seconds and 0.3531 GiB. These are saved-mesh diagnostic costs, not live
workflow timings. Status remains `repair-incomplete`, `qualified=false`, and
no FE or Tet10 accuracy claim follows.

Stop this two-sweep recipe. More sweeps would repeat an unchanged state.
Before a coupled relocation or topology change, diagnose one remaining worst
star with an independent convex-feasibility/slack check and supporting-plane
bound. Separate failure to find a position from a checked fixed-neighbor
obstruction. Do not weaken quality or select another MMG grading value from
these results. The admitted actual Gyroid comparator is still missing.

Records are in `runs/g3-gyroid-interior-star-feasibility-v1` and
`runs/g3-gyroid-firstpass-fem43-interior-relocation-v1`. Maintained drivers are
`interior_star_probe.py` and `interior_relocation_trial.py`.

## Declared single-star convex check

The next bounded diagnostic freezes the partially repaired mesh and selects
its worst tet's single free interior vertex. For each incident tet, let `D`
be its signed determinant and `S` its sum of squared edges. The existing
quality is `12*sqrt(3)*D/S**1.5` for positive orientation. At target `q`, the
constraint `D - q*S**1.5/(12*sqrt(3)) >= 0` is concave in this vertex: `D`
is affine, while `S` is a nonnegative convex quadratic. Each tet uses a
positive fixed normalization to avoid scale-dependent solver tolerances.

First minimize a common nonnegative slack, with analytic derivatives. Then
form concave upper tangents at that point. A separate small linear program
selects nonnegative weights summing to one. Recompute their weighted tangent's
maximum over the unchanged comparison box. A value below -1e-8 is the declared
floating-point obstruction screen. It does not rely on the slack optimizer's
success flag. Keep weights, expansion point, coordinate normalization, box,
bound and residuals for review. If no bound or feasible point resolves the
question, report it as unresolved rather than infer impossibility.

Controls cover derivatives, direct MMG-quality equivalence, global tangent
inequalities, a feasible star, and individually good faces with incompatible
joint constraints. No source mesh changes, another repair sweep, topology
change, FEA or reference promotion occurs in this diagnostic. A checked
fixed-neighbor obstruction would not rule out coupled relocation or remeshing,
and floating-point checks are not interval or continuum certification.

## Declared oriented-facet consistency audit

The local 2-to-3 preflight stops rather than treating a non-oriented source
cavity as valid. Its first offending pair, cells 1297 and 1298, shares face
`[704, 9126, 9369]`. Both corner determinants are positive, but their opposite
vertices have same-signed scalar triple products against that face. This is
a necessary source-consistency issue, not a successful topology candidate.
Preserve the failed preflight and do not bypass its cavity guard.

Before another repair, audit the original VTK43 seed, MMG stage and final
output, partially repaired output, and representative admitted Schwarz,
lattice and native controls. Check oriented incidence, nonmanifold and
duplicate faces/cells, and agreement with the saved boundary. Synthetic
controls include overlapping positive-volume tets that must fail. This is
one supplemental validity study, not a mesh-parameter sweep or a certificate
against every possible curved or nonlocal intersection. Frozen qualification
flags remain unchanged in place; any affected earlier claim must be stated
explicitly after reading the audit.

## Checked obstruction and invalid source lineage

Job 336 passes six controls and selects the remaining worst node 704 on the
partially repaired one-pass FEM source. Its 22-tet slack problem converges
with normalized slack 0.00178742. The independently checked weighted tangent
upper bound is -0.0017874161, well below the declared -1e-8 screen. A second
implementation uses scalar triple products, all eight box corners, and a
64-bit mantissa. It agrees within 1.3e-16. Three nonzero weights support
cells 1297, 1298 and 6378; every opposite vertex in that support is a boundary
vertex. Moving unrelated interior vertices cannot alter those three constraints
while boundary coordinates and these tets remain fixed.

This is an algebraic fixed-neighbor obstruction, not a general remeshing or
continuum impossibility claim. The subsequent source audit changes its physical
interpretation. Job 337 passes three synthetic cavity controls but stops at an
incoherent original cavity. Cells 1297 and 1298 share face `[704, 9126, 9369]`.
Both have positive determinants and their opposite vertices lie on the same
physical side. No flips are applied. Do not bypass that source guard.

Job 338 passes three incidence controls and audits twelve saved meshes:

| Saved mesh group | Same-orientation shared faces |
|---|---:|
| Gyroid VTK43 seed, before MMG | 684 |
| Gyroid one-pass FEM MMG stage, final output, partial repair | 146 each |
| Gyroid two-pass output | 0 |
| Admitted Schwarz linear and stronger Tet10 parent | 0 |
| Representative native Gyroid, Schwarz, spinodoid and lattice | 0 |
| Admitted tuned lattice Microgen parent | 0 |

All twelve have positive corner determinants, no duplicate tets, and matching
saved/extracted boundaries. Those checks alone therefore do not detect the
Gyroid conflicts. The one-pass source already has this problem before the
interior-relocation trial. Retire repair work on that source. The two-pass
Gyroid remains coherent under this necessary check but still fails quality.
The VTK seed locates the issue upstream of MMG; it does not yet isolate which
clipping or tetrahedralization operation introduces it.

Job 339 checks every archived Schwarz timing arm, not just representatives.
All 24 arm meshes exactly match their audited controls. All six manual Tet10
lifts retain the audited parent corners and surfaces. Earlier Schwarz timing
conclusions are unaffected. These are necessary corner-incidence checks, not
certificates against curved or nonlocal intersections.

## Fail-closed future admission and next extraction check

The common `geometry` helper now records `facet_consistency` and requires
oriented-facet consistency before admitting any newly measured mesh. The
0.002 distance, 1% volume, 0.1 quality, positive sampled Jacobian and projection
conditions are unchanged. A high-quality numerical row with incoherent
incidence is red before the predicate fix and green afterward.

Queued check 340 passes sixteen controls. Fresh invalid-Gyroid and admitted-
Schwarz geometry replays reproduce every legacy metric while rejecting the
former and retaining the latter. Frozen rows and timings are not rewritten.
Its old-study verifier still checks all 96 mixed-load comparisons and 205
hashes. Exact measured driver archives and explicitly declared Git source
`2f0363f` account for thirteen source records after code changes. Mesh or field
data cannot use that fallback. Sixty-eight focused checks pass afterward.

The selected next extraction preflight compares unmodified clipping against
triangulating the same sampled background before Microgen's two scalar clips.
Use only resolution 31, not a grid sweep. Preserve the same analytic Gyroid,
offset 1.1, sample coordinates and scalar arrays; label the variant explicitly
non-default. Bound background and returned material tets separately by the
existing 250,000 cap. Check oriented incidence, cap feasibility, periodicity
and independent geometry before any MMG call or FEA. This is a hypothesis
about the upstream extraction, not a completed fix or an admitted comparator.

Records are in `runs/g3-gyroid-convex-star-feasibility-v1`,
`runs/g3-gyroid-convex-star-witness-audit-v1.json`,
`runs/g3-gyroid-convex-star-witness-audit-v2.json`,
`runs/g3-gyroid-two-to-three-preflight-v1`,
`runs/g3-facet-incidence-lineage-controls-v1`,
`runs/g3-schwarz-admitted-facet-replay-v1.json`, and
`runs/g3-facet-admission-check-v1`. The failed cavity preflight and earlier
source snapshots retain their exact recorded hashes.

Queued check 341 also verifies the retired repair entrypoint. An overlapping-
positive-tet source regression is red before its guard and green afterward.
The actual one-pass source now returns `source-incidence-rejected` before
relocation. Seven controls pass; no new mesh, repair or FEA is attempted.
The earlier 108-move diagnostic remains frozen as evidence, not a supported
repair route.

## Spinodoid seed41 feasibility, 19 October 2026

The earlier actual VTK seed41 has 261,650 tets and was skipped solely by the
250,000 seed cap. Its measured generation peak was 0.24 GiB. Job 243 therefore
checks only cap feasibility in an explicit 300,000-tet seed stage, with the
same 8 GiB worker ceiling and 120-second timeout. The new runner option is
restricted to cap-only mode. Default seed, final mesh, and FE caps do not change.

The fixed FFT16 realization and all 52 coefficients remain bitwise unchanged.
Seed41 has 190 cap triangles with adjacent-tet quality ceilings below 0.1;
the minimum ceiling is 0.00015428. Its 2.56-second screen peaks at 0.30 GiB.
Retaining those particular caps cannot satisfy the quality gate, so no MMG or
FE solve follows. This extends the tested obstruction, not a universal
impossibility result. Records are in `runs/g3-spinodoid-seed41-cap-only300k-v1`.

## Required-entity grading diagnostic, 19 October 2026

Saved-output audit 212 finds 189 Gyroid two-pass tets and 222 one-pass FEM tets
below quality 0.1. Every bad tet touches a cap vertex, yet no cap or wall face
has a quality ceiling below 0.1. This supports inspecting near-required-entity
gradation, not a guarantee that sizing alone repairs these meshes. Records
are in `runs/g3-gyroid-mmg-quality-localization-v1`.

The installed MMG CLI supports `-hgradreq`, which the installed Microgen wrapper
does not expose. Job 215 adds only `-hgradreq 1.1` to the actual recorded
VTK43/one-pass FEM recipe, retaining first-pass hausd 0.0015, hmax 0.06,
hgrad 1.3 and unchanged admission gates. The full actual CLI is archived.
This is an explicit override, not Microgen defaults or a Meshers seed.

The result has 32,194 nodes and 150,873 tets. Sampled distance 0.0014971,
volume error +0.3494%, and periodicity pass, but minimum quality deteriorates
from 0.0079758 to 0.00018695. No FEA runs. Stop this variant. Its unsuccessful
option is removed from the maintained runner; restore the frozen driver in
`runs/g3-gyroid-first-pass-fem43-required11-v1` to reproduce it. The admitted
Gyroid Microgen comparator remains missing.

The opening tables and failure list record the initial 17 October experiment.
Later sections include native fixes, admitted spinodoid and lattice meshes,
and tuned Microgen controls. For current conclusions, start with the
"Current G3 continuation" section in `PLAN.md`; the four-fixture goal remains
unmet.

The real Microgen VTK and CLI MMG workflows ran on all four fixed
specimens. This is not the earlier MMG-on-a-Meshers-seed comparison.
There is no demonstrated four-fixture win for Meshers.

## Gyroid first-pass FEM diagnostic, 18 October 2026

Job 167 removes `nofem=True` from the actual Microgen first MMG call. The
installed MMG help defines this flag as disabling forced finite-element mesh
creation. This is an explicitly modified one-pass diagnostic, not a default
workflow. Keep VTK resolution 43, first-pass `hausd=0.0015`, `hmax=0.06` and
`hgrad=1.3`. Omit the second pass exactly as recorded in the frozen driver.

The returned mesh has 22,482 nodes and 95,904 tets. Two-sided sampled distance
is 0.0015123 and volume error is +0.3676%, but minimum quality is 0.0079758,
far below 0.1. No stiffness solve runs. Removing the flag does not establish
an admitted comparator. Intermediate audits are charged, so this diagnostic's
28.56 seconds is not a fair live timing comparison. See
`runs/g3-gyroid-first-pass-fem43-v1`.

The useful negative results are:

- Meshers rejects all tested periodic spinodoids at quality 0.1.
- Meshers rejects the smooth lattice's raw-field periodicity, although
  Microgen produces matching discrete caps and completes diagnostic FEA.
- Curved Tet10 is much more accurate than linear Tet4 on the two admitted
  fixtures, but its strict single-thread CG solve costs more at the tested
  node counts. No matched-accuracy workflow speed comparison was run.
- Default Microgen MMG coarsens the default VTK meshes and worsens boundary
  distance on these fixtures. Every default control fails admission.
  This does not establish superiority over a tuned Microgen workflow.
- Spinodoid and smooth-lattice stiffness errors remain unproven. Neither
  has a justified, admitted fine reference in this experiment.

Native meshing code is unchanged. The branch started at `2e35c7e441b5ccaaa42ea324502ca1a2afabef1b`.
The support checkouts remain at `2b1ea80b36d4ba03842a7d9087638abf3868fe3e`
and `602dc900cd028d205b36e29dcdb5b4a9cf4f581a`.
The exact specimens, defaults, resource limits and metric definitions are
in `MICROGEN.md`. Versions and source/native hashes are in the run manifest.

## Counts, geometry and full-tensor disagreement

H is the independently sampled two-sided maximum in a unit box. Tet10
reverse queries use eight subdivisions per curved edge. The subsequent
verification also checks sixteen. These are finite samples, not continuous
Hausdorff bounds. Volume error is signed. Quality is minimum MMG corner
quality, with separate sampled Jacobian checks for Tet10.

The fixed gates are H <= 0.002, absolute relative volume error <= 1%,
quality >= 0.1, positive sampled Jacobians, and matching periodic nodes
and cap triangles. Only the Meshers Gyroid and Schwarz-solid rows below
pass. The reported Microgen stiffness disagreements are diagnostics.

| Specimen and route | Nodes | Tets | Sampled H | Volume error | q_min | Stiffness disagreement |
|---|---:|---:|---:|---:|---:|---:|
| Gyroid, VTK default | 7,808 | 30,642 | 0.00561 | +1.5190% | 0.0000270 | 14.7909% |
| Gyroid, MMG default | 3,527 | 12,282 | 0.00867 | +2.1768% | 0.000483 | 16.9218% |
| Gyroid, Meshers linear 24 | 17,647 | 72,757 | 0.00158 | +0.1891% | 0.1503 | 3.7152% |
| Gyroid, Meshers curved 12 | 19,931 | 10,832 | 0.000573 | -0.0280% | 0.3714 | 0.2951% |
| Schwarz solid, VTK default | 6,288 | 26,636 | 0.00219 | -0.0029% | 0.006997 | 3.6002% |
| Schwarz solid, MMG default | 3,151 | 10,829 | 0.00846 | -0.0551% | 0.05484 | 5.9304% |
| Schwarz solid, Meshers linear 24 | 19,004 | 90,792 | 0.00189 | +0.0086% | 0.3997 | 1.3458% |
| Schwarz solid, Meshers curved 12 | 19,628 | 11,976 | 0.000242 | -0.0002% | 0.4657 | 0.0738% |
| Spinodoid, VTK default | 5,268 | 20,702 | 0.0301* | -0.0053% | 0.00000000118 | No justified reference |
| Spinodoid, MMG default | 5,509 | 22,907 | 0.0354* | +0.0160% | 0.001835 | Solve withheld |
| Smooth lattice, VTK default | 48,353 | 233,532 | 0.00265* | -0.4017% | 0.000000232 | No justified reference |
| Smooth lattice, MMG default | 6,709 | 27,717 | 0.0106 | -2.8544% | 0.009052 | No justified reference |

An asterisk marks unresolved forward projections. There are 99 for
spinodoid VTK, 23 for spinodoid MMG and 288 for lattice VTK. Their H values
are diagnostics, not trustworthy complete maxima. The geometry gate fails.
Four Sobol volume scrambles have full spreads of 0.000191 for Gyroid,
0.000427 for spinodoid and 0.000156 for the lattice. The Schwarz-solid
reference volume is exactly 0.5. These spreads are not confidence intervals.
In particular, the tiny reported fine Gyroid volume bias is not resolved
against the sampling uncertainty.

Meshers spinodoid rejects linear 24 and curved 12/16/20 before returning
meshes. Reported minimum qualities are 0.03438, 0.03344, 0.03983 and
0.01583. The unchanged floor remains 0.1. Microgen's MMG result passes its
own node-periodicity check, but the independent cap-triangle check fails
on axis 1, so no homogenized tensor is admitted or solved for that mesh.

All four Meshers lattice configurations reject with
`field/grading is incompatible with periodic axis 0`. Sequential cubic
smooth unions are not associative. Opposite-face scalar differences reach
about 0.0000764 in void. A separate 131,072-probe face check found zero
material-membership differences and zero field differences within 0.03
of the material boundary. The Microgen meshes pass the actual periodic
node and cap-triangle checks. No field rewrite or gate relaxation was used
to make Meshers accept this specimen.

## Same-host wall time

Seconds are single live runs on G3, with one numerical thread. Generation
includes specimen setup and actual VTK/MMG or native generation.
Meshing-only time is `generation_s - specimen_setup_s` in each raw record.
Total time includes geometry audits, serialization, conversion, assembly
and six strict CG solves. It excludes interpreter imports. Fine reference
solves are not charged to a candidate's workflow. These are not repeated
performance estimates and must not be compared with OMEN seconds.

| Specimen and route | Generation, s | Six-load FEA, s | Total workflow, s |
|---|---:|---:|---:|
| Gyroid, VTK default | 0.128 | 12.387 | 13.478 |
| Gyroid, MMG default | 4.600 | 1.977 | 7.302 |
| Gyroid, Meshers linear 24 | 1.245 | 14.836 | 17.480 |
| Gyroid, Meshers curved 12 | 0.282 | 45.558 | 47.080 |
| Schwarz solid, VTK default | 0.068 | 4.365 | 5.236 |
| Schwarz solid, MMG default | 3.050 | 1.096 | 4.819 |
| Schwarz solid, Meshers linear 24 | 0.649 | 13.525 | 15.210 |
| Schwarz solid, Meshers curved 12 | 0.139 | 36.688 | 38.073 |
| Spinodoid, VTK default | 0.080 | 11.466 | 14.846 |
| Spinodoid, MMG default | 4.035 | Withheld | Incomplete |
| Smooth lattice, VTK default | 0.449 | 190.207 | 196.387 |
| Smooth lattice, MMG default | 12.623 | 5.086 | 18.669 |

The short default MMG workflows use fewer nodes and fail the accuracy
gates. Their times are not matched-accuracy wins. Conversely, the curved
Meshers results do not prove a faster end-to-end workflow than Microgen.
Only Jacobi-preconditioned CG was tested, not a direct or multigrid solver.

## References and geometry/order controls

The score is the maximum absolute generalized eigenvalue of
`C_candidate - C_reference` against `C_reference`, across all strain
combinations. It is not a maximum of six diagonal tests.

The reference extrapolates curved 16/20 tensors with assumed order four.
For Gyroid, the fine-to-extrapolated disagreement is 0.04273%, the
12/16 versus 16/20 extrapolation difference is 0.01665%, and the archived
linear-64/80 extrapolation differs by 0.00716%. For Schwarz solid, these
first two values are 0.01300% and 0.01296%. Its archived quadratic tensor
agrees to numerical precision, which checks reproduction but is not an
independent-method reference. Both fine inputs pass the geometry gates.
No value is a certified continuum error bound.

Saved 12-cell skeletons provide separate approximation-order and geometry
controls. They have identical connectivity and material for the linear
versus straight-quadratic comparison. Curved versus straight quadratic
then changes only the geometric mapping at fixed polynomial order.

| Specimen | Linear corners | Straight Tet10 | Curved Tet10 |
|---|---:|---:|---:|
| Gyroid, stiffness disagreement | 13.2473% | 1.6762% | 0.2951% |
| Schwarz solid, stiffness disagreement | 5.4935% | 0.2645% | 0.0738% |

Relative to the fine reference, the worst magnitudes of the tensor-change
eigenvalues are 11.5777% for Gyroid's order change and 1.3987% for its
geometry change. For Schwarz solid they are 5.3065% and 0.2256%.
The signed tensor changes add exactly. Their worst scalar magnitudes do
not generally add, and neither is a certified discretization-only error.
The straight controls fail the common geometry gate and remain diagnostic.
Their replay timings are not presented as live generation times.

## Evidence and checks

Jobs 4 and 5 finished sequentially and both passed their progress audits.
The first full pilot used at most 2.37 GiB sampled worker RSS. No worker
hit its 8 GiB address-space or aggregate RSS ceiling, 1,200-second timeout,
250,000-tet cap, or 100,000-node solve cap. No unrelated jobs were changed.

There are 28 final configurations, including controls and fine references:
eight admitted solves, eleven unqualified diagnostic solves, eight native
rejections, and one MMG topology rejection. Nineteen complete six-load
tensors are saved. Public native Tet10 ordering also passes a real solid
cube affine-patch solve. All numerical solver residuals are below 2e-9.

Job 7 rechecked all twenty saved meshes and nineteen tensors. It verified
frozen/inherited source hashes, positive linear orientation, exact counts,
periodic pairs for admitted meshes, mixed-load eigenvalues and solver
records. Independent five-point Tet10 volume integration agrees with the
original fifteen-point result within 7.8e-15. Sixteen-edge boundary
tessellations retain admission for all six fine/coarse curved meshes.
Every original artifact hash is unchanged. Historical saved displacement
archives were not transferred and were not revalidated.

The first verifier attempt, job 6, stopped after two controls because its
support import path was not initialized. It changed no experiment files.
The corrected read-only verifier passed as job 7.

Persistent local artifacts are at:

- `/home/kmarchais/dev/meshers/experiments/adaptive/runs/g3-microgen-defaults-v3`
- `/home/kmarchais/dev/meshers/experiments/adaptive/runs/g3-microgen-defaults-v4`
- `/home/kmarchais/dev/meshers/experiments/adaptive/g3-verification.json`
- `/home/kmarchais/dev/meshers/.audit/g3-microgen.tsv`

Text records and hashes are committed. Mesh archives and installed-source
snapshots remain local. No historical artifact or admission flag was edited.

To repeat the read-only verification on G3, keep the local mesh archives,
source snapshots and support checkouts, then run from the repository root:

```sh
qrun -l 'Verify saved G3 geometry tensors and hashes' --progress report \
  --interrupt none --interrupt-reason 'Read-only verification with one final report' \
  --check-start --wait env \
  LD_LIBRARY_PATH=/home/kmarchais/dev/meshers/.local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu \
  OPENBLAS_NUM_THREADS=1 timeout 300 \
  .venv/bin/python experiments/adaptive/verify_microgen_workflow.py \
  --run experiments/adaptive/runs/g3-microgen-defaults-v4 \
  --output /tmp/g3-microgen-verification.json
```

Set the library path before Python starts. A clone with only the committed
text evidence cannot replay the mesh checks without those local archives.

## What remains open

This is a reviewable diagnostic result, not completion of a four-fixture
matched-accuracy study. Spinodoid and smooth-lattice references are missing.
Projection failures also leave some geometry distances unresolved. A tuned
Microgen accuracy/cost comparison, repeated workflow timings and local
stress recovery have not been measured. The experiment does not justify
batching changes to the mesher or claiming continuum certification.

Kevin's later requested order remains unchanged. Direct fedoo Tet10 now
has G3 execution evidence, but public Microgen Tet10 conversion and its
periodic curved-node integration still need a separate integration result.

## Retained cap triangles limit quality

A subsequent read-only audit, job 10, explains part of the default Microgen
quality failure. The VTK seeds contain thin box-cap triangles that the MMG
workflow retains. Their shape alone limits the quality of any adjacent
tetrahedron, regardless of its fourth vertex or the requested MMG size.

For triangle area A and sum of its three squared edge lengths S, the quality
definition used in this study satisfies `q_tet <= 4*sqrt(3)*A/S`. To derive
the bound, place the fourth vertex at height h above the face and at planar
distance r from its centroid. The tetrahedron volume is A*h/3, and the sum
of all six squared edge lengths is `4*S/3 + 3*h*h + 3*r*r`. Quality is
maximized at r = 0 and `h*h = 2*S/9`, which gives the stated bound. Additional
constraints on the fourth vertex can only lower that maximum.

| Specimen | Retained cap triangles below the 0.1 quality ceiling | Smallest ceiling |
|---|---:|---:|
| Gyroid | 60 | 0.009126 |
| Schwarz solid | 0 | 0.138143 |
| Spinodoid | 46 | 0.009166 |
| Smooth lattice | 96 | 0.017177 |

The audit matched every seed cap triangle in the saved MMG mesh at coordinate
tolerance 1e-10. Spinodoid MMG also has one extra cap triangle, consistent with
its earlier topology rejection. The analytic bound was checked numerically
on regular and thin triangles, including a fourth vertex that attains it.
All eight input mesh hashes are unchanged. The job finished before the
progress-start audit could observe it; its five result records are available.

Thus tuning only MMG sizes while retaining these particular seed caps cannot
pass the existing quality gate for Gyroid, spinodoid or lattice. Schwarz has
no obstruction from this bound, though its actual mesh still fails quality.
The bound does not rule out a different VTK extraction grid or a separately
identified boundary-repair workflow. It says nothing about stiffness error.

The executable audit is `cap_quality_bound.py`; its output and input hashes
are in `g3-cap-quality-bound.json`. Run it through the shared queue with
`--source experiments/adaptive/runs/g3-microgen-defaults-v4` and a fresh
`--output` JSON path.

## Tuned Schwarz baseline pilot

Jobs 11 through 14 tested one candidate at a time using the actual Microgen
VTK generation and two-pass MMG function. Each completed result was read before
selecting the next candidate. The specimen remains Schwarz P solid at offset
zero. These are explicitly tuned runs, with `hausd=0.0005`, `hmax=0.06` and
`hgrad=1.3`, and extraction resolutions 28, 32 and 20 points per axis.
Microgen applies the size arguments to its second pass. The first pass still
receives only `nofem=True`, as recorded in every call list.

At resolution 28, all five built-in first-pass attempts failed. Job 11 exposed
the failure, but the Microgen wrapper hid the captured MMG text. Job 12 repeated
that candidate with exception-output capture. MMG reports a worst element of
zero quality and stops in `MMG5_minQualCheck`. No cells were deleted or repaired.

| Resolution | Output nodes | Output tets | Sampled H | Signed volume error | Minimum quality | Outcome |
|---|---:|---:|---:|---:|---:|---|
| 28 | None | None | Unavailable | Unavailable | Unavailable | MMG first-pass failure |
| 32 | 28,472 | 131,516 | 0.005800 | -0.01289% | 0.02092 | Geometry and quality rejection |
| 20 | 22,928 | 108,979 | 0.005347 | -0.07307% | 0.03994 | Geometry and quality rejection |

Both completed meshes pass periodic node and cap matching and have no unresolved
projection samples. Resolution 32 has 48 cap triangles whose adjacent-tet quality
ceiling is below 0.1. Resolution 20 has none, but its actual element quality still
fails. Their complete geometry-audit workflows take 12.98 and 8.80 seconds on G3.
FEA is withheld because admission fails, so these are not total meshing-plus-FEA
times and no new stiffness errors are reported.

This bounded pilot does not find an admissible tuned competitor. It does not
rule out other extraction grids, controls applied to the first MMG pass, or
separately identified periodic boundary repair. The actual first-pass geometry
has not been archived or audited, so its contribution to the final boundary
error is still a hypothesis.

The runner is `microgen_tuned.py`. It uses one thread, an 8 GiB worker ceiling,
a 600-second worker timeout, a 250,000-tet cap and a 100,000-node solve cap.
Its manifests link the frozen reference inputs by hash. Records are under
`runs/g3-microgen-tuned-schwarz28-v1`, `runs/g3-microgen-tuned-schwarz28-v2`,
`runs/g3-microgen-tuned-schwarz32-v1` and `runs/g3-microgen-tuned-schwarz20-v1`.
All runs preserve the same independent admission gates.

## Audited intermediate MMG meshes

Job 18 repeats the resolution-20 tuned candidate without changing either MMG
call. It archives and independently measures the VTK seed and each MMG output.
The final mesh reproduces job 14. Intermediate audits add work, so this run is
explicitly excluded from comparable workflow timings.

| Stage | Nodes | Tets | Sampled two-sided H | Minimum quality |
|---|---:|---:|---:|---:|
| VTK seed | 6,288 | 26,636 | 0.002190 | 0.00700 |
| First MMG pass, nofem | 3,118 | 10,519 | 0.007183 | 0.06860 |
| Second MMG pass, tuned | 22,928 | 108,979 | 0.005347 | 0.03994 |

This confirms first-pass boundary degradation for this candidate. The seed
already narrowly misses the 0.002 distance gate, and the second pass does not
recover the lost geometry. It does not prove that every extraction grid behaves
the same way. Records are in `runs/g3-microgen-stage-audit-v1`.

The next diagnostic uses resolution 24 and overrides the first call with
`hausd=0.0001`, `hmax=0.06` and `hgrad=1.3`, keeping `nofem=True`. It still calls
Microgen's VTK generation and both MMG commands, but this override is not exposed
by its public remeshing function. The runner labels this modified route
`microgen-vtk-mmg-first-pass-override`, never as the default pipeline.

Job 19 completed that first pass with 346,003 vertices and 1,730,349 tetrahedra.
The stage audit stopped at the unchanged 250,000-tet gate before the second pass
or FEA. The VTK seed itself has H=0.001500 but minimum quality 0.000344.
This is an over-refinement rejection, not an MMG executable failure. The frozen
`result.json` uses the overly broad status `microgen-mmg-failed`; its error text
records the actual intermediate cap. The runner now distinguishes that status
as `intermediate-element-cap`. Raw first-pass meshes stay local and ignored.
Records are in `runs/g3-microgen-firstpass24-v1`.

Job 20 relaxes only the first-pass `hausd` override to 0.0005. The first pass
now has 129,222 tets and the second 168,621, both inside the cap. Final two-sided
distance is 0.001440 and signed volume error is -0.01052%, so the geometry loss
is resolved for this specimen and extraction grid. Minimum quality remains
0.02402. There are 72 retained cap triangles whose adjacent-tet quality ceiling
is below 0.1, making further size tuning insufficient on this seed. Periodic
checks pass. FEA remains withheld. Records are in
`runs/g3-microgen-firstpass24-v2`. These instrumented timings include intermediate
audits and are not used for a workflow speed comparison.

## Nearby VTK seed grids

Jobs 22 and 23 screen Schwarz resolutions 21, 22, 23, 25, 26 and 27 before
spending time in MMG. The first screen failed while serializing nonfinite quality
from zero-volume seed cells. The corrected screen records those cells as a
rejection instead of deleting them or admitting nonfinite metrics.

| Resolution | Zero-volume tets | Caps with quality ceiling below 0.1 |
|---|---:|---:|
| 21 | 538 | 48 |
| 22 | 245 | 33 |
| 23 | 0 | 48 |
| 25 | 459 | 81 |
| 26 | 0 | 96 |
| 27 | 0 | 48 |

None offers a feasible set of unchanged retained caps under the fixed quality
gate. This screen does not rule out other grids or explicit boundary repair.
It does rule out obtaining an admitted competitor merely by tuning MMG on
these six seeds. Records and the runnable driver are in
`runs/g3-microgen-seed-screen-v1` and `runs/g3-microgen-seed-screen-v2`.

## First-pass override on the default Schwarz seed

Job 61 applies the explicit first-pass hausd=0.0005 override to the original
20-point VTK seed, keeping the actual two Microgen MMG calls. It produces 35,084
nodes and 172,800 tets in 15.25 seconds including audit. Distance improves from
0.005347 to 0.002098, but still fails 0.002. Minimum quality is 0.02980 and
41 elements fail 0.1. Volume error is -0.02849%, periodic checks pass, and no
retained cap triangle has an adjacent-tet quality ceiling below 0.1. No FEA runs.
Records are in `runs/g3-microgen-tuned-schwarz20-firstpass-v1`.

Job 62 localizes the worst probes to box-cap curves and the ten worst elements
to tets adjoining a box cap. Its optional detail metrics reproduce every prior
aggregate geometry metric exactly. Jobs 63 and 64 match ten rejected curve
probes to edges retained from the default VTK seed. Job 64 recomputes the probes
from saved mesh connectivity and barycentric weights rather than trusting the
diagnostic coordinates. Endpoint matching errors are below 4.2e-16.

For a worst probe at (0.2479023, 0, 0.5), the Schwarz field stays positive
throughout its coordinate box of radius 0.002. The analytic lower estimate is
0.0005350, from the minimum cosine in each coordinate interval. The Euclidean
ball is contained in that coordinate box. Retaining this chord therefore leaves
an obstruction to the distance gate, regardless of interior MMG sizing. This is
a geometry argument, not FE or local stress certification. Records are in
`runs/g3-microgen-firstpass20-diagnostic-v1`, `runs/g3-microgen-fixed-curve-v1`
and `runs/g3-microgen-fixed-curve-v2`.

Before adding a boundary-repair algorithm, the next bounded experiment screens
the actual VTK seeds at resolutions 28, 30, 34, 36, 38 and 40 for degeneracy,
retained-cap quality ceilings and independent geometry. The specimen and fixed
admission gates remain unchanged. These are tuned extraction grids, not defaults.

Job 65 finds cap-feasible seeds at 30 and 38 points. Neither has a zero-volume
tetrahedron or degenerate cap. Their seed distances are 0.0009615 and 0.0005868,
and cap quality ceilings are at least 0.1448 and 0.1507. Raw seed tet quality
still fails 0.1, so neither is an admitted comparator before remeshing. The
other four grids have cap obstructions, and 28, 34 and 40 also contain
zero-volume tets. Records are in `runs/g3-microgen-seed-screen-large-v1`.

Select the smaller 30-point seed for the next actual two-pass Microgen/MMG
experiment, with the explicit first-pass 0.0005 override. Intermediate audits
enforce the 250,000-tet cap. Their extra cost makes this a diagnostic run, not
a comparable workflow timing.

Job 66's actual two-pass workflow retains the 30-point seed distance, H=0.0009615.
Its first pass has 139,788 tets and q_min=0.10623, passing geometry and quality.
The second pass has 163,887 tets and reduces q_min to 0.08494. The final mesh
therefore fails admission and receives no FEA. Intermediate diagnostics are in
`runs/g3-microgen-firstpass30-v1`.

The next diagnostic retains the measured first-pass output and omits the second
native call. This is explicitly a one-pass Microgen VTK/MMG override, not the
default or a public Microgen API option. The API's output cleanup and periodic
checks remain, followed by the fixed independent admission tests. If admitted,
its full six-load stiffness and live workflow time can expose a stronger
competitor than the failing two-pass controls.

Job 67's explicitly modified one-pass route passes all independent geometry,
quality and periodic gates with 29,713 nodes and 139,788 tets. H=0.0009615 and
q_min=0.10623. It completes six converged loads in 34.98 seconds and a full live
workflow in 44.35 seconds. Worst mixed-load disagreement is 1.13403% against the
justified Schwarz reference, so it still misses the 1% stiffness target. This
is an admitted geometry competitor, not a matched-accuracy winner. Records are
in `runs/g3-microgen-singlepass30-v1`.

Select the other screened cap-feasible seed, 38 points, for the next one-pass
candidate with unchanged MMG settings. Keep the standard two-pass failures and
the explicit omission of its second native call visible in the records.

Job 68's 38-point one-pass alternative rejects at minimum quality 0.09785,
despite H=0.0005868 and matching periodic caps. Its 31,172 nodes and 141,591
tets receive no FEA. Records are in `runs/g3-microgen-singlepass38-v1`.
The next candidate returns to the geometry-admitted 30-point route and reduces
only first-pass hmax from 0.06 to 0.055, seeking the missing stiffness accuracy
without changing the seed, quality floor or other MMG settings.

Job 69's hmax=0.055 candidate rejects quality at 0.08943. It has 31,949 nodes,
153,277 tets and unchanged H=0.0009615. No FEA runs. Records are in
`runs/g3-microgen-singlepass30-h055-v1`.

The next diagnostic replaces the original second level-set call with MMG's
interior optimization, `optim=True, nosurf=True, noinsert=True`. The installed
CLI documents those switches as mesh optimization, no surface modifications,
and no point insertion or deletion. It is an explicit modified two-pass recipe,
not a public Microgen default. The purpose is to repair measured quality defects
without coarsening the refined first-pass mesh or changing its surface. All
output admission gates and the full stiffness target remain fixed.

Job 70 fails all five API retries before optimization. MMG reports that its
optimization option cannot be used with an input metric. The first pass leaves
an automatically discovered sibling `.sol` file. This is an unsupported call
combination, not evidence about optimized mesh quality. Records, including raw
native errors, are in `runs/g3-microgen-interior-optim30-h055-v1`.
The corrected explicit optimization variant hashes and removes only that fresh
temporary metric before calling MMG. No frozen input or saved output is removed.

Job 71's corrected optimization runs but worsens minimum quality from 0.08943
to 0.08610 at the same 31,949 nodes. It remains rejected and receives no FEA.
Records are in `runs/g3-microgen-interior-optim30-h055-v2`. Discard this modified
optimization recipe rather than weaken quality admission. The next candidate
returns to the original one-pass variant at hmax=0.05, with intermediate audits
to enforce the unchanged 250,000-tet cap before any solve.

Job 72's hmax=0.05 one-pass variant passes geometry, quality and periodic gates
with 32,070 nodes and 154,012 tets. H=0.0009615, volume error=-0.01016%,
q_min=0.10701, and all six loads converge. Worst mixed-load disagreement is
0.92806%, below 1% against the justified Schwarz reference. This is the first
tested Microgen-derived variant to meet both admission and the stiffness target.
It is explicitly a modified one-pass recipe, not the default two-pass pipeline.
Records are in `runs/g3-microgen-singlepass30-h050-v1`.

Its 48.45-second total includes intermediate audits, so it is not a comparable
workflow timing. Before repeated timing, test a coarser eight-cell curved Meshers
candidate. The existing 12-cell candidate is much more accurate than the target,
and its FE solve dominates total cost. Independent geometry and stiffness
admission must still pass. Paired live timing must use the same six-load
formulation and solver, not mix the stress-controlled Microgen screen with the
strain-controlled Meshers admission runner.

## Schwarz matched-target timing protocol

Job 73's eight-cell curved Meshers candidate passes independent admission with
6,864 nodes, 3,852 Tet10 elements, H=0.0009637, volume error=-0.0008914% and
q_min=0.47840. All sampled Jacobians are positive and periodic nodes and caps
match. Its complete six-load tensor differs from the justified reference by
0.25870% in the worst mixed load. Its single pilot takes 7.84 seconds total.
Job 74 repeats reverse distance at eight and sixteen subdivisions, retaining
the same admitted maximum. The subdivision difference drops to 0.0001922.
Records are in `runs/g3-schwarz-curved8-pilot-v1` and
`runs/g3-schwarz-curved8-distance-v1`.

The paired experiment declares one warm-up per arm and five live repetitions,
alternating arm order between pairs. It regenerates both specimens and meshes,
audits geometry and periodicity, serializes meshes and displacement fields,
and solves six prescribed engineering strains with the same strict CG helper.
There is no saved-mesh timing replay. Fresh worker processes use one thread,
the hard 8 GiB ceiling, 250,000-tet and 100,000-node caps and 600-second timeout.
Imports and fine-reference computation are excluded consistently from reported
workflow time. Every repetition must meet all admission gates and 1% stiffness
disagreement before any timing summary is accepted. Failed repetitions stop the
experiment rather than disappear from the sample.

The Microgen arm is the explicit one-pass 30-point VTK recipe, hausd=0.0005,
hmax=0.05 and hgrad=1.3. The Meshers arm is eight-cell nonadaptive curved Tet10,
with one optimization pass and native corner tolerance 1.0. Independent final
distance admission remains 0.002. `schwarz_matched.py` records pinned sources,
the native extension hash, counts, per-load residuals, full tensors and paired
time ratios. A result applies only to Schwarz against this modified recipe,
not to the default two-pass pipeline or all four specimens.

Job 75 completes both warm-ups and all five interleaved pairs. Every run meets
the fixed geometry, quality, periodicity, positive sampled Jacobian and 1%
stiffness gates. Meshers has 6,864 nodes and 3,852 Tet10 elements; the explicit
one-pass Microgen variant has 32,070 nodes and 154,012 Tet4 elements. Their
worst mixed-load disagreements are 0.25870% and 0.92806%, respectively.

Median total workflow times are 7.7916 seconds for Meshers and 40.5036 seconds
for Microgen. Paired speedups range from 5.1155 to 5.2908. Median generation
times are 0.06008 and 7.72084 seconds; median six-load solve times are 6.70432
and 31.04057 seconds. Full workflow times range from 7.6838 to 7.8803 seconds
for Meshers and 40.2588 to 41.2876 seconds for Microgen. Both arms use the same
prescribed-strain solver. These are same-G3 repeated measurements, not OMEN
comparisons. Records are in `runs/g3-schwarz-matched-v1`.

This supports a Schwarz matched-target workflow advantage against the explicitly
modified one-pass comparator. It does not establish a four-fixture victory,
repair the default two-pass Microgen pipeline, or establish lattice stiffness
accuracy. Representative measured displacement fields still require the next
independent stress, energy and reaction audit.

Job 76 independently integrates measured fields from both arms. Stress, energy
and reaction tensors agree within 1.32e-11 relative. Periodic displacement jumps,
including curved mid-edge pairs, are at most 5.56e-17. The representative fields
therefore pass the independent FE consistency check. Records are in
`runs/g3-schwarz-matched-field-audit-v1`. This establishes the stated one-case
advantage, not local stress accuracy or a continuum certification.

## Gyroid retained-cap screening

Jobs 103 and 104 screen actual VTK Gyroid sheet seeds at resolutions 24, 28,
30, 32, 36, 38 and 31, 33, 35, 37, 39, 41. The specimen remains offset 1.1,
the original sheet, and no Meshers seed is involved. All twelve seeds contain
fewer than 250,000 tets and no zero-volume tetrahedra. Eleven have retained
cap triangles whose analytic tetrahedral quality ceiling is below 0.1. MMG
interior sizing cannot admit those fixed cap triangulations.

Resolution 31 has no such cap obstruction, with minimum ceiling 0.15055,
22,199 nodes and 95,776 seed tets. Its seed distance is 0.002496, so it is
not admitted as generated. Select one actual two-pass Microgen/MMG run,
explicitly overriding first-pass hausd=0.0005, hmax=0.05, hgrad=1.3. The
second pass uses the existing tuned hausd=0.0005, hmax=0.06, hgrad=1.3.
Stage audits distinguish a retained-curve obstruction from remeshing losses.
This is an instrumented tuned diagnostic, not default Microgen or a comparable
workflow timing. Solve only if all unchanged admission gates pass. Records are
in `runs/g3-gyroid-seed-screen-v1`, `runs/g3-gyroid-seed-screen-odd-v1` and
`runs/g3-gyroid-tuned31-v1`.

Job 105 stops after the first native MMG pass produces 288,870 tetrahedra,
above the fixed 250,000 intermediate cap. No second pass or FEA runs. Peak
RSS is 0.73 GiB, so this is an element-cap rejection, not a memory failure.
Select hmax=0.06 for one successor while leaving all other settings and gates
unchanged. Do not enlarge the cap or report the aborted route as an admitted
competitor.

Job 106 also exceeds the intermediate cap, at 289,317 tets after the first
pass. Increasing hmax from 0.05 to 0.06 does not reduce this count. The next
single-parameter test changes only first-pass hausd from 0.0005 to 0.001.
The independent 0.002 distance gate, second-pass settings and resource caps
remain unchanged. This tests sizing sensitivity, not a relaxed admission gate.

Job 107 reaches the 600-second worker timeout during its first native MMG
call. The worker peaks at 0.72 GiB; its entire process tree is terminated.
`runs/g3-gyroid-tuned31-v3/worker.json` records return code -9 and
`limit_failure="timeout"`. The child result remains an incomplete `started`
checkpoint, not a returned mesh, geometry rejection or completed workflow.
No second pass or FEA occurs.

Job 108 independently recomputes the frozen VTK seed's geometry and confirms
every prior metric exactly. Its 32 worst forward probes are not on required
cap vertices or edges. This does not prove every cap curve passes, but gives
no retained-cap explanation for the observed maximum. Record this finite
diagnostic separately from the analytic cap-quality obstruction. The next
seed-only screen covers remaining denser resolutions 34, 40, 42, 43 and 44,
seeking both seed geometry and cap feasibility before more native tuning.

Job 109 finds one additional cap-feasible seed at resolution 43. It has 50,883
nodes, 233,408 tets, H=0.001226 and minimum cap-quality ceiling 0.13254. Its
raw minimum tetrahedral quality still fails, so MMG remains necessary. The
other four denser seeds have cap-quality obstructions. Select the genuine
public two-pass route on resolution 43 with its first native call unchanged,
`nofem=True` only, and the public tuned second-pass sizing parameters already
recorded. Capture both stages and solve only if all gates pass. No first-pass
override or omitted pass is used in this candidate.

Job 110 captures the first-pass geometry loss directly. The unmodified native
call coarsens the seed from 233,408 to 25,869 tets, with H worsening from
0.001226 to 0.007353 and minimum quality only 0.005458. The second pass
produces 299,836 tets and exceeds the fixed intermediate cap. No FEA runs.
Select one first-pass sizing override, hausd=0.0015 and hmax=0.06, retaining
the actual two-pass path, stage captures and all independent gates. This
tests the now-observed first-pass loss rather than attributing it by guess.

Job 111 keeps first-pass distance within the gate at H=0.001461, with 85,640
tets, but minimum quality is only 0.02611. The real second pass produces
261,249 tets, again above the unchanged 250,000 cap. No solve runs. Select
one change to the public second-pass hausd, from 0.0005 to 0.0006, keeping
first-pass settings, hmax, hgrad and admission gates unchanged. The runner's
new `--hausd` argument preserves the previous default of 0.0005.

Job 112 completes both actual passes under the cap, producing 47,262 nodes
and 210,503 tets. H=0.001192, volume error=+0.24209%, sampled Jacobians are
positive and periodic nodes and cap triangles match. Minimum quality is
0.04072, below 0.1, despite a first-percentile quality of 0.37106. The route
therefore remains unqualified and receives no FEA. Select a saved-mesh
worst-element audit before another parameter or algorithm change. The
instrumented total of 34.62 seconds is not an admitted timing comparison.

Job 113 identifies 189 tets below 0.1. The worst have three boundary vertices
and one interior apex, making interior-only optimization a testable hypothesis.
Select one explicitly modified third MMG call
with `optim=True`, `nosurf=True`, `noinsert=True`, retaining the same actual
two-pass precursor. Microgen's own serializer marks cap triangles required.
This is a three-pass diagnostic variant, not the public or default two-pass
pipeline. Generation and all audits remain charged if it reaches FEA.

Job 114's third pass changes connectivity and raises the tet count to 210,669,
but minimum quality remains exactly 0.04072. Distance remains 0.001192 and
the mesh still rejects. Standard interior optimization is not a demonstrated
repair. The installed MMG 5.8 help and Microgen wrapper also expose
`optimLES`, explicitly for skewness improvement. Select that single method
change for one further third-pass diagnostic, retaining `nosurf` and
`noinsert`, the same precursor, caps and all admission gates.

Job 115's skewness mode worsens minimum quality to 0.007873 while increasing
the count to 216,087 tets. Distance remains 0.001192. No FEA runs. Reject
both third-pass hypotheses and remove their new controls from the maintained
runner; preserve the frozen experimental drivers and output records. The
useful public `--hausd` parameter remains. Next change only seed resolution
from 43 to the other cap-feasible seed, 31, in the two-pass recipe at
first-pass hausd=0.0015 and second-pass hausd=0.0006. This retains all gates
and tests whether the fixed cap triangulation affects attainable quality.

Job 116's alternate seed also yields a quality-only rejection. Gyroid31 now
returns 48,000 nodes and 218,723 tets, H=0.001528 and minimum quality 0.04059.
Changing between the two cap-feasible grids does not establish an admitted
competitor. No solve is added to the normal gated runner.

Select a separate, explicitly unranked diagnostic on these two returned
meshes to determine whether quality is the only remaining obstacle to the
1% stiffness target. `microgen_stiffness_diagnostic.py` requires positive
sampled Jacobians, matching periodic caps, the original node and tet caps,
and all geometry gates except the already failed quality floor. It retains
`qualified=false` in every output and excludes these cases from admitted
rankings. Six prescribed-strain loads and independent stress, energy and
periodic-jump checks quantify the discrete tensors. Saved-mesh solve times
are not full live workflow timings or evidence of a matched-accuracy win.
This diagnostic does not relax the quality gate for either workflow.

Job 117 completes both six-load diagnostics. Gyroid31 has 2.58544% worst
mixed-load stiffness disagreement and Gyroid43 has 2.29279%, so quality repair
alone would not satisfy the 1% target at these discretizations. Their solve
times are 77.82 and 75.38 seconds, with peak RSS below 1.95 GiB. Independent
stress and energy agree with reactions within 7.26e-12 relative and periodic
jumps are zero. Both outputs remain unqualified. No full-workflow speedup or
matched-accuracy victory is inferred. Records are in
`runs/g3-gyroid-unqualified-stiffness-v1`.

Before further MMG tuning, test whether the admitted curved Meshers Gyroid
control can reduce its initial density from 12 to 10 cells while retaining
all admission gates and six-load accuracy. Use the verified native build and
charge fresh generation, geometry and FEA. This is a Meshers cost/accuracy
pilot, not a comparison against an admitted Microgen competitor.

Job 118's curved10 Meshers pilot passes admission with 13,953 nodes and 7,398
Tet10 elements. H=0.001163, volume error=-0.02437%, q_min=0.33550 and sampled
Jacobians are positive. Periodic coordinates and cap triangles match.
The six-load tensor differs from the justified Gyroid reference by 0.45586%,
below the 1% target. Solve time is 23.24 seconds and live total workflow is
24.71 seconds. The preceding curved12 control had 0.2951% disagreement and
47.08 seconds total, but used the earlier stress-controlled formulation and
native revision. Do not attribute that timing difference solely to density.
These are individual G3 pilots, not repeated paired timings. Independently
audit the curved10 fields before another candidate. Microgen remains unqualified.

Job 119's curved10 field audit passes within 6.79e-12 relative reaction error
and has zero periodic or mid-edge jumps. Select a fresh curved12 control with
the same native build, prescribed-strain helper and other settings as curved10.
This supplies a valid same-helper cost baseline before further coarsening.

Job 120's fresh curved12 control passes admission at 19,931 nodes, H=0.000573
and 0.29515% stiffness disagreement. Its prescribed-strain solve takes 37.47
seconds and live total is 39.17 seconds. Curved10's observed total is 36.9%
lower under this same-helper protocol, while both meet the 1% target. These
are single runs, not a repeated timing result. Select one eight-cell curved
candidate next, retaining every gate and withholding FEA if admission fails.

Job 121's eight-cell Meshers candidate is rejected natively at minimum quality
0.07832, with eight elements below 0.1. No mesh is admitted or solved. This
is a Meshers limitation, not a Microgen-only problem. Select four existing
optimization passes instead of one for a single successor, leaving density,
geometry and the required quality minimum unchanged. No optimizer code change
or diagnostic quality bypass is introduced.

Job 122's four-pass eight-cell mesh passes corner quality at 0.12588 but
fails boundary distance at H=0.007253. Six mid-edge corrections are partially
curved to preserve sampled validity; there are no fully straightened edges.
Positive sampled Jacobians and a good volume error do not compensate for the
distance failure. No FEA runs. Select the nine-cell density midpoint with the
ordinary one-pass setting before proposing a new meshing algorithm.

Job 123's nine-cell candidate passes admission with 11,857 nodes and 6,008
Tet10 elements. H=0.001189, volume error=-0.02633%, q_min=0.14469 and all
sampled Jacobians are positive. All mid-edge corrections are fully curved.
Its primary-reference disagreement is 0.72921%, solve time 19.16 seconds and
live total workflow 20.56 seconds. Job 124 verifies the nine-cell and fresh
12-cell fields within 8.20e-12 relative reaction error, with zero periodic
jumps. Job 125 retains admission after refining reverse tessellation from
8 to 16; H is unchanged and reverse-distance delta is 0.000172.

The nine-cell pilot also stays below 1% across assumed orders 1, 2, 3 and 4
for the existing fine16/fine20 window, using both seed-spacing and element-
density proxies. The full envelope is 0.72921% to 0.95990%. These are declared
empirical reference models, not continuum error bounds. The nine-cell result
is a cheaper admitted Meshers pilot, not a matched Microgen victory.

## Fixed-realization Spinodoid cap feasibility

Job 126 screens actual Spinodoid VTK seeds at sampling resolutions 20, 24,
28, 32, 36 and 40 cells. Crucially, constructor FFT resolution remains 16.
Changing that constructor resolution would change the random specimen even
with the same seed. Only the sampling-grid attribute changes after field
construction. All 52 coefficients match the frozen realization exactly and
4,096 field evaluations remain bitwise unchanged for every stage.

All six seeds have cap-quality obstructions, with 68 to 122 required cap
triangles below the analytic 0.1 quality ceiling. All remain within 250,000
tets, with the largest at 243,508. These are necessary-feasibility screens,
not geometry audits, MMG timings or stiffness results. Do not run an interior
sizing sweep on these fixed cap triangulations. Select one reproduction of
the default 16-cell arrays and five odd sampling grids, still with the exact
same FFT16 field, before choosing any native successor.

Job 127 reproduces the default 16-cell point, Tet4 and surface arrays exactly,
including its 46 cap obstructions and minimum quality ceiling 0.009166.
The five odd sampling grids also have cap obstructions. Select one final
bounded screen to complete every integer sampling resolution from 16 through
41, including the first expected over-cap seed. Keep the exact realization
and all caps. End this seed search after that result rather than treat a
finite sweep as proof about arbitrary resolutions or modified cap meshes.

Job 128 completes that bounded range. Every integer sampling grid from 16
through 40 has cap obstructions; 41 produces 261,650 tets and exceeds the
250,000 cap before quality analysis. The original field is unchanged at every
stage. End the seed search. This establishes no feasible fixed-cap route in
that tested range, not impossibility for every extraction or changed cap mesh.

The same fine-window reference sensitivity check is also applied to the
representative repeatedly timed Schwarz meshes. Meshers spans 0.25842% to
0.32050%, and the explicitly one-pass Microgen comparator spans 0.92778% to
0.99022%, across orders 1 through 4 and both length proxies. Both remain
below 1% under those empirical models. This strengthens the narrow supported
Schwarz result without changing its route or claiming continuum certification.

## Lattice cap feasibility

Select one bounded cap-only screen of actual smooth-union lattice VTK seeds.
Use the unchanged original ordered eight-strut field, verify shape and case
values at 4,096 points, and reproduce the default 50-point mesh arrays exactly.
Test grids 50, 33, 37, 41, 45, 49, 51 and 52 under the existing tet cap before
choosing any MMG successor. This does not use the Meshers symmetric-gradient
expression or void clamp as a Microgen seed.

Job 129 reproduces the default50 arrays exactly, including 96 cap
obstructions. Grids 33, 41, 45 and 51 have no cap-quality obstruction and
contain 75,774, 126,648, 175,888 and 249,432 tets respectively. Grids 37 and
49 remain obstructed, and 52 exceeds the cap at 256,848 tets. Field values
remain unchanged at every stage. These four feasible grids warrant independent
geometry audits before selecting one MMG successor; cap feasibility alone is
not admission or evidence of stiffness accuracy.

Job 130 measures all four feasible seeds. Grid51 alone passes seed distance,
with H=0.001778 and no failed projections; its minimum cap-quality ceiling is
0.17813. Grids33 and 41 exceed distance, while 45 also has 48 failed
projections. Select actual two-pass Microgen/MMG on grid51 with the original
field, first-pass hausd=0.0015, second-pass hausd=0.0006, hmax=0.06 and
hgrad=1.3. Capture both stages and solve only if geometry, quality and
periodicity all pass. The runner now preserves any admitted lattice tensor
as `solved-needs-reference`, rather than pretending the old summary contains
a justified lattice reference or assigning a stiffness error silently.

Job 131's first pass returns 115,619 tets with H=0.001551 but minimum quality
0.09904, still below 0.1. The second pass produces 282,630 tets and exceeds
the cap. No solve runs. Do not round the first-pass value into admission.
Select one second-pass hausd change from 0.0006 to 0.00075, keeping the
first pass, hmax, hgrad and every independent gate unchanged.

Job 132 completes the actual two-pass route and all six prescribed-strain
loads. It is admitted with 46,219 nodes and 218,394 Tet4 elements,
H=0.001570, volume error=-0.38358%, q_min=0.11491, positive sampled
Jacobians and matching periodic caps. This is a genuine Microgen VTK seed
with explicit first-pass sizing overrides, not the default pipeline or a
Meshers-seed hybrid. The tensor differs by 3.89090% from the current finer
adaptive Meshers tensor. That is pairwise disagreement, not a certified error.

The instrumented solve takes 95.12 seconds and total is 141.87 seconds.
Intermediate stage audits are included, so do not compare this total as a
matched timing against uninstrumented Meshers runs. Select independent field
audits of the admitted comparator and both reference families before assessing
model-dependent accuracy or choosing a fresh timing experiment. The lattice
1% reference question remains separate.

Job 133 passes all six independent field audits. Maximum reaction disagreement
is 4.02e-11 and periodic jumps are zero; the admitted Microgen mesh is within
1.06e-11. Across both existing reference families, Microgen's disagreement
spans 3.35444% to 5.20673%, while the historical Meshers linear32 candidate
spans 3.37702% to 5.23059%. Their actual-fine differences are 3.89090% and
3.91273%. These are closely matched accuracy levels, but neither supports
the original 1% target. The curved24 Meshers candidate spans 0.86871% to
2.66674%, at a substantially higher solve cost.

Select a fresh Meshers linear32 control with the current verified native
build, the same prescribed-strain solver and the original explicit void
clamp. This checks admission, accuracy and live cost before choosing a
proper uninstrumented comparison. Do not compare instrumented Microgen
time directly with the historical Meshers timing or claim the 1% goal.

Jobs 134 and 135 complete that control and its independent field audit.
Current-native linear32 has 48,259 nodes and 219,507 Tet4 elements,
H=0.001802, q_min=0.16532 and all admission gates pass. Six prescribed-strain
loads take 77.64 seconds; the live workflow takes 94.43 seconds. Reaction
disagreement is 1.36e-11, Hill-Mandel disagreement is 4.98e-12 and periodic
jumps are zero. The tensor differs by 3.95566% from the fine adaptive tensor
and by 0.30894% from the admitted tuned Microgen tensor. The latter is a
pairwise difference, not a reference error. Microgen's stage-instrumented
141.87 seconds is still unsuitable for a fair timing ratio.

Select a saved-field recovery-localization diagnostic before another solve.
It measures where the existing periodic recovery indicator concentrates and
whether its proposed size field has useful variation. It does not establish
observed stiffness error, local stress accuracy or a live adaptive speedup.

Job 136 completes in 8.19 seconds at 0.772 GiB peak RSS. The highest-density
20% of elements cover only 0.758% of volume and contain 15.70% of recovery
mass. The median proposed target length is 0.065, with the lower and upper
deciles clipped to 0.03 and 0.1. This does not show that a recovery-driven
solve would repay its coarse-solve cost. Defer that solve and measure the
admitted Microgen recipe without intermediate stage audits first. The
localization diagnostic reuses saved fields and is not a workflow timing.

Job 137 repeats the admitted grid51 Microgen recipe without intermediate
stage audits. Its mesh counts, quality and complete stiffness tensor reproduce
job 132 exactly. Generation takes 24.60 seconds, the final audit 4.87 seconds,
six-load solve 101.25 seconds and full live workflow 131.63 seconds. The fresh
Meshers linear32 total is 94.43 seconds. These are single runs at roughly 4%
fine-tensor disagreement, not a repeated matched-1% result.

Job 138 tests a possible assembly-memory reduction before extending the
reference study. Four-point reintegration of the saved admitted curved24
fields shifts worst mixed-load energy by 0.0102316%, exceeding the declared
0.0100000% threshold. The affine Tet10 patch passes and the relative volume
shift is only 1.83e-10, but the mixed-load gate still rejects the shortcut.
No curved re-solve runs. Keep the fifteen-point default; do not round this
result into a pass.

Select one finer adaptive lattice level at tolerance 0.0028, with the same
24-cell seed, explicit symmetric-union expression and void clamp as the
existing adaptive reference family. Keep the 200,000 solve-node stage,
250,000-tet cap, 8 GiB ceiling and all admission gates. This tests another
refinement increment, not an assumed resolution of cross-family uncertainty.

Jobs 139 and 140 complete this level and seven independent field audits.
The finer mesh is admitted with 193,242 nodes and 124,141 Tet10 elements.
The six-load workflow takes 1,178.53 seconds with 6.12 GiB peak worker RSS.
Its tensor differs by 0.03654% from the tolerance-0.003 level. Reaction
disagreement is 3.46e-11 and periodic jumps are zero. All seven audits pass.

The small successive change does not settle the reference. The finer mesh
differs by 0.24757% to 1.13038% across the newly declared fine-window models,
and by 0.11349% to 1.26988% across the retained earlier families. No candidate
passes every declared model at 1%. The fresh linear candidate reaches
5.25237% and the admitted Microgen candidate reaches 5.20673% across the
combined families. Do not discard the first-order models or the older family
to obtain a pass. Results are in `runs/g3-lattice-finer-reference-comparison-v1.json`.

The next diagnostic addresses the requested volume-correction priority.
The admitted linear lattice underfills the sampled volume target by about
0.43%. At fixed connectivity, test half and full correction using averaged
periodic wall-normal displacements, with box-normal components locked.
Use the original sampled target and record its scramble spread rather than
claim exact volume knowledge. Independently audit distance, quality, volume,
Jacobians and periodic caps before choosing any stiffness solve. These are
saved-mesh geometry probes, not changes to native meshing or live timings.

Job 141 returns one admitted offset. Half correction moves wall nodes by at
most 0.0001558, halves the sampled volume bias from -0.43146% to -0.21573%,
and retains H=0.001947 and q_min=0.15617. Full correction meets the sampled
volume target but has H=0.002096, so it is rejected without a stiffness solve.
The target's full scramble spread remains 0.0001564 in absolute volume.

Job 142 solves only the admitted half correction. Independent reaction and
energy checks pass, with reaction disagreement 1.25e-11 and zero periodic
jump. Against the new fine adaptive tensor, worst mixed-load disagreement
worsens from the original 3.99223% to 4.49776%. Better volume alone is not a
better stiffness mesh. Do not promote uniform volume offsets into native
meshing from this experiment. Its 76.74-second worker time uses a saved mesh
and is not a live generation-plus-correction workflow benchmark.

Relocation is not the current workflow bottleneck. Native diagnostics charge
0.816 seconds to the relocation phase of the finer 1,178.53-second workflow,
and 1.399 seconds in the 94.43-second linear control. These phases also include
settling. Even removing that entire phase would save less than 0.1% and 1.5%
respectively in these runs. Defer a parallel rewrite until a measured workload
justifies it. The unresolved stiffness/reference cost has much larger impact.

Job 143 refreshes the older adaptive20 rejection with the verified native
build. The quality fix now returns 134,109 nodes and 85,306 Tet10 elements,
q_min=0.11744, but H=0.002284 still rejects it. There are 84 partially curved
edges. It also uses more nodes than admitted uniform24, so it is not an
established cheaper route. No solve runs.

Job 144 reproduces the geometry metrics and locates the worst 32 probes.
None is on a required cap edge or vertex. The worst faces include corners
and mid-edge nodes still off the implicit surface, away from box planes.
This does not justify a box-curve repair as the explanation of this failure.

Job 145 tests enrolling inherited off-surface wall vertices in the existing
pending-projection/star-splitting routine. The rebuild, 25 core tests and
36 Python tests pass. However, the target's distance, quality and minimum
sampled Jacobian are unchanged; its volume differs only by floating-point
roundoff. Reject this change for the present problem and restore the verified
`5914097` native build. The proposal and test remain in the frozen
`native-proposal.patch`. The coarse candidate used adaptive tolerance 0.004,
so its failure of the stricter independent 0.002 distance gate does not itself
prove a native tolerance-contract violation.

Job 146 tightens the coarse20 adaptive budget to 0.002 without changing the
independent gates. Geometry now passes at H=0.00022394, but the mesh has
288,800 nodes and 186,619 Tet10 elements. No solve runs above the existing
200,000-node stage. The coarser seed is not a cheaper qualified solution.

Before considering any larger solve stage, job 147 tests four-point sensitivity
on the already-solved tolerance-0.003 reference, within the existing 200,000-node
and 8 GiB stage. Keep the same 0.0001 relative tensor threshold that rejected
the coarse24 shortcut. A coarse-mesh rejection does not establish sensitivity
on a finer mesh, but no reduction is accepted without the actual comparison.

Job 147 passes the fine same-mesh comparison. Four-point stiffness differs
from the solved fifteen-point control by 0.00203264%, below the unchanged
0.01% sensitivity limit. Own-rule stress/energy/reaction disagreement is
4.12e-11 and periodic jumps are zero. Peak worker RSS is 4.358 GiB, versus
the recorded 5.64 GiB fifteen-point control. The solve takes 1,001.12 seconds
versus the earlier 977.30 seconds, so this is not a demonstrated speedup.
It supports a memory-limited reference-study option, not replacing the
default integration rule or assuming coarse meshes have the same sensitivity.

Job 148 extends only the reference study to one uniform34 mesh. The measured
fine-control memory supports a 250,000-node stage, projected near 6.2 GiB at
the cap. The 8 GiB worker ceiling, 250,000-tet cap and existing 1,800-second
timeout remain unchanged. This mesh passes geometry with 233,964 nodes and
150,384 Tet10 elements. Explicit four-point equilibrium must pass independent
own-rule field checks and fifteen-point fixed-field reintegration at the
same 0.0001 relative sensitivity limit. Fifteen-point reintegration is not
a fifteen-point equilibrium re-solve. No reference-model exclusion or 1%
accuracy claim follows from queuing this level.

Jobs 148 and 149 complete this level and eight independent field audits.
The audited reference workflow takes 1,481.66 seconds with 5.92 GiB peak
RSS. Own-rule reaction disagreement is 1.63e-11, Hill-Mandel disagreement is
8.55e-12 and periodic jumps are zero. Fifteen-point fixed-field energy differs
by 0.00521056%, below 0.01%. The new tensor is 0.101736% from uniform32 but
still 0.973372% from the finer adaptive tensor. It spans 1.21347% to 2.08205%
across the retained fine-adaptive models and 0.43492% to 2.22983% across the
other retained families. No candidate passes every retained model at 1%.
These exploratory ranges do not propagate quadrature-model uncertainty and
must not be treated as continuum error bounds or a new stopping certificate.

Comparison helpers now reject failed solver rows and four-point rows missing
the declared sensitivity pass, or using a weaker sensitivity limit. Five
admission regression cases pass against actual saved controls. The default
solver quadrature remains fifteen points.

The next comparison question is straight Tet10 on a genuine Microgen/MMG
mesh through direct fedoo conversion, not Microgen's unsupported public Tet10
converter. Exact edge counts show that the admitted two-pass candidate would
need 326,716 Tet10 nodes, outside the current 250,000-node stage.
Job 150 coarsens only the second-pass hausd from 0.00075 to 0.0009. It has
H=0.001640 and a 267,550-node lift count, but q_min=0.097405 rejects it. No
Tet10 lift or solve runs, and the quality value is not rounded into admission.

Select one explicitly modified actual MMG variant on the smaller first-pass
mesh: replace the second level-set call with `optim`, `nosurf`, `noinsert`.
This is not the default public Microgen pipeline. Preserve the original VTK
seed, first-pass parameters and all independent admission gates. Read its
geometry and prospective lift count before choosing any direct-fedoo solve.

Job 151 rejects that variant. Interior optimization returns 24,923 nodes
with H=0.001551, but worsens minimum quality to 0.07025. No Tet10 solve runs.
This is another unsupported admitted route under the fixed gates, not a
Microgen default failure inferred from a Meshers seed.

The admitted Schwarz Microgen parent provides a stronger supported integration
control without another mesh-parameter search. Its exact straight Tet10 lift
has 230,337 nodes and 154,012 elements, inside the current 250,000-node stage.
Job 152 independently audits its unchanged affine geometry, converts the
actual Microgen-derived arrays through `Mesh.to_fedoo()`, and uses the converted
arrays in six prescribed-strain solves. No native Meshers generation is called.
Four-point integration is exact for affine Tet10 linear elasticity; also check
fifteen-point fixed-field energy and independent reactions/periodic jumps.
This is a saved-mesh diagnostic. It excludes the live Microgen extraction/MMG
cost and must not be ranked as an end-to-end timing or a newly supported public
Microgen Tet10 converter.

Job 152 completes that stronger Schwarz control. Its declared full mixed-load
reference-model range is 0.054006% to 0.062702%. All six loads and independent
field checks pass, including fifteen-point reintegration of the affine
Tet10 fields. The saved-mesh solve takes 1,128.92 seconds and peaks at 6.05 GiB.
This strengthens the accuracy control, not a new end-to-end timing claim.

At a supplementary 0.1% target, current Meshers curved14 and curved15 reach
0.101239% and 0.103035%, so neither passes. Independent field audits confirm
these results. The known curved16 level passes after a fresh current-native
run and independent audit: 42,928 nodes, 0.031469% to 0.093778% model
disagreement and an 88.01-second single-run live workflow. Its tensor reproduces
the archived sixteen-cell reference input within 2.67e-11. This checks that
input's numerical reproduction, not independent continuum convergence.

Before comparing strict-target workflow costs, job 159 doubles only the actual
Microgen first-pass hmax from 0.05 to 0.1. The original VTK30 seed, hausd and
hgrad, explicit one-pass route and admission gates remain unchanged. The result
is admitted and would lift to 207,200 Tet10 nodes rather than 230,337. Job 160
measures its six-load accuracy as a saved-mesh control. Do not assume that the
denser control was the cheapest eligible Microgen route. `TET10_RESULTS.md`
contains the exact protocol and limits. The original four-fixture 1% objective
remains unmet.
