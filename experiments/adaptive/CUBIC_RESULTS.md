# Reference-only cubic enrichment on G3

## Measured edge refinement improves lattice error but still fails, 25 October 2026

Jobs 369–374 replay all 24 reference formulas, verify exact parent provenance,
and evaluate one residual-ranked periodic edge-star Tet10 refinement on the
unchanged current native domain. At 249,998 nodes / 163,802 tets, geometry,
quality, periodicity and the 8 GiB resource limit pass. Six fresh loads and
independent field/quadrature/Galerkin checks pass. Worst all-model disagreement
improves from 1.283822914% to 1.216217882%, still failing 1%. No production
default, native candidate, reference or timing victory is promoted. Do not
repeat this h recipe automatically. See [edge residual investigation](EDGE_RESIDUAL_RESULTS.md).
All four cases and all 24 models remain in scope. Earlier next-step descriptions
below are historical; a separately reproduced physical-side gradient control
must pass before any native evaluation change is retained.


## Fixed native and actual wall geometry inspected, 25 October 2026

Three existing surface-projection controls pass in job 366, including exact
planar area/distance and explicit rejection reporting for feet outside the box.
Read-only jobs 367 and 368 compare the saved native lattice26 and admitted
actual lattice meshes against the same ordered implicit wall. The frozen
specimen JSON and sample NPZ hashes match between routes. No mesh coordinate,
connectivity, sample or material field is regenerated or changed.

Positive 9/16-point sampling initially changes native RMS by about 10%, so
job 368 adds 36/64 points per wall face. Native full-wall RMS settles at
1.579072e-5, changing 0.00222% from 36 to 64 points. All 987,776 native wall
probes project within the box. The curve-incident wall geometry is especially
close to the analytic wall, although this does not certify local stress.

Actual wall faces without curve vertices have RMS 4.306095e-4, stable across
the denser rules, with zero failed projections in that group. The two routes'
excluded curve strips differ, so these group values are not a comparison on
identical wall patches. Actual curve-incident probes include 4,659 projections
outside the box at 64 points per face. Their raw RMS and signed integrals are
saved with explicit failures and must not be treated as fully valid bounded
wall-projection metrics. This diagnostic does not replace the original
admission method or silently reclassify the admitted actual competitor.

The native mesh already has better original sampled distance, smaller absolute
volume error, higher minimum corner quality and fewer nodes/tets than this
actual route. It still fails the retained-model stiffness target. These new
geometric observations do not prove the cause of that disagreement, justify
fitting a geometry offset, or certify a converged reference. Worker 368 takes
13.35 seconds at 0.302 GiB, solely a saved-domain diagnostic, not new workflow
timing. See [fixed-domain diagnostic](../../handoff/macbook/lattice-fixed-domain-comparison-20261025.json).

Next audit the original lattice reference extrapolation assumptions and their
fixed-domain evidence before choosing another mesh allocation change. Retain
all 24 model tensors and all four cases. The residual screen remains below
its threshold for another solve, and no native/reference victory is declared.
Earlier selected-next descriptions below are historical.


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


## Completed curve-priority nested diagnostic on Mac

Mac job 11 passes sixteen controls and exact-domain preflight at 249,999
nodes and 158,967 tets. It accepts 544 of 549 preferred wall/cap curve edge
orbits. Minimum quality is 0.10271976. Domain and old-field interpolation
errors are below 3.34e-16, and periodic caps, coherent facets, distance,
volume, and sampled positive Jacobians pass. Job 12 confirms actual reduced
assembly at 6.230 GiB peak RSS under the unchanged 8 GiB watchdog.

Job 13 completes six fresh loads with verified Accelerate single-thread
mode. Worker time is 216.33 seconds at 5.472 GiB. Independent fields and
positive quadrature pass. The largest parent-to-refined FE correction is
-0.0132094%, with Galerkin balance error 4.82986e-10. Worst retained-model
disagreement remains 1.2556246%, failing 1%. All 24 models remain unchanged.

This gives little extra gain over longest-edge refinement at 1.2565145%.
Stop these one-stage bisection recipes without a native feature, replacement
reference, or workflow claim. Independent checkpoint verification rechecks
all transferred evidence, exact measured source versions, and every full
model score for the three Mac nested studies. The objective remains active.
Next return to missing admitted actual-comparator coverage, using preserved
failure evidence before selecting a new bounded experiment.

## Lattice26 fixed-domain negative control and Mac thread correction

Mac job 6 adds 50,028 nodes to the best saved native lattice26 mesh, reaching
exactly 250,000 nodes and 160,291 tets. It preserves the original curved
physical domain and six old fields within 3.33e-16. Every admission gate
passes. Job 7 confirms actual reduced assembly at 5.270 GiB RSS.

Job 8 completes six fresh loads and independent field, positive-rule, and
Galerkin checks. The largest generalized FE correction is -0.0122671%, and
the full balance agrees within 6.07082e-10. Worst disagreement across all
24 retained models is 1.2565145%, still above 1%. Stop this bulk longest-edge
recipe. This negative test does not refute targeted refinement, promote a
reference, or establish a live workflow timing.

The installed Mac NumPy BLAS is Apple Accelerate. Initial OpenBLAS/OMP
settings did not explicitly control its native threading mode. Numerical
timings from Mac jobs 1–8 are excluded from verified one-thread comparisons;
their tensor and field evidence remains frozen. The shared runner now uses
`VECLIB_MAXIMUM_THREADS=1` and verifies `BLASSetThreading(1)` with
`BLASGetThreading()`. The SDK defines mode 1 as single threaded on the
calling thread. Jobs 9–10 pass 55 focused tests and twelve fresh AMG loads.
The separate clarification preserves original manifests and source hashes.

Job 11 then tests one geometry-only ordering change on original lattice26:
prioritize actual wall/cap intersection edge orbits before descending length.
The earlier full correction is concentrated near curve vertices, while the
bulk test resolves little within its budget. No field-derived sizing strength,
geometry projection, specimen change, or cap relaxation is introduced.
All parent-map, field-nesting, quality, and periodic guards remain. Read the
preflight before selecting any solve. No native API change is promoted.

## Second nested correction on native Mac, 24 October 2026

Mac jobs 0–2 establish the native environment and queue, pass 54 focused
tests, and pass twelve fresh cube/Schwarz solver-control loads. Exact G3
package versions have arm64 wheels. macOS rejects Linux's `RLIMIT_AS` setup;
the shared bounded runner retains the 8 GiB process-tree RSS watchdog and
records its enforcement explicitly. Small timeout and memory-stop controls
pass. Neither mesh cap changes.

Job 3 refines the completed G3 level1 mesh on its exact quadratic physical
domain. Fourteen controls pass. The new mesh has 218,135 nodes and 141,144
tets, from 64,197 parent bisections. Geometry and six-field interpolation
errors are 2.78e-16 and 3.33e-16. Relative volume change is -2.52e-14.
Distance 0.0006492808, quality 0.10010709, positive sampled Jacobians,
periodic caps, and coherent oriented facets pass. Preflight worker time is
87.72 seconds at 1.047 GiB peak RSS.

Job 4 captures only actual reduced-system metrics before any solve. The
646,029-square matrix has 51,030,214 stored entries and 614,946,688 CSR bytes.
It fits at 4.464 GiB peak RSS in 9.61 seconds. No large matrix archive or
fake solution is produced.

After reading that resource check, job 5 completes all six fresh loads.
Solve worker time is 191.45 seconds at 4.508 GiB. Independent fields and
positive quadrature pass. All six generalized stiffness corrections from
level1 are negative, with maximum magnitude 0.3282719%, versus 0.4995130%
for level1. The full prolonged-parent/fine/correction energy balance agrees
within 3.78355e-8. Worst retained-model disagreement is 2.1351167%, failing
1%. All 24 models remain unchanged. No reference is promoted.

This isolates another FE correction on one fixed domain. It does not
establish an asymptotic rate or supply a live workflow timing. Records are
`runs/macbook-lattice-nested23-level2-preflight-v1`,
`runs/macbook-lattice-nested23-level2-assembly-v1`, and
`runs/macbook-lattice-nested23-level2-v1`.

The controlled FE softening justifies one different candidate investigation:
the best native lattice26 tensor has positive worst-model error 1.2686985%.
Job 6 tests its saved mesh with the unchanged nested bisection rule, first
checking geometry, periodicity, and actual counts within 250,000 nodes/tets.
Read before a solve. This tests FE resolution without changing geometry or
strengthening a sizing recipe. It is not yet an improved native candidate,
matched-target workflow, or substitute reference.

## Nested fixed-domain correction, 23 October 2026

Job 473 completes lattice26 saved-field and positive-quadrature checks.
Reaction error is 1.0611e-11, Hill-Mandel error 4.1723e-12, periodic jumps
zero, and 125-point fixed-field energy shift 4.0667e-8. Its retained-model
score still fails 1% at 1.2686985%.

Job 475 tests complete periodic edge-star bisection with disjoint incident
stars, exact parent Tet10 geometry restriction and field prolongation.
Ten controls pass. Parent23 grows from 68,095 nodes and 43,844 tets to
123,514 nodes and 76,947 tets. Interpolation errors are 3.33e-16 and volume
change 1.40e-14. Independent distance 0.0006492808, quality 0.10385108,
positive Jacobians, periodic caps and oriented facets pass.

Job 476 completes six fresh loads in this nested space. Solve time is
245.83 seconds, worker 247.63 seconds at 3,304,448,000 bytes peak RSS.
Independent fields and positive rules pass. Every generalized correction
from the parent is negative; largest magnitude is 0.499513%. The full
prolonged-parent/fine/correction energy balance agrees within 6.6756e-8.
Worst disagreement across all 24 models is 2.4480846%, above 1%. This is
a controlled correction on a coarser fixed domain, not a new reference or
live meshing workflow timing. One level does not establish convergence.

Records are `runs/g3-lattice-adapt26-field-audit-v1`,
`runs/g3-lattice-nested23-preflight-v1`, and
`runs/g3-lattice-nested23-level1-v1`. Kevin requests MacBook primary
continuation on 24 October. Read `handoff/macbook/HANDOFF.md`; next prepare
one second-level nested geometry/resource preflight before choosing a solve.

## Final lattice26 workflow before pause, 20 October 2026

Job 379 completes the fresh six-load workflow selected after preflight 378.
It retains all 24 models and the fixed specimen. The geometry-admitted curved
Tet10 mesh has 199,972 nodes and 128,300 elements. Cached elastic-block AMG
uses 108, 108, 108, 109, 108, and 110 iterations, with maximum reported
relative residual 9.41e-13. Generation is 9.03 seconds, solve 716.94 seconds,
and charged total workflow 736.51 seconds. Worker wall is 738.82 seconds at
5,781,733,376 bytes peak RSS, below 8 GiB.

A lightweight read-only comparison of its symmetric tensor to all 24 models
gives worst generalized energy disagreement 1.2686985%. Three models exceed
1%; the worst is index 0 of `g3-lattice-reference-sensitivity-v2.json`.
This is a preliminary rejection, not an independently audited stiffness
result. Saved-field, reaction, Hill-Mandel, and higher-positive-rule checks
remain pending. The run status is `solved-needs-reference`, and `qualified`
refers to geometry admission only. No controlled lattice24 refinement or
timing ratio follows because native binaries and solver configurations differ.

The user requested a pause. No successor or audit was queued. The completed
record is `runs/g3-lattice-adapt26-tol0028-amg-solve-v1`. Mesh and all six
displacement archives remain local and ignored. Resume instructions and the
project-wide unresolved objective are recorded in `PROJECT_TRACKING.md`.

## Retained lattice target is algebraically feasible, 20 October 2026

Job 377 checks the unchanged 24 empirical stiffness models before another
native mesh or reference solve. It does not promote a new reference.

If a symmetric candidate is within relative energy error `e` of every model,
then for any strain direction its energy must be at least `(1-e)*A` and at
most `(1+e)*B` for every pair of model directional energies. Therefore an
ordered generalized ratio `r` requires `e >= (r-1)/(r+1)`. All 552 ordered
pairs of the full 6-by-6 matrices are checked. The worst direction and both
energies are saved and independently recomputed in the assessment.

The worst ratio is 1.0182028025, below the necessary one-percent limit
1.0202020202. Its lower bound is 0.9019312871%. The pair consists of models
7 and 0 in `g3-lattice-reference-sensitivity-v2.json`, using element-density
order 4 and uniform-seed order 1 extrapolations. The bound is necessary, not
alone sufficient for a common tensor across 24 matrices.

The pair's harmonic trial matrix is then checked against every retained
model using the full generalized-eigenvalue metric. Its worst disagreement
is 0.9084715371%, below 1%. This establishes algebraic feasibility, so an
inconsistent target is not an explanation for the existing near misses.
The trial is not a physical mesh, equilibrium solve, continuum estimate,
replacement reference, or permission to fit geometry to its coefficients.
No model or admission threshold changes.

Six queued controls pass, including exact scaled-SPD bounds, mixed-direction
and basis-congruence checks, identical models, invalid stiffness rejection,
and standard progress controls. Inputs are hashed before and after the run.
The record is `runs/g3-lattice-reference-envelope-feasibility-v1`, with
`reference_envelope.py` and `reference_envelope_job.py` as maintained code.

The selected next candidate is current-native lattice26 at the existing
adaptive tolerance 0.0028. Geometry-only preflight 378 passes at 199,972 nodes
and 128,300 tets, below the fixed 250,000-node solve stage. Sampled distance
is 0.0003236401 and minimum quality 0.149329. A fresh regenerated six-load
AMG workflow is selected before any accuracy statement. The archived lattice24
binary differs, so node-count scaling is only a budget estimate, not a
controlled convergence or timing claim.

## Cubic geometry sensitivity, 20 October 2026

Job 309 separates a change in boundary mapping from another equilibrium
enrichment. The existing cubic displacement space embeds a quadratic physical
domain. The probe retains its 222,331 nodes, 44,976 elements, connectivity,
corner coordinates, and all six saved displacement coefficient arrays. It
projects only 30,972 higher-order wall nodes onto the fixed implicit boundary,
locking box-plane coordinates. It introduces no native Meshers Tet20 API.

Three focused controls pass: cubic basis identities, independent cubic surface
interpolation against the registered reference element, and rejection of an
implicit linear-geometry fallback for nonstandard elements. The experimental
geometry auditor now accepts an explicitly registered cubic element, cubic
boundary interpolation, and positive-rule volume quadrature. Default Tet4 and
Tet10 paths remain unchanged. Baseline distance and volume reproduce the saved
quadratic-domain report within 1e-10.

Both old and projected domains pass the unchanged sampled admission gates:

| Measurement | Embedded quadratic domain | Projected cubic domain |
|---|---:|---:|
| Two-sided sampled distance | 0.00054354 | 0.00013101 |
| Volume | 0.24823394598 | 0.24823372775 |
| Minimum corner quality | 0.34102 | 0.34102 |
| Minimum sampled Jacobian | 3.6811e-6 | 2.2608e-6 |
| Failed forward projections | 0 | 0 |

Periodic cubic nodes and full ten-node cap faces match. Maximum coordinate
mismatch is below 1e-15. The largest nodal movement is 0.00039995; median
movement is 3.24e-7. This is not an equilibrium or accuracy result.

The saved fields are then integrated without changing their nodal coefficients.
Positive 64-point and 125-point rules pass affine elasticity controls. The
125-point mixed-load energy change between domains is 0.0105922%, below the
predeclared 0.05% geometry-signal screen. All six signed changes are negative.
The largest 64-to-125-point change is 9.965e-7 relative. Baseline energy
reproduces the earlier saved-equilibrium quadrature difference, 6.869e-6.
Projected-field stress and energy disagree by 1.11e-5 relative; those fields
are not newly equilibrated.

Geometry takes 18.64 seconds at 0.475 GiB, and reintegration takes 105.57
seconds at 0.531 GiB. Do not run a costly projected-domain equilibrium solve
from this sub-threshold signal alone. This does not bound re-equilibrated
geometry error or continuum error. The original cubic result and all 24
empirical models remain unchanged. Before another reference solve, inspect
one finer native parent and its exact cubic node budget for an independent
higher-order refinement level.

Records are in `runs/g3-lattice-cubic-geometry-sensitivity-v1`; the maintained
runner is `cubic_geometry_probe.py`.

## Finer reference resource preflight, 20 October 2026

Job 310 generates one nonadaptive curved25 parent with the same symmetric
lattice representation and positive void clamp used by the original cubic24
parent. It passes admission at 75,113 nodes and 48,696 elements: sampled
distance 0.00038198, minimum corner quality 0.32248, positive sampled
Jacobians, no unresolved projections, and matching periodic caps.

Exact topological counting predicts 241,555 cubic nodes, inside the existing
250,000-node stage. The count is corner vertices plus twice the unique edge
count plus the unique triangular face count. It reproduces the old cubic24
count of 222,331 and passes single-tet and shared-face checks. The assembly
runner now checks and records this budget before constructing the enrichment.

Job 311 builds the actual cubic25 mesh and matrix without solving. Its
241,555 nodes match the prediction. Matrix dimensions are 724,665 squared,
with 99,084,987 nonzeros and 1,191,918,508 CSR bytes. The worker takes 58.35
seconds and peaks at 5.818 GiB, below the same 8 GiB ceiling. This is resource
feasibility, not convergence or a workflow performance claim.

One fresh six-load parent25 control is selected next. Its fields are needed
to verify exact quadratic-to-cubic embedding and the later enrichment energy
identity. Read and independently audit that parent before selecting a cubic25
equilibrium solve. No reference model is removed or replaced.

Records are in `runs/g3-lattice-cubic-parent25-preflight-v1`,
`runs/g3-lattice-cubic25-assembly-memory-v1`, and
`runs/g3-lattice-cubic-parent25-fields-v1`.

Job 312 completes all six fresh parent loads in a 128.64-second workflow.
Its arrays reproduce the preflight exactly and reduced residuals are below
9.35e-13. Independent audit 313 reproduces reactions within 7.68e-12 relative
with zero periodic jumps. The source is ready for an enrichment check.

The first fresh cubic25 solve, job 314, fails before completing any load.
Exact geometry embedding agrees within 1.67e-16 and embedded parent energy
within 9.08e-12 relative, but the full `A - A.T` symmetry check requests a
1.47 GiB temporary and hits the unchanged 8 GiB address-space ceiling. Peak
RSS is 6.068 GiB. Assembly fitting in memory did not prove solve feasibility.
The failed record remains in `runs/g3-lattice-cubic25-reference-v1`.

The research AMG helper now computes the same complete relative Frobenius
symmetry norm in blocks of at most 32,768 rows. The threshold stays 1e-12,
every row is checked, and complete matrix hashes still guard reuse. Ten
focused tests compare against the full norm, detect off-diagonal asymmetry,
and prohibit whole-matrix subtraction on the forced-block path. Job 315
passes twelve integrated cube and Schwarz loads with the same equations and
solver tolerances. Cube disagreement is 3.45e-12; Schwarz reproduces its saved
tensor within 4.42e-11 worst mixed-load disagreement. Neither is a large-case
memory proof.

Only after reading those controls, one fresh cubic25 retry is selected in
`runs/g3-lattice-cubic25-reference-blocked-symmetry-v2`. No prefix is replayed.
The 8 GiB ceiling, 250,000-node stage, 240-second first-load deadline, and
1,800-second worker timeout stay fixed. Do not infer convergence or reference
accuracy until all six loads, own-rule checks, and positive-rule sensitivity
are complete. Control records are in `runs/g3-blocked-symmetry-amg-controls-v1`.

## Completed cubic25 refinement and unresolved models

Job 316 completes all six fresh cubic25 loads, without prefix replay. The
first recorded solver call takes 185.43 seconds, and the full first-load
continuation gate passes within 240 seconds.
Final solve time is 1,020.87 seconds; worker time is 1,071.18 seconds and
peak RSS is 6.753 GiB. The original 8 GiB ceiling is unchanged. This confirms
large-case feasibility of blocked symmetry verification, not a workflow-speed
comparison with the cubic24 run that replayed four verified fields.

Own-rule stress and energy reproduce reactions within 1.03e-11 relative.
Periodic displacement jumps are below 2.23e-16. The 125-point energy shift
is 2.6354e-5 relative, below the unchanged 1e-4 limit. The minimum eigenvalue
of the parent-minus-cubic tensor is +0.000140792, consistent with same-domain
variational enrichment. Status remains `solved-needs-reference`.

Fresh-process audit 317 passes in 62.11 seconds at 0.555 GiB. Its own-rule
reaction disagreement is 9.12e-12, and the 64-to-125-point energy change is
6.8313e-6 relative. All sampled Jacobians and affine controls pass. Every
earlier empirical model is retained.

| Compared tensor | Worst mixed-load difference relative to cubic25 |
|---|---:|
| Cubic24 | 0.248177% |
| Finest adaptive Meshers Tet10 | 0.651535% |
| Eight-round curve-graded Meshers Tet10 | 0.483742% |
| Strongest Microgen-derived straight Tet10 | 0.917416% |

Cubic24 is softer than cubic25 in all six generalized modes. Two neighboring
levels do not establish an asymptotic reference. Candidate ranking also depends
on reference choice: the curve pilot is closer to cubic25 than the finest
adaptive mesh, despite being worse across the original model envelope.
Do not select cubic25 alone to manufacture a target-qualified timing comparison.

The separate model-normalized score uses the original benchmark denominator,
not the reverse model-versus-cubic diagnostic. Cubic25 itself spans 0.114924%
to 1.920819% across all 24 models and fails their 1% envelope. The earlier
candidate model envelopes remain exactly unchanged. No model is deleted or
reclassified as a certified continuum tensor.

Before a third higher-order solve, job 318 preflights one neighboring parent23.
It admits 68,095 nodes and 43,844 Tet10 elements, with distance 0.00059463,
minimum quality 0.15825, positive sampled Jacobians, no unresolved projections,
and matching periodic caps. Exact topology counting gives 218,508 cubic nodes,
inside the existing stage. A fresh parent23 six-load control is selected next.
This permits a third refinement point; it does not guarantee a reliable
three-point extrapolation or settle sensitivity to nonnested meshes.

Records are in `runs/g3-lattice-cubic25-reference-blocked-symmetry-v2`,
`runs/g3-lattice-cubic25-assessment-v1`,
`runs/g3-lattice-cubic25-reference-envelope-v1`, and
`runs/g3-lattice-cubic-parent23-preflight-v1`.

The fresh parent23 control, job 319, completes all six loads in a 119.43-second
workflow and exactly reproduces the preflight arrays. Independent audit 320
has reaction disagreement 6.57e-12 and zero periodic jumps. Assembly-only
job 321 confirms 218,508 cubic nodes, 43,844 elements, and 89,362,584 matrix
nonzeros. Its worker takes 52.91 seconds and peaks at 5.014 GiB. Only after
these gates pass, one fresh cubic23 solve is selected as job 322 in
`runs/g3-lattice-cubic23-reference-v1`. No prefix is replayed and no reference
model is removed. Read its actual completion and fresh-process audit before
interpreting a three-level sequence.

## Declared three-level consistency screen

Before the cubic23 tensor is complete, `cubic_refinement.py` declares a
diagnostic for `C(h) = C_limit + A*h**p`. It uses one common symmetric leading
matrix, not independent fits to whichever scalar entries look favorable.
Differences are whitened by the finest tensor. Unequal mesh spacings enter
the correction-ratio equation explicitly.

The assumed positive-order range is 0.25 to 6. Corrections must have positive
proportionality, best common-matrix relative Frobenius residual at most 10%,
and mixed-load signal at least ten times the largest observed fixed-field
quadrature shift. An inferred limit must be positive definite. Seed spacing
and element-density proxies are reported separately, with neither silently
chosen after the result. Exact audited model tensors must match the retained
24-model set, not merely have the same count.

Synthetic controls test known integer and fractional orders, correction
reversals, nonproportional matrices, and quadrature-unresolved corrections
before the dataset is analyzed. Compatibility only supports that stated local
model. Failure does not prove absence of convergence under another model;
success does not certify a continuum limit or replace benchmark admission.
No existing model is deleted and no geometry or stiffness threshold changes.

## Completed third level and rejected leading-order fit

Job 322 completes all six fresh cubic23 loads. Solve time is 972.61 seconds,
worker time is 1,018.27 seconds, and peak RSS is 5.976 GiB. Reduced residuals
are below 9.76e-13. Own-rule error is 1.29e-11 relative and periodic jumps are
below 2.23e-16. The 125-point energy shift is 7.6754e-6, below the unchanged
1e-4 limit; parent-minus-cubic has minimum eigenvalue +0.000137181.
Independent audit 323 passes in 56.75 seconds at 0.521 GiB, including
64-to-125-point energy change 1.1212e-6.

Job 324 first passes all seven synthetic tests, then reads the audited
23/24/25 tensors. Every field hash, parent provenance, specimen, and exact
24-model tensor set is verified. Both declared length proxies reject the
common leading-order model:

- The 23-to-24 correction is nominally positive definite, while the 24-to-25
  correction is negative definite in the common cubic25 energy metric.
- Correction-direction cosine is -0.959707. The best proportionality factor
  is -1.153881, inconsistent with a positive-order common leading term.
- Relative matrix-model residual is 28.1003%, above the declared 10% screen.
- The dominant correction signal is over 94 times the largest measured
  fixed-field quadrature shift. That sensitivity does not explain the largest
  observed step, though it is not a rigorous quadrature-error bound.

No observed order or extrapolated tensor is promoted. Cubic25 differs from
cubic23 by 0.08948% with cubic23 as denominator, but fitting only those two
odd levels would discard the valid middle result. Do not do that. This rejects
the stated local model, not every possible convergence mechanism.

Against cubic23 alone, the curve pilot differs by 0.53964%, the finest adaptive
native tensor by 0.70541%, and the strongest Microgen-derived Tet10 by
0.96909%. Those values do not override the unchanged retained-model failures
or justify a target-qualified timing comparison. Reference uncertainty remains
the blocker to interpreting another native sizing change.

Records are in `runs/g3-lattice-cubic23-reference-v1`,
`runs/g3-lattice-cubic23-assessment-v1`, and
`runs/g3-lattice-cubic-three-level-screen-v1`.

## Quartic edge-witness protocol

Do not select another nonnested grid or a favorable extrapolation from this
oscillation. The next proposed diagnostic uses the saved cubic24 fields on
their unchanged quadratic physical domain. It adds one quartic edge mode
`16*lambda_i*lambda_j*(lambda_i-lambda_j)**2` per periodic wall/box-curve edge
orbit, with a shared three-component displacement coefficient. Each mode is
supported on the complete incident-tet star. Geometry, topology, and the
existing fields remain frozen inputs.

Assemble the small coupled mode stiffness `K` and six-load residual `R`,
including overlap terms between edge stars. Minimizing only these extra
coefficients gives `z = solve(K, R)` and a trial-energy reduction `R.T @ z`.
Do not sum isolated overlapping patch corrections as if they were independent.

Before a lattice measurement, controls must verify the quartic gradient,
constant-strain periodic patch, coupled stiffness symmetry and positivity,
and independently integrated trial-energy identity. Check continuity and
periodic traces, including quartic face samples, and retain the source's
geometry admission. Positive 64-point and 125-point rules must agree on the
energy reduction within 1%. A worst mixed-load reduction of at least 0.05%
relative to the saved cubic tensor is the predeclared signal for a useful
curve-edge witness. Preserve the 8 GiB ceiling and bounded stages.

This is a proposed residual/enrichment diagnostic, not a new Meshers API or
another recovery-strength trial. Its trial fields are not globally equilibrated;
their energy is not a new homogenized reference. A small signal does not rule
out other unresolved FE modes, and a large signal does not certify continuum
or local-stress accuracy. The measured experiments following this protocol
are recorded below.

## Coupled quartic edge-witness results

The prototype is `quartic_edge_probe.py`; it changes no native meshing API.
Three unit tests check the quartic derivatives, support and endpoint symmetry,
and engineering-strain mapping. Jobs 325 and 326 then validate affine and
manufactured periodic cube controls. The strengthened control has eight edge
orbits, 24 scalar displacement DOFs, and 29 support tets. Affine residual is
1.26e-16, trial-energy identity error is 1.00e-15, and positive 125-point versus
216-point mode matrices differ by 1.79e-15 relative. Nonzero face and box-edge
periodic traces agree at all 15 quartic face nodes. Off-orbit coupling is
nonzero and tested explicitly. A nonzero pin-gauge control is exercised;
uniform translations can restore the source pin without changing energy.

Job 327 uses the saved cubic24 domain and fields without modifying them.
It assembles all overlap terms for 102 periodic curve-edge orbits, 306 mode
DOFs, and 332 support tets. Independent trial-energy checks and continuity
pass. However, the 64-to-125-point gain changes by 6.0323%, failing the
predeclared 1% precision gate. That run remains precision-rejected.

One explicitly higher-rule check, job 328, preserves the geometry, fields,
modes, overlap matrix, 1% precision gate, and 0.05% signal screen. Its
125-to-216-point gain change is 0.017934%, so precision passes. The measured
worst mixed-load trial-energy reduction is only 6.4791e-7 relative, or
0.000064791%. This is 0.648 ppm versus the required 500 ppm, about 772 times
below the signal screen. All six reduction eigenvalues are positive.
The local mode solve residual is below 1.05e-16, direct energy-identity error is
2.04e-13, continuity and periodic quartic-trace mismatches are zero, and no
pin translation is needed on the actual lattice source.

The higher-rule lattice worker takes 26.57 seconds and peaks at 0.3404 GiB.
These are diagnostic costs, not live meshing-workflow timings. Existing
geometry and benchmark references remain unchanged. Do not report the
reduction as a stiffness error, a new equilibrium tensor, or a continuum
bound. The resolved small witness does not justify another native sizing
strength trial or a new meshing feature.

The original cubic coefficients were fixed. A small mode-only reduction does
not bound the benefit of re-equilibrating the original FE space with the new
modes, or exclude errors in untested modes. The earlier full Tet10-to-cubic
audit covers 1,408 curve-vertex-incident tets and 75.02% of its correction
trace; the edge-incident subset has 332 tets and 35.55%. Those are a different
enrichment diagnostic, not a claim about the remaining cubic error.

Records are in `runs/g3-quartic-edge-controls-v1`,
`runs/g3-quartic-edge-controls-v2`,
`runs/g3-lattice-quartic-edge-witness-v1`, and
`runs/g3-lattice-quartic-edge-witness-higher-v2`.

## One remaining coarse-coupling check

Before closing this mode family, test only one specific omission: allow the
existing cubic corner basis functions at periodic curve vertices to relax
jointly with the same quartic edge modes. Include their complete vertex stars
and all cross terms, without changing geometry, adding another degree family,
or selecting a favorable grid. Retain the source pin through a uniform gauge
translation where needed. This is a proposed diagnostic, not an implemented
native feature or a selected global reference solve.

Report edge-only, coarse-corner-only, and joint gains separately. The signal
must use the incremental quartic gain after coarse relaxation, not credit
coarse-only integration or re-equilibration effects. Check it independently
against the Schur complement and require nonnegative nested-space gain
increments. The original edge-only result must reproduce before comparison.
Keep the validated 125/216-point pair, 1% precision gate, 0.05% signal screen,
8 GiB ceiling, and existing mode/support caps. If that increment remains
small, stop this polynomial-mode route rather than escalating degree or
tuning its strength. The completed check is recorded below.

This experiment addresses lattice reference uncertainty. It does not add a
Meshers Tet20 product API, change the Microgen pipeline or establish a timing
advantage. All existing uniform and adaptive reference models remain retained.

## Completed corner-coupling check, 20 October 2026

Job 329 passes five algebraic tests and manufactured periodic controls for
the combined space. The controls cover registered cubic corner functions,
complete periodic corner copies, nonzero face traces, affine residuals,
nonzero pin gauges, coupled trial energy, and the Schur-complement identity.
The manufactured affine residual is 1.88e-16 and energy identity error is
2.70e-15. Positive 125-point and 216-point matrices agree within 4.65e-15.

Job 330 measures the unchanged cubic24 domain and saved six-load fields.
The space has 102 quartic edge orbits and 96 cubic corner orbits, totaling
594 displacement DOFs on 1,408 tets. All cross terms are assembled. Original
edge-only gains reproduce exactly at both rules. The incremental gain after
corner relaxation agrees with the independent Schur calculation within
5.37e-16 relative. Nested-space positivity checks pass.

| Positive-rule measurement | Result |
|---|---:|
| Edge-only maximum relative reduction | 0.647912 ppm |
| Corner-only maximum relative reduction | 0.000873 ppm |
| Joint maximum relative reduction | 0.676767 ppm |
| Incremental quartic maximum after corner relaxation | 0.676219 ppm |
| Incremental gain change, 125 to 216 points | 0.018333% |
| Direct trial-energy identity relative error | 3.86e-13 |
| Maximum continuity and periodic-trace mismatch | 3.39e-21 |
| Lattice worker time and peak RSS | 263.45 s, 0.6863 GiB |

Each maximum uses its own worst mixed load. Differences of these scalar
maxima need not equal the maximum of the matrix difference. The signal uses
the full joint-minus-corner matrix, not subtraction of scalar maxima.

Precision passes the unchanged 1% screen. The incremental signal fails the
500 ppm screen by about 739 times. No source-pin translation is needed.
The 300-second worker stage and 8 GiB ceiling are unchanged and satisfied.
Geometry, source fields, and all 24 model tensors remain unchanged.

Stop this polynomial-mode route. Do not escalate degree, tune mode strength,
or promote the trial energy as a new homogenized reference. This closes the
specified corner-coupling omission, not all possible FE or geometry errors.
The four-fixture benchmark objective remains unachieved.

Records are in `runs/g3-quartic-corner-controls-v1` and
`runs/g3-lattice-quartic-corner-witness-v1`. The maintained driver and tests
are `quartic_edge_probe.py` and `test_quartic_edge_probe.py`.

## Analytic feasibility

Job 178 registers a degree-three tetrahedral Lagrange element through fedoo's
existing element interface. It checks all cubic monomials and their derivatives,
nodal interpolation, partition of unity and degree-four quadrature moments.
The cubic basis reproduces a nonaffine quadratic geometry mapping within
5.56e-16, and its derivatives within 1.34e-15.

A periodic solid patch contains 768 tets. The existing Tet10 control has
1,241 nodes and the enriched mesh 3,925. Both solve six prescribed strains.
Worst mixed-load stiffness errors against the analytic tensor are respectively
1.65e-12 and 4.75e-12. Independent positive 125-point reintegration agrees with
reactions and energy within 1.46e-11 absolute. Cubic periodic jumps are below
2.23e-16. This validates the experimental element and the shared helper's
unchanged default Tet10 route, not a lattice reference.

Records are in `runs/g3-cubic-reference-patch-v1`. Nine existing comparison
admission regression checks also pass after the helper extension. Unsupported
node counts no longer silently select Tet4 or Tet10 in these helpers; an
experimental element must be explicitly supplied.

## Measured lattice assembly stage

The admitted `g3-lattice-curved24-symmetric-v1` mesh has 9,907 corner nodes,
59,118 edges and 94,188 faces. Cubic enrichment therefore needs 222,331 nodes
on the same 44,976 tets, inside the existing 250,000-node stage.

Job 179 constructs this mesh, verifies periodicity and assembles its matrix
without solving. Peak RSS is 5.893 GiB under the unchanged 8 GiB ceiling.
The 666,993-square CSR matrix has 91,401,219 stored entries and occupies
1,099,482,604 bytes. This resource stage takes 54.81 seconds. It is not an
end-to-end meshing benchmark. Records are in
`runs/g3-lattice-cubic-assembly-memory-v1`.

## Initial reference solve stopped

Job 180 uses that prepared mesh. Before solving, it verifies the source
quadratic mapping at every cubic node and embeds the saved parent fields.
Mapping disagreement is 1.67e-16; integrated parent stress and energy reproduce
the archived tensor within 2.17e-11 relative. Periodic checks pass.

The run has an 8 GiB ceiling and a hard 1,800-second worker timeout. Its
240-second post-first-load continuation gate becomes unreachable while the
first load is still running. At 328 seconds after the preceding status write,
no load has completed. Stop the run rather than wait for a result that cannot
pass that gate. Final worker wall time is 381.19 seconds and peak RSS is
5.742 GiB. No lattice reference tensor or displacement archive is obtained.

Pueue cancellation initially leaves this known controller and worker alive.
Their exact argv and process group are checked, the worker group is terminated
and both known PIDs are confirmed gone before the next heavy job starts. No
other jobs are touched. Frozen partial results remain unchanged; the separate
stop record is `runs/g3-lattice-cubic-reference-jacobi-stop-v1.json`.

## Cached AMG validation and retry

PyAMG 5.3.0 is already installed and pinned in the G3 runtime requirements.
Job 181 tests smoothed-aggregation V-cycle preconditioning with three reduced
translational candidates. The original CG tolerances remain `rtol=1e-12` and
`atol=1e-14`, with independently recomputed residuals below 2e-9. Reuse is
allowed only when the complete reduced CSR matrix hashes identically.

Six-load cubic cube and saved curved Schwarz checks both pass. Cube analytic
stiffness disagreement is 4.42e-12; Schwarz reproduces its archived tensor
within 4.74e-11. Independent own-rule stress and energy disagreement remains
below 3.33e-12, and periodic jumps below 2.23e-16. Each case builds one
hierarchy and reuses it for the remaining five loads. Records are in
`runs/g3-amg-reference-validation-v1`.

The small Schwarz solve takes 6.80 seconds, versus the archived diagonal-CG
control's 6.74 seconds. This is not evidence of a speed advantage. The purpose
is to test whether the large reference becomes feasible. Benchmark defaults
are unchanged, and any future workflow comparison must use the same solver
policy on both meshers and charge hierarchy setup.

Job 182 tests the same admitted cubic lattice mesh with this validated solver.
Its interruptible 240-second first-load deadline fires inside a multigrid
coarse correction. Peak RSS is 6.976 GiB. No load completes, so there is still
no new lattice tensor or accuracy claim. The worker exits normally through
the deadline exception, without an orphaned process. Records are in
`runs/g3-lattice-cubic-reference-amg-v1`.

## Isolated system and node-block diagnosis

PyAMG's installed documentation states that its default block Gauss-Seidel
smoother becomes point-wise Gauss-Seidel on CSR input. The prior solver feeds
CSR in component-major order. Test actual node-block coupling before changing
tolerances or increasing the reference budget.

Job 183 captures the actual first-load reduced matrix and RHS through fedoo's
solver callback. It raises a dedicated capture exception rather than returning
a fake solution. The 661,029-square matrix has 91,117,159 entries. Its free
unknowns form complete three-component blocks on 220,343 nodes, verified from
fedoo's free-DOF metadata. Capture peaks at 5.850 GiB and takes 65.76 seconds.
The matrix and vectors remain local ignored archives, with hashes in
`runs/g3-lattice-cubic-reduced-system-v1/result.json`.

Job 184 permutes this unchanged system to node-major order, uses BSR blocks of
size three and retains the three translational candidates. It reaches 185 CG
iterations but does not converge before the same 240-second algorithm cap.
The hierarchy has levels of 661,029, 2,160 and 30 unknowns; setup takes 9.90
seconds. Peak RSS is only 3.048 GiB outside the FE assembly process. Reordering
alone does not establish reference feasibility. Records are in
`runs/g3-lattice-cubic-block-amg-v1`.

Job 185 adds the three rotational candidate fields used by PyAMG's own
`gallery.linear_elasticity_p1` example, retaining the original matrix, RHS,
CG tolerances and deadline. These are algebraic near-nullspace candidates,
not exact periodic rigid modes or modified physical constraints. The original
CSR residual passes after undoing the permutation. The first system converges
in 117 iterations and 165.05 seconds of algorithm wall time. Its original-system
relative residual is 7.72e-13. Hierarchy setup takes 16.19 seconds, CG 146.85
seconds and peak RSS is 2.909 GiB. This is still one load, not a full mixed-load
reference or a live workflow speed result. Records are in
`runs/g3-lattice-cubic-block-amg-rigid-v1`.

## Six-load block integration and full reference retry

Job 186 verifies the integrated block solver on the same cubic cube and curved
Schwarz controls. It passes all twelve solves and reuses the hierarchy for five
loads in each case, conditional on an identical reduced matrix hash. Cube
analytic disagreement is 3.45e-12 and Schwarz reproduces its archived tensor
within 4.42e-11. Independent own-rule stress and energy error remains below
3.16e-12, with periodic jumps below 2.23e-16. No small-case speed advantage is
established. Records are in `runs/g3-elastic-block-amg-validation-v1`.

That validation's manifest correctly records `elastic_blocks=true`, as do its
solver records, but inherits a stale scalar-method description. The frozen
description is not silently rewritten. Its correction is recorded separately
in `runs/g3-elastic-block-amg-validation-v1-method-clarification.json`, and the
maintained controller now writes the selected method explicitly.

Job 187 retries all six lattice loads through the validated integrated solver.
The actual mesh, coefficients, boundary conditions, residual tolerances,
240-second first-load deadline, 1,800-second worker cap and 8 GiB ceiling remain
unchanged. Completed fields are saved before stress/energy and 125-point audits.
The parent embedding check passes again. The first load also passes, with the
same 117 iterations and 7.72e-13 residual as the isolated system. The second
load fails before solving, while allocating a repeated `A - A.T` symmetry
check under the 8 GiB address-space cap. The rejected allocation is 695 MiB
for sparse indices; observed peak RSS is 7.105 GiB. There is no full tensor or
reference admission. Records are in
`runs/g3-lattice-cubic-reference-elastic-block-v1`.

## Hash-guarded symmetry reuse

The matrix is already hashed completely before hierarchy reuse. Check symmetry
once for the first matrix; on subsequent loads, require an identical shape,
CSR structure and data hash before reusing that verified property. A changed
matrix still fails rather than silently reusing the certificate. Independent
residual and stress/energy checks remain per load. This avoids a redundant
large sparse subtraction without changing any physical or numerical gate.

Job 188 passes all twelve integrated control solves with five verified symmetry
reuses in each case. The tensors and iterations are unchanged. Three focused
reference tests plus the nine comparison-admission regressions also pass.
These cover polynomial identities, exact quadratic-map reproduction, explicit
element selection and preserving all six solved field checkpoints in order.
Records are in `runs/g3-elastic-block-symmetry-cache-validation-v1`.

Job 189 retries the full reference with this allocation fix. Each completed
displacement field is now archived immediately, so a later allocation or audit
failure does not destroy earlier solved fields. The complete reference still
requires all six loads, own-rule stress/energy, periodic jumps, Galerkin
monotonicity and the retained 125-point sensitivity check. It solves and saves
two displacement fields, then fails during load-two reaction extraction. Fedoo
evaluates `M @ A @ X` left-to-right and attempts another large sparse matrix
product; a 349 MiB index allocation fails under the same address-space cap.
Peak RSS is 7.216 GiB. Only one force column reaches the partial result file,
so there is no complete tensor. Both field archives are intact and finite.
Records are in
`runs/g3-lattice-cubic-reference-elastic-cache-v1`.

## Vector-first reactions and verified prefix replay

The reference path now computes the equivalent `M @ (A @ X - D)` using vector
operations. It reproduces fedoo's stored MPC transpose and free/blocked identity
contributions, without forming `M @ A`. The installed dependency and benchmark
default reaction path are not patched. An explicit experimental reaction reader
selects this path for the cubic reference.

Job 190 compares native and vector-first reactions on all six cube and Schwarz
loads. Maximum relative differences are 5.62e-15 and 7.34e-14. Its two saved
prefix fields per case are accepted only after checking the current reduced
equations and reconstructing the same boundary-constrained displacements. The
remaining four loads are solved normally. Tensors, independent stress/energy
and periodic checks pass. These are twelve load checks with four verified
replays, not twelve fresh solves or a live-workflow timing. Records are in
`runs/g3-vector-reaction-prefix-validation-v1`.

Fourteen focused and comparison tests pass, including rejection of a wrong RHS
and rejection of a changed matrix even when the saved vector solves it. The
controller also now hashes the source field archive whenever a saved prefix is
used. Frozen earlier validation manifests are not rewritten.

Job 191 resumes the lattice reference from the two retained fields rather than
pay for their solves again. It verifies the identical prepared-mesh hash, each
current reduced residual and the reconstructed boundary values before accepting
them. Both retained fields pass, including a 9.56e-13 residual for load two.
The first new solve, load three, passes in 116 iterations at 9.37e-13 residual.
The 240-second deadline applies to this first newly solved load, after verified
prefix preparation. The full worker remains capped at 1,800 seconds and 8 GiB.
No field replay time is presented as a fresh meshing or workflow speed result.

It also solves and saves load four. Before load five, another reduced-matrix
allocation fails under the unchanged address-space cap. Four force columns and
four intact fields remain, not a complete reference. Records are in
`runs/g3-lattice-cubic-reference-vector-resume-v1`.

## Completed-problem cache lifetime

The replay wrapper updates its delegated problem only when the solver callback
starts, which is too late to release the previous problem before the next
reduction. Forwarding that binding immediately passes a weak-reference test
and the four-prefix control cases in job 192. It is a real retention fix, but
job 193 shows it is insufficient: replay reaches four fields and the fifth
reduction still fails. Do not claim that this binding fix alone solves the
large-case memory problem.

Inspection of pinned fedoo source also finds a problem/boundary-condition cycle,
including `problem.bc._problem = problem`. Completed problems can therefore
retain their private reduced-matrix cache until cyclic garbage collection runs.
Each cached reduction is about 1.1 GB. The reference helper now explicitly drops
that unused cache after copying the field and extracting reactions. It still
constructs a new problem for the next load, and never clears a matrix while
that problem is solving. Benchmark helper defaults remain unchanged.

A regression using actual fedoo solves and weak references verifies release of
all six completed CSR matrix objects. Sixteen focused and admission tests pass.
Job 194 also passes the two six-load control cases with four verified prefix
fields each. Numeric tensors remain unchanged. Records are in
`runs/g3-reduced-cache-release-validation-v1`.

Job 195 resumes from the four retained fields with explicit per-load cache
release. All four are revalidated against current equations and boundary values
before fresh solves. The 8 GiB ceiling, 1,800-second worker cap, first-new-load
240-second deadline and final independent audits remain unchanged. There is no
complete reference or four-fixture accuracy claim while it is pending. Records
are in `runs/g3-lattice-cubic-reference-cache-release-v1`.

## Complete same-domain cubic tensor

Job 195 completes all six loads. It revalidates four retained fields and solves
the final two in 119 iterations each. All current reduced residuals are below
9.57e-13. Explicit completed-cache release brings peak RSS to 6.341 GiB under
the unchanged 8 GiB ceiling. Worker wall time is 369.81 seconds, including
verified replay. This is not a fresh workflow timing.

Independent own-rule stress and energy reproduce reactions within 7.10e-12
relative. Periodic displacement jumps are below 2.23e-16. Positive 125-point
reintegration changes worst mixed-load energy by 6.87e-6 relative, below the
unchanged 1e-4 quadrature-sensitivity limit. The smallest eigenvalue of the
parent-minus-cubic tensor is +0.000140585, consistent with variational enrichment
on the identical quadratic physical domain. All six fields remain archived.

The tensor passes FE consistency, geometry and quadrature gates. It still has
status `solved-needs-reference`: one higher-order same-domain result does not
establish convergence to the fixed implicit specimen. All 24 earlier empirical
reference models remain retained. The next bounded stage audits these saved
fields in a fresh process and compares every retained model without re-solving.

Job 196 passes that fresh-process check at 0.514 GiB peak RSS in 53.94 seconds.
It uses explicit Tet20 registration rather than a Tet4 fallback. Affine analytic
patches pass for every rule, all sampled Jacobians are positive, and 64-to-125
point fixed-field energy changes by only 9.54e-7 relative. Records are in
`runs/g3-lattice-cubic-assessment-v1`.

Against cubic24, uniform Tet10 meshes at 24, 32 and 34 cells differ by 0.87761%,
0.56789% and 0.49251%. The finest adaptive tensor differs by 0.52086%, in the
opposite energy direction. The actual tuned Microgen Tet4 tensor differs by
3.49233%. These pairwise values do not resolve specimen accuracy. The cubic
tensor normalizes retained-model disagreements of 0.15415% to 1.75050% across
all 24 models. No tested Tet10 candidate passes every retained model at 1%.
Neither reference selection
nor a four-fixture timing victory follows from this result.

Before another reference solve, test a stronger actual Microgen control.
Its previously admitted lattice parent needs 326,716 straight Tet10 nodes.
Second-pass hausd 0.0009 reduces the lift to 267,550 but fails quality at
0.097405. One geometry-only stage changes that second-pass value to 0.001,
retaining VTK51, first-pass 0.0015 and all admission gates. Check exact lift
size and admission before any solve. This is explicitly tuned Microgen,
not its recorded defaults and not MMG on a Meshers seed.

## Worst-mode discrete error localization

Job 202 uses the worst retained-model strain of the finest adaptive Meshers
candidate, dominated by deviatoric normal strain. On the unchanged curved24
physical domain it embeds the parent Tet10 field exactly into the cubic space
and integrates the energy of its difference from the cubic solution. The
own-rule Galerkin identity reproduces the corresponding parent-minus-cubic
quadratic form within 9.50e-9 relative. This is a discrete enrichment check,
not a local stress or continuum certificate.

The highest-energy 20% of coarse tets contain 98.394% of the correction energy.
The largest-volume 20% contain only 48.515%. Elements incident on cap vertices
are 6.990% of the mesh but contain 77.461% of the correction energy. Wall-vertex
incident elements contain 95.416%. Positive 125-point integration changes total
correction energy by 0.05995% relative to the own rule; the 64-to-125 change is
0.00800%. These percentages concern the small correction, not the full tensor's
separate quadrature admission check. Records are in
`runs/g3-lattice-cubic-error-localization-v1`.

This supports one cap-focused sizing pilot, not an inference about the error
location on every refined mesh. The current finest adaptive cap edges still
reach 0.06241 despite median length 0.01080. Start from the admitted 24-cell,
adapt-tolerance 0.004 recipe and add a periodic size callback with target 0.015
on box planes. A smoothstep transitions to maximum length 0.1 over distance
0.05. Keep native budgets, geometry, quality and periodicity gates unchanged.
Check geometry and counts before selecting any solve. No native API change or
oracle field is added to the workflow.

Job 203 passes geometry at 159,147 nodes and 99,604 tets. Job 204 completes
six fresh 15-point loads with the validated explicit AMG preset. Fresh audit
205 reproduces reactions within 7.08e-12 and Hill-Mandel energy within
2.73e-12; periodic jumps are zero. Positive 125-point reintegration changes
full mixed-load energy by 4.25e-7 relative. The cap-sizing callback is positive
and has bitwise-equal opposite-plane traces in its recorded preflight.

Its worst retained-model error is 1.411382%, versus 1.427185% for the
123,570-node untargeted baseline and 1.269882% for the finest Meshers result.
Adding 35,577 nodes produces little accuracy gain and does not meet 1%.
Do not promote this as a workflow improvement or keep tuning cap lengths.
Records are in `runs/g3-lattice-cap-size015-band005-comparison-v1`.

The host suspends for 56 minutes during job 204, confirmed by the system
journal. Its monotonic timer excludes sleep, and qrun records a matching long
reporting gap. All timing claims for that run are excluded by the separate
`runs/g3-lattice-cap-size015-band005-timing-exclusion-v1.json`. Frozen numerical
results remain unchanged. The standard Job reporter is not replaced or exempted
to hide the host interruption.
