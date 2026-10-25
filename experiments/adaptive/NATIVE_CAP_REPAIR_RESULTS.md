# Native periodic cap repair investigation

The fixed linear spinodoid failure now passes independent strict mesh admission
in job 470. The unchanged 52-mode specimen and fixed samples are hash-identical
to the baseline. Distance improves from 0.004834255 to 0.001079175. The repaired
mesh has 30,336 nodes and 131,011 tets, versus 30,499 nodes and 131,598 tets.
Volume, coherent facets, positive determinants, periodic nodes and caps pass.
Minimum quality is 0.100652964, down from 0.152883712 and just above the
unchanged 0.1 gate. This fixes the geometry failure; it is not dominance on
all quality measures or a matched-accuracy performance victory.

The final source passes 104 native tests, one ignored research test, Clippy,
and 11 awake-worker/frozen-recipe controls. Rebuilt arm64 binary SHA256 is
`0f94afb012903f917b1d2c18ee7e7f5330a048cb8f2394fb3a446da74bf0f33e`.
Complete Python suite 471 passes 182 tests with one skip. All four curved
routes pass in 472 through 475. Spinodoid and lattice meshes change and need
fresh mechanical evidence. Lattice preflight 476 passes at 65,721 nodes and
41,340 tets. Job 477 completes six fresh AMG loads and audit 478 passes reactions,
energy, periodicity and positive 64/125-point reintegration. All 24 lattice
models remain retained. Assessment 479 gives 2.412306220% worst disagreement,
versus 2.427923243% before repair, still failing 1%.

Fresh curved-spinodoid loads 481 and field audit 482 pass. Assessment 483
checks every entry in four archived reports. They contain six distinct model
tensors repeated across 24 entries, not 24 independent reference tensors.
Worst disagreement is 0.183380673%, versus 0.184393513% before repair, passing
1%. Both old and new meshes fit the smaller historical 100,000-node solve
guard despite different recorded guard values. Study admission stays 250,000.
The shared assessor replay 484 preserves every lattice tensor and comparison
exactly. Python-feature Clippy 485 passes. Gyroid and Schwarz have exactly
identical point, Tet10 and Tri6 arrays to the accepted periodic-normal build.
No new solve is needed for their unchanged approximation spaces.

No stiffness solve is performed for the new linear spinodoid route, and no
new reference or four-case comparative winner is promoted.

## Fixed failure and negative guard ablation

Job 420 reproduces the failure on native binary
`eab1823f1234fd39223a89ec103df7d87db90c2418f59dceb7cfeebf59ceabd7`.
The unchanged linear recipe has 30,499 nodes and 131,598 tets. Distance is
0.004834255, versus the fixed 0.002 gate. Minimum corner quality is
0.152883712. The adapter aims for 0.001 and reports maximum error 0.006709921.
The fixed 52-mode field and specimen/sample hashes are preserved.

Job 421 changes only the internal operation guard from 0.3 to 0.1. The hard
admission quality remains 0.1. This ablation still fails, at distance
0.004631774 with 29,736 nodes, 127,707 tets, and minimum quality 0.123016421.
Maximum adapter error remains exactly 0.006709921. Do not promote a lower
guard default or run another guard sweep from this result.

Read-only job 422 finds that every baseline wall vertex already lies on the
field to about 9.9e-14 first-order distance. Unfinished baseline vertex
projection does not explain the miss. The worst probe is the midpoint of a
curved wall edge on a periodic cap. Its endpoints are on the exact field.

## Coupled periodic edge obstruction

Job 423 extracts the exact incident stars on both cap copies. There are two
parents on one side and one on the other. Full projection of a new midpoint
would leave minimum quality 0.165118 on the first side but 0.077945 on the
second. Accepting the second would violate the fixed 0.1 gate. A quarter
projection leaves 0.115143 on the obstructed side and is valid under the
original split floor. The extracted local neighborhood has 48 nodes, 68 tets,
and the unchanged 52 Fourier modes. It is a local incidence control, not a
new admitted study mesh.

Native control 424 inserts that valid partial step but stops above target
with a pending wall vertex. Diagnostic 425 shows that another call advances
the same geometry without any parameter changes. Accounting-only candidate
426 reduces error further but still stalls. Trying all facing edges in 427
does not help and is reverted. No package promotion follows from these
negative controls.

In diagnostic 428, three existing local quality relocations around the blocked
star let projection finish. First-order cap error becomes 0.001754939, below
the predeclared public 0.002 geometry gate. The exploratory control still
fails its stronger internal-aim assertion of 0.001. That archive stays failed.

The final regression checks the actual fixed public 0.002 gate and 0.1 quality
floor while keeping the adapter aim at 0.001. Job 429 runs this final regression
against the unchanged production implementation and fails. Job 430 runs the
identical assertions against the candidate and passes in 0.03 seconds.
This is an explicit regression correction, not a rewrite or reclassification
of the earlier exploratory controls.

## Candidate implementation and remaining checks

`settle` reports whether pending projection actually moved a vertex or changed
a cavity. Error enforcement counts that progress rather than stopping solely
because no additional error edge was split. When insertion is blocked and
neither projection nor splitting makes progress, enforcement tries the existing
guarded relocation on vertices of the blocked edge stars. Every
existing projection, topology, minimum-quality, length, periodic, and tetrahedron
budget check remains active. Successful relocations are counted in diagnostics.

No new sizing value, guard default, specimen, smooth-union order, distance
threshold, quality floor, or node/tet cap is introduced. The unsuccessful
alternate-facing-edge implementation is absent from the candidate. Debug
prints remain only in archived diagnostic source snapshots.

Initial helper candidate 431 passes 102 native tests with one ignored
research test. Native build 432 has SHA256
`659d52436795c10ad74dc4c53672e5104b3deee8db627298fa5998c34166c6b3`.
Whole-specimen job 433 still fails at distance 0.004631774. Counts remain
30,499 nodes and 131,598 tets, while minimum quality improves to 0.161066270.
The local regression and complete suite do not promote that candidate as a
whole-case fix.

A fresh-state control in 434 also fails without manual midpoint insertion.
Diagnostic 435 identifies the ordering obstruction. Enforcement handles an
edge of length 0.046458 with error 0.001748 before the shorter cap edge of
length 0.040479 with error 0.004854. The first split spends quality margin
needed by the worse defect. The bad edge remains present afterward.

The next candidate canonicalizes edge identities, retains the worst duplicate
error report, and prioritizes geometry error before length. This changes no
sizing value or admission threshold. Fresh-state control 436 then reaches
0.001961840 with two error splits, 14 guarded relocations, and no pending
vertex. Clean final control 437 checks both partial and fresh states, pending
completion, and every live tet's 0.1 quality floor. It passes with diagnostic
prints removed. A pipeline-order change was considered but not implemented.

Native build 438 has SHA256
`625f1435250b2ad2fa5314a941fcfd7cd2f0ce420b3b224d1821abf5ec2e350a`.
Whole-field replay 439 still fails at distance 0.004631774. This candidate
is preserved for reproducibility and is not promoted.

Strict public native control 440 fails at sampled distance 0.004648402
against tolerance 0.002. Analytic control 441 exposes another shared guard
bug. When a wall nearly touches a periodic cap, its 3D normal turning can
hide a much larger error of the wall's trace within the cap plane. The
control has sag below 0.0001 and cap-trace error 0.00625. Guard correction
442 detects that error but prevents an already invalid local star from
improving. Candidate 443 permits a changed edge only up to the larger of
its old error and the existing target. All four local cap controls pass.
Good edges remain protected; the public geometry gate stays 0.002.

The strict whole-field test 444 still fails at 0.004609891. Preserved accepted
binary `eab1823f1234fd39223a89ec103df7d87db90c2418f59dceb7cfeebf59ceabd7`
also fails the same strict call in 445, at 0.004853506. These controls establish
an existing failure, not an admitted repaired case.

Job 446 adds a geometry-only closing pass after final quality improvement.
The strict whole-field error remains exactly 0.004609890725. The extra pass
is removed; its source and hashes are archived in that run. Trace 447 shows the same worst cap error in all 16 repair rounds.
Even the quarter-projected midpoint would have minimum quality below
0.069, versus the hard 0.1 floor. Relocation converges without making that
insertion feasible. Job 448 tests existing quality-increasing interior edge
removal within blocked stars. These operations preserve the surface. Job 448 leaves the full-field error unchanged and its interior-edge attempt is
reverted. Diagnostic 449 measures insertion positions on that exact cavity.
It finds an unprojected position 0.2 along the edge with minimum child quality
0.133069101. It also tests one hard-floor relocation, which improves the
whole error to 0.004112444 but still fails. That diagnostic mutation is removed.

Candidate 450 preserves midpoint insertion as the first choice. If every
midpoint child fails the existing quality floor, it tries quarter-edge positions
with the existing full and partial projection choices. Every child must meet
the same quality and orientation checks; partial projection uses the existing
pending-node mechanism. Whole-case 450 reduces sampled error to 0.003089486 but still fails the
0.002 gate. Diagnostic 451 checks whether repair is still advancing at its
bound or has reached a new obstruction. The whole-specimen regression remains
the acceptance check. This is a local cavity operation, not a specimen or sizing sweep.

Trace 451 leaves two pending wall nodes through the late repair rounds with
no further error splits. Candidate 452 tests existing relocation at the hard
user quality floor only inside blocked repair stars and only after ordinary
relocation fails. This is a diagnostic quality trade, not a new default guard.
Job 452 still fails at 0.003132254, slightly worse than quarter insertion
alone. Its hard-floor relocation is reverted. Build 453 preserves the
quarter-insertion candidate as SHA256
`cfae5059346e6cc855f9c251d30c1324f21a34a6c226e8149b771bda3288d263`.
This package still fails the strict native whole-field control. Job 454
replays the unchanged loose native call with independent external gates so
its remaining bad patch can be extracted without changing admission.
Independent whole-field audit 454 rejects that package at distance
0.003085223, with 30,320 nodes and 130,980 tets. Compared with baseline 420,
minimum quality improves from 0.152883712 to 0.196677363. The periodic pairs,
cap triangles, positive determinants, facets and volume gate pass. Geometry
still fails 0.002. Audit 455 finds only two off-field wall vertices, both below
0.000338 first-order wall distance. The worst distance remains a cap chord.

Read-only extraction 456 finds that the final quality-improved cap star admits
a valid partially projected midpoint. This differs from the earlier pre-quality
obstruction. Candidate 457 closes geometry after the final quality pass using
zero additional sizing rounds. The strict public whole-field regression passes
for the first time, including its original 0.002 distance, 0.1 quality and
250,000-node/tet checks. Full native suite 458 finds a regression in the local cap-star control. The
whole spinodoid test passes, but the small control leaves a pending wall
vertex. The suite stops with 62 passed, one failed and one ignored unit
test before the integration tests. Diagnostic 459 inspects that pending
projection without changing assertions. The candidate remains held.
Scheduling-only candidate 460 leaves the local pending regression and is
rejected. Restricting quarter insertion to cap curves in 461 restores the
local control but misses strict full geometry at 0.002048200. Concurrent
quality repair in 462 improves that to 0.002031349 and still fails. Extra
interior edge and face swaps in 463 change it by less than 5e-10; those
ineffective additions are removed.

Candidate 464 distinguishes the operation phases. During early sizing, the
midpoint-first fallback can use quarter positions before the full quality
stage. During late error repair, quarter insertion is restricted to wall-cap
curves. Existing quality-improving relocation acts on blocked stars. All ten
adaptation controls pass together, including the unchanged pending-empty
local control, strict whole fixed spinodoid, periodic image controls, sphere,
and Gyroid band. No threshold or sizing-round count changes. Complete native
suite 465 passes 104 tests, with one ignored research test. Clippy 467 also
passes. The failed archives remain failed. The queue reports a 116-second
progress gap in 465; macOS logs confirm idle sleep during that verification
run. Its duration is not a comparative performance measurement.

Bounded Mac workers now hold a per-command `caffeinate -i` assertion and record
monotonic and real elapsed time so suspension or clock differences remain
visible. This changes no persistent power setting, timeout or memory ceiling.
Controls 466 check real success, exit-code propagation, timeout, RSS stop and
one-thread BLAS. Expanded controls 468 also observe the actual macOS assertion
and verify frozen specimen/recipe identities. All 11 pass. Native Python
build 469 is next; independent full-mesh and four-case checks remain pending.

No rebuilt package, independent admission or performance result is promoted
from the native regression alone.
No sizing sweep or gate relaxation follows from these failures.

The complete native/Python checks and all four curved geometry replays are
completed in 465 through 475. Changed lattice and spinodoid meshes receive
fresh six-load fields and independent audits in 477 through 484. This
validates the native bug fix while preserving the lattice accuracy failure.

The last complete lattice reference investigation remains the fixed-domain
quadratic level5 result at 1.325163420% worst disagreement, still failing 1%.
The next level is below its existing trial-gain screen and is not solved.
All four original cases and all 24 lattice models remain in scope.
