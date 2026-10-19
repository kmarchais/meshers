# G3 study checkpoint, 19 October 2026

The four-fixture matched-accuracy goal remains unachieved. This checkpoint
preserves one verified native admission fix and a smaller audited spinodoid
configuration. It does not schedule shutdown, change admission, or declare
the active goal complete.

## New result

Native pending-wall settling now tries bounded quality-preserving partial
steps before splitting facing edges. The red regression fails before the
change; all 25 core and 37 focused Python tests pass afterward. The native
change and first validation are pushed in commit `56e1b59`.

The old fixed spinodoid curved24 recipe changes from rejected to admitted.
At adaptation tolerance 0.002 it has 105,618 nodes, distance 0.00173079,
zero failed projections, and worst model disagreement 0.13529%. Its slower
640.08-second single workflow is not a performance improvement.

The useful tuned recipe retains cells24, one optimizer pass, the compiled
fixed FFT16 field, minimum length 0.002, and curved Tet10, but sets adaptation
tolerance 0.003. It has 73,951 nodes and 45,994 elements. Distance is
0.00133298, volume error +0.02586%, corner quality 0.26858, sampled Jacobians
are positive, and periodic caps match. Six prescribed-strain loads take
351.27 seconds, total workflow is 364.07 seconds, and peak RSS is 2.61 GiB.

Independent stress, energy, and periodic-jump audits pass. Worst mixed-load
disagreement is 0.10333% to 0.18439% across all six unchanged reference models.
Positive 125-point fixed-field energy differs from the solved tensor by
3.01e-8 relative. These are empirical stiffness results, not continuum or
local-stress certification. No interleaved timing claim is made.

The archived pre-patch binary rejects the exact 0.003 recipe at distance
0.00237022. Its binary hash and original source are independently checked.
The observer workspace source hash in its raw manifest is not the archived
binary's source hash; `attribution.json` records that distinction.

## Controls and failures

- All four existing geometry controls remain admitted. Schwarz curved8 and
  Gyroid curved9 are bitwise unchanged. Lattice adaptive24 matches coordinates
  and canonical cell node sets after renumbering. Spinodoid32 remains admitted
  but its new tensor is not measured; retain the older audited references.
- The 20-cell seed fails quality. Two optimizer passes fix that gate but leave
  two failed independent boundary projections. Neither receives an FE solve.
- Actual Microgen seed41 has 190 retained-cap quality obstructions. Its
  separate 300,000-tet seed-only screen uses 0.30 GiB. Default seed, final mesh,
  solve caps, and the 8 GiB ceiling are unchanged. No MMG or FE follows.

## Remaining conclusion and next work

Schwarz retains the only narrow demonstrated matched-target workflow win.
Gyroid and spinodoid still lack admitted actual Microgen/MMG comparators.
Lattice's strongest Microgen straight Tet10 control reaches 1.004354% worst
retained-model disagreement, versus 1.269882% for the best tested Meshers
candidate. Neither passes 1%. Keep all 24 lattice models and their fields.

Stop spinodoid coarsening here. Return to missing comparator coverage and
lattice reference disagreement before another expensive spinodoid solve.
Read `MICROGEN_RESULTS.md`, `PERIODIC_RESULTS.md`, `PLAN.md`, `RESULTS.md`,
and `.audit/g3-microgen.tsv`. Do not repeat failed required-entity grading,
selective midpoint restoration, affine corner postprocessing, or cap-only
lattice sizing without new evidence.

## Local runtime

Current native SHA256 is
`908e326211792142f83f854ad2696cd879d927bd81deacf38eb6a0cb63bc4dc1`.
The archived pre-patch SHA256 is
`d41ccb980fbf3339a86e8f9a195e675abaa3910be67de8589dca138fd401b700`.
Its package is in `.local-work/g3-builds/baseline-090a71c`.

Preserve ignored NPZ meshes, fields, and native binaries locally. Use Rust
1.95.0, two build threads, one numerical thread, the pinned `.venv`, and
`LD_LIBRARY_PATH=/home/kmarchais/dev/meshers/.local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu`.
All study jobs through 243 finish before this checkpoint. No successor is
queued here. The shared queue belongs to other sessions too; do not touch
unrelated jobs or timers.
