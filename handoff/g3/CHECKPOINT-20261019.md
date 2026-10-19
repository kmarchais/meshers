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

Jobs 250 through 255 subsequently close one all-load energy-sizing study.
A non-oracle, non-recovery stress-energy indicator captures 94.02% of full
same-domain cubic correction trace and at least 92.24% for every mixed load.
Its localization screen passes without reversing the earlier recovery failure.
The first live half-length target exceeds the unchanged tet budget. A single
geometry-only floor correction admits 204,956 nodes, above its 200,000-node
solve stage. A fresh charged workflow then uses the previously validated
250,000-node stage under 8 GiB, without reclassifying the old rejection.

All twelve fresh coarse and final loads complete. Total charged workflow is
812.17 seconds and peak RSS is 5.573 GiB. Independent field and positive-rule
audits pass, but worst disagreement across all 24 unchanged models is
1.352008%, worse than the previous best Meshers result. Higher-rule energy
changes by only 2.49e-7 relative. Stop that energy-sizing recipe. Read
`STRESS_RESULTS.md` and `runs/g3-lattice-live-energy-floor025-v1`.

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

Before volume correction, improve the independent volume estimate. The
lattice reference still uses four scrambled Sobol sets of 2^18 points. Their
absolute spread 0.00015640 exceeds the apparent native curved-mesh deficit,
about 0.000078. Do not bias geometry to this uncertain mean or to the stiffness
models. Preserve the frozen estimates and fixed gates when auditing them.

Job 258 subsequently performs the volume audit: four nested 2^24-point
sequences reproduce the original prefixes exactly. Mean volume is
0.2482374460, with absolute spread 0.0000137091 and final nested mean shift
0.0000007004. Both fixed precision screens pass. The finest Meshers volume
error is -0.001898% at that mean, and its empirical replicate range includes
zero. Do not add a volume offset or reinterpret stiffness agreement from this.
Old outputs and all 24 stiffness models remain unchanged. Read
`VOLUME_RESULTS.md` and `runs/g3-lattice-volume-reference-audit-v1`.

Next inspect the raw smooth-union gradient mismatch at exact branch ties.
`expression.py` differentiates min and max ties by selecting one operand,
but differentiates abs at zero as zero. A focused smooth-composition
regression should prove the mismatch before changing that policy. Existing
manual symmetric-minimum and void-clamp representations are still not raw
default support. Do not weaken periodicity or claim a workflow victory from
a compiler fix alone.

## Local runtime

Current native SHA256 is
`908e326211792142f83f854ad2696cd879d927bd81deacf38eb6a0cb63bc4dc1`.
The archived pre-patch SHA256 is
`d41ccb980fbf3339a86e8f9a195e675abaa3910be67de8589dca138fd401b700`.
Its package is in `.local-work/g3-builds/baseline-090a71c`.

Preserve ignored NPZ meshes, fields, and native binaries locally. Use Rust
1.95.0, two build threads, one numerical thread, the pinned `.venv`, and
`LD_LIBRARY_PATH=/home/kmarchais/dev/meshers/.local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu`.
All study jobs through 258 finish before this checkpoint. No successor is
queued here. The shared queue belongs to other sessions too; do not touch
unrelated jobs or timers.
