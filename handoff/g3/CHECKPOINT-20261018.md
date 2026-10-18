# Shutdown checkpoint, 18 October 2026

Kevin requested stopping G3 within one hour of the message received around
20:53 UTC. Finish the bounded study, preserve its evidence and stop work for
shutdown. The four-case Meshers superiority goal is not achieved. Do not mark
it complete or resume work automatically across the requested shutdown.

## Conclusion

- Schwarz has one supported 1% matched-target win against an explicitly
  modified actual Microgen one-pass recipe. Five interleaved live pairs give
  median workflow times 7.79 versus 40.50 seconds. This is not Microgen defaults
  or a result covering all four specimens.
- Gyroid and spinodoid still lack admitted actual Microgen/MMG comparators
  within the tested recipes and bounds. Meshers has admitted candidates, but
  no matched-accuracy workflow victory is established for these cases.
- Lattice remains unproven. All 24 empirical reference models are retained.
  The strongest admitted actual Microgen-derived straight Tet10 control has
  worst model disagreement 1.004354%, versus 1.269882% for the best tested
  Meshers candidate. Neither passes the unchanged 1% target. Do not round the
  Microgen miss into admission or discard this stronger comparator.
- The cap-focused Meshers pilot is negative. It uses 159,147 nodes and has
  worst model disagreement 1.411382%, barely better than its 123,570-node
  untargeted baseline at 1.427185%. It does not beat the best Meshers result
  or the stronger Microgen control. Stop blind cap-length tuning.

## Latest evidence

Read these files in `/home/kmarchais/dev/meshers/experiments/adaptive`:

- `TET10_RESULTS.md`: actual Microgen lattice Tet10 control, five-load CG
  timeout, explicit validated AMG retry and fresh-process field audit.
- `CUBIC_RESULTS.md`: complete same-domain cubic reference, independent
  positive quadrature and worst-mode error localization.
- `RESULTS.md`, `PLAN.md` and `.audit/g3-microgen.tsv` at the repository root
  relative audit directory contain conclusions and the decision trail.

Important run directories under
`/home/kmarchais/dev/meshers/experiments/adaptive/runs`:

- `g3-lattice-cubic-reference-cache-release-v1`: six-load reference-only cubic
  enrichment, 222,331 nodes. FE consistency passes; not a continuum certificate
  or a fresh workflow timing. Four fields are equation-verified replays.
- `g3-lattice-microgen-straight-tet10-hausd001-amg-v2`: six fresh loads on an
  actual tuned VTK51/two-pass MMG parent, manually elevated and converted
  directly to fedoo. 238,411 nodes, 158,252 tets, peak RSS 4.378 GiB. Public
  Microgen Tet10 conversion remains unsupported. This is a saved-mesh study.
- `g3-lattice-tet10-retained-comparison-v1`: all model families, audited input
  hashes and worst mixed-load witnesses. The worst model for every tested
  candidate is the retained uniform-seed first-order extrapolation.
- `g3-lattice-cubic-error-localization-v1`: on the coarse unchanged domain,
  cap-vertex incident 6.99% of tets carry 77.46% of discrete enrichment energy.
  This does not establish error locations on every refined mesh.
- `g3-lattice-cap-size015-band005-solve-v1`: geometry-admitted live pilot with
  cap target 0.015 and smooth transition band 0.05. Six loads use explicit AMG
  and default 15-point quadrature. All fields are retained.
- `g3-lattice-cap-size015-band005-final-audit-v1`: fresh field audit and positive
  64-point and 125-point reintegration. Reaction error 7.08e-12, Hill-Mandel
  error 2.73e-12 and periodic jumps zero. Higher-rule energy shift 4.25e-7.
- `g3-lattice-cap-size015-band005-comparison-v1`: negative all-model comparison.

The timing exclusion is the separate frozen file
`g3-lattice-cap-size015-band005-timing-exclusion-v1.json` in that runs directory.
The system journal confirms suspend from 19:51:26 to 20:47:41 UTC. Qrun records
a 3,378-second reporting gap. Monotonic worker durations exclude the sleep
interval, so the recorded 367.23-second pilot is not admitted timing evidence.
Do not edit the original result or replace the standard Job reporter to hide
the host interruption. The final fresh audit runs after resume.

## Restart safely

Use `/home/kmarchais/dev/meshers`, branch `handoff/g3-microgen-20261017`.
The source state is the commit containing this checkpoint, following `a26d2ec`.
Latest changes are research helpers and documentation, not new native meshing
features. The existing Linux extension SHA256 remains
`d41ccb980fbf3339a86e8f9a195e675abaa3910be67de8589dca138fd401b700`.

Support repositories remain at the documented commits. NPZ meshes, fields and
captured sparse matrices are ignored but retained locally; selected small text
records and drivers are committed. Do not delete ignored archives or rerun
expensive solves to recover fields that already exist.

Use `.venv` from that repository and set:

```sh
export LD_LIBRARY_PATH=/home/kmarchais/dev/meshers/.local-work/g3-runtime/root/usr/lib/x86_64-linux-gnu
export OPENBLAS_NUM_THREADS=1
```

Read global instructions, the computers inventory and local queue/progress
skills again. Check shutdown state and the shared queue before starting. No
study job should remain active at this checkpoint. Preserve the 8 GiB worker
ceiling, conservative threads, staged limits, one queued experiment at a time,
fixed admission gates and no delegates. Use Rust 1.95.0 for native rebuilds.
Before GitHub operations switch `gh` to `kmarchais` and verify auth status.

## Next decision, not another queued sweep

The oracle same-domain enrichment localization does not predict a useful live
cap-sizing gain. More global refinement or another cap-length sweep is not
justified by this result. A possible next diagnostic is validating a cheap
Tet10 stress-recovery indicator against the archived Tet10-to-cubic correction
in the actual worst mixed mode. Keep that discrete energy study separate from
local stress accuracy. Only select a live sizing trial after reading such a
diagnostic; charge all predictor solves and preparation to workflow time.

Reference uncertainty and missing admitted Gyroid/spinodoid comparators remain
the main obstacles. No reference model is removed, no limit is raised and no
new timing victory is inferred from this shutdown checkpoint.
