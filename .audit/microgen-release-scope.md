# Microgen release scope audit

Source: study commit `4eea25af3d6f0aa6f8eb15a60030458263224183`.
Baseline: freshly fetched main `d1324b7e1e94d4ced47a781e4c118d617b0589d5`.
The local `main` branch was stale; it was not used as the release baseline.

This is a scope and dependency audit, not a proof of universal mesh admission
or a replacement for the cross-platform release checks.

## Release decision

Prepare this branch from current main, carrying the library changes, regression
tests and three relevant documentation updates. Do not merge the full study
branch. Its 13,076 changed paths include 12,650 experiment paths and 363 handoff
paths. Those results stay preserved on the study branch.

Every changed library source and native fixture path is classified in
`microgen-release-scope.tsv`. All crate files in this branch were checked byte
for byte against the source used for the locally tested default-feature wheel.
The crate manifests add feature switches but no new Rust dependencies relative
to current main. Root packaging and publishing workflows already exist on main.

| Change group | Needed for this integration | Decision |
| --- | --- | --- |
| Scalar field API, owned arrays, compiled batch evaluation | Yes | Keep the existing API; add no Microgen-specific SDK entry point |
| Geometry, determinant, quality and resource admission | Yes | Keep; removing failed experiments must not remove their guards |
| Adaptive Tet4 and pending projection restoration | Yes for the measured spinodoid recipes | Keep, including `ddc376d`; this is a correctness fix, not a new sampling policy |
| Periodic cap repair, physical image evaluation and material trace matching | Yes for periodic specimens | Keep; do not require fabricated images without material support |
| Exact-tie derivatives and affine/nonlinear constraint preservation | Yes for the broader primitive and composed implicit scope | Keep; TPMS-only tests cannot justify deleting these paths |
| Tet10 elevation and periodic midside correspondence | Yes for the branch's optional order 2 | Keep; corner quality is not a curved-volume or FEA certificate |
| Static interpolation tracing | Not required by current linear callable grading | Keep tested support; introduce no grading table dependency or default |
| Curved wall fitting and interior-refinement budget | No | Leave opt-in and off by default; no automatic fitting or reference promotion |
| Surface extraction/remeshing and wide constraint masks | No | Leave feature-gated and disabled in the release wheel |
| Linear observations and supplied-background PMY experiments | No | Keep diagnostic hooks gated; do not expose or enable them through Microgen |
| Volume correction probes and periodic feature research | No | Keep under `cfg(test)`; do not turn investigations into production corrections |
| Direct fedoo conversion | No | Retain the small optional export; no solver dependency is added to Microgen |
| Experiment meshes, solver fields, dashboards, handoffs and logs | No | Exclude from this branch; preserve on the study branch |

Large shared modules are not safe to delete wholesale: scalar, periodic,
polyhedral and adaptive paths share projection and topology helpers. The
release narrows activation and repository contents without discarding measured
correctness fixes. More aggressive removal needs separate before/after evidence.

## Microgen branch

Keep an opt-in selector on the existing volume methods and one shared adapter
for three callers. That adapter imports the optional SDK, forwards options,
converts connectivity, applies local placement and preserves periodic metadata.
There is no new geometry generator, fallback or automatic MMG pass.

Remove the integration Markdown file, its README pointer and diagnostic/quality
arrays computed solely for the integration experiment. Independent quality
checks belong in the tests. Keep the shared Spinodoid evaluator, TPMS raw-offset
semantics, sampled/mapped-domain rejections and unchanged default VTK path.

No additional Meshers runtime code change is necessary for the current adapter.
Its required API is `generate`, the returned mesh arrays and `compile_field`.
It needs neither `generate_surface` nor linear-observation hooks. There is no
reason to add a PyVista dependency to Meshers for this integration.

## Verification and remaining release gates

The Mac arm64 wheel was rebuilt from the candidate with no optional Rust
features. Its surface API and linear observations were confirmed disabled.
The combined SDK and Microgen selection passed 340 tests, skipped 15 and
deselected the two previously reproduced upstream CAD spinodoid failures.
The first combined test attempt stopped at collection because the new test
environment lacked h5py. Installing the recorded test dependency `h5py==3.16.0`
resolved collection. It is not a Microgen runtime dependency.

The clean release checkout also passed `cargo +1.95.0 test --release --locked`:
93 passed, zero failed, two ignored diagnostic tests. The bounded native run
used 1.484 GiB peak RSS. Its receipt is
`/Users/kmarchais/dev/meshers/.local-work/microgen-release-scoped-rust-tests-v1`.

Local receipts remain on the study checkout:
`/Users/kmarchais/dev/meshers/.local-work/microgen-release-default-build-v1`,
`/Users/kmarchais/dev/meshers/.local-work/microgen-release-default-tests-v1`, and
`/Users/kmarchais/dev/meshers/.local-work/microgen-release-default-tests-v2`.

Before publication, run the existing main-commit CI, installed-wheel checks
on Linux, Windows and both macOS architectures, source archive installation,
and performance gates on the exact final source. No package has been published.
Scientific comparison and lattice FEA convergence gaps are not closed by this
integration or by passing software tests.
