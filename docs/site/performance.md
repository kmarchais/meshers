# Performance and regression checks

Use `threads=1` as a starting point. Native workers share memory and release the
Python GIL; extra workers do not parallelize Python callbacks. Compiled fields
avoid callback overhead. `batch_size` limits callback evaluation batches,
including optional intersection surface verification.

Measure representative geometries before changing thread count, resolution,
optimization passes or evaluator. Keep those settings identical when comparing.

## Every pull request

The Performance workflow builds the base and proposed source into separate
wheels and installs them in separate environments with the same Python and
NumPy versions. It runs both on the same Linux runner, alternates their order,
and measures five fresh processes per case after one discarded warm-up pair.
Medians are compared, and every raw sample is uploaded with the source hashes.

Cases cover compiled and callback periodic Gyroids, a graded sheet, a mapped
named intersection, Schwarz P, a fixed Fourier spinodoid and an ordered lattice. Each uses 24 lattice cells and two
optimization passes. Timing includes field preparation and meshing, but excludes
imports, independent validation and export. Peak resident memory is the Linux
process high-water mark through mesh generation, including imports and JIT.

| Metric | Fails when |
| --- | --- |
| Median elapsed time | More than 25% plus 20 ms above base |
| Median peak RSS | More than 15% plus 8 MiB above base |
| Minimum, maximum, mean, median, p05, p10 or p90 MMG quality | More than 1% below base |
| Sampled surface error | More than 5% plus 1e-10 above base |
| Mesh volume | Changes by more than 0.5% |
| Tetrahedron count | Changes by more than 1% |

Every sample must also have finite coordinates, positive cell volumes, closed
boundary edge incidence and the requested periodic correspondence. Its minimum
quality must be at least 0.01 and sampled surface error at most 0.1.

A changed mesh count or volume may be intentional, but it needs review rather
than being counted silently as a performance improvement. Do not weaken a gate
just to make a change pass. Inspect raw measurements and reproduce an unexpected
failure before changing thresholds.

Shared runners are noisy. These allowances detect substantial regressions;
they cannot promise zero slowdown, cover every geometry, or prevent cumulative
small changes. Review the uploaded trends, and compare against the release tag
periodically with the workflow's manual baseline input. The initial commit
compares against itself to establish measurements and exercise the gate.

## Run locally on Linux, WSL or macOS

Build and install the two source versions into separate environments, then run:

```sh
python tools/performance_gate.py --base-python /path/base/bin/python \
  --head-python /path/head/bin/python --output performance.json
```

The script writes JSON samples and a Markdown comparison. Its comparison logic
has tests that deliberately inject time, memory, quality and geometry regressions.

## Release suite and public history

The [benchmark history](benchmark-history.md) keeps immutable JSON measurements,
within-run timing plots and full quality tables for each archived run. It records
source commits, the benchmark-script hash, loaded native hashes, Python/NumPy
versions, host and CPU information, every measured sample, counts and failures.
Existing labels cannot be overwritten. A failed run remains marked **FAIL**;
missing samples never become a successful timing. A comparison of one source
commit against itself is marked **REFERENCE**, not a release-to-release speedup.

The automatic PR suite now also covers Schwarz P solid, a fixed six-mode periodic
Fourier spinodoid fixture and a fixed ordered smooth-union cylinder lattice.
The manual **release** suite adds adaptive Gyroid, spinodoid and lattice fixtures,
and a curved Tet10 Gyroid. These small deterministic regression fixtures are
separate from the full Microgen accuracy study. They do not certify geometry,
FEA convergence, or superiority over Microgen+MMG.

All cases retain minimum, maximum, mean, median, p05, p10 and p90 quality. All
seven statistics have the same 1% regression allowance. Node counts have the
same 1% change guard as tetrahedron counts. Tet10 quality and sampled surface
error describe the linear parent; its curved volume is independently integrated
and Jacobians are checked at the reference nodes and quadrature samples. These
are sampled validity checks, not a proof of global invertibility.

For every release:

1. Run **Performance** manually on the candidate commit. Choose `suite=release`,
   provide the previous compatible release tag in `baseline`, and use a unique
   `history_label`, such as `v0.2.0-rc.1`.
2. Inspect the CI summary and download `performance-comparison`. The artifact
   includes raw JSON and a generated documentation history. Failures must be
   explained or fixed; do not weaken thresholds to publish a passing result.
3. Archive the downloaded `performance.json` locally and commit the generated
   documentation with its evidence:

   ```sh
   python -m tools.benchmark_history --input /path/performance.json \
     --label v0.2.0-rc.1
   ```

4. For the final release record, run against the immutable release tag and archive
   under that version. Keep the candidate records as well. Documentation can be
   updated after the release without changing the measured source identity.

Version 0.1 cannot execute adaptive/Tet10 fixtures. For the first 0.2 release,
use the reviewed compatible reference commit
`afc49e583a0450454b02440557f39b89bbc39ca6` for the full suite. The automatic
seven-case PR suite still compares against the actual PR base. Later releases
can compare the complete suite against the previous 0.2-or-newer release.
The initial history entry is an M4 reference run of that approved commit on
both sides, not a published-release comparison.

Local full-suite runs use the same command with `--suite release`. Linux and
macOS normalize their different OS RSS units. Numerical workers and BLAS use
one thread. Run on an idle host; absolute timings from different hosts or fixture
versions must not be treated as a performance trend. The CI comparison remains
within one runner and one run, with alternating base/candidate order.

Timing and RSS gates use medians to tolerate host noise. Quality also checks
**every candidate repetition** against the base median, so an occasional bad
mesh cannot hide behind four good meshes. Each sample is capped at 250,000
nodes and tetrahedra. Changing fixture formulas or settings requires a new
fixture version and a reviewed comparison protocol, not a silent history reset.

To archive a local run, pass full commit SHAs with `--base-ref` and `--head-ref`;
the default `local` labels deliberately cannot enter the public history.
