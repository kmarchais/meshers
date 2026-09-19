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

Cases cover a compiled periodic gyroid, a callback periodic gyroid, a graded
sheet and a mapped named intersection. Each uses 24 lattice cells and two
optimization passes. Timing includes field preparation and meshing, but excludes
imports, independent validation and export. Peak resident memory is the Linux
process high-water mark through mesh generation, including imports and JIT.

| Metric | Fails when |
| --- | --- |
| Median elapsed time | More than 25% plus 20 ms above base |
| Median peak RSS | More than 15% plus 8 MiB above base |
| Minimum, fifth-percentile or median MMG quality | More than 1% below base |
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

## Run locally on Linux or WSL

Build and install the two source versions into separate environments, then run:

```sh
python tools/performance_gate.py --base-python /path/base/bin/python \
  --head-python /path/head/bin/python --output performance.json
```

The script writes JSON samples and a Markdown comparison. Its comparison logic
has tests that deliberately inject time, memory, quality and geometry regressions.
