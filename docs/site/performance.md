# Performance and regression checks

Both Schwarz P meshes meet the 1% stiffness target. The native Tet4 route took less time on all five configurations. [Comparison protocol and raw runs](https://github.com/kmarchais/meshers/blob/main/demos/solver/release-wheel/schwarz-matched-cost-20261118/README.md)

![Median time for three qualifying runs per route on each CPU and OS. Native Tet4 was faster than the repaired route on all five configurations for this Schwarz P specimen.](assets/schwarz-cost.svg)

Median of three runs per setting. Stiffness discrepancy was 0.992% for native Tet4 and 0.927–0.928% after Microgen repair.

## Pull request gates

On one Linux runner with matching Python and NumPy, discard one warm-up pair. Alternate base and candidate across five fresh-process pairs. Time preparation and meshing. Peak memory includes imports and JIT.

| Metric | Regression limit against base |
| --- | --- |
| Median time | +25% plus 20 ms |
| Median peak memory | +15% plus 8 MiB |
| MMG quality (min, max, mean, median, p05, p10, p90) | 1% decrease for each statistic. Each candidate repeat is also checked against the base median. |
| Sampled wall error and mapped cap error | +5% plus 10⁻¹⁰ for each |
| Volume | ±0.5% |
| Nodes and tetrahedra | ±1% each |

- **Seven cases** — compiled and callback gyroids, graded sheet, mapped intersection, Schwarz P, Fourier spinodoid and ordered lattice. The fixtures use 24 cells and two optimization passes.
- **Absolute gates** — finite coordinates, positive volumes, closed boundary edge incidence, matching requested periodic faces, minimum quality ≥ 0.01 and sampled error ≤ 0.1. Cap each sample at 250,000 nodes and 250,000 tetrahedra.
- **Geometry comparisons** — fixture version 4 measures walls and mapped caps separately on both wheels. The [implicit coverage study](https://github.com/kmarchais/meshers/blob/main/demos/solver/implicit-coverage/README.md) also records admitted and rejected geometries.

## Release suite and public history

The release suite adds adaptive gyroid, spinodoid, lattice and curved Tet10. For Tet10, quality and surface error describe the parent. Volume integration and Jacobian checks use the curved element. Subdivision must resolve inconclusive determinant bounds.

1. Run **Performance** on the candidate with `suite=release`, a compatible release `baseline` and a unique `history_label`. For the first 0.2 suite the compatible reference was `afc49e583a0450454b02440557f39b89bbc39ca6`.
2. Inspect the `performance-comparison` artifact, including failed or missing samples. Keep the same fixtures and thresholds when investigating a regression.
3. Archive the JSON under the candidate label. Repeat on the immutable final tag and retain both records.

```sh
python -m tools.benchmark_history --input /path/performance.json --label v0.2.0-rc.1
```

## Run a local comparison

Use separate environments, one numerical thread and an idle host. On Windows, both workers load the `psutil` memory probe before timing. Its memory is included in the process peak. Memory is peak working set on Windows and peak RSS on Linux/macOS, in bytes.

```sh
python tools/performance_gate.py --base-python /path/base/bin/python \
  --head-python /path/head/bin/python --output performance.json
```

Use `--suite release` for all fixtures. Archive with full `--base-ref` and `--head-ref` SHAs. Each new record fingerprints the driver (`script_sha256`) and the numerical worker (`worker_sha256`). Compare within one configuration and fixture version. Bump the fixture version when formulas or settings change. [Benchmark history](benchmark-history.md)

## Candidate wheel FEA gate

Python 3.11 CI reopens the Schwarz P VTKHDF export before solving.

- **Geometry** — two-sided sampled wall/cap distance ≤ 0.002, volume within 1%, MMG quality ≥ 0.1. Check topology, periodic faces and material probes.
- **Stiffness** — all six strain and six stress solves must pass. Worst relative strain-energy discrepancy ≤ 1% against the reference tensor.

```sh
uv run --with fedoo==1.0.0b2 python demos/solver/schwarz_workflow.py native-candidate candidate-acceptance
```
