# Resource gate calibration

The gate compares five adjacent base/candidate pairs in fresh processes. It discards one warm-up pair and alternates which side runs first. Version 5 allows an additional 15% plus 10 ms for generation and 5% plus 4 MiB for peak RSS.

## Choosing the limits

The [five configuration snapshots](../../demos/solver/release-wheel/configurations-20261118/README.md) contain three repeats of each of the 11 release fixtures using one wheel. The widest individual time spread (maximum divided by minimum, minus one) was 14.76% on the Core i7-10750H. The widest RSS spread was 1.41% on the Core Ultra 9 285K. These ranges guided the choice of allowances. They do not estimate the false-alarm rate of a five-pair median.

The [earlier same-commit release comparison](../../docs/site/benchmarks/v0.2.0-rc-reference.json) supplied a paired check on Apple M4. Across its 11 fixtures, swapping the two sides changed median time by at most 5.7% and median RSS by less than 0.2%. We chose the limits before reading the fresh comparisons below. The absolute allowances cover short jobs and small working sets.

## Fresh comparisons

All three configurations passed all 11 cases with the proposed policy in either direction. Each cell below shows the median paired change in generation time / RSS, in percent. Positive values mean the candidate side was higher.

| Fixture | [AMD EPYC 7763](linux-epyc-7763/same-commit.json) | [AMD Ryzen 7 7730U](linux-ryzen-7-7730u/same-commit.json) | [Intel Core Ultra 9 285K](windows-core-ultra-9-285k/same-commit.json) |
| --- | ---: | ---: | ---: |
| compiled_gyroid | -0.1 / 0.0 | -0.7 / -0.2 | +0.1 / +0.2 |
| callback_gyroid | -0.1 / -0.1 | +0.2 / -0.1 | +0.4 / 0.0 |
| graded | +0.2 / -2.8 | -0.7 / -0.2 | +0.1 / +0.2 |
| mapped_intersection | -0.5 / 0.0 | +0.6 / 0.0 | +0.2 / +0.4 |
| schwarz_p | -0.7 / 0.0 | +1.3 / 0.0 | +0.3 / +0.2 |
| spinodoid | -0.3 / 0.0 | -1.4 / 0.0 | -1.7 / +0.3 |
| smooth_lattice | -0.3 / 0.0 | +1.5 / +0.4 | +0.2 / +0.1 |
| adaptive_gyroid | +2.1 / +1.0 | -2.1 / 0.0 | +0.2 / -0.1 |
| adaptive_spinodoid | +0.6 / +0.8 | +1.5 / 0.0 | +0.2 / -0.1 |
| adaptive_lattice | +0.1 / +1.3 | -7.7 / +0.4 | 0.0 / +0.1 |
| gyroid_tet10 | -0.2 / 0.0 | +3.0 / +0.1 | -1.3 / +0.2 |

The hosted EPYC runner built commit `fc54032b` twice. Its native binary hashes differ. The Ryzen system on Linux and Core Ultra 9 on Windows each used one installed wheel from `eb7cb9d4` for both sides. Each linked JSON records the CPU, usable RAM, OS, architecture, worker hashes and all sample metrics.

The Core i7-10750H and i7-12700H systems were offline. The collecting system had no configured route to the Mac. Their earlier samples are in the five-configuration archive.

## Known regression

On the Ryzen compiled Gyroid, the candidate slept for 427 ms and retained 12 MiB of additional memory during generation. Baseline time shifted from about 2.1 to 2.6 s during the run. The ratio of group medians showed a 12.4% time increase and missed the delay even with the tighter allowance. This prompted the switch to a paired comparison.

The median paired increases were 18.3% for time and 14.8% for RSS. Version 5 flagged both, with median excesses of 76 ms and 4.0 MiB above the allowances. Version 4 accepted both. All geometry and quality checks passed. [Injected runs](linux-ryzen-7-7730u/injection.json)

On a held-out Schwarz P mesh, we repeated the experiment after fixing the paired rule. The added delay was 173 ms. Median paired time rose 18.2% and RSS 15.4%. Version 5 rejected both resources and passed every mesh check. Version 4 accepted the run. [Held-out measurements](linux-ryzen-7-7730u/schwarz-injection.json)

On the held-out Windows configuration, the compiled Gyroid control added 184 ms and retained 12 MiB. Median paired generation time rose 20.6% and peak working set rose 21.3%. Version 5 flagged both resources. Version 4 accepted the run. [Windows measurements](windows-core-ultra-9-285k/injection.json)

With five pairs, three must exceed `base × factor + allowance` to fail. The EPYC run was held out when the paired rule was chosen. These controls establish sensitivity on these configurations. The number of independent runners is too small to estimate a CI-wide false-alarm rate.

## Reproduce

For the same-wheel comparison, install the wheel once and pass its interpreter for both sides. The hosted comparison instead builds the same commit twice using the [Performance workflow](https://github.com/kmarchais/meshers/actions/runs/44893152648).

```sh
python tools/performance_gate.py --suite release --repeats 5 \
  --base-python "$WHEEL_PYTHON" --head-python "$WHEEL_PYTHON" \
  --base-ref "$WHEEL_COMMIT" --head-ref "$WHEEL_COMMIT" --output same-commit.json
```

For the positive control, copy the version 4 driver and worker from commit `fc54032b`. In the candidate copy, insert these lines immediately after `start = time.perf_counter()` in `_measure_generation`. Substitute the recorded `delay_seconds` from the injection JSON. Keeping the array alive through generation makes its pages resident when peak memory is measured.

```python
time.sleep(delay_seconds)
calibration_pages = np.ones(12 * 1024**2, dtype=np.uint8)
```

Run each copy with the same interpreter and `tools/performance_gate.py --sample compiled_gyroid` or `--sample schwarz_p` for the held-out case. Discard one base/candidate warm-up pair. Measure five pairs starting with candidate/base and alternating pair order. The delay is 20% of the earlier same-wheel base median. The replay test recomputes medians and checks both versions of the gate against the measurements.

```sh
python -m pytest -q tools/tests/test_performance_calibration.py
```
