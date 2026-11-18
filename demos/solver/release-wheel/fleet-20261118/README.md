# Cross-platform Tet4 replay, 18 November 2026

All four hosts ran the same three cases with the released Meshers 0.2.1 wheel. Each host completed three repeats per case and 12 linear systems per repeat. The largest recomputed relative residual across the 432 systems was 1.07e-12.

Each timing entry gives median generation seconds / median homogenization seconds across the three repeats. Homogenization includes the six prescribed strains, six stress-controlled checks and solver audits. The last column gives the worst mixed-load strain-energy difference between that host's adaptive and uniform Gyroid solutions, summarized by the median over paired repeats.

| Raw records | CPU | Uniform Gyroid (s) | Adaptive Gyroid (s) | Adaptive Schwarz P (s) | Gyroid energy difference |
| --- | --- | ---: | ---: | ---: | ---: |
| [Runs](intel-i7-12700h.json) | Core i7-12700H | 1.91 / 41.54 | 3.74 / 17.21 | 2.53 / 13.27 | 2.93861% |
| [Runs](ryzen-7-7730u.json) | Ryzen 7 7730U | 1.05 / 17.07 | 1.79 / 6.83 | 1.14 / 5.42 | 2.93861% |
| [Runs](apple-m4.json) | Apple M4 | 0.42 / 10.37 | 0.81 / 4.35 | 0.55 / 3.45 | 2.93033% |
| [Runs](intel-i7-10750h.json) | Core i7-10750H | 1.16 / 27.39 | 2.29 / 10.90 | 1.33 / 8.90 | 2.93861% |

The x86-64 hosts reproduced the element counts from the original Windows run. On Apple M4, uniform Gyroid also had 72,757 tetrahedra. Adaptive Gyroid had 34,854 versus 34,863 on the original Windows run, and adaptive Schwarz P had 35,994 versus 35,989. Comparing the same recipe's stiffness tensor against the original Windows run, the largest observed mixed-load energy discrepancy was 0.0103%. This measures cross-platform reproducibility. Continuum error needs an independent reference. These runs exercised the PR #9 solver checks. PR #10 reports a separate two-sided geometry audit on the original Windows run.

## Replay conditions

Use the [pinned environment and three release-wheel commands](../README.md#replay). We transferred the two Python driver files and frozen baseline from PR #9, installed native platform wheels in isolated Python 3.12 environments, and used one shared queue slot per host. Each host used one Python process for all nine runs. Repeat 1 used uniform Gyroid, adaptive Gyroid, Schwarz P. Each subsequent repeat rotated the first case to the end. No warmup was discarded. The original Windows baseline has one observation.

The runs requested one numerical thread through the mesher and the BLAS, OpenMP and Accelerate environment settings. The JSON files retain package versions, input and extension hashes, platform, detected thread pools, actual order, all per-load checks and individual timings. Threadpoolctl 3.6.0 detected no pool on macOS, so the Accelerate thread count was not independently observed. Core i7-12700H ran on AC with its Balanced plan, the Core i7-10750H on AC with the powersave governor, and the Apple M4 on AC with low-power mode off. Timing includes small progress-reporting overhead. Compare mesh recipes within each configuration. OS, libraries and power management also vary across the systems.
