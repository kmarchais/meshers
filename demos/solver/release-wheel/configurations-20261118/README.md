# Release-wheel measurements by configuration

Each row identifies the configuration observed at the start of a benchmark run. RAM is physical memory visible to the OS, in GiB. The JSON retains the exact byte count, CPU model, logical CPU count, architecture, platform version and UTC observation time. If hardware or software changes, save the new run under its own date and configuration. Earlier observations retain their original specifications.

| CPU | Logical CPUs | RAM (GiB) | OS / architecture | Tet4 solver records | Release meshing records |
| --- | ---: | ---: | --- | --- | --- |
| Apple M4 | 10 | 24.00 | macOS / arm64 | [9 replays](apple-m4-091916z-a3fb9327/tet4.json) | [33/33 samples](apple-m4-091916z-a3fb9327/release-fixtures.json) |
| Core Ultra 9 285K | 24 | 63.39 | Windows / AMD64 | [9 replays](core-ultra-9-285k-095358z-610e8c8f/tet4.json) | [33/33 samples](core-ultra-9-285k-095358z-610e8c8f/release-fixtures.json) |
| Core i7-10750H | 12 | 14.88 | Linux / x86_64 | [9 replays](core-i7-10750h-091913z-4f4ba086/tet4.json) | [33/33 samples](core-i7-10750h-091913z-4f4ba086/release-fixtures.json) |
| Core i7-12700H | 20 | 31.73 | Windows / AMD64 | [9 replays](core-i7-12700h-091916z-2ebf2eee/tet4.json) | [33/33 samples](core-i7-12700h-091916z-2ebf2eee/release-fixtures.json) |
| Ryzen 7 7730U | 16 | 30.25 | Linux / x86_64 | [9 replays](ryzen-7-7730u-091912z-fb87984d/tet4.json) | [33/33 samples](ryzen-7-7730u-091912z-fb87984d/release-fixtures.json) |

## What ran

- Each configuration ran uniform Gyroid, adaptive Gyroid and adaptive Schwarz P three times on Meshers 0.2.1. Each replay solved 12 systems. All 540 solves passed, with a maximum recomputed relative residual of 1.07e-12. Across paired repeats, the worst mixed-load energy difference between adaptive and uniform Gyroid ranged from 2.930% to 2.939%. See the [solver replay protocol](../README.md#replay).
- Each configuration also ran the 11 release fixtures from `tools/performance_gate.py`. Each case had one discarded warmup and three measured fresh processes, with one numerical thread requested. 165 of 165 measured samples passed the per-mesh checks. These are single-wheel observations. A regression verdict requires a paired base and candidate on the same configuration.

For a solver replay, use the frozen baseline and the two solver files at commit `3ca1f726`. Use one process, cycle the three cases in the table order from the [original replay](../fleet-20261118/README.md), and rotate the starting case after each repeat. The records include the input and native extension hashes and detected thread pools. For a meshing sample run `python tools/performance_gate.py --sample compiled_gyroid`, substituting each entry in `RELEASE_CASES`. The study used the gate at `e11810a7`. Windows reported the peak process working set through psutil 7.2.2. Linux and macOS reported peak resident memory through the OS. Use within-configuration timings for speed comparisons.
