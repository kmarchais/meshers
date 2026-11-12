# Benchmark history

Versioned measurements from fixed regression fixtures. Failed records remain visible. Release candidates and same-binary references are explicitly labelled.

Compare base/candidate measurements within each run. Different machines, suite versions or settings cannot establish an absolute performance trend. Raw JSON retains every measured sample and source/native identity.

| Record | Suite | Result | Platform |
| --- | --- | --- | --- |
| [v0.2.0-rc-reference](benchmarks/v0.2.0-rc-reference.md) | release | REFERENCE | macOS-26.5.1-arm64-arm-64bit / arm64 |
| [pr6-22cd611-quality-regression](benchmarks/pr6-22cd611-quality-regression.md) | pr | FAIL | Linux-7.0.0-1019-azure-x86_64-with-glibc2.43 / x86_64 |
| [pr6-0c3c4c1-quality-recovery](benchmarks/pr6-0c3c4c1-quality-recovery.md) | pr | PASS | macOS-26.5.1-arm64-arm-64bit / arm64 |
| [v0.2.0](benchmarks/v0.2.0.md) | release | PASS | Linux-7.0.0-1019-azure-x86_64-with-glibc2.43 / x86_64 |

See [regression protocol and release instructions](performance.md) for thresholds and how to add a release record.
